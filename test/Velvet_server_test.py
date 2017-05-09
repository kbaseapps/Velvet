# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import os.path
import json  # noqa: F401
import time
import requests
import Queue

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

from Velvet.VelvetImpl import Velvet
from Velvet.VelvetServer import MethodContext
from Velvet.authclient import KBaseAuth as _KBaseAuth
from Workspace.WorkspaceClient import Workspace as workspaceService
from ReadsUtils.baseclient import ServerError
from ReadsUtils.ReadsUtilsClient import ReadsUtils

class VelvetTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        cls.token = token
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('Velvet'):
            print(nameval[0] + '=' + nameval[1])
            cls.cfg[nameval[0]] = nameval[1]

        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)

        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'Velvet',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL, token=token)
        cls.serviceImpl = Velvet(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

        cls.shockURL = cls.cfg['shock-url']
        cls.handleURL = cls.cfg['handle-service-url']
        cls.queue = Queue.Queue()

        cls.readUtilsImpl = ReadsUtils(cls.callback_url, token=cls.token)
        cls.staged = {}
        cls.nodes_to_delete = []
        cls.handles_to_delete = []
        #cls.setupTestData()

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # Helper script borrowed from the transform service, logger removed
    @classmethod
    def upload_file_to_shock(cls, file_path):
        """
        Use HTTP multi-part POST to save a file to a SHOCK instance.
        """

        header = dict()
        header["Authorization"] = "Oauth {0}".format(cls.token)

        if file_path is None:
            raise Exception("No file given for upload to SHOCK!")

        with open(os.path.abspath(file_path), 'rb') as dataFile:
            files = {'upload': dataFile}
            print('POSTing data')
            response = requests.post(
                cls.shockURL + '/node', headers=header, files=files,
                stream=True, allow_redirects=True)
            print('got response')

        if not response.ok:
            response.raise_for_status()

        result = response.json()

        if result['error']:
            raise Exception(result['error'][0])
        else:
            return result["data"]

    @classmethod
    def upload_file_to_shock_and_get_handle(cls, test_file):
        '''
        Uploads the file in test_file to shock and returns the node and a
        handle to the node.
        '''
        print('loading file to shock: ' + test_file)
        node = cls.upload_file_to_shock(test_file)
        pprint(node)
        cls.nodes_to_delete.append(node['id'])

        print('creating handle for shock id ' + node['id'])
        handle_id = cls.hs.persist_handle({'id': node['id'],
                                           'type': 'shock',
                                           'url': cls.shockURL
                                           })
        cls.handles_to_delete.append(handle_id)

        md5 = node['file']['checksum']['md5']
        return node['id'], handle_id, md5, node['file']['size']

    @classmethod
    def upload_reads(cls, wsobjname, object_body, fwd_reads,
                     rev_reads=None, single_end=False, sequencing_tech='Illumina',
                     single_genome='1'):

        ob = dict(object_body)  # copy
        ob['sequencing_tech'] = sequencing_tech
#        ob['single_genome'] = single_genome
        ob['wsname'] = cls.getWsName()
        ob['name'] = wsobjname
        if single_end or rev_reads:
            ob['interleaved']= 0
        else:
            ob['interleaved']= 1
        print('\n===============staging data for object ' + wsobjname +
              '================')
        print('uploading forward reads file ' + fwd_reads['file'])
        fwd_id, fwd_handle_id, fwd_md5, fwd_size = \
            cls.upload_file_to_shock_and_get_handle(fwd_reads['file'])

        ob['fwd_id']= fwd_id
        rev_id = None
        rev_handle_id = None
        if rev_reads:
            print('uploading reverse reads file ' + rev_reads['file'])
            rev_id, rev_handle_id, rev_md5, rev_size = \
                cls.upload_file_to_shock_and_get_handle(rev_reads['file'])
            ob['rev_id']= rev_id
        obj_ref = cls.readUtilsImpl.upload_reads(ob)
        objdata = cls.wsClient.get_object_info_new({
            'objects': [{'ref': obj_ref['obj_ref']}]
            })[0]
        cls.staged[wsobjname] = {'info': objdata,
                                 'ref': cls.make_ref(objdata),
                                 'fwd_node_id': fwd_id,
                                 'rev_node_id': rev_id,
                                 'fwd_handle_id': fwd_handle_id,
                                 'rev_handle_id': rev_handle_id
                                 }

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_Velvet_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    # Uncomment to skip this test
    @unittest.skip("skipped test_run_velveth")
    def test_run_velveth(self):
        # run velveth
        rc = {
            'read_type': 'short',
            'file_format': 'sam',
            'file_layout': 'interleaved',
            'read_file_info' : {
                                'read_file': 'test_reads.sam',
                                'reference_file': 'test_reference.fa',
                                'left_file': '',
                                'right_file': ''
                               }
        }
        params = {
            'workspace_name': self.getWsName(),
            'out_folder': 'velvet_outfolder',
            'hash_length': 21,
            'reads_channels': [rc]
        }

        result = self.getImpl().run_velveth(self.getContext(), params)
        self.assertTrue(os.path.isfile(os.path.join(self.scratch, params['out_folder'] + '/Roadmaps')))
        self.assertTrue(os.path.isfile(os.path.join(self.scratch, params['out_folder'] + '/Sequences')))
        print('RESULT from velveth is saved in:\n' + os.path.join(self.scratch, params['out_folder']))
        pprint('Returned value by Velveth is: ' + str(result))
        return result

        # check the output

    # Uncomment to skip this test
    @unittest.skip("skipped test_run_velvetg")
    def test_run_velvetg(self):
        # run velvetg
        #work_folder = self.velveth()[0]
        #print "Returned work folder from velveth call: " + work_folder
        params = {
            'workspace_name': self.getWsName(),
            'wk_folder': 'velvet_outfolder',
            'output_contigset_name': 'test_contigset', 
            'min_contig_length': 100,
            'cov_cutoff': 5.2
        }

        result = self.getImpl().run_velvetg(self.getContext(), params)
        self.assertTrue(os.path.isfile(os.path.join(self.scratch, params['wk_folder'] + '/LastGraph')))
        self.assertTrue(os.path.isfile(os.path.join(self.scratch, params['wk_folder'] + '/Log')))
        print('RESULT from velvetg is saved in:\n' + os.path.join(self.scratch, params['wk_folder']))
        pprint('Returned value by Velvetg is: ' + str(result))
        return result

    # Uncomment to skip this test
    #@unittest.skip("skipped test_run_velvet")
    def test_run_velvet(self):
        # velveth parameters
        rc = {
            'read_type': 'short',
            'file_format': 'sam',
            'file_layout': 'interleaved',
            'read_file_info' : {
                                'read_file': 'test_reads.sam',
                                'reference_file': 'test_reference.fa',
                                'left_file': '',
                                'right_file': ''
                               }
        }
        h_params = {
            'workspace_name': self.getWsName(),
            'out_folder': 'velvet_outfolder',
            'hash_length': 21,
            'reads_channels': [rc]
        }

        # velvetg parameters
        g_params = {
            'workspace_name': self.getWsName(),
            'wk_folder': 'velvet_outfolder',
            'output_contigset_name': 'test_contigset'#, 
            #'cov_cutoff': 5.2
            #'min_contig_length': 100#,
        }

        params = {'h_params': h_params, 'g_params': g_params}

        result = self.getImpl().run_velvet(self.getContext(), params)

        # check the report. We assume that kb_quast and KBaseReport do what they're supposed to do
        rep = self.wsClient.get_objects2({'objects': [{'ref': result[0]['report_ref']}]})['data'][0]
        print('REPORT object:')
        pprint(rep)

        self.assertEqual(rep['info'][1].rsplit('_', 1)[0], 'kb_velvet_report')
        self.assertEqual(rep['info'][2].split('-', 1)[0], 'KBaseReport.Report')
