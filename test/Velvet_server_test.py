# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import os.path
import json  # noqa: F401
import time
import requests
import shutil

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

    # borrowed from Megahit - call this method to get the WS object info of a Paired End Library (will
    # upload the example data if this is the first time the method is called during tests)
    def getPairedEndLibInfo(self):
        if hasattr(self.__class__, 'pairedEndLibInfo'):
            return self.__class__.pairedEndLibInfo
        # 1) upload files to shock
        shared_dir = "/kb/module/work/tmp"
        forward_data_file = 'data/small.forward.fq'
        forward_file = os.path.join(shared_dir, os.path.basename(forward_data_file))
        shutil.copy(forward_data_file, forward_file)
        reverse_data_file = 'data/small.reverse.fq'
        reverse_file = os.path.join(shared_dir, os.path.basename(reverse_data_file))
        shutil.copy(reverse_data_file, reverse_file)

        ru = ReadsUtils(os.environ['SDK_CALLBACK_URL'])
        paired_end_ref = ru.upload_reads({'fwd_file': forward_file, 'rev_file': reverse_file,
                                          'sequencing_tech': 'artificial reads',
                                          'interleaved': 0, 'wsname': self.getWsName(),
                                          'name': 'test.pe.reads'})['obj_ref']

        new_obj_info = self.ws.get_object_info_new({'objects': [{'ref': paired_end_ref}]})
        self.__class__.pairedEndLibInfo = new_obj_info[0]
        print ('paired reads uploaded:\n')
        pprint (pformat(new_obj_info))

        return new_obj_info[0]

    @classmethod
    def make_ref(self, object_info):
        return str(object_info[6]) + '/' + str(object_info[0]) + \
            '/' + str(object_info[4])

    # Uncomment to skip this test
    #@unittest.skip("skipped test_run_velveth")
    def test_run_velveth(self):
        # get the test data
        pe_lib_info = self.getPairedEndLibInfo()
        pprint(pe_lib_info)

        # run velveth
        rc1 = {
            'read_type': 'short',
            'file_format': 'sam',
            'file_layout': 'interleaved',
            'read_file_info' : {
                                'read_file_name': 'test_reads.sam',
                                'reference_file': 'test_reference.fa',
                                'left_file': '',
                                'right_file': ''
                               }
        }
        rc = {
            'read_type': 'short',
            'file_format': 'fastq',
            'file_layout': 'separated',
            'read_file_info' : {
                                'read_library_ref': pe_lib_info[7] + '/' + pe_lib_info[1],
                                'read_file': '',
                                'reference_file': '',
                                'left_file': '',
                                'right_file': ''
                               }
        }
        params = {
            'workspace_name': pe_lib_info[7],
            #'workspace_name': self.getWsName(),
            'out_folder': 'velvet_outfolder',
            'hash_length': 21,
            'reads_channels': [rc1, rc]
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
    @unittest.skip("skipped test_run_velvet")
    def test_run_velvet(self):
        # velveth parameters
        rc = {
            'read_type': 'short',
            'file_format': 'sam',
            'file_layout': 'interleaved',
            'read_file_info' : {
                                'read_file_name': 'test_reads.sam',
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
            'output_contigset_name': 'Velvet_test_contigset'#, 
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
