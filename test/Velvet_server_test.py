# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import os.path
import json  # noqa: F401
import time
import requests
import shutil

from pprint import pprint, pformat
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
        #forward_data_file = '../work/small.forward.fq'
        forward_data_file = 'data/small.forward.fq'
        #forward_data_file = '../work/GW456A_trim_reads_unpaired_rev.single.fastq'
        forward_file = os.path.join(shared_dir, os.path.basename(forward_data_file))
        shutil.copy(forward_data_file, forward_file)
        #reverse_data_file = '../work/small.reverse.fq'
        reverse_data_file = 'data/small.reverse.fq'
        reverse_file = os.path.join(shared_dir, os.path.basename(reverse_data_file))
        shutil.copy(reverse_data_file, reverse_file)

        ru = ReadsUtils(os.environ['SDK_CALLBACK_URL'])
        paired_end_ref = ru.upload_reads({'fwd_file': forward_file,
                                          #'rev_file': reverse_file,
                                          'sequencing_tech': 'artificial reads',
                                          'interleaved': 0, 'wsname': self.getWsName(),
                                          'name': 'test.pe.reads'})['obj_ref']

        new_obj_info = self.wsClient.get_object_info_new({'objects': [{'ref': paired_end_ref}]})
        self.__class__.pairedEndLibInfo = new_obj_info[0]
        print ('paired reads uploaded:\n')
        pprint (pformat(new_obj_info))

        return new_obj_info[0]

    @classmethod
    def make_ref(self, object_info):
        return str(object_info[6]) + '/' + str(object_info[0]) + \
            '/' + str(object_info[4])

    # Uncomment to skip this test
    @unittest.skip("skipped test_run_velveth")
    def test_velveth(self):
        # get the test data
        out_folder = os.path.join(self.scratch, 'velvet_output_dir')
        if not os.path.exists(out_folder):
            os.makedirs(out_folder)

        rc1 = {
            'read_type': 'long',
            'file_format': 'fastq.gz',
            'file_layout': 'interleaved',
            'read_file_info' : {
                                'read_file_name': 'ecoli_ref-5m-trim.fastq.gz'
                               }
        }
        rc2 = {
            'read_type': 'longPaired',
            'file_format': 'fasta.gz',
            'file_layout': 'interleaved',
            'read_file_info' : {
                                'read_file_name': 'ecoli-reads-5m-dn-paired.fa.gz'
                               }
        }
        rc3 = {
            'read_type': 'shortPaired',
            'file_format': 'fastq',
            'file_layout': 'separate',
            'read_file_info' : {
                                'read_file_name': 'small.reverse.fq',
                                'left_file': 'small.forward.fq',
                                'right_file': 'small.reverse.fq',
                               }
        }

        pe_lib_info = self.getPairedEndLibInfo()
        print(pe_lib_info)

        obj_ids = [{'ref':pe_lib_info[7] + '/' + pe_lib_info[1]}]
 
        ws_info = self.wsClient.get_object_info_new({'objects': obj_ids})
        reads_params = []

        reftoname = {}
        for wsi, oid in zip(ws_info, obj_ids):
            ref = oid['ref']
            reads_params.append(ref)
            obj_name = wsi[1]
            reftoname[ref] = wsi[7] + '/' + obj_name

        readcli = ReadsUtils(self.callback_url, token=self.token)

        typeerr = ('Supported types: KBaseFile.SingleEndLibrary ' +
                   'KBaseFile.PairedEndLibrary ' +
                   'KBaseAssembly.SingleEndLibrary ' +
                   'KBaseAssembly.PairedEndLibrary')
        try:
            reads = readcli.download_reads({'read_libraries': reads_params,
                                            'interleaved': 'false',
                                            'gzipped': None
                                            })['files']
        except ServerError as se:
            print('logging stacktrace from dynamic client error')
            print(se.data)
            if typeerr in se.message:
                prefix = se.message.split('.')[0]
                raise ValueError(
                    prefix + '. Only the types ' +
                    'KBaseAssembly.PairedEndLibrary ' +
                    'and KBaseFile.PairedEndLibrary are supported')
            else:
                raise

        print ('Got reads data from converter:\n' + pformat(reads))

        reads_data = []
        for ref in reads:
            reads_name = reftoname[ref]
            f = reads[ref]['files']
            seq_tech = reads[ref]["sequencing_tech"]
            if f['type'] == 'interleaved':
                reads_data.append({'fwd_file': f['fwd'], 'type':'interleaved',
                                   'seq_tech': seq_tech})
            elif f['type'] == 'paired':
                reads_data.append({'fwd_file': f['fwd'], 'rev_file': f['rev'],
                                   'type':'separated', 'seq_tech': seq_tech})
            elif f['type'] == 'single':
                reads_data.append({'fwd_file': f['fwd'], 'type':'single',
                                   'seq_tech': seq_tech})
            else:
                raise ValueError('Something is very wrong with read lib' + reads_name)

        params = {
            'workspace_name': pe_lib_info[7],
            'out_folder': out_folder,
            'hash_length': 21,
            'reads_channels': [rc1, rc2, rc3]#tests passed
            #'reads_files': reads_data
        }

        result = self.getImpl().exec_velveth(params)
        self.assertTrue(os.path.isfile(os.path.join(self.scratch, params['out_folder'] + '/Roadmaps')))
        self.assertTrue(os.path.isfile(os.path.join(self.scratch, params['out_folder'] + '/Sequences')))
        print('RESULT from velveth is saved in:\n' + os.path.join(self.scratch, params['out_folder']))
        pprint('Returned value by Velveth is: ' + str(result))
        return result

    # Uncomment to skip this test
    @unittest.skip("skipped test_run_velvetg")
    def test_velvetg(self):
        # run velvetg
        #work_folder = self.velveth()[0]
        #print "Returned work folder from velveth call: " + work_folder
        params = {
            'workspace_name': self.getWsName(),
            'output_contigset_name': 'test_contigset', 
            'min_contig_length': 500,
            'cov_cutoff': 5.2
        }

        result = self.getImpl().run_velvetg(self.getContext(), params)
        self.assertTrue(os.path.isfile(os.path.join(self.scratch, params['wk_folder'] + '/LastGraph')))
        self.assertTrue(os.path.isfile(os.path.join(self.scratch, params['wk_folder'] + '/Log')))
        print('RESULT from velvetg is saved in:\n' + os.path.join(self.scratch, params['wk_folder']))
        pprint('Returned value by Velvetg is: ' + str(result))
        return result

    # Uncomment to skip this test
    # HIDE @unittest.skip("skipped test_run_velvet")
    def test_run_velvet(self):
        # get the test data
        pe_lib_info = self.getPairedEndLibInfo()
        pprint(pe_lib_info)

        # velvet parameters
        params = {
            'workspace_name': self.getWsName(),
            'output_contigset_name': 'Velvet_test_contigset',
            'hash_length': 21,
            'read_libraries':[self.make_ref(pe_lib_info)],
            'min_contig_length': 300,
            'cov_cutoff': 5.2,
            'read_trkg': '',
            'amos_file': 'yes',
            'exp_cov': 21.3,
            'ins_length': 400
        }

        result = self.getImpl().run_velvet(self.getContext(), params)

        if not result[0]['report_ref'] is None:
                rep = self.wsClient.get_objects2({'objects': [{'ref': result[0]['report_ref']}]})['data'][0]
                print('REPORT object:')
                pprint(rep)

                self.assertEqual(rep['info'][1].rsplit('_', 1)[0], 'kb_velvet_report')
                self.assertEqual(rep['info'][2].split('-', 1)[0], 'KBaseReport.Report')
        else:
                print('Velvet failed!')
