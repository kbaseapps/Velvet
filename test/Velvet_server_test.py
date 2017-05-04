# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401
from biokbase.AbstractHandle.Client import AbstractHandle as HandleService

#from biokbase.workspace.client import Workspace as workspaceService
from Velvet.VelvetImpl import Velvet
from Velvet.VelvetServer import MethodContext
from Velvet.authclient import KBaseAuth as _KBaseAuth
from Workspace.WorkspaceClient import Workspace as workspaceService


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

    def velveth(self):
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
        print('RESULT from velveth:\n')
        pprint(result)
        return result

        # check the output


    def velvetg(self, work_folder, output_contigset):
        # run velvetg
        params = {
            'workspace_name': self.getWsName(),
            'wk_folder': work_folder,
            'output_contigset_name': output_contigset 
            #'min_contig_length': 20
        }

        result = self.getImpl().run_velvetg(self.getContext(), params)
        print('RESULT from velvetg:\n')
        pprint(result)

        # check the report. We assume that kb_quast and KBaseReport do what they're supposed to do
        rep = self.wsClient.get_objects2({'objects': [{'ref': result[0]['report_ref']}]})['data'][0]
        print('REPORT object:')
        pprint(rep)

        self.assertEqual(rep['info'][1].rsplit('_', 1)[0], 'kb_velvet_report')
        self.assertEqual(rep['info'][2].split('-', 1)[0], 'KBaseReport.Report')


    def test_run_velvet(self):
        # run velveth
        work_folder = self.velveth()
        work_folder = 'velvet_outfolder'#self.velveth()
       
        # run velvetg
        if(not work_folder == ""):
            self.velvetg(work_folder, 'test_contigset')
        else:
            print('velvet failed!')
