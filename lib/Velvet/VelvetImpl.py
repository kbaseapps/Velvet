# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
import re
import shutil
import subprocess
import numpy as np
from Bio import SeqIO
from pprint import pprint, pformat
from datetime import datetime
import time
import uuid

from KBaseReport.KBaseReportClient import KBaseReport
from KBaseReport.baseclient import ServerError as _RepError
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from kb_quast.kb_quastClient import kb_quast
from kb_quast.baseclient import ServerError as QUASTError
from ReadsUtils.ReadsUtilsClient import ReadsUtils 
from ReadsUtils.baseclient import ServerError
from Workspace.WorkspaceClient import Workspace as workspaceService
#END_HEADER


class Velvet:
    '''
    Module Name:
    Velvet

    Module Description:
    Name of module: Velvet

   This is a KBase module that wraps the open source package "Short read de novo assembler using de Bruijn graphs"
   Velvet_1.2.10

   References:
   https://github.com/dzerbino/velvet
   https://github.com/dzerbino/velvet/blob/master/Columbus_manual.pdf
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/kbaseapps/kb_Velvet/"
    GIT_COMMIT_HASH = "6bf7071a438f1284b6169522760ff06dad8c987e"

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block
    VELVETH = '/kb/module/velvet/velveth'
    VELVETG = '/kb/module/velvet/velvetg'
    VELVET_DATA = '/kb/module/work/tmp'
    #VELVET_DATA = '/kb/module/test/data'
    PARAM_IN_WS = 'workspace_name'
    PARAM_IN_CS_NAME = 'output_contigset_name'
    PARAM_IN_MIN_CONTIG_LENGTH = 'min_contig_length'
    PARAM_IN_HASH_LENGTH = 'hash_length'

    INVALID_WS_OBJ_NAME_RE = re.compile('[^\\w\\|._-]')
    INVALID_WS_NAME_RE = re.compile('[^\\w:._-]')

    PARAM_IN_LIB = 'read_libraries'
    PARAM_IN_SEQ = 'sequence_files'


    def log(self, message, prefix_newline=False):
            print(('\n' if prefix_newline else '') +
                str(time.time()) + ': ' + str(message))

    def process_params(self, params):
        if (self.PARAM_IN_WS not in params or
                not params[self.PARAM_IN_WS]):
            raise ValueError(self.PARAM_IN_WS + ' parameter is required')
        if self.PARAM_IN_LIB not in params:
            raise ValueError(self.PARAM_IN_LIB + ' parameter is required')
        if type(params[self.PARAM_IN_LIB]) != list:
            raise ValueError(self.PARAM_IN_LIB + ' must be a list')
        if not params[self.PARAM_IN_LIB]:
            raise ValueError('At least one reads library must be provided')
        if self.PARAM_IN_HASH_LENGTH not in params or params[self.PARAM_IN_HASH_LENGTH] is None:
            raise ValueError(self.PARAM_IN_HASH_LENGTH + ' parameter is required')
        if self.PARAM_IN_HASH_LENGTH in params:
            if not isinstance(params[self.PARAM_IN_HASH_LENGTH], int):
                raise ValueError(self.PARAM_IN_HASH_LENGTH + ' must be of type int')
        if (self.PARAM_IN_CS_NAME not in params or
                not params[self.PARAM_IN_CS_NAME]):
            raise ValueError(self.PARAM_IN_CS_NAME + ' parameter is required')
        if self.INVALID_WS_OBJ_NAME_RE.search(params[self.PARAM_IN_CS_NAME]):
            raise ValueError('Invalid workspace object name ' +
                             params[self.PARAM_IN_CS_NAME])
        if self.PARAM_IN_MIN_CONTIG_LENGTH in params:
            if not isinstance(params[self.PARAM_IN_MIN_CONTIG_LENGTH], int):
                raise ValueError(self.PARAM_IN_MIN_CONTIG_LENGTH + ' must be of type int')

    def construct_velveth_cmd(self, params):
        if params.get('reads_channels', None) is not None:
                reads_channels = params['reads_channels']
        else:
                reads_channels = []
        # STEP 1: fetch the reads files and build the reads channel
        if params.get('reads_files', None) is not None:
            #print('Input reads files:' + pformat(params['reads_files']))
            for reads in params['reads_files']:
                ftype = reads['type']
                fwd = reads['fwd_file']
                pprint('forward: ' + fwd)
                fext = (os.path.splitext(fwd)[1]).replace(".", "")
                if ftype == 'single':
                    file_info = {'read_file_name': fwd}
                    reads_channels.append({
                        'read_type': 'short',
                        'file_format': fext,
                        'read_file_info': file_info,
                        'file_layout': ''
                    })
                elif ftype == 'paired':
                    rev = reads['rev_file']
                    file_info = {
                        'read_file_name': fwd,
                        'left_file': fwd,
                        'right_file': rev
                    }
                    reads_channels.append({
                        'read_type': 'shortPaired',
                        'file_format': fext,
                        'read_file_info': file_info,
                        'file_layout': 'separate'
                    })
                elif ftype == 'interleaved':
                    file_info = {
                        'read_file_name': fwd
                    }
                    reads_channels.append({
                        'read_type': 'shortPaired',
                        'file_format': fext,
                        'read_file_info': file_info,
                        'file_layout': ''
                        })

        # STEP 2: build the reads channels from the sequence files
        if 'sequence_files' in params:
                sq_files = ' '.join(params['sequence_files'])
                if( sq_files != ''):
                        fext = (os.path.splitext(params['sequence_files'][0])[1]).replace(".", "")
                        file_info = {
                                'read_file_name': sq_files
                        }
                        reads_channels.append({
                                'read_type': 'short',
                                'file_format': fext,
                                'read_file_info': file_info
                        })

        # STEP 3: construct the command for running velveth
        out_folder = params['out_folder']
        hash_length = params[self.PARAM_IN_HASH_LENGTH]
        wsname = params[self.PARAM_IN_WS]
        vh_cmd = [self.VELVETH]
        vh_cmd.append(out_folder)
        vh_cmd.append(str(hash_length))

        for rc in reads_channels:
            vh_cmd.append('-' + rc['file_format'])
            vh_cmd.append('-' + rc['read_type'])
            if 'read_reference' in rc and rc['read_reference'] == 1:
                vh_cmd.append(os.path.join(self.VELVET_DATA, rc['read_file_info']['reference_file']))

            if 'file_layout' in rc and rc['file_layout'] == 'separate':
                vh_cmd.append('-' + rc['file_layout'])
                vh_cmd.append(os.path.join(self.VELVET_DATA, rc['read_file_info']['left_file']))
                vh_cmd.append(os.path.join(self.VELVET_DATA, rc['read_file_info']['right_file']))
            else:
                vh_cmd.append(os.path.join(self.VELVET_DATA, rc['read_file_info']['read_file_name']))

        # STEP 3 return vh_cmd
        print ('Velveth CMD:')
        print ' '.join(vh_cmd)
        return vh_cmd

    def construct_velvetg_cmd(self, params):
        # STEP 1: get the working folder housing the velveth results as well as the reads info
        out_folder = params['out_folder']
        wsname = params[self.PARAM_IN_WS]

        # STEP 2: construct the command for running velvetg
        vg_cmd = [self.VELVETG]
        vg_cmd.append(out_folder)
        #appending the standard optional inputs
        if (params.get(self.PARAM_IN_MIN_CONTIG_LENGTH, None) is not None and
                isinstance(params[self.PARAM_IN_MIN_CONTIG_LENGTH], int)):
            vg_cmd.append('-min_contig_lgth')
            vg_cmd.append(str(params[self.PARAM_IN_MIN_CONTIG_LENGTH]))
        if (params.get('cov_cutoff', None) is not None and
                isinstance(params['cov_cutoff'], float)):
            vg_cmd.append('-cov_cutoff')
            vg_cmd.append(str(params['cov_cutoff']).lower())
        if (params.get('ins_length', None) is not None and
                isinstance(params['ins_length'], int)):
            vg_cmd.append('-ins_length')
            vg_cmd.append(str(params['ins_length']))
        if params.get('read_trkg', None) is not None:
            vg_cmd.append('-read_trkg')
            vg_cmd.append('yes' if (params['read_trkg']==1 or str(params['read_trkg']).lower()=='yes') else 'no')
        if params.get('amos_file', None) is not None:
            vg_cmd.append('-amos_file')
            vg_cmd.append('yes' if (params['amos_file']==1  or str(params['amos_file']).lower()=='yes') else 'no')
        if params.get('exp_cov', None) is not None:
            vg_cmd.append('-exp_cov')
            if isinstance(params['exp_cov'], float):
                vg_cmd.append(str(params['exp_cov']).lower())
            else:
                vg_cmd.append('auto')
        if (params.get('long_cov_cutoff', None) is not None and 
                isinstance(params['long_cov_cutoff'], float)):
            vg_cmd.append('-long_cov_cutoff')
            vg_cmd.append(str(params['long_cov_cutoff']))

        # appending the advanced optional inputs--TODO

        # STEP 3 return vg_cmd
        print ('Velvetg CMD:')
        print ' '.join(vg_cmd)
        return vg_cmd

    def exec_velveth(self, params):
        self.log('Running run_velveth with params:\n' + pformat(params))
        velveth_cmd = self.construct_velveth_cmd(params)

        p = subprocess.Popen(velveth_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()

        self.log('Return code: ' + str(retcode))
        if p.returncode != 0:
            raise ValueError('Error running VELVETH, return code: ' + str(retcode) + '\n')

        return p.returncode

    def exec_velvetg(self, params):
        self.log('Running run_velvetg with params:\n' + pformat(params))
        velvetg_cmd = self.construct_velvetg_cmd(params)
        p = subprocess.Popen(velvetg_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()
        self.log('Return code: ' + str(p.returncode))
        if p.returncode != 0:
            raise ValueError('Error running VELVETG, return code: ' + str(p.returncode) + '\n')

        return p.returncode


    def exec_velvet(self, params, reads_data):
        outdir = os.path.join(self.scratch, 'velvet_output_dir')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        tmpdir = os.path.join(self.scratch, 'velvet_tmp_dir')
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)

        # build the parameters
        params_h = {
                'workspace_name': params[self.PARAM_IN_WS],
                'hash_length': params[self.PARAM_IN_HASH_LENGTH],
                'reads_files': reads_data,
                'out_folder': outdir
        }
        params_g = {
                'workspace_name': params[self.PARAM_IN_WS],
                'output_contigset_name': params[self.PARAM_IN_CS_NAME],
                'out_folder': outdir
        }
        if self.PARAM_IN_MIN_CONTIG_LENGTH in params and not (params[self.PARAM_IN_MIN_CONTIG_LENGTH] is None):
            params_g[self.PARAM_IN_MIN_CONTIG_LENGTH] = params.get(self.PARAM_IN_MIN_CONTIG_LENGTH, 1)
        if 'cov_cutoff' in params and not (params['cov_cutoff'] is None):
            params_g['cov_cutoff'] = params['cov_cutoff']
        if 'ins_length' in params and not (params['ins_length'] is None) :
            params_g['ins_length'] = params['ins_length']
        if 'read_trkg' in params:
            params_g['read_trkg'] = params['read_trkg']
        if 'amos_file' in params:
            params_g['amos_file'] = params['amos_file']
        if 'exp_cov' in params and not (params['exp_cov'] is None):
            params_g['exp_cov'] = params['exp_cov']
        if 'long_cov_cutoff' in params and not (params['long_cov_cutoff'] is None):
            params_g['long_cov_cutoff'] = params['long_cov_cutoff']

        ret = 1
        try:
            ret = self.exec_velveth(params_h)
            while( ret != 0 ):
                time.sleep(1)
        except ValueError as eh:
            self.log('Velveth raised error:\n')
            print(eh)
        else:#no exception raised by Velveth and Velveth returns 0, then run Velvetg
            ret = 1 
            try:
                ret = self.exec_velvetg(params_g)
                while( ret != 0 ):
                    time.sleep(1)
            except ValueError as eg:
                self.log('Velvetg raised error:\n')
                print(eg)
            else:#no exception raised by Velvetg and Velvetg returns 0, then move to saving and reporting  
                ret = outdir
        return ret

    # adapted from
    # https://github.com/kbaseapps/kb_SPAdes/blob/master/lib/kb_SPAdes/kb_SPAdesImpl.py
    # which was adapted from
    # https://github.com/kbase/transform/blob/master/plugins/scripts/convert/trns_transform_KBaseFile_AssemblyFile_to_KBaseGenomes_ContigSet.py
    def load_stats(self, input_file_name):
        self.log('Starting conversion of FASTA to KBaseGenomeAnnotations.Assembly')
        self.log('Building Object.')
        if not os.path.isfile(input_file_name):
            raise Exception('The input file name {0} is not a file!'.format(input_file_name))
        with open(input_file_name, 'r') as input_file_handle:
            contig_id = None
            sequence_len = 0
            fasta_dict = dict()
            first_header_found = False
            # Pattern for replacing white space
            pattern = re.compile(r'\s+')
            for current_line in input_file_handle:
                if (current_line[0] == '>'):
                    # found a header line
                    # Wrap up previous fasta sequence
                    if not first_header_found:
                        first_header_found = True
                    else:
                        fasta_dict[contig_id] = sequence_len
                        sequence_len = 0
                    fasta_header = current_line.replace('>', '').strip()
                    try:
                        contig_id = fasta_header.strip().split(' ', 1)[0]
                    except:
                        contig_id = fasta_header.strip()
                else:
                    sequence_len += len(re.sub(pattern, '', current_line))
        # wrap up last fasta sequence
        if not first_header_found:
            raise Exception("There are no contigs in this file")
        else:
            fasta_dict[contig_id] = sequence_len
        return fasta_dict

    def generate_report(self, input_file_name, params, out_folder, wsname):
        self.log('Generating and saving report')

        fasta_stats = self.load_stats(input_file_name)
        lengths = [fasta_stats[contig_id] for contig_id in fasta_stats]

        assembly_ref = params[self.PARAM_IN_WS] + '/' + params[self.PARAM_IN_CS_NAME]

        report = ''
        report += 'Velvet results saved to: ' + wsname + '/' + out_folder + '\n'
        report += 'Assembly saved to: ' + assembly_ref + '\n'
        report += 'Assembled into ' + str(len(lengths)) + ' contigs.\n'
        report += 'Avg Length: ' + str(sum(lengths) / float(len(lengths))) + ' bp.\n'

        # compute a simple contig length distribution
        bins = 10
        counts, edges = np.histogram(lengths, bins)  # @UndefinedVariable
        report += 'Contig Length Distribution (# of contigs -- min to max ' + 'basepairs):\n'
        for c in range(bins):
            report += '   ' + str(counts[c]) + '\t--\t' + str(edges[c]) + ' to ' + str(edges[c + 1]) + ' bp\n'
        print('Running QUAST')
        kbq = kb_quast(self.callbackURL)
        quastret = kbq.run_QUAST({'files': [{'path': input_file_name,
                                             'label': params[self.PARAM_IN_CS_NAME]}]})
        print('Saving report')
        kbr = KBaseReport(self.callbackURL)
        report_info = kbr.create_extended_report(
            {'message': report,
             'objects_created': [{'ref': assembly_ref, 'description': 'Assembled contigs'}],
             'direct_html_link_index': 0,
             'html_links': [{'shock_id': quastret['shock_id'],
                             'name': 'report.html',
                             'label': 'QUAST report'}
                            ],
             'report_object_name': 'kb_velvet_report_' + str(uuid.uuid4()),
             'workspace_name': params[self.PARAM_IN_WS]
            })
        reportName = report_info['name']
        reportRef = report_info['ref']
        return reportName, reportRef

    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR

        # Any configuration parameters that are important should be parsed and
        # saved in the constructor.
        self.cfg = config
        self.callbackURL = os.environ['SDK_CALLBACK_URL']
        self.log('Callback URL: ' + self.callbackURL)
        self.workspaceURL = config['workspace-url']
        self.scratch = os.path.abspath(config['scratch'])
        if not os.path.exists(self.scratch):
            os.makedirs(self.scratch)

        #END_CONSTRUCTOR
        pass


    def run_velvet(self, ctx, params):
        """
        Definition of run_velvet
        :param params: instance of type "VelvetParams" (Arguments for
           run_velvet string workspace_name - the name of the workspace from
           which to take input and store output. int hash_length - an odd
           integer (if even, it will be decremented) <= 31 string
           output_contigset_name - the name of the output contigset
           list<paired_end_lib> read_libraries - Illumina PairedEndLibrary
           files to assemble min_contig_length - integer to filter out
           contigs with length < min_contig_length from the Velvet output.
           Default value is 500 (where 0 implies no filter). @optional
           min_contig_length @optional cov_cutoff @optional ins_length
           @optional read_trkg @optional amos_file @optional exp_cov
           @optional long_cov_cutoff) -> structure: parameter
           "workspace_name" of String, parameter "hash_length" of Long,
           parameter "read_libraries" of list of type "read_lib" (The
           workspace object name of a SingleEndLibrary or PairedEndLibrary
           file, whether of the KBaseAssembly or KBaseFile type.), parameter
           "output_contigset_name" of String, parameter "min_contig_length"
           of Long, parameter "cov_cutoff" of Double, parameter "ins_length"
           of Long, parameter "read_trkg" of type "bool" (A boolean - 0 for
           false, 1 for true. @range (0, 1)), parameter "amos_file" of type
           "bool" (A boolean - 0 for false, 1 for true. @range (0, 1)),
           parameter "exp_cov" of Double, parameter "long_cov_cutoff" of
           Double
        :returns: instance of type "VelvetResults" (Output parameter items
           for run_velvet report_name - the name of the KBaseReport.Report
           workspace object. report_ref - the workspace reference of the
           report.) -> structure: parameter "report_name" of String,
           parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_velvet
        self.log('Running run_velvet with params:\n' + pformat(params))

        token = ctx['token']
        wsname = params[self.PARAM_IN_WS]
        self.process_params(params)

        # STEP 0: preprocess the reads in KBase way
        obj_ids = []
        for r in params[self.PARAM_IN_LIB]:
            obj_ids.append({'ref': r if '/' in r else (wsname + '/' + r)})
        ws = workspaceService(self.workspaceURL, token=token)
        ws_info = ws.get_object_info_new({'objects': obj_ids})
        reads_params = []

        reftoname = {}
        for wsi, oid in zip(ws_info, obj_ids):
            ref = oid['ref']
            reads_params.append(ref)
            obj_name = wsi[1]
            reftoname[ref] = wsi[7] + '/' + obj_name

        readcli = ReadsUtils(self.callbackURL, token=ctx['token'])

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
            self.log('logging stacktrace from dynamic client error')
            self.log(se.data)
            if typeerr in se.message:
                prefix = se.message.split('.')[0]
                raise ValueError(
                    prefix + '. Only the types ' +
                    'KBaseAssembly.SingleEndLibrary ' +
                    'KBaseAssembly.PairedEndLibrary ' +
                    'KBaseFile.SingleEndLibrary ' +
                    'and KBaseFile.PairedEndLibrary are supported')
            else:
                raise

        self.log('Got reads data from converter:\n' + pformat(reads))

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
                                   'type':'paired', 'seq_tech': seq_tech})
            elif f['type'] == 'single':
                reads_data.append({'fwd_file': f['fwd'], 'type':'single',
                                   'seq_tech': seq_tech})
            else:
                raise ValueError('Something is very wrong with read lib' + reads_name)

        # STEP 1: run velveth and velvetg sequentially
        velvet_out = self.exec_velvet(params, reads_data)
        self.log('Velvet final return: ' + str(velvet_out))

        # STEP 2: parse the output and save back to KBase, create report in the same time
        if isinstance(velvet_out, str) and velvet_out != '':
                output_contigs = os.path.join(velvet_out, 'contigs.fa')

                self.log('Uploading FASTA file to Assembly')

                assemblyUtil = AssemblyUtil(self.callbackURL, token=ctx['token'], service_ver='release')

                min_contig_len = params.get(self.PARAM_IN_MIN_CONTIG_LENGTH, 0)
                if min_contig_len > 0:
                        assemblyUtil.save_assembly_from_fasta(
                                {'file': {'path': output_contigs},
                                'workspace_name': wsname,
                                'assembly_name': params[self.PARAM_IN_CS_NAME],
                                'min_contig_length': min_contig_len
                                })
                else:
                        assemblyUtil.save_assembly_from_fasta(
                        {'file': {'path': output_contigs},
                        'workspace_name': wsname,
                        'assembly_name': params[self.PARAM_IN_CS_NAME]
                        })
                # generate report from contigs.fa
                report_name, report_ref = self.generate_report(output_contigs, params, velvet_out, wsname)

                # STEP 3: contruct the output to send back
                output = {'report_name': report_name, 'report_ref': report_ref}
        else:
            output = {'report_name': 'Velvet failed', 'report_ref': None}

        #END run_velvet

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_velvet return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
