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
from pprint import pformat, pprint
import time
import uuid

from KBaseReport.KBaseReportClient import KBaseReport
from KBaseReport.baseclient import ServerError as _RepError
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from kb_quast.kb_quastClient import kb_quast
from kb_quast.baseclient import ServerError as QUASTError
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
    GIT_URL = "https://github.com/kbaseapps/kb_Velvet"
    GIT_COMMIT_HASH = "6efca462c5b4e221751c6d448d386829d8e776f1"

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block
    VELVETH = '/kb/module/velvet/velveth'
    VELVETG = '/kb/module/velvet/velvetg'
    VELVET_DATA = '/kb/module/velvet/data/'
    PARAM_IN_WS = 'workspace_name'
    PARAM_IN_CS_NAME = 'output_contigset_name'

    def log(self, message, prefix_newline=False):
            print(('\n' if prefix_newline else '') +
                          str(time.time()) + ': ' + str(message))

    def process_params_h(self, params):
        if self.PARAM_IN_WS not in params:
            raise ValueError('a string reprsenting workspace_name parameter is required')
        if 'out_folder' not in params:
            raise ValueError('a string reprsenting out_folder parameter is required')
        if 'hash_length' not in params:
            raise ValueError('an integerreprsenting  hash_length parameter is required')
        if 'hash_length' in params:
            if not isinstance(params['hash_length'], int):
                raise ValueError('hash_length must be of type int')
        if 'reads_channels' not in params:
            raise ValueError('key-val pairs representing the reads_channels parameter is required')
        if type(params['reads_channels']) != list:
            raise ValueError('reads_channels must be a list')
        for rc in params['reads_channels']:
            self.log('Read file channel info :\n' + pformat(rc))
            if 'read_file_info' not in rc:
                raise ValueError('a read_channel must have a read_file_info dictionary')
            else:
                if 'read_file' not in rc['read_file_info'] or rc['read_file_info']['read_file'] == '':
                    raise ValueError('a non-blank read file name is required')
            if 'read_type' not in rc:
                raise ValueError('a read_type is required')
            if 'file_format' not in rc:
                raise ValueError('a file_format is required')

    def process_params_g(self, params):
        if self.PARAM_IN_WS not in params:
            raise ValueError('a string reprsenting workspace_name parameter is required')
        if 'wk_folder' not in params:
            raise ValueError('a string reprsenting wk_folder parameter is required')
        if self.PARAM_IN_CS_NAME not in params:
                raise ValueError('output_contigset_name parameter is required')
        rm_dir = os.path.join(self.scratch, params['wk_folder'] + '/Roadmaps')
        sq_dir = os.path.join(self.scratch, params['wk_folder'] + '/Sequences')
        if not os.path.exists(rm_dir) or not os.path.exists(sq_dir):
            raise ValueError('no valid subfolders named %s and %s in the working directory for running velvetg' % (rm_dir, sq_dir))

    def construct_velveth_cmd(self, params):
        # STEP 1: get the reads channels as reads file info
        out_folder = params['out_folder']
        hash_length = params['hash_length']
        wsname = params[self.PARAM_IN_WS]
        reads_channels = params['reads_channels']

        # STEP 2: construct the command for run_velveth
        vh_cmd = [self.VELVETH]

        # set the output location
        output_dir = os.path.join(self.scratch, out_folder)
        vh_cmd.append(output_dir)
        vh_cmd.append(str(hash_length))

        for rc in reads_channels:
            vh_cmd.append('-' + rc['file_format'])
            read_type = rc['read_type']
            vh_cmd.append('-' + read_type)
            if 'read_reference' in rc and rc['read_reference'] == 1:
                vh_cmd.append(self.VELVET_DATA + rc['read_file_info']['reference_file'])

            if 'file_layout' in rc and rc['file_layout'] == 'separate':
                vh_cmd.append('-' + rc['file_layout'])
                vh_cmd.append(self.VELVET_DATA + rc['read_file_info']['left_file'])
                vh_cmd.append(self.VELVET_DATA + rc['read_file_info']['right_file'])
            else:
                vh_cmd.append(self.VELVET_DATA + rc['read_file_info']['read_file'])

        # STEP 3 return vh_cmd
        return wsname, vh_cmd

    def construct_velvetg_cmd(self, params):
        # STEP 1: get the working folder housing the velveth results as well as the reads info
        wk_folder = params['wk_folder']
        wsname = params[self.PARAM_IN_WS]
        # set the output location
        work_dir = os.path.join(self.scratch, wk_folder)

        # STEP 2: construct the command for run_velvetg
        vg_cmd = [self.VELVETG]
        vg_cmd.append(work_dir)
        #appending the standard optional inputs
        if 'cov_cutoff' in params:
            vg_cmd.append('-cov_cutoff ' + str(params['cov_cutoff']))
        if 'ins_length' in params:
            vg_cmd.append('-cov_cutoff ' + str(params['cov_cutoff']))
        if 'ins_length' in params:
            vg_cmd.append('-ins_length ' + str(params['ins_length']))
        if 'read_trkg' in params:
            vg_cmd.append('-read_trkg ' + str(params['read_trkg']))
        if 'min_contig_length' in params:
            vg_cmd.append('-min_contig_lgth ' + str(params['min_contig_length']))
        if 'amos_file' in params:
            vg_cmd.append('-amos_file ' + str(params['amos_file']))
        if 'exp_cov' in params:
            vg_cmd.append('-exp_cov ' + str(params['exp_cov']))
        if 'long_cov_cutoff' in params:
            vg_cmd.append('-long_cov_cutoff ' + str(params['long_cov_cutoff']))

        # appending the advanced optional inputs--TODO

        # STEP 3 return vg_cmd
        return wsname, vg_cmd

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

    def generate_report(self, input_file_name, params, wsname):
        self.log('Generating and saving report')

        fasta_stats = self.load_stats(input_file_name)
        lengths = [fasta_stats[contig_id] for contig_id in fasta_stats]

        assembly_ref = params[self.PARAM_IN_WS] + '/' + params[self.PARAM_IN_CS_NAME]

        report = ''
        report += 'Velvet results saved to: ' + wsname + '/' + params['wk_folder'] + '\n'
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


    def run_velveth(self, ctx, params):
        """
        Definition of run_velveth
        :param params: instance of type "VelvethParams" (Arguments for
           velveth input string workspace_name - the name of the workspace
           for input/output string out_folder - the folder name for output
           files int hash_length - EITHER an odd integer (if even, it will be
           decremented) <= 31 (if above, will be reduced)L) -> structure:
           parameter "out_folder" of String, parameter "workspace_name" of
           String, parameter "hash_length" of Long, parameter
           "reads_channels" of list of type "ReadsChannel" (Define a
           structure that mimics the concept of "channel" used by the Velvet
           program. string read_type - the read type, e.g., -short,
           -shortPaired, short2, shortPaired2, -long, or -longPaired string
           file_format - the format of the input file, e.g., -fasta, -fastq,
           -raw,-fasta.gz, -fastq.gz, -raw.gz, -sam, -bam, -fmtAuto string
           read_file_info - the hash that holds the details about the read
           file string file_layout - the layout of the file, e.g.,
           -interleaved or -separate bool read_reference - indicating if a
           reference file is used) -> structure: parameter "read_type" of
           String, parameter "file_format" of String, parameter
           "read_file_info" of type "ReadFileInfo" (Define a structure that
           holds the read file name and its use. Note: only read_file_name is
           required, the rest are optional. e.g., {"reference_file" =>
           "test_reference.fa", "read_file_name" => "mySortedReads.sam",
           "left_file" => "left.fa", "right_file" => "right.fa"}) ->
           structure: parameter "read_file" of String, parameter
           "reference_file" of String, parameter "left_file" of String,
           parameter "right_file" of String, parameter "file_layout" of
           String, parameter "read_reference" of type "bool" (A boolean - 0
           for false, 1 for true. @range (0, 1))
        :returns: instance of Long
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_velveth
        self.log('Running run_velveth with params:\n' + pformat(params))

        token = ctx['token']

        # STEP 1: basic parameter checks + parsing
        self.process_params_h(params)

        # STEP 2: construct the command for run_velveth
        wsname, velveth_cmd = self.construct_velveth_cmd(params)

        # STEP 3: run velveth
        self.log('Running velveth with command:\n' + pformat(velveth_cmd))
        #p = subprocess.Popen(velveth_cmd, cwd=self.scratch, shell=False)
        p = subprocess.Popen(velveth_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()

        self.log('Return code: ' + str(retcode))
        if p.returncode != 0:
            raise ValueError('Error running VELVETH, return code: ' + str(retcode) + '\n')

        output = p.returncode 

        #END run_velveth

        # At some point might do deeper type checking...
        if not isinstance(output, int):
            raise ValueError('Method run_velveth return value ' +
                             'output is not type int as required.')
        # return the results
        return [output]

    def run_velvetg(self, ctx, params):
        """
        Definition of run_velvetg
        :param params: instance of type "VelvetgParams" (Arguments for
           run_velvetg string workspace_name - the name of the workspace from
           which to take input and store output. string wk_folder - the name
           of the folder where the velvet results are created and saved
           output_contigset_name - the name of the output contigset
           list<paired_end_lib> float cov_cutoff - the removal of low
           coverage nodes AFTER tour bus or allow the system to infer it
           (default: no removal) int ins_length - expected distance between
           two paired end reads (default: no read pairing) int read_trkg; - 
           (1=yes|0=no) tracking of short read positions in assembly
           (default:0) int min_contig_length - minimum contig length exported
           to contigs.fa file (default: hash length * 2) int amos_file -
           (1=yes|0=no) #export assembly to AMOS file (default: 0) float
           exp_cov - <floating point|auto>, expected coverage of unique
           regions or allow the system to infer it (default: no long or
           paired-end read resolution) float long_cov_cutoff - removal of
           nodes with low long-read coverage AFTER tour bus(default: no
           removal)) -> structure: parameter "workspace_name" of String,
           parameter "wk_folder" of String, parameter "output_contigset_name"
           of String, parameter "cov_cutoff" of Double, parameter
           "ins_length" of Long, parameter "read_trkg" of Long, parameter
           "min_contig_length" of Long, parameter "amos_file" of Long,
           parameter "exp_cov" of Double, parameter "long_cov_cutoff" of
           Double
        :returns: instance of Long
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_velvetg
        self.log('Running run_velvetg with params:\n' + pformat(params))

        token = ctx['token']

        # STEP 1: basic parameter checks + parsing
        self.process_params_g(params)

        # STEP 2: construct the command for run_velvetg
        wsname, velvetg_cmd = self.construct_velvetg_cmd(params)

        # STEP 3: run velvetg
        self.log('running velvetg with command:\n' + pformat(velvetg_cmd))
        p = subprocess.Popen(velvetg_cmd, cwd=self.scratch, shell=False)
        retcode = p.wait()

        self.log('Return code: ' + str(retcode))
        if p.returncode != 0:
            raise ValueError('Error running VELVETG, return code: ' + str(retcode) + '\n')

        output = p.returncode 
        #END run_velvetg

        # At some point might do deeper type checking...
        if not isinstance(output, int):
            raise ValueError('Method run_velvetg return value ' +
                             'output is not type int as required.')
        # return the results
        return [output]

    def run_velvet(self, ctx, params):
        """
        Definition of run_velvet
        :param params: instance of type "VelvetParams" (Arguments for
           run_velvet) -> structure: parameter "h_params" of type
           "VelvethParams" (Arguments for velveth input string workspace_name
           - the name of the workspace for input/output string out_folder -
           the folder name for output files int hash_length - EITHER an odd
           integer (if even, it will be decremented) <= 31 (if above, will be
           reduced)L) -> structure: parameter "out_folder" of String,
           parameter "workspace_name" of String, parameter "hash_length" of
           Long, parameter "reads_channels" of list of type "ReadsChannel"
           (Define a structure that mimics the concept of "channel" used by
           the Velvet program. string read_type - the read type, e.g.,
           -short, -shortPaired, short2, shortPaired2, -long, or -longPaired
           string file_format - the format of the input file, e.g., -fasta,
           -fastq, -raw,-fasta.gz, -fastq.gz, -raw.gz, -sam, -bam, -fmtAuto
           string read_file_info - the hash that holds the details about the
           read file string file_layout - the layout of the file, e.g.,
           -interleaved or -separate bool read_reference - indicating if a
           reference file is used) -> structure: parameter "read_type" of
           String, parameter "file_format" of String, parameter
           "read_file_info" of type "ReadFileInfo" (Define a structure that
           holds the read file name and its use. Note: only read_file_name is
           required, the rest are optional. e.g., {"reference_file" =>
           "test_reference.fa", "read_file_name" => "mySortedReads.sam",
           "left_file" => "left.fa", "right_file" => "right.fa"}) ->
           structure: parameter "read_file" of String, parameter
           "reference_file" of String, parameter "left_file" of String,
           parameter "right_file" of String, parameter "file_layout" of
           String, parameter "read_reference" of type "bool" (A boolean - 0
           for false, 1 for true. @range (0, 1)), parameter "g_params" of
           type "VelvetgParams" (Arguments for run_velvetg string
           workspace_name - the name of the workspace from which to take
           input and store output. string wk_folder - the name of the folder
           where the velvet results are created and saved
           output_contigset_name - the name of the output contigset
           list<paired_end_lib> float cov_cutoff - the removal of low
           coverage nodes AFTER tour bus or allow the system to infer it
           (default: no removal) int ins_length - expected distance between
           two paired end reads (default: no read pairing) int read_trkg; - 
           (1=yes|0=no) tracking of short read positions in assembly
           (default:0) int min_contig_length - minimum contig length exported
           to contigs.fa file (default: hash length * 2) int amos_file -
           (1=yes|0=no) #export assembly to AMOS file (default: 0) float
           exp_cov - <floating point|auto>, expected coverage of unique
           regions or allow the system to infer it (default: no long or
           paired-end read resolution) float long_cov_cutoff - removal of
           nodes with low long-read coverage AFTER tour bus(default: no
           removal)) -> structure: parameter "workspace_name" of String,
           parameter "wk_folder" of String, parameter "output_contigset_name"
           of String, parameter "cov_cutoff" of Double, parameter
           "ins_length" of Long, parameter "read_trkg" of Long, parameter
           "min_contig_length" of Long, parameter "amos_file" of Long,
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
        wsname = params['g_params']['workspace_name']

        # STEP 1: run velveth and velvetg sequentially
        ret = []
        try:
            ret = self.run_velveth(ctx, params['h_params'])
            while( ret[0] != 0 ):
                time.sleep(1)
        except ValueError as eh:
            self.log('Velveth raised error:\n')
            print(eh)
        else:#no exception raised by Velveth and Velveth returns 0, then run Velvetg
            try:
                ret = self.run_velvetg(ctx, params['g_params'])
                while( ret[0] != 0 ):
                    time.sleep(1)
            except ValueError as eg:
                self.log('Velvetg raised error:\n')
                print(eg)
            else:#no exception raised by Velvetg and Velvetg returns 0, then move to STEP 3  
                ret[0] = 0

        # STEP 2: parse the output and save back to KBase, create report in the same time
        if( ret[0] == 0 ):
                work_dir = os.path.join(self.scratch, params['g_params']['wk_folder'])
                self.log('Velvet output folder: ' + work_dir)

                output_contigs = os.path.join(work_dir, 'contigs.fa')

                self.log('Uploading FASTA file to Assembly')

                assemblyUtil = AssemblyUtil(self.callbackURL, token=ctx['token'], service_ver='release')

                min_contig_length = (params['g_params']).get('min_contig_length', 0)
                if min_contig_length > 0:
                        assemblyUtil.save_assembly_from_fasta(
                                {'file': {'path': output_contigs},
                                'workspace_name': wsname,
                                'assembly_name': params['g_params'][self.PARAM_IN_CS_NAME],
                                'min_contig_length': params['g_params']['min_contig_length']
                                })
                else:
                        assemblyUtil.save_assembly_from_fasta(
                        {'file': {'path': output_contigs},
                        'workspace_name': wsname,
                        'assembly_name': params['g_params'][self.PARAM_IN_CS_NAME]
                        })
                # generate report from contigs.fa
                report_name, report_ref = self.generate_report(output_contigs, params['g_params'], wsname)

                # STEP 3: contruct the output to send back
                output = {'report_name': report_name, 'report_ref': report_ref}
        else:
            output = {'report_name': 'Velvet failed', 'report_ref': null}

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
