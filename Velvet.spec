/*

   Name of module: Velvet

   This is a KBase module that wraps the open source package "Short read de novo assembler using de Bruijn graphs"
   Velvet_1.2.10

   References:
   https://github.com/dzerbino/velvet
   https://github.com/dzerbino/velvet/blob/master/Columbus_manual.pdf

*/

module Velvet {
    /* 
        A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int bool;

    /* The workspace object name of a read library file of the KBaseFile type.
    */
    typedef string seq_file_name;

    /* 
        Arguments for run_velveth
     */

    /*
        Define a structure that holds the read file name and its use.
        Note: only read_file_name is required, the rest are optional.
        e.g., 
        {"reference_file" => "test_reference.fa", "read_file_name" => "mySortedReads.sam", 
        "left_file" => "left.fa", "right_file" => "right.fa"}
    */ 
    typedef structure { 
        string read_file;
        string reference_file;
        string left_file;
        string right_file;
    } ReadFileInfo;

    /* 
        Define a structure that mimics the concept of "channel" used by the Velvet program.
        string read_type - the read type, e.g., -short, -shortPaired, short2, shortPaired2, -long, or -longPaired
        string file_format - the format of the input file, e.g., -fasta, -fastq, -raw,-fasta.gz, -fastq.gz, -raw.gz, -sam, -bam, -fmtAuto
        string read_file_info - the hash that holds the details about the read file
        string file_layout - the layout of the file, e.g., -interleaved or -separate 
        bool read_reference - indicating if a reference file is used
    */
    typedef structure {
        string read_type; 
        string file_format; 
        ReadFileInfo read_file_info;
        string file_layout; 
        bool read_reference; 
    } ReadsChannel;

    /* 
        Arguments for velveth input
        string workspace_name - the name of the workspace for input/output
        string out_folder - the folder name for output files
        int hash_length - an odd integer (if even, it will be decremented) <= 31
        list<seq_file_name> sequence_files - sequence files to assemble
        list<ReadsChannel> reads_channels - a list/an array of ReadsChannel defining {read_type, file_format, {read_file[,...]}[, file_layout, read_reference]}
   */
    typedef structure {
        string out_folder; 
        string workspace_name;
        int hash_length; 
        list<seq_file_name> sequence_files;
        list<ReadsChannel> reads_channels; 
    } VelvethParams;
    
    /* 
        Definition of run_velveth
     */
    funcdef run_velveth(VelvethParams params) returns (int output) authentication required;


    /* 
        Arguments for run_velvetg

        string workspace_name - the name of the workspace from which to take input and store output.
        string wk_folder - the name of the folder where the velvet results are created and saved
        string output_contigset_name - the name of the output contigset
        float cov_cutoff - the removal of low coverage nodes AFTER tour bus or allow the system to infer it (default: no removal)
        int ins_length - expected distance between two paired end reads (default: no read pairing)
        int read_trkg; -  (1=yes|0=no) tracking of short read positions in assembly (default:0)
        int min_contig_length - minimum contig length exported to contigs.fa file (default: hash length * 2)
        int amos_file - (1=yes|0=no) #export assembly to AMOS file (default: 0)
        float exp_cov - <floating point|auto>, expected coverage of unique regions or allow the system to infer it (default: no long or paired-end read resolution)
        float long_cov_cutoff - removal of nodes with low long-read coverage AFTER tour bus(default: no removal)
     */

    typedef structure {
        string workspace_name;
        string wk_folder;
        string output_contigset_name; 
        float cov_cutoff; 
        int ins_length; 
        int read_trkg; 
        int min_contig_length; 
        int amos_file; 
        float exp_cov; 
        float long_cov_cutoff; 
    } VelvetgParams;
   
    /* 
        Definition of run_velvetg
     */
    funcdef run_velvetg(VelvetgParams params) returns (int output) authentication required;


    /*
        Arguments for run_velvet
     */
    typedef structure {
        VelvethParams h_params;
        VelvetgParams g_params;
    } VelvetParams;
   
    /* Output parameter items for run_velvet

    report_name - the name of the KBaseReport.Report workspace object.
    report_ref - the workspace reference of the report.

    */
    typedef structure {
        string report_name;
        string report_ref;
    } VelvetResults;
    
    /* 
        Definition of run_velvet
     */
    funcdef run_velvet(VelvetParams params) returns (VelvetResults output) authentication required;

};
