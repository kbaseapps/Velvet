/*

        A KBase module: Velvet

*/

module Velvet {
    /* 
        A boolean - 0 for false, 1 for true.
        @range (0, 1)
    */
    typedef int bool;


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
    */
    typedef structure {
        string out_folder; 
        string workspace_name;
        int hash_length; 
        list<ReadsChannel> reads_channels; 
    } VelvethParams;
    
    /* 
        Definition of run_velveth
     */
    funcdef run_velveth(VelvethParams params) returns (string output) authentication required;


    /* 
        Arguments for run_velvetg
     */

    typedef structure {
        string workspace_name;
        string wk_folder; 
        float cov_cutoff; 
        int ins_length; 
        int read_trkg; 
        int min_contig_length; 
        int amos_file; 
        float exp_cov; 
        float long_cov_cutoff; 
    } VelvetgParams;
   
    /* Output parameter items for run_velvetg

    report_name - the name of the KBaseReport.Report workspace object.
    report_ref - the workspace reference of the report.

    */
    typedef structure {
        string report_name;
        string report_ref;
    } VelvetResults;
    
    /* 
        Definition of run_velvetg
     */
    funcdef run_velvetg(VelvetgParams params) returns (VelvetResults output) authentication required;

};
