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

    /* The workspace object name of a SingleEndLibrary or PairedEndLibrary file, whether of the
       KBaseAssembly or KBaseFile type.
    */
    typedef string read_lib;

    /* 
        Arguments for run_velvet

        string workspace_name - the name of the workspace from which to take input and store output.
        int hash_length - an odd integer (if even, it will be decremented) <= 31
        string output_contigset_name - the name of the output contigset
        list<paired_end_lib> read_libraries - Illumina PairedEndLibrary files to assemble
        min_contig_length - integer to filter out contigs with length < min_contig_length
                     from the Velvet output. Default value is 500 (where 0 implies no filter).

        @optional min_contig_length
        @optional cov_cutoff
        @optional ins_length
        @optional read_trkg
        @optional amos_file
        @optional exp_cov
        @optional long_cov_cutoff
     */

    typedef structure {
        string workspace_name;
        int hash_length;
        list<read_lib> read_libraries; 
        string output_contigset_name;
 
        int min_contig_length;
        float cov_cutoff;
        int ins_length;
        bool read_trkg;
        bool amos_file;
        float exp_cov;
        float long_cov_cutoff;
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
