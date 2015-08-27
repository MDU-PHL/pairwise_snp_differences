################################################################################
#### Summarise pairwise SNP differences                                    #####
####                                                                       #####
####  An R script to summarise pairwise SNP differences from nullabor      #####
####  output.                                                              #####
####                                                                       #####
####  Input:                                                               #####
####    - cat: a CSV or tab-delimited file with a column of sequence IDs   #####
####        that match the IDs in the FASTA file, and one or more columns  #####
####        indicating to which group the sequence belongs to.             #####
####    - seqs: a FASTA alignment file with all the SNP differences        #####
####        outputted by nullabor <optional>                               #####
####    - diff: a CSV file containing the pairwise SNP difference count    #####
####        outputted by nullabor <optional>                               #####
####                                                                       #####
####  Output:                                                              #####  
####    - tab: a table with summary information for eqch possible          #####
####      pairwise category comparison, including:                         #####
####          -mu: mean differences                                        #####
####          -sd: standard deviation                                      #####
####          -min: minimum difference observed                            #####
####          -max: maximum difference observed                            #####
####    - fig: a figure depicting the information in the table             #####
####                                                                       #####
####  Running from the command line:                                       #####
####                                                                       #####
####    Rscript summarise_snp_differences.R cat <seqs | diff>              #####
####                                                                       #####
####  Running from within R:                                               #####
####                                                                       #####
####    Change the appropriate parameters below to point to the            #####  
####    necessary files.                                                   #####
####                                                                       #####
####  History:                                                             #####
####    - version 0.1: 27 August 2015                                      #####
################################################################################



