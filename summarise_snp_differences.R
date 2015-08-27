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


################################################################################
# Load some necessary libraries
# 
require(ape)

################################################################################



################################################################################
# If not running off the command line, change these parameters to point the 
# approriate files, e.g.:
#   cat = '/home/user/cat.csv'

cat_file = NULL

# Only one of these files needs to be specified. If both are specified, the
# diff file will have precedence
seq_file = NULL
diff_file = NULL

################################################################################



################################################################################
# The main function
# 
summ_distances <- function(categories, dist_obj){
  #dist_obj is a distance object produced by using the dist.dna() function of ape
  #categories is a data.frame with two columns:
  #   - seq_id: that matches the sequence ids in dist_obj
  #   - groups: that assigns the individual seq_ids to a group
  
  # some sanity checks
  if(class(dist_obj) != 'dist' & class(dist_obj) != 'matrix') {
    stop("dist_obj is not an object of type dist or matrix! 
         Please use dist.dna() to create a distance object first OR
         input a CSV file with count of differences produced by 
         nullabor")
  }
  
  if(!is.data.frame(categories)){
    stop("categories must be a data.frame! 
         Please create a data.frame with the metadata first.")
  }
  
  if(ncol(categories) > 2) {
    warning("Number of colums in categories is >2, 
            taking the first two columns only")
    categories <- categories[,c(1,2)]
  }
  
  if(!all(sort(names(categories)) == sort(c("seq_id", "groups")))) {
    warning("The columns of categories do not have names this function 
            recognizes. It will assume that the first column contains seq_ids, 
            and the second column the relevant categories.")
    names(categories) <- c("seq_id", "categories")
  }
  
  dat <- as.matrix(dist_obj)
  taxa <- unique(as.character(categories[,'groups']))
  n_taxa <- length(taxa)
  total_comp <- (n_taxa^2 + n_taxa)/2
  out <- data.frame(grp1 = character(total_comp), 
                    grp2 = character(total_comp),
                    comp = character(total_comp),
                    type = rep("inter-group", total_comp),
                    mu = numeric(total_comp), 
                    sd = numeric(total_comp), 
                    min_dist = numeric(total_comp),
                    max_dist = numeric(total_comp), stringsAsFactors = F)
  n_comp = 1
  for(i in 1:n_taxa){
    g1 <- taxa[i]
    seq_1 <- as.character(categories[categories$groups == g1, 'seq_id'])
    for(j in i:n_taxa){
      g2 <- taxa[j]
      seq_2 <- as.character(categories[categories$groups == g2, 'seq_id'])
      tmp_dat <- dat[seq_1, seq_2]
      if(i == j) {
        tmp_dat <- tmp_dat[lower.tri(tmp_dat)]
        out[n_comp, 'type'] <- 'intra-group'
      }
      out[n_comp, "grp1"] <- g1
      out[n_comp, "grp2"] <- g2
      out[n_comp, "comp"] <- paste(g1, g2, sep='_')
      out[n_comp, "mu"] <- mean(tmp_dat)
      out[n_comp, "sd"] <- sd(tmp_dat)
      out[n_comp, "min_dist"] <- min(tmp_dat)
      out[n_comp, "max_dist"] <- max(tmp_dat)
      n_comp = n_comp + 1
    }
  }
  return(out)
}

################################################################################


################################################################################
# if running off a FASTA file, it is necessary to calculate the pairwise distance
# matrix.

calc_pairwise_distance <- function(seq_file, model = 'raw') {
  # this function takes as input a string defining a path to a FASTA file
  # and a string to pass on to the function dna.dist() specifying the distance
  # model to use. 
  # The model is normally assumed to be 'raw', which means it takes the 
  # proportional pairwise differences. In most cases at MDU, this is a reasonable
  # measure if, however, finite mutation models are required to account for 
  # mutation saturation, other methods can be used. Type ?dna.dist to read the 
  # manual page.
   
  if(!file.exits(seq_data)) {
    stop(paste("Could not find file:", seq_file, "\n"))  
    }
  seq_data <- read.FASTA(file = seq_file)
  raw_dist <- dist.dna(x = seq_data, model = "raw")
  return(raw_dist)
}

################################################################################
# if running off a CSV/TSV file with counts of SNP differences
#

read_diff_file <- function(diff_file) {
  if(!file.exists(diff_file)) {
    stop(paste("Could not find file:", diff_file, "\n"))
  }
  file_sep = ','
  if(!grepl(pattern = 'csv', x = tolower(diff_file))) {
    file_sep = '\t'
  }
  diff_mat <- as.matrix(
                read.table(file = diff_file, 
                         header = TRUE, 
                         row.names = 1, 
                         check.names = F, 
                         sep = file_sep)
                )
  return(diff_mat)
}

################################################################################

################################################################################
# output the results to a pretty table
#

################################################################################

################################################################################
# output the results to a pretty figure
#

plot_figure <- function(summ_table, outfile = NULL, file_type = 'png') {
  # here, we take the output from the summarise_snp_differences function, and 
  # make a plot that includes the mean, the sd, and the min/max for each of 
  # the possible comparisons
  # If outfile is specified, the figure is outputted as a png or pdf.
  require(ggplot2)
  p1 <- ggplot(summ_table, aes(x = comp, y = mu, colour = type)) + 
          geom_point(size = 4) + 
          geom_errorbar(aes(ymax = mu + sd, ymin = mu - sd, width = 0.05)) + 
          geom_point(aes(x = comp, min_dist), size = 3) +
          geom_point(aes(x = comp, max_dist), size = 3) +
          xlab("Pairwise comparisons") + 
          ylab("Mean proportional SNP differences\n(errorbars: sd; points: min/max)") +
          scale_colour_discrete(name = "Comparison type") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if(is.null(outfile)) {
    print(p1)
  } else {
    if (file_type == 'png') {
      png(filename = outfile, width = 2048, height = 1536, res = 300)
    } else {
      pdf(file = outfile, width = 7, height = 5.5)
    }
    print(p1)
    dev.off()
  }
}

################################################################################


################################################################################
# If running from the command line:
# 

if(!interactive()) {
  # collect the arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # check the argument length
  if(length(args) < 1 | length(args) > 3) {
    args <- c("--help")
  }
  
  ## Help section
  if("--help" %in% args) {
    cat("
        Summarise SNP differences
        
        Arguments:
        filename    - a string defining the path to the categories file
        --seq=filename    - a string defining the path to
                            the FASTA file <optional>
        --diff=filename   - a string defining the path to 
                            the SNP differences matrix <optional>
        --help            - print this text
        
        Example:
        ./summarise_snp_differences.R \"/home/user/cat.csv\" --seq=\"/home/user/seq.fa\"\n\n")
    
    q(save="no")
  }
  
  cat(args, "\n")
  
  cat_file = args[1]
  
  args <- args[2:length(args)]
  
  parse_args <- function(x) strsplit(sub("^--", "", x), "=")
  args_df <- as.data.frame(do.call("rbind", parse_args(args)), stringsAsFactors = FALSE)
  names(args_df) <- c("key", "value")
  
  if('diff' %in% args_df$key) {
    diff_file = args_df[args_df$key == 'diff', 'value']
  } else {
    seq_file = args_df[args_df$key == 'seq', 'value']
  }
}

cat(cat_file, seq_file, diff_file, "\n")

