################################################################################
#### Summarise pairwise SNP differences                                    #####
####                                                                       #####
####  An R script to summarise pairwise SNP differences from nullabor      #####
####  output.                                                              #####
####                                                                       #####
####  Minimum input:                                                       #####
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
# If not running off the command line, change these parameters to point to the 
# approriate files, e.g.:
#   cat = '/home/user/cat.csv'
#   
#   To test the script substitute below as follows:
#   
#     cat_file = 'test/woodm_grouping.csv'
#     seq_file = 'test/woodm.fa'

cat_file = NULL

# Only one of these files needs to be specified. If both are specified, the
# diff file will have precedence
seq_file = NULL
diff_file = NULL

# Options:
#   Change the following options to set output
#   
out_basename = 'snp_diff'
tab_fmt = "csv" # options are "csv" or "md"
tab_type = "raw" # options "pretty" or "raw" --- "pretty" formats numbers in
                      # scientific format (e.g., 1.05e-9), while "raw" gives
                      # raw value outputs
fig_fmt = "png"     # options are "png" or "pdf"
fig_type = "median" # option to print out mean +/- sd and range ("mean") or
                  # median +/- interquartile range ("median")
fig_filter = "both" # produce plot with either "inter" or "intra" or "both"
                    # comparisons
                    # it will also accept "specific". If this flag is set
                    # then comparison below must also be set. This flag 
                    # will limit to plotting only comparisons that include
                    # the group specified in "comp" below.
fig_comp = NULL     # A string that specifies a unique cat:group pair found in your
                    # categories file.
                    # here cat refers to a column name heading, and group refers
                    # to one of the grouping units within that column.
                    # (e.g., 'clade:clade_a' in the test file). Note it must be
                    # specified with the colon mark.
exclude_ids = NULL # a string to a path to a file with one sequence ID per
                   # line. these sequences will be excluded from the 
                   # analysis.
save_long = FALSE  # 

################################################################################

################################################################################
# for testing purposes
#      cat_file = 'test/woodm_grouping.csv'
#      seq_file = 'test/woodm.fa'
#      fig_filter = 'specific'
#      fig_comp = 'clade:clade_a'
#      save_long = TRUE

################################################################################

################################################################################
# The function that summarises the data
#  
summ_distances <- function(categories, dist_obj, save_long, outfile){
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
    names(categories) <- c("seq_id", "groups")
  }
  
  # calculations
  dat <- as.matrix(dist_obj)
  taxa <- unique(as.character(categories[,'groups']))
  n_taxa <- length(taxa)
  total_comp <- (n_taxa^2 + n_taxa)/2
  out <- data.frame(grp1 = character(total_comp), 
                    grp2 = character(total_comp),
                    comp = character(total_comp),
                    N = numeric(total_comp),
                    type = rep("inter-group", total_comp),
                    mu = numeric(total_comp), # store the mean
                    sd = numeric(total_comp), # store the standard deviation
                    med = numeric(total_comp), # store the median
                    lqr = numeric(total_comp), # store the lower quartile
                    uqr = numeric(total_comp), # store the upper quartile
                    min_dist = numeric(total_comp), # store the min value
                    max_dist = numeric(total_comp), # store the maximum value
                    stringsAsFactors = F)
  if (save_long) {
    n_inds <- nrow(dat)
    n_entries = n_inds^2
    out_long <- data.frame(taxa1 = character(n_entries),
                           taxa2 = character(n_entries),
                           grp1 = character(n_entries),
                           grp2 = character(n_entries),
                           type = rep("inter-group", n_entries),
                           count = numeric(n_entries),
                           stringsAsFactors = F)
    outlong_fn = paste(outfile, "_long.csv", sep = "")
    count_entries = 1
  }
  n_comp = 1
  for(i in 1:n_taxa){
    g1 <- taxa[i]
    seq_1 <- as.character(categories[categories$groups == g1, 'seq_id'])
    for(j in i:n_taxa){
      g2 <- taxa[j]
      seq_2 <- as.character(categories[categories$groups == g2, 'seq_id'])
      tmp_dat <- dat[seq_1, seq_2]
      if (save_long) {
        for (ii in 1:length(seq_1)) {
          for (jj in 1:length(seq_2)) {
            out_long[count_entries,'taxa1'] = seq_1[ii]
            out_long[count_entries,'taxa2'] = seq_2[jj]
            out_long[count_entries,'grp1'] = g1
            out_long[count_entries,'grp2'] = g2
            if (i == j) {
              out_long[count_entries,'type'] = 'intra-group'
            }
            out_long[count_entries,'count'] = dat[seq_1[ii], seq_2[jj]]
            count_entries = count_entries + 1
          }
        }
        if (n_entries - count_entries < 100) {
          # in case we start to run out of space
          tmp_long <- data.frame(taxa1 = character(n_entries),
                                 taxa2 = character(n_entries),
                                 grp1 = character(n_entries),
                                 grp2 = character(n_entries),
                                 type = rep("inter-group", n_entries),
                                 count = numeric(n_entries),
                                 stringsAsFactors = F)
          out_long <- rbind(out_long, tmp_long)
        }
      }
      out[n_comp, "grp1"] <- g1
      out[n_comp, "grp2"] <- g2
      out[n_comp, "N"] <- length(tmp_dat)
      out[n_comp, "comp"] <- paste(g1, g2, sep='_')
      if(i == j) {
        if(length(tmp_dat) > 1){
          #if length is one, this results in a empty set. 
          #so, added this condition to fix the problem
          tmp_dat <- tmp_dat[lower.tri(tmp_dat)]
        }
        out[n_comp, 'type'] <- 'intra-group'
        out[n_comp, "comp"] <- g1
      }
      if (length(tmp_dat) > 1 & max(tmp_dat) > min(tmp_dat)) {
        out[n_comp, "mu"] <- mean(tmp_dat)
        out[n_comp, "sd"] <- sd(tmp_dat)
        out[n_comp, "med"] <- quantile(tmp_dat, p = 0.50)
        out[n_comp, "lqr"] <- quantile(tmp_dat, p = 0.25)
        out[n_comp, "uqr"] <- quantile(tmp_dat, p = 0.75)
        out[n_comp, "min_dist"] <- min(tmp_dat)
        out[n_comp, "max_dist"] <- max(tmp_dat)
      } else {
          out[n_comp, "mu"] <- mean(tmp_dat)
          out[n_comp, "sd"] <- 0
          out[n_comp, "med"] <- quantile(tmp_dat, p = 0.50)
          out[n_comp, "lqr"] <- quantile(tmp_dat, p = 0.25)
          out[n_comp, "uqr"] <- quantile(tmp_dat, p = 0.75)
          out[n_comp, "min_dist"] <- min(tmp_dat)
          out[n_comp, "max_dist"] <- max(tmp_dat)
        }
      n_comp = n_comp + 1
      }
  }
  out$type <- factor(out$type, levels = c("intra-group", "inter-group"))
  out$comp <- factor(out$comp, levels = out$comp[order(out$type, out$comp)])
  if (save_long) {
    out_long <- out_long[1:(count_entries-1),]
    write.table(x = out_long, file = outlong_fn, quote = F, sep = ",", row.names = F)
  }
  return(out)
}

################################################################################

################################################################################
# read the categories table
#

read_cat_file <- function(categories, exclude_ids = NULL) {
  if(!file.exists(categories)) {
    stop(paste("Could not find file:", categories, "\n"))
  }
  if(!is.null(exclude_ids) && !file.exists(exclude_ids)) {
    stop(paste("Could not find file:", exclude_ids, "\n"))
  }
  file_sep = ','
  if(!grepl(pattern = 'csv', x = tolower(categories))) {
    file_sep = '\t'
  }
  cat_df <- read.table(file = categories, 
               header = TRUE, 
               check.names = F, 
               sep = file_sep,
               stringsAsFactors = F)
  if(!is.null(exclude_ids)) {
    exclude_ids <- read.table(file = exclude_ids, 
                              header = F, 
                              stringsAsFactors = F)
  }
  #if file is an mlst.tab output from nullabor
  if(all(c("FILE", "ST") %in% names(cat_df))) {
    seq_id <- gsub(pattern = "\\/contigs\\.fa", replacement = "", cat_df$FILE)
    ST <- cat_df$ST
    ix_calls <- which(ST != "-")
    cat_df <- data.frame(seq_id = seq_id[ix_calls], ST = ST[ix_calls], stringsAsFactors = F)
    if(!is.null(exclude_ids)) {
      cat_df <- cat_df[!(cat_df[,1] %in% exclude_ids), ]
    }
    cat_list <- list(ST = cat_df)
    return(cat_list)
  }
  #otherwise, treat as a file prepared by the user
  if(!is.null(exclude_ids)) {
    cat_df <- cat_df[!(cat_df[,1] %in% exclude_ids), ]
  }
  if(ncol(cat_df) == 2) {
    cat_list <- list(cat_df)
    names(cat_list) <- names(cat_df)[2]
  } else if(ncol(cat_df) > 2) {
    warning("Number of identified columns in category file is >2, 
            assuming that the first column contains the sequence IDs")
    cats <- names(cat_df)[-1]
    seqid <- names(cat_df)[1]
    cat_list <- lapply(cats, function(cat) subset(cat_df, select = c(seqid, cat)))
    names(cat_list) <- cats
  }
  return(cat_list)
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
   
  if(!file.exists(seq_file)) {
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

write_summ_table <- function(summ_table, 
                             outfile, 
                             file_type = 'csv', 
                             method = "pretty") {
  # here, we take the output from the summarise_snp_differences function and 
  # make a CSV table, which can be read into EXCEL, or someother spreadsheet
  # application. 
  # If the file_type is 'md', the table will be outputted in markdown format.
  
  outfile = paste(outfile, file_type, sep = ".")
  if(method == 'pretty') {
    pretty_column_names <- c("Group 1", "Group 2", "N", "Comparison", "Mean (±SD)", "Median", "Inter-Quartile Range", "Range")
    pretty_mean <- format(summ_table[,'mu'], scientific = T, digits = 3)
    pretty_sd <- format(summ_table[,'sd'], scientific = T, digits = 3)
    pretty_musd <- paste(pretty_mean, " (±", pretty_sd,")", sep = "")
    pretty_median <- format(summ_table[,'med'], scientific = T, digits = 3)
    pretty_lqr <- format(summ_table[,'lqr'], scientific = T, digits = 3)
    pretty_uqr <- format(summ_table[,'uqr'], scientific = T, digits = 3)
    pretty_iqr <- paste(pretty_lqr, pretty_uqr, sep = "; ")
    pretty_min <- format(summ_table[,'min_dist'], scientific = T, digits = 3)
    pretty_max <- format(summ_table[,'max_dist'], scientific = T, digits = 3)
    pretty_range <- paste(pretty_min, pretty_max, sep = "; ")
    tab <- data.frame(summ_table$grp1, summ_table$grp2, summ_table$type, summ_table$N)
    tab$musd <- pretty_musd
    tab$median <- pretty_median
    tab$iqr <- pretty_iqr
    tab$range <- pretty_range
    names(tab) <- pretty_column_names
  } else {
    column_names <- c("Group 1", "Group 2", "Comparison", "N", "Mean", "SD", "Media", "Lower IQR", "Upper IQR", "Min", "Max")
    tab <- summ_table[,c(1,2,4,5,6,7,8,9,10,11,12)]
    names(tab) <- column_names
  }
  if(file_type == 'md') {
    require(pander)
    capture.output(pander(tab, split.tables = Inf), file = outfile)
  } else {
    write.table(x = tab, file = outfile, sep = ",", quote = F, row.names = F)
  }
}

################################################################################

################################################################################
# output the results to a pretty figure
#

plot_figure <- function(summ_table, outfile = NULL, 
                        file_type = 'png', 
                        fig_type = "mean", 
                        fig_filter = "both",
                        fig_comp = NULL) {
  # here, we take the output from the summarise_snp_differences function and 
  # make a plot that includes the mean, the sd, and the min/max for each of 
  # the possible comparisons
  # If outfile is specified, the figure is outputted as a png or pdf.
  require(ggplot2)
  tmp_tab <- summ_table
  if(fig_filter == 'inter') {
    tmp_tab <- summ_table[summ_table$type == "inter-group",]
    outfile = paste(outfile, "inter", sep = "_")
  } else if (fig_filter == 'intra') {
    tmp_tab <- summ_table[summ_table$type == "intra-group",]
    outfile = paste(outfile, "intra", sep = "_")
  } else if (fig_filter == "specific" && !is.null(fig_comp)) {
    tmp_tab <- summ_table[grepl(pattern = fig_comp, x = summ_table$comp),]
    outfile = paste(outfile, fig_comp, sep = "_")
  }
  
  if(fig_type == 'mean'){
    p1 <- ggplot(tmp_tab, aes(x = comp, y = mu, colour = type)) + 
            geom_point(size = 3 ,shape = 18) + 
            geom_errorbar(aes(ymax = mu + sd, ymin = mu - sd, width = 0.05)) + 
            geom_point(aes(x = comp, min_dist), size = 2) +
            geom_point(aes(x = comp, max_dist), size = 2) +
            xlab("Pairwise comparisons") + 
            ylab("Mean SNP distance\n(errorbars: sd; points: min/max)") +
            scale_colour_discrete(name = "Comparison type") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3), 
                  legend.position="bottom")
  } else if(fig_type == 'median') {
    p1 <- ggplot(tmp_tab, aes(x = comp, y = med, colour = type)) + 
      geom_point(size = 3 ,shape = 18) + 
      geom_errorbar(aes(ymax = uqr, ymin = lqr, width = 0.05)) + 
      geom_point(aes(x = comp, min_dist), size = 2) +
      geom_point(aes(x = comp, max_dist), size = 2) +
      xlab("Pairwise comparisons") + 
      ylab("Median SNP distance\n(errorbars: inter-quartile range; points: min/max)") +
      scale_colour_discrete(name = "Comparison type") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 3), 
            legend.position="bottom")
  }
  
  if(is.null(outfile)) {
    print(p1)
  } else {
    outfile = paste(outfile, file_type, sep = ".")
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
# the main() function
# 
main <- function(categories, 
                 seq_file = NULL, 
                 diff_file = NULL, 
                 out_base = NULL, 
                 tab_fmt = NULL,
                 tab_type = NULL,
                 fig_filter = NULL,
                 fig_type = NULL,
                 fig_fmt = NULL,
                 fig_comp = NULL,
                 exclude_ids = NULL,
                 save_long = NULL) {
  ################################################################################
  ## check dependencies
  
  miss_packages = c()
  
  if(!require(ggplot2, quietly = T)) {
    miss_packages = c(miss_packages, "ggplot2")
  }
  
  if(!require(ape, quietly = T)){
    miss_packages = c(miss_packages, "ape")
  }
  
#   if(!require(DT, quietly = T)) {
#     miss_packages = c(miss_packages, "DT")
#   }
  
  if (length(miss_packages) >= 1) {
    mp <- paste(miss_packages, sep = ",", collapse = "")
    stop(paste("It seems we are missing some dependencies: ", mp, ". To install type the following:\n
               install.packages(c(\'",mp,"\'))", sep = ""))
  }
  # Load some necessary libraries
  library(ape)
  library(ggplot2)
  
  #load the categories
  cats_list <- read_cat_file(categories = categories, exclude_ids = exclude_ids)
  
  #load the data
  if(is.null(diff_file)) {
    dist_obj <- calc_pairwise_distance(seq_file = seq_file, model = "raw")
  } else {
    dist_obj <- read_diff_file(diff_file = diff_file)
  }
  
  #summarise the information
  cats <- names(cats_list)
  #check fig.comp, and parse
  if (!is.null(fig_comp)) {
    tmp <- strsplit(x = fig_comp, ":")
    cat_keep <- tmp[[1]][1]
    fig_comp <- tmp[[1]][2]
    cats <- cats[which(cats %in% cat_keep)]
  }
  for(cat in cats) {
    outf_b <- paste(out_base, cat, sep = "_")
    results <- summ_distances(categories = cats_list[[cat]], 
                              dist_obj = dist_obj, 
                              save_long = save_long,
                              outfile = outf_b)
  
    #output table
    write_summ_table(summ_table = results, 
                     file_type = tab_fmt, 
                     method = tab_type,
                     outfile = outf_b)
  
    #output figure
    plot_figure(summ_table = results,
                fig_filter = fig_filter,
                fig_type = fig_type,
                file_type = fig_fmt,
                fig_comp = fig_comp,
                outfile = outf_b)
  }
}

################################################################################
# If running from the command line:
# 

if(!interactive()) {
  # collect the arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  # check the argument length
  if(length(args) < 2 | length(args) > 7) {
    args <- c("--help")
  }
  
  ## Help section
  if("--help" %in% args) {
    cat("
        Summarise SNP differences
        
        Necessary arguments:
        filename    - a string defining the path to the categories file

        Optional arguments:
        --seq=filename    - a string defining the path to
                            the FASTA file (ignored if --diff is defined)
        --diff=filename   - a string defining the path to 
                            the SNP differences matrix (must be defined if no
                            --seq is defined)
        --out_basename    - basename for output files (default: snp_diff)
        --tab_fmt         - output format for table (default: \"csv\")
                            \"csv\" or \"md\" for CSV or Markdown, respectively
        --tab_type        - make table \"pretty\" by formatting numbers or 
                              or output \"raw\" numbers (default: \"pretty\")
        --fig_filter      - one of \"both\", \"intra\", \"inter\", or \"specific\".
                            use \"both\", \"intra\" or \"inter\" if wanting to plot 
                            both inter and intra distance comparisons on 
                            the same plot, or only intra or inter, respectively.
                            if wanting to plot just comparisons that include a 
                            single category, use \"specific\", and then specify
                            the group name in --fig_comp.
                            (default: \"both\")
        --fig_comp        - when specifying --fig_filter=\"specific\", this must
                            be specified. A string identifying one category in the 
                            comp file to plot (e.g., \"clade:clade_a\"). To avoid saving
                            over previous analyses, the group name is added as an
                            extension. Note the use of the colon to specify the 
                            category (i.e., column in the cat file) and group 
                            (i.e., name of the grouping unit within that column).
        --fig_type        - output figure of mean +/- sd and range (\"mean\") or
                            median +/- inter-quartile range (\"median\")
        --fig_fmt         - output format for figure (default: \"png\")
                            \"png\" or \"pdf\" for PNG or PDF, respectively
        --save_long       - \"TRUE\" or \"FALSE\" if wanting to save the raw data
                            in long format, along with with all the additional
                            metadata (default: \"FALSE\")
        --exclude_ids     - string defining the path to a file with sequence
                            ids to be excluded, one per line (default: None)
        --help            - print this text
        
        Example:
        ./summarise_snp_differences.R \"/home/user/cat.csv\" --seq=\"/home/user/seq.fa\"\n\n")
    
    q(save="no")
  }
  
  #parse arguments
  cat_file = args[1]
  args <- args[2:length(args)]
  parse_args <- function(x) strsplit(sub("^--", "", x), "=")
  args_df <- as.data.frame(
                            do.call("rbind", parse_args(args)), 
                            stringsAsFactors = F)
  names(args_df) <- c("key", "value")
  if('diff' %in% args_df$key) {
    diff_file = args_df[args_df$key == 'diff', 'value']
  } else {
    seq_file = args_df[args_df$key == 'seq', 'value']
  }
  for(arg in c("out_basename", "tab_fmt", "tab_type", "fig_fmt", "exclude_ids")) {
    if (arg %in% args_df$key) {
      assign(arg, args_df[args_df$key == arg, 'value'])
    }
  }
}

#run the main function
main(categories = cat_file, 
     seq_file = seq_file, 
     diff_file = diff_file, 
     out_base = out_basename, 
     tab_fmt = tab_fmt, 
     tab_type = tab_type,
     fig_filter = fig_filter,
     fig_type = fig_type,
     fig_fmt = fig_fmt,
     fig_comp = fig_comp,
     save_long = save_long,
     exclude_ids = exclude_ids)

cat("The script has ended successfully!\n")
################################################################################