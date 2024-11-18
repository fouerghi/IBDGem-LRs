#5/14/24
#Code to simulate and analyze diploid individuals
#under the IBDGem model under complete LD.

set.seed(203484)

# generate all possible joint genotypes based on the number of loci 
# store in a matrix where every row is a locus and pairs of columns represent each joint genotype
generate_genotype_matrix <- function(num_loci) {
  # Generate haplotypes in the specified format
  haplotypes <- character(num_loci + 1)
  
  for (i in 0:num_loci) {
    # Construct the haplotype
    haplotypes[i + 1] <- paste0(c(rep(0, num_loci - i), rep(1, i)), collapse = "")
  }
  
  # Generate unique combinations of haplotypes
  joint_genotypes <- character()
  
  for (i in 1:length(haplotypes)) {
    for (j in i:length(haplotypes)) {
      # Combine haplotypes to create joint genotypes
      joint_genotype <- paste(haplotypes[i], haplotypes[j], sep = "/")
      joint_genotypes <- c(joint_genotypes, joint_genotype)
    }
  }
  
  # Create a matrix to store the individual loci of the joint genotypes
  genotype_matrix <- matrix(0, nrow = num_loci, ncol = length(joint_genotypes) * 2)
  
  # Fill the matrix with haplotype elements for each joint genotype
  for (col in 1:length(joint_genotypes)) {
    hap1 <- strsplit(joint_genotypes[col], "/")[[1]][1]  # First haplotype
    hap2 <- strsplit(joint_genotypes[col], "/")[[1]][2]  # Second haplotype
    
    for (locus in 1:num_loci) {
      # Fill the matrix with individual loci from haplotypes as integers
      genotype_matrix[locus, col * 2 - 1] <- as.integer(substr(hap1, locus, locus))
      genotype_matrix[locus, col * 2] <- as.integer(substr(hap2, locus, locus))
    }
  }
  
  return(genotype_matrix)
}

# now based on whether the genotype is heterozygous or homozygous, make two copies or one copy respectively
# store in a new matrix
process_genotype_matrix <- function(genotype_matrix) {
  num_loci <- nrow(genotype_matrix)  # Number of loci
  num_genotypes <- ncol(genotype_matrix) / 2  # Number of unique joint genotypes
  
  # Create an empty list to hold processed columns
  processed_genotypes <- list()
  
  for (i in 1:num_genotypes) {
    # Extract the two columns (haplotypes) to compare
    hap1 <- genotype_matrix[, (i * 2 - 1)]
    hap2 <- genotype_matrix[, (i * 2)]
    
    # Check if the haplotypes are equal (homozygous) or not (heterozygous)
    if (all(hap1 == hap2)) {
      # Homozygous: keep one copy (hap1 and hap2 are the same)
      processed_genotypes[[length(processed_genotypes) + 1]] <- hap1
      processed_genotypes[[length(processed_genotypes) + 1]] <- hap2
    } else {
      # Heterozygous: keep two copies of the full set (hap1 hap2, hap1 hap2)
      processed_genotypes[[length(processed_genotypes) + 1]] <- hap1
      processed_genotypes[[length(processed_genotypes) + 1]] <- hap2
      processed_genotypes[[length(processed_genotypes) + 1]] <- hap1
      processed_genotypes[[length(processed_genotypes) + 1]] <- hap2
    }
  }
  
  # Combine the list into a matrix
  result_matrix <- do.call(cbind, processed_genotypes)
  
  # Set the correct number of rows
  result_matrix <- result_matrix[1:num_loci, , drop = FALSE]
  
  return(result_matrix)
}

# function to generate the number of reads according to a Poisson distribution depending on coverage depth
generate_reads_poisson <- function(lambda) {
  while(TRUE) {
    value <- rpois(1, lambda)
    if (value != 0) {
      return(value)
    }
  }
}

# function to generate reads according to the individual's genotype
generate_reads <- function(target_genotype, readnumber ,epsilon,k) {
  reads <- matrix(0, nrow = k, ncol = readnumber) 
  for (i in 1:k) {
    allele_0 <- target_genotype[i,1]
    allele_1 <- target_genotype[i, 2]
    
    if (allele_0 == allele_1) {
      # Homozygous case (0/0 or 1/1)
      if (allele_0 == 0) {
        # 0/0 case
        prob_1 <- epsilon
      } else {
        # 1/1 case
        prob_1 <- 1 - epsilon
      }
    } else {
      # Heterozygous case (0/1 or 1/0)
      prob_1 <- 0.5
    }
    
    # Generate reads for the SNP
    reads_output <- rbinom(readnumber,1,prob_1)
    for (j in 1:readnumber){
      reads[i,j] <- reads_output[j]
    }
  }
  return(reads)
}

# function to simulate the reads for the target
simulate_reads_for_target <- function(target_genotype,eps,number_of_reads,k){
  reads <- generate_reads(target_genotype, readnumber = number_of_reads, epsilon = eps,k)
  return(reads)
}

# function to simulate .legend file (necessary for IBDGem input)
simulate_legend_file <- function(k){
  # Generate IDs 
  IDs <- paste("SNP", 1:k, sep = "")
  # Generate positions: 100, 200, 300, ...
  positions <- seq(100, length.out = k, by = 100)
  # Define alleles
  alleles <- c("A", "T", "C", "G")
  # Initialize empty vectors for alleles
  allele0 <- character(k)
  allele1 <- character(k)
  # Loop through each SNP
  for (i in 1:k) {
    # Randomly select allele0
    allele0[i] <- sample(alleles, 1)
    # Select allele1 from the remaining alleles
    allele1[i] <- sample(alleles[alleles != allele0[i]], 1)
  }
  snp_matrix <- matrix(data = c(IDs, positions, allele0, allele1), nrow = k, ncol = 4, byrow = FALSE)
  colnames(snp_matrix) <- c("ID", "pos", "allele0", "allele1")
  return(snp_matrix)
}

# function to simulate the input files for IBDGem
getsims <- function(n, k, p,eps, nsim,number_of_reads,order_number, reference){
  
  new_directory_name <- paste0("n", n, "_k", k,"_r",number_of_reads,"_order",order_number)
  new_directory <- file.path(directory, new_directory_name)
  dir.create(new_directory)
  setwd(new_directory)
  
  num_individuals <- ncol(reference)/2
  new_sub_directory_name <- paste0("n", n, "_k", k,"_s",num_individuals,"_entire")
  new_sub_directory <- file.path(new_directory, new_sub_directory_name)
  dir.create(new_sub_directory)
  setwd(new_sub_directory)
  
  for (m in 1:nsim){
    
    # create appropriate directories
    another_sub_directory_name <- paste0("n", n, "_k", k,"_s",num_individuals,"_entire_sim",m)
    another_sub_directory <- file.path(new_sub_directory, another_sub_directory_name)
    dir.create(another_sub_directory)
    
    # Set the new directory as the working directory
    setwd(another_sub_directory)
    
    #create directory to put in IBDGem results LD
    dir.create(file.path(another_sub_directory, "output_LD"))
    #create directory to put in IBDGem results non-LD
    dir.create(file.path(another_sub_directory, "output_non_LD"))
    
  
    haplo_numbers <- sample(1:ncol(reference), 2)
    target_genotype_k <-cbind(reference[,haplo_numbers[1]],reference[,haplo_numbers[2]])
    snp_matrix <- simulate_legend_file(k=k)
    reads <- simulate_reads_for_target(target_genotype_k, eps = eps, number_of_reads = number_of_reads, k = k)
    genotypes <- cbind(reference,target_genotype_k)
    #to generate .indv file
    samples <- paste0("sample", 1:(num_individuals+1))
    
    # Sample data for CHR (replace with your actual data)
    chr <- rep("chr1", k)
    
    # Extract position information from snp_matrix (replace snp_matrix$position with your actual data)
    pos <- as.numeric(snp_matrix[,"pos"])
    # Set frequency value
    freq <- rep(p, k)
    
    # Create a data frame with CHR, POS, and FREQ
    my_table <- data.frame(CHR = chr, POS = pos, FREQ = freq)
    
    write.table(genotypes, file = "test.hap", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(snp_matrix, file = "test.legend", sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(samples, "test.indv")
    write.table(my_table, file = "output.frq", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    snp_positions <- snp_matrix[, "pos"]
    reference_alleles <- snp_matrix[, "allele0"] 
    bases <- matrix(nrow = nrow(reads), ncol = number_of_reads)
    
    for (i in 1:nrow(reads)){
      for (j in 1:ncol(reads)){
        if (reads[i,j]==0){
          bases[i,j] <- "."
        }
        else{
          bases[i,j] <- snp_matrix[i, "allele1"]
        }
      }
    }
    
    if (eps == 0.0001){
      symbol = "I"
    }
    if (eps == 0.001){
      symbol = "?"
    }
    if (eps == 0.01){
      symbol = 5
    }
    if (eps == 0.02){
      symbol = 2
    }
    if (eps == 0.05){
      symbol = "."
    }
    if (eps == 0.1){
      symbol = "+"
    }
    if (eps == 0.2){
      symbol = "("
    }
    
    concatenated_column <- apply(bases, 1, paste, collapse = "")
    basequals <- paste(rep(symbol, number_of_reads), collapse = "")
    alignmapping = paste(rep("!", number_of_reads), collapse = "")
    pileup_data <- data.frame(
      CHROM = rep("chr1", k),
      POS = snp_positions,
      REF = reference_alleles,
      DEPTH = rep(number_of_reads, k),
      BASES = concatenated_column, 
      BASEQUALS = rep(basequals, k),
      ALIGNMAPPING = rep(alignmapping, k)
    )
    
    write.table(pileup_data, file = "unknowndna.pileup", sep = "\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
    setwd(new_sub_directory)
  }

  # make a new subdirectory for n,k,downsampled (not the entire dataset)
  setwd(new_directory)
  new_sub_directory_name <- paste0("n", n, "_k", k,"_downsampled")
  new_sub_directory <- file.path(new_directory, new_sub_directory_name)
  dir.create(new_sub_directory)
  setwd(new_sub_directory)
  
  for (m in 1:nsim){
    
    # create appropriate directories
    another_sub_directory_name <- paste0("n", n, "_k", k,"_downsampled_sim",m)
    another_sub_directory <- file.path(new_sub_directory, another_sub_directory_name)
    dir.create(another_sub_directory)
    
    # Set the new directory as the working directory
    setwd(another_sub_directory)
    
    #create directory to put in IBDGem results LD
    dir.create(file.path(another_sub_directory, "output_LD"))
    #create directory to put in IBDGem results non-LD
    dir.create(file.path(another_sub_directory, "output_non_LD"))
    
    #col1 <- c(rep(0, 3/4*k),rep(1, 1/4*k))
    #col2 <- c(rep(0, 1/4*k),rep(1, 3/4*k))
    
    #target_genotype_k <- matrix(c(col1, col2), nrow = k, ncol = 2)
    
    haplo_numbers <- sample(1:ncol(reference), 2)
    target_genotype_k <-cbind(reference[,haplo_numbers[1]],reference[,haplo_numbers[2]])
    snp_matrix <- simulate_legend_file(k=k)
    reads <- simulate_reads_for_target(target_genotype_k, eps = eps, number_of_reads = number_of_reads, k = k)
    
    #to generate .indv file
    samples <- paste0("sample", 1:n)
    
    # Sample data for CHR (replace with your actual data)
    chr <- rep("chr1", k)
    
    # Extract position information from snp_matrix (replace snp_matrix$position with your actual data)
    pos <- as.numeric(snp_matrix[,"pos"])
    # Set frequency value
    freq <- rep(p, k)
    
    # Create a data frame with CHR, POS, and FREQ
    my_table <- data.frame(CHR = chr, POS = pos, FREQ = freq)
    
    #to generate .hap file
    # downsample with replacement to extract n genotypes from the reference dataset
    num_individuals <- ncol(reference)/2
    sampled_individuals <- sample(1:num_individuals, n-1, replace = TRUE)
    sampled_dataset <- reference[, c(rbind(sampled_individuals * 2 - 1, sampled_individuals * 2))]
    genotypes <- cbind(sampled_dataset,target_genotype_k)
    
    #save the .hap file
    write.table(genotypes, file = "test.hap", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    #save the .legend file
    write.table(snp_matrix, file = "test.legend", sep = "\t", quote = FALSE, row.names = FALSE)
    
    #save the .indv file
    writeLines(samples, "test.indv")
    
    # save the allele frequencies file
    write.table(my_table, file = "output.frq", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    # to generate the .pileup file
    
    snp_positions <- snp_matrix[, "pos"]
    reference_alleles <- snp_matrix[, "allele0"] 
    
    #convert the reads to the bases used in pileup format
    #if 0 (i.e. matches reference base), then "."
    #else (i.e. mismatches the reference base), then either A,C,G, or T
    
    # Create a new matrix with the same dimensions as reads_matrix
    bases <- matrix(nrow = nrow(reads), ncol = number_of_reads)
    
    for (i in 1:nrow(reads)){
      for (j in 1:ncol(reads)){
        if (reads[i,j]==0){
          bases[i,j] <- "."
        }
        else{
          bases[i,j] <- snp_matrix[i, "allele1"]
        }
      }
    }
    
    if (eps == 0.0001){
      symbol = "I"
    }
    if (eps == 0.001){
      symbol = "?"
    }
    if (eps == 0.01){
      symbol = 5
    }
    if (eps == 0.02){
      symbol = 2
    }
    if (eps == 0.05){
      symbol = "."
    }
    if (eps == 0.1){
      symbol = "+"
    }
    if (eps == 0.2){
      symbol = "("
    }
    
    concatenated_column <- apply(bases, 1, paste, collapse = "")
    basequals <- paste(rep(symbol, number_of_reads), collapse = "")
    alignmapping = paste(rep("!", number_of_reads), collapse = "")
    pileup_data <- data.frame(
      CHROM = rep("chr1", k),
      POS = snp_positions,
      REF = reference_alleles,
      DEPTH = rep(number_of_reads, k),
      BASES = concatenated_column, 
      BASEQUALS = rep(basequals, k),
      ALIGNMAPPING = rep(alignmapping, k)
    )
    
    write.table(pileup_data, file = "unknowndna.pileup", sep = "\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
    
    setwd(new_sub_directory)
    
  }
}

# take the entire reference dataset and run IBDGem LD mode on it to get the standard LR
# take the downsampled (to n=20 or n=100) dataset with replacement from the dataset above and run IBDGem in LD mode on it to get the IBDGem LR
# take the downsampled (to n=20 or n=100) dataset with replacement from the dataset above and run IBDGem non-LD mode to get LR assuming LR.

kvals <- c((1:10)*10)
#kvals <- c(12,24,36,48,60,72,84)
nsim <- 100
reference_datasets <- list()

# Loop through each k in kvals to generate the appropriate reference dataset for each, and the target genotype -- depends on k
for (k in kvals) {
  genotype_matrix <- generate_genotype_matrix(k)
  reference_dataset <- process_genotype_matrix(genotype_matrix)
  reference_datasets[[as.character(k)]] <- reference_dataset
}


#run the simulations
#for n=20
current_directory <- getwd()
directory_name <- "complete_LD_n20"
directory <- file.path(current_directory, directory_name)
dir.create(directory)
setwd(directory)

for (i in 1:length(kvals)) {
  k <- kvals[i]
  getsims(n = 20, k = k, p=1/2, eps = 0.02, nsim = nsim, number_of_reads = 2,
          order_number = i,reference = reference_datasets[[as.character(k)]])
}

#for n=100
current_directory <- getwd()
directory_name <- "complete_LD_n100"
directory <- file.path(current_directory, directory_name)
dir.create(directory)
setwd(directory)

for (i in 1:length(kvals)) {
  k <- kvals[i]
  getsims(n = 100, k = k, p=1/2, eps = 0.02, nsim = nsim, number_of_reads = 2,
          order_number = i,reference = reference_datasets[[as.character(k)]])
}

############################### Run IBDGem on inputs from above scenario (i.e. complete LD)
system("complete_LD.sh")

############################### Analyze the output of IBDGem
generate_main_directory <- function(list_name, kval, base_directory, i,x,nval) {
  main_directories <- file.path(base_directory, paste0("n",nval,"_k", kval, "_r2_order", i))
  return(main_directories)
}

run_summary <- function(base_directory,xvals,nval){
  i <- 0
  list_of_dfs <- c()
  
  for (x in xvals){
    i <- i+1
    kval <- xvals[[i]]
    main_directory <- generate_main_directory(list_name = xvals,kval,base_directory,i,x,nval)
    # List all subdirectories within the main directory
    subdirectories <- list.dirs(main_directory, recursive = FALSE)
    
    # Function to find the "output_LD" sub-subdirectory within a subdirectory
    find_output_LD <- function(directory) {
      output_LD_path <- file.path(directory, "output_LD")
      if (file.exists(output_LD_path) && file.info(output_LD_path)$isdir) {
        return(output_LD_path)
      } else {
        return(NULL)
      }
    }
    
    # Function to find the "output_non_LD" sub-sub-directory within a sub-directory 
    find_output_non_LD <- function(directory) {
      output_non_LD_path <- file.path(directory, "output_non_LD")
      if (file.exists(output_non_LD_path) && file.info(output_non_LD_path)$isdir) {
        return(output_non_LD_path)
      } else {
        return(NULL)
      }
    }
    read_summary_file <- function(directory,n) {
      summary_file <- file.path(directory, paste0("sample", n, ".","sample",n,".summary.txt"))
      if (file.exists(summary_file)) {
        # Read the file line by line
        lines <- readLines(summary_file)
        
        # Initialize an empty list to store data
        data <- list()
        
        # Process each line
        for (line in lines) {
          # Check if the line starts with a digit
          if (grepl("^\\d+", line)) {
            # Split the line by whitespace
            parts <- strsplit(line, "\\s+")[[1]]
            # Create a named list from parts
            row_data <- setNames(as.list(parts), c("SEGMENT", "START", "END", "LIBD0", "LIBD1", "LIBD2", "NUM_SITES"))
            # Append the row data to the list
            data <- c(data, list(row_data))
          }
        }
        
        # Convert list to dataframe
        result_df <- do.call(rbind, data)
        
        # Convert data types
        result_df <- type.convert(result_df, as.is = TRUE)
        
        return(result_df)
      } else {
        return(NULL)
      }
    }
    
    ####### Applying to LD and non-LD outputs
    
    summary_data_LD_entire <- list()
    summary_data_LD_downsampled <- list()
    summary_data_non_LD_downsampled <- list()
    
    # Loop through each subdirectory
    for (subdir in subdirectories) {
      if(grepl("entire", subdir)){
        sub_subdirectories <- list.dirs(subdir, recursive = FALSE)
        for (sub_subdir in sub_subdirectories) {
          output_LD_entire <- find_output_LD(sub_subdir)
          if (!is.null(output_LD_entire)) {
            n <- as.integer(gsub(".*_s(\\d+)_.*", "\\1", basename(sub_subdir))) + 1
            summary_data_LD_entire[[sub_subdir]] <- read_summary_file(output_LD_entire,n=n)
          }
        }
      }
      else{
        sub_subdirectories <- list.dirs(subdir, recursive = FALSE)
        for (sub_subdir in sub_subdirectories) {
          output_LD_downsampled <- find_output_LD(sub_subdir)
          output_non_LD_downsampled <- find_output_non_LD(sub_subdir)
          if (!is.null(output_LD_downsampled)){
            summary_data_LD_downsampled[[sub_subdir]] <- read_summary_file(output_LD_downsampled,nval)
          }
          if (!is.null(output_non_LD_downsampled)){
            summary_data_non_LD_downsampled[[sub_subdir]] <- read_summary_file(output_non_LD_downsampled,nval)
          }
        }
      }
    }
    
    result_df <- data.frame(log10LR_LD_entire = numeric(nsim),log10LR_LD_downsampled = numeric(nsim), log10LR_non_LD_downsampled = numeric(nsim))
    
    for (s in 1:length(summary_data_non_LD_downsampled)){
      entire_LD_LIBD0 = as.numeric(summary_data_LD_entire[[s]][1,][4])
      entire_LD_LIBD2 = as.numeric(summary_data_LD_entire[[s]][1,][6])
      downsampled_LD_LIBD0 = as.numeric(summary_data_LD_downsampled[[s]][1,][4])
      downsampled_LD_LIBD2 = as.numeric(summary_data_LD_downsampled[[s]][1,][6])
      downsampled_nonLD_LIBD0 = as.numeric(summary_data_non_LD_downsampled[[s]][1,][4])
      downsampled_nonLD_LIBD2 = as.numeric(summary_data_non_LD_downsampled[[s]][1,][6])
      result_df[s,1] = log10(entire_LD_LIBD2/entire_LD_LIBD0)
      result_df[s,2] = log10(downsampled_LD_LIBD2/downsampled_LD_LIBD0)
      result_df[s,3] = log10(downsampled_nonLD_LIBD2/downsampled_nonLD_LIBD0)
    }
    df_name <- paste0("LR_", x)
    assign(df_name, result_df)
    list_of_dfs[[df_name]] <- result_df 
    print('finished')
  }
  return(list_of_dfs)
}

compress_and_modify <- function(df) {
  df_compressed <- rbind(df, df, df)
  df_compressed$log10LR_LD_entire[101:200] <- df_compressed$log10LR_LD_downsampled[101:200]
  df_compressed$log10LR_LD_entire[201:nrow(df_compressed)] <- df_compressed$log10LR_non_LD_downsampled[201:nrow(df_compressed)]
  df_compressed <- df_compressed[, -c(2,3), drop = FALSE]
  return(df_compressed)
}

plot_results <- function(xvals, df, xlab, ylab, main, xlim=NULL) {
  plot(xvals, type="n", ylim = range(c(0, df)), bty = "n", pch = 19, xlim = xlim,
       xlab = xlab, ylab = ylab, las=1)
  mtext(main, adj = 0, cex = 1.2)
  
  # Keeping purple dots in place (no adjustment)
  for(i in 1:length(xvals)){
    points(rep(xvals[i], nsim), df[101:200,i], col = "#7570b3", pch = 19)
  }
  
  # Move black dots further to the left (increase offset)
  for(i in 1:length(xvals)){
    points(rep(xvals[i] - 1, nsim), df[1:100,i], col = "black", pch = 19)
  }
  
  # Move orange dots further to the right (increase offset)
  for(i in 1:length(xvals)){
    points(rep(xvals[i] + 1, nsim), df[201:nrow(df),i], col = "#d95f02", pch = 19)
  }
  
  # The lines remain unchanged
  lines(xvals, log10(colMeans(10^df[101:200,])), col = "#7570b3", lty = 2)
  lines(xvals, colMeans(df[101:200,]), col = "#7570b3", lty = 1)
  
  lines(xvals, log10(colMeans(10^df[1:100,])), col = "black", lty = 2)
  lines(xvals, colMeans(df[1:100,]), col = "black", lty = 1)
  
  lines(xvals, log10(colMeans(10^df[201:nrow(df),])), col = "#d95f02", lty = 2)
  lines(xvals, colMeans(df[201:nrow(df),]), col = "#d95f02", lty = 1)
}

################## plotting the output of complete_LD.sh script

pdf("complete_LD.pdf", width = 8, height = 6)
par(mfrow = c(2,2), mar = c(4.1, 4.1, 1.1, 1.1), mgp = c(2.2, 0.8, 0))

base_directory <- "complete_LD_n20"
df_list <- run_summary(base_directory,kvals,nval=20)
compressed_dataframes <- lapply(df_list, compress_and_modify)
combined_df<- do.call(cbind, compressed_dataframes)
plot_results(kvals, combined_df,
             "number of linked loci", "log likelihood ratio (base 10)", "C",xlim= range(c(0,kvals)))

base_directory <- "complete_LD_n100"
df_list <- run_summary(base_directory,kvals,nval=100)
compressed_dataframes <- lapply(df_list, compress_and_modify)
combined_df <- do.call(cbind, compressed_dataframes)
plot_results(kvals, combined_df,
             "number of linked loci", "log likelihood ratio (base 10)", "C",xlim= range(c(0,kvals)))

dev.off()
