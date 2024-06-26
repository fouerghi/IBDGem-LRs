#5/14/24
#Code to simulate and analyze diploid individuals
#under the IBDGem model.

set.seed(8675309)

############################### simulate reference panels

getsims <- function(n, k, eps, p = p, nsim,number_of_reads,snp_matrix,order_number,target_genotype_k,reads){
  
  new_directory_name <- paste0("n", n, "_k", k,"_r",number_of_reads,"_order",order_number)
  new_directory <- file.path(directory, new_directory_name)
  dir.create(new_directory)
  setwd(new_directory)
  
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
  
  for (m in 1:nsim){
    # create appropriate directories
    new_sub_directory_name <- paste0("n", n, "_k", k,"_sim",m)
    new_sub_directory <- file.path(new_directory, new_sub_directory_name)
    dir.create(new_sub_directory)
    
    # Set the new directory as the working directory
    setwd(new_sub_directory)
    
    #create directory to put in IBDGem results LD
    dir.create(file.path(new_sub_directory, "output_LD"))
    #create directory to put in IBDGem results non-LD
    dir.create(file.path(new_sub_directory, "output_non_LD"))
    
    #to generate .hap file
    genotypes <- cbind(do.call(cbind, lapply(1:(n-1), function(x) generate_genotype(k))), target_genotype_k)
    
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
    
    setwd(new_directory)
    
  }
}

############################### functions to simulate genotypes and reads

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

# function to generate the genotype
generate_genotype <- function(k) {
  if (k == 1) {
    # For k = 1, just generate a single genotype randomly
    matrix( sample( c(0,1), 2, replace = TRUE ), nrow = 1, ncol = 2 )
  } else {
    # For k > 1
    matrix <- matrix( NA, nrow = k, ncol = 2 )
    if (k %% 2 == 0) {
      # If k is even, proceed with the original logic
      half_rows <- sample( 1:k, k/2, replace = FALSE )
      matrix[half_rows, ] <- matrix( sample( c(0, 1), 2 * (k/2), replace = TRUE ), ncol = 2 )
      matrix[-half_rows, ] <- matrix( replicate( 2, sample( c(0, 1), k/2, replace = TRUE ) ), ncol = 2 )
    } else {
      # If k is odd, randomly select (k-1)/2 rows and fill them randomly
      half_rows <- sample( 1:k, (k-1)/2, replace = FALSE )
      matrix[half_rows, ] <- matrix( sample( c(0, 1), 2 * ((k-1)/2), replace = TRUE ), ncol = 2 )
      # The remaining row will be filled with 0 for the first column and 1 for the second column
      remaining_row <- setdiff(1:k, half_rows)
      matrix[remaining_row, ] <- c(0, 1)
    }
    return(matrix)
  }
}

# function to calculate the number of each type of genotype (00,10/01,11)
calculate_counts <- function(target_genotype) {
  count_01_10 <- sum(target_genotype[, 1] != target_genotype[, 2])
  count_00 <- sum(target_genotype[, 1] == 0 & target_genotype[, 2] == 0)
  count_11 <- sum(target_genotype[, 1] == 1 & target_genotype[, 2] == 1)
  
  return(list(count_01_10 = count_01_10, count_00 = count_00, count_11 = count_11))
}

# function to get the indices for k=50
sample_indices <- function(target_genotype, num_samples) {
  indices_00 <- which(target_genotype[, 1] == 0 & target_genotype[, 2] == 0)
  indices_11 <- which(target_genotype[, 1] == 1 & target_genotype[, 2] == 1)
  indices_10_01 <- which(target_genotype[, 1] != target_genotype[, 2])
  
  sampled_indices_00 <- sample(indices_00, num_samples[1])
  sampled_indices_11 <- sample(indices_11, num_samples[2])
  sampled_indices_10_01 <- sample(indices_10_01, num_samples[3])
  
  sampled_indices <- list(indices_00 = sampled_indices_00,
                          indices_11 = sampled_indices_11,
                          indices_10_01 = sampled_indices_10_01)
  
  return(sampled_indices)
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

############################### simulate the genotype and reads for the target individuals with the default parameters being k = 100 and number of reads = 2.
# Generate target genotypes
target_genotype_100_loci <- generate_genotype(k=100)

# Define the number of samples for each genotype
num_samples <- c(12, 13, 25)

# Sample indices
sampled_indices <- sample_indices(target_genotype_100_loci, num_samples)

target_genotype_50_loci <- rbind(target_genotype_100_loci[sampled_indices$indices_00, ],
                                 target_genotype_100_loci[sampled_indices$indices_11, ],
                                 target_genotype_100_loci[sampled_indices$indices_10_01, ])


# generate the snp matrix (.legend file)
snp_matrix_100_loci <- simulate_legend_file(k=100)

#indices we extracted 
indices <- sort(unlist(sampled_indices,use.names=FALSE))
snp_matrix_sampled_50_loci <- snp_matrix_100_loci[indices,]

# generate the reads for the target genotype when k = 50
reads_50 <- simulate_reads_for_target(target_genotype_50_loci,eps=0.02,number_of_reads=2,k=50)

############################### simulate different scenarios
nvals <- c(10, 50, 100, 200, 400, 600, 800, 1000)
kvals <- c(1,25, 50, 75, 100)
readvals <- c(1,2,4,6,8,10)
nsim <- 100
epsvals <- c(1e-4, 1e-3, .01, .02, .05, 0.1, 0.2)

################## Starting with the cases where k is set to 50 (varying n, varying number of reads)
current_directory <- getwd()


################## Vary the reference panel size according to nvals while k is set to 50 with number of reads exactly 2
directory_name <- "reference_panel"
directory <- file.path(current_directory, directory_name)
dir.create(directory)
setwd(directory)

for(i in 1:length(nvals)){
  n <- nvals[i]
  getsims(n, k = 50, eps = 0.02, p = 0.5, nsim = nsim,number_of_reads=2,snp_matrix_sampled_50_loci,order_number=i,target_genotype_50_loci,reads = reads_50)
}

##################  Vary the number of reads for the target genotype while n = 100 and k = 50

directory_name <- "read_numbers"
directory <- file.path(current_directory, directory_name)
dir.create(directory)
setwd(directory)

for (i in 1:length(readvals)){
  r <- readvals[i]
  generated_reads <- simulate_reads_for_target(target_genotype_50_loci,eps=0.02,number_of_reads=r,k=50)
  getsims(n=100, k = 50, eps = 0.02, p = 0.5, nsim = nsim,number_of_reads=r,snp_matrix_sampled_50_loci,order_number=i,target_genotype_50_loci,reads = generated_reads)
}

##################  Vary the number of loci according to kvals while n is set to 100, p=0.5, eps=0.02, number of reads = 2

#distinguish heterozygous loci from homozygous loci
hets_indices <- which((target_genotype_100_loci[, 1] == 0 & target_genotype_100_loci[, 2] == 1) |
                        (target_genotype_100_loci[, 1] == 1 & target_genotype_100_loci[, 2] == 0))

homs_indices <- which(target_genotype_100_loci[, 1] == target_genotype_100_loci[, 2])


loci_to_sample <- c(25, 75)
# Initialize lists to store the results
target_genotypes <- list()
snp_matrices <- list()

#sample 1 locus for target genotype
sampled_index <- sample(nrow(target_genotype_100_loci), 1)
target_genotype_1_loci <- matrix(unlist(target_genotype_100_loci[sampled_index, ]), nrow = 1, ncol = 2)
snp_matrix_sampled_1_loci <- snp_matrix_100_loci[sampled_index,,drop=FALSE]

target_genotypes[["target_genotype_1_loci"]] <- target_genotype_1_loci
snp_matrices[["snp_matrix_sampled_1_loci"]] <- snp_matrix_sampled_1_loci

# we already have 50 loci defined
target_genotypes[["target_genotype_50_loci"]] <- target_genotype_50_loci
snp_matrices[["snp_matrix_sampled_50_loci"]] <- snp_matrix_sampled_50_loci

# we already have 100 loci defined
target_genotypes[["target_genotype_100_loci"]] <- target_genotype_100_loci
snp_matrices[["snp_matrix_sampled_100_loci"]] <- snp_matrix_100_loci

# Loop through each number of loci to sample
for (num_loci in loci_to_sample) {
  # Sample indices for hets and homs
  sampled_index_hets <- sample(hets_indices, num_loci %/% 2)
  sampled_index_homs <- sample(homs_indices, num_loci - (num_loci %/% 2))
  
  # Extract sampled hets and homs
  sampled_hets <- target_genotype_100_loci[sampled_index_hets, ]
  sampled_homs <- target_genotype_100_loci[sampled_index_homs, ]
  
  # Combine sampled hets and homs into target genotype
  target_genotypes[[paste0("target_genotype_", num_loci, "_loci")]] <- rbind(sampled_hets, sampled_homs)
  
  # Combine sampled indices
  sampled_indices <- sort(c(sampled_index_hets, sampled_index_homs))
  
  # Extract corresponding snp matrix
  snp_matrices[[paste0("snp_matrix_sampled_", num_loci, "_loci")]] <- snp_matrix_100_loci[sampled_indices, , drop=FALSE]
}


#simulate reads for k=1,k=25,k=50,k=75,k=100
# Define the list of k values
k_values <- c(1, 25, 75, 100)

# Initialize an empty list to store the reads
reads_list <- list()

#add reads generated for 50 loci earlier 
reads_list[[as.character(50)]] <- reads_50

# Loop through each k value and simulate reads
for (k in k_values) {
  target_name <- paste0("target_genotype_", k, "_loci")
  reads_list[[as.character(k)]] <- simulate_reads_for_target(target_genotypes[[target_name]], eps = 0.02, number_of_reads = 2, k = k)
}

directory_name <- "loci_number"
directory <- file.path(current_directory, directory_name)
dir.create(directory)
setwd(directory)

for (i in 1:5) {
  k <- kvals[i]
  target_name <- paste0("target_genotype_", k, "_loci")
  snp_matrix_name <- paste0("snp_matrix_sampled_", k, "_loci")
  getsims(n = 100, k = k, eps = 0.02, p = 0.5, nsim = nsim, number_of_reads = 2,
          snp_matrix = snp_matrices[[snp_matrix_name]], order_number = i, target_genotype_k = target_genotypes[[target_name]],
          reads = reads_list[[as.character(k)]])
}

##################  Vary the sequencing error rates using epsvals while n is set to 100, k=50, p=1/2, number of reads = 2
directory_name <- "error_rate"
directory <- file.path(current_directory, directory_name)
dir.create(directory)
setwd(directory)

for (i in 1:length(epsvals)){
  eps <- epsvals[i]
  reads <- reads_50
  getsims(n = 100, k = 50, eps = eps, p = 0.5, nsim = nsim,number_of_reads=2,snp_matrix=snp_matrix_sampled_50_loci,order_number=i,target_genotype_k = target_genotype_50_loci,reads = reads_50)
}

############################### Run IBDGem on inputs from above scenarios
system("run_ibdgem.sh")

############################### Analyze the output of IBDGem


################## Functions to get the likelihood ratios under both LD and non-LD modes and plot results
generate_main_directory <- function(list_name, list_element, base_directory, i,type) {
  if (type == "e") {
    main_directories <- file.path(base_directory, paste0("n100_k50_r2_order", i))
  } else if (type == "k") {
    main_directories <- file.path(base_directory, paste0("n100_k", list_element, "_r2_order", i))
  } else if (type == "n") {
    main_directories <- file.path(base_directory, paste0("n", list_element, "_k50_r2_order", i))
  } else if (type == "r") {
    main_directories <- file.path(base_directory, paste0("n100_k50_r", list_element, "_order", i))
  }
  return(main_directories)
}

run_summary <- function(base_directory,xvals, type){
  i=0
  list_of_dfs <- c()
  for (x in xvals){
    i=i+1
    main_directory <- generate_main_directory(list_name = xvals,list_element = x,base_directory,i,type)
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
    
    # Function to read the "sample100.sample100.summary.txt" file
    if (type == "n"){
      n <- x
    }else{
      n <- 100
    }
    read_summary_file <- function(directory) {
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
    
    summary_data_LD <- list()
    summary_data_non_LD <- list()
    
    # Loop through each subdirectory
    for (subdir in subdirectories) {
      output_LD <- find_output_LD(subdir)
      output_non_LD <- find_output_non_LD(subdir)
      if (!is.null(output_LD)) {
        summary_data_LD[[subdir]] <- read_summary_file(output_LD)
      }
      if (!is.null(output_non_LD)) {
        summary_data_non_LD[[subdir]] <- read_summary_file(output_non_LD)
      }
    }
    
    result_df <- data.frame(log10LR_LD = numeric(nsim),log10LR_non_LD = numeric(nsim))
    
    for (s in 1:length(summary_data_LD)){
      LD_LIBD0 = as.numeric(summary_data_LD[[s]][1,][4])
      LD_LIBD2 = as.numeric(summary_data_LD[[s]][1,][6])
      nonLD_LIBD0 = as.numeric(summary_data_non_LD[[s]][1,][4])
      nonLD_LIBD2 = as.numeric(summary_data_non_LD[[s]][1,][6])
      result_df[s,1] = log10(LD_LIBD2/LD_LIBD0)
      result_df[s,2] = log10(nonLD_LIBD2/nonLD_LIBD0)
    }
    
    df_name <- paste0("LR_", x)
    assign(df_name, result_df)
    list_of_dfs[[df_name]] <- result_df 
  }
  return(list_of_dfs)
}

compress_and_modify <- function(df) {
  df_compressed <- rbind(df, df)
  df_compressed$log10LR_LD[101:nrow(df_compressed)] <- df_compressed$log10LR_non_LD[101:nrow(df_compressed)]
  df_compressed <- df_compressed[, -2, drop = FALSE]
  return(df_compressed)
}

plot_results <- function(xvals, df, xlab, ylab, main,xlim=NULL) {
  plot(xvals, type="n",ylim = range(c(0,df)), bty = "n", pch = 19, xlim = xlim,
       xlab = xlab, ylab = ylab,las=1)
  mtext(main, adj = 0, cex = 1.2)
  for(i in 1:length(xvals)){
    points(rep(xvals[i], nsim), df[1:100,i], col = "#7570b3", pch = 19)
  }
  for(i in 1:length(nvals)){
    points(rep(xvals[i], nsim), df[101:200,i], col = "black", pch = 19)
  }
  lines(xvals, log10(colMeans(10^df[1:100,])), col = "#7570b3", lty = 2  )
  lines(xvals, colMeans(df[1:100,]), col = "#7570b3", lty = 1  )
  
}

################## When we vary the panel size

base_directory <- "/Users/ouerghi1/reference_panel/"
df_list <- run_summary(base_directory,nvals,"n")
compressed_dataframes <- lapply(df_list, compress_and_modify)
combined_df <- do.call(cbind, compressed_dataframes)

pdf("Fig4.pdf", width = 8, height = 6)
par(mfrow = c(2,2), mar = c(4.1, 4.1, 1.1, 1.1), mgp = c(2.2, 0.8, 0))

plot_results(nvals, combined_df,
             "reference panel size", "log likelihood ratio (base 10)", "A",xlim=c(0,max(nvals)))

################## When we vary the number of reads

base_directory <- "read_numbers/"
df_list <- run_summary(base_directory,readvals,"r")
compressed_dataframes <- lapply(df_list, compress_and_modify)
combined_df <- do.call(cbind, compressed_dataframes)
plot_results(readvals, combined_df,
             "number of reads per site", "log likelihood ratio (base 10)", "B",xlim= range(c(0,readvals)))


################## When we vary the number of loci

base_directory <- "loci_number/"
df_list <- run_summary(base_directory,kvals,"k")
compressed_dataframes <- lapply(df_list, compress_and_modify)
combined_df <- do.call(cbind, compressed_dataframes)
plot_results(kvals, combined_df,
             "number of loci", "log likelihood ratio (base 10)", "C",xlim= range(c(0,kvals)))

legend("topleft", pch = c(19,19,26, 26), col = c("black", "#7570b3", "#7570b3", "#7570b3"), lty = c(0,0,1,2), legend = c("Standard LR", "IBDGem LR (simulated)", "mean log(IBDGem LR)", "mean IBDGem LR"), bty = "n")

################## When we vary the error rate

base_directory <- "error_rate/"
df_list <- run_summary(base_directory,epsvals,"e")
compressed_dataframes <- lapply(df_list, compress_and_modify)
combined_df <- do.call(cbind, compressed_dataframes)
plot_results(log10(epsvals), combined_df,
             "log sequencing error rate (base 10)", "log likelihood ratio (base 10)", "D",xlim= range(c(0,log10(epsvals))))


dev.off()





