#10/15/24
#Code to simulate and analyze diploid individuals
#under the IBDGem model under perfect LD.

set.seed(304804)

############################### simulate reference panels

getsims <- function(n, k, eps, nsim,number_of_reads,snp_matrix,order_number,target_genotype_k,reads,xvals,loci_frequency){
  
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
  freq <- loci_frequency
  
  # Create a data frame with CHR, POS, and FREQ
  my_table <- data.frame(CHR = chr, POS = pos, FREQ = freq)
  
  duplicate_snp_matrix <- function(snp_matrix, copy) {
    new_rows <- list()  # List to store the original and new rows
    
    for (j in 1:nrow(snp_matrix)) {
      # Add the original SNP row
      original_row <- snp_matrix[j, ]
      new_rows[[length(new_rows) + 1]] <- original_row
      
      # If it's not the last row, create the additional SNPs between current and next row
      if (j < nrow(snp_matrix)) {
        current_id <- as.numeric(sub("SNP", "", original_row[1]))
        current_pos <- as.numeric(original_row[2])
        next_row <- snp_matrix[j + 1, ]
        next_id <- as.numeric(sub("SNP", "", next_row[1]))
        next_pos <- as.numeric(next_row[2])
        
        # Generate the new SNP rows based on the 'copy' value
        for (k in 1:(copy - 1)) {
          # Interpolate the new SNP ID and position
          new_snp_id <- current_id + k * (next_id - current_id) / copy
          new_position <- current_pos + k * (next_pos - current_pos) / copy
          
          # Create a new SNP row
          new_snp_row <- c(paste0("SNP", round(new_snp_id)), round(new_position), original_row[3], original_row[4])
          new_rows[[length(new_rows) + 1]] <- new_snp_row
        }
      }
      
      # If it's the last row, create additional SNPs beyond the last original SNP
      if (j == nrow(snp_matrix)) {
        current_id <- as.numeric(sub("SNP", "", original_row[1]))
        current_pos <- as.numeric(original_row[2])
        
        # Generate the new SNP rows beyond the last SNP based on the 'copy' value
        for (k in 1:(copy - 1)) {
          # Generate SNP ID and position with fixed intervals
          new_snp_id <- current_id + k
          new_position <- current_pos + k * 100  # Assuming position increments by 100
          
          # Create a new SNP row
          new_snp_row <- c(paste0("SNP", round(new_snp_id)), round(new_position), original_row[3], original_row[4])
          new_rows[[length(new_rows) + 1]] <- new_snp_row
        }
      }
    }
    
    # Combine the list of rows into a matrix and return
    duplicated_matrix <- do.call(rbind, new_rows)
    return(duplicated_matrix)
  }
  
  for (copy in mvals){
    new_sub_directory_name <- paste0("n", n, "_k", k,"_r",number_of_reads,"_order",order_number,"_m",copy)
    new_sub_directory <- file.path(new_directory,new_sub_directory_name)
    dir.create(new_sub_directory)
    setwd(new_sub_directory)
    for (msim in 1:nsim){
      # create appropriate directories
      new_sub_sub_directory_name <- paste0("n", n, "_k", k,"_r",number_of_reads,"_order",order_number,"_m",copy,"_sim",msim)
      new_sub_sub_directory <- file.path(new_sub_directory, new_sub_sub_directory_name)
      dir.create(new_sub_sub_directory)
      
      # Set the new directory as the working directory
      setwd(new_sub_sub_directory)
      
      no_redundant_directory_name <- paste0("n", n, "_k", k,"_r",number_of_reads,"_order",order_number,"_m",copy,"_sim",msim,"_no_redundant")
      no_redundant_direct <- file.path(new_sub_sub_directory, no_redundant_directory_name)
      
      redundant_directory_name <- paste0("n", n, "_k", k,"_r",number_of_reads,"_order",order_number,"_m",copy,"_sim",msim,"_redundant")
      redundant_direct <- file.path(new_sub_sub_directory, redundant_directory_name)
      
      dir.create(no_redundant_direct)
      dir.create(redundant_direct)
      
      #create directory to put in IBDGem results LD
      dir.create(file.path(no_redundant_direct, "output_LD"))
      #create directory to put in IBDGem results non-LD
      dir.create(file.path(no_redundant_direct, "output_non_LD"))
      
      #create directory to put in IBDGem results LD
      dir.create(file.path(redundant_direct, "output_LD"))
      #create directory to put in IBDGem results non-LD
      dir.create(file.path(redundant_direct, "output_non_LD"))
      
      #to generate .hap file
      #generate one without the redundant markers
      genotypes <- cbind(do.call(cbind, lapply(1:(n-1), function(x) generate_genotype_panel(k,loci_frequency)$ref.test)), target_genotype_k)
      
      setwd(no_redundant_direct)
      # Do this for the non-redundant matrix first 
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
      bases <- matrix(nrow = nrow(reads), ncol = copy*2)
      
      for (i in 1:nrow(reads)){
        for (j in 1:(2*copy)){
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
      
      basequals <- paste(rep(symbol, 2*copy), collapse = "")
      alignmapping = paste(rep("!", 2*copy), collapse = "")
      pileup_data <- data.frame(
        CHROM = rep("chr1", k),
        POS = snp_positions,
        REF = reference_alleles,
        DEPTH = rep(2*copy, k),
        BASES = concatenated_column, 
        BASEQUALS = rep(basequals, k),
        ALIGNMAPPING = rep(alignmapping, k)
      )
      
      write.table(pileup_data, file = "unknowndna.pileup", sep = "\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
      
      # then do this for the genotypes matrix with the redundant markers 
      
      setwd(redundant_direct)
      if (copy >1){
        
        # Do this for the non-redundant matrix first 
        #generate a hap file with the m-1 copies of the markers
        genotypes_redundant <- genotypes[rep(1:nrow(genotypes), each = copy), ]
        write.table(genotypes_redundant, file = "test.hap", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        #save the .legend file
        duplicated_snp_matrix <- duplicate_snp_matrix(snp_matrix, copy)
        write.table(duplicated_snp_matrix, file = "test.legend", sep = "\t", quote = FALSE, row.names = FALSE)
        
        #save the .indv file
        writeLines(samples, "test.indv")
        
        # save the allele frequencies file after updating it to match the duplication
        duplicated_positions <- as.numeric(duplicated_snp_matrix[, 2])
        duplicated_my_table <- do.call(rbind, replicate(copy, my_table, simplify = FALSE))
        duplicated_my_table$POS <- duplicated_positions
        write.table(duplicated_my_table, file = "output.frq", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
        
        # to generate the .pileup file
        snp_positions <- duplicated_snp_matrix[, "pos"]
        reference_alleles <- duplicated_snp_matrix[, "allele0"] 
        
        #convert the reads to the bases used in pileup format
        #if 0 (i.e. matches reference base), then "."
        #else (i.e. mismatches the reference base), then either A,C,G, or T
        
        # Create a new matrix with the same dimensions as reads_matrix
        # duplicate reads first
        bases <- matrix(nrow = copy*nrow(reads), ncol = 2)
        row_idx <- 1
        for (i in 1:nrow(reads)) {
          for (j in seq(1, (2*copy), by=2)) {
            bases[row_idx, ] <- reads[i, j:(j+1)]
            row_idx <- row_idx + 1
          }
        }
        
        for (i in 1:nrow(bases)) {
          for (j in 1:ncol(bases)) {
            if (bases[i, j] == 0) {
              bases[i, j] <- "."
            } else if (bases[i, j] == 1) {
              bases[i, j] <- duplicated_snp_matrix[i, "allele1"]
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
        basequals <- paste(rep(symbol, 2), collapse = "")
        alignmapping = paste(rep("!", 2), collapse = "")
        pileup_data <- data.frame(
          CHROM = rep("chr1", k*copy),
          POS = snp_positions,
          REF = reference_alleles,
          DEPTH = rep(2, k*copy),
          BASES = concatenated_column, 
          BASEQUALS = rep(basequals, k*copy),
          ALIGNMAPPING = rep(alignmapping, k*copy)
        )
        
        write.table(pileup_data, file = "unknowndna.pileup", sep = "\t", row.names = FALSE, col.names = FALSE,quote = FALSE)
      }
      else{
        files_to_copy <- list.files(no_redundant_direct, full.names = TRUE)
        file.copy(files_to_copy, redundant_direct, overwrite = TRUE)
      }
      setwd(new_sub_directory)
    }
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

# generating genotype(s) that obey(s) allele frequencies derived from neutral SFS

#generate n.loci allele frequencies from neutral SFS that
#would arise for a sample of size n.step
gen.neut.sfs <- function(n.loci, n.step = 100){
  probs.uw <- 1/(1:(n.step-1))
  probs <- probs.uw/sum(probs.uw)
  cdf <- cumsum(probs)
  percs <- runif(n.loci, 0 , 1)
  get.cdfind <- function(perc){
    (sum(cdf <= perc) + 1)/n.step
  }
  sapply(percs, get.cdfind)
}

#simulates a gentoype with alleles at frequencies
#in the vector ps, all in linkage equilibrium
sim.ref.varfreq.diploid <- function(n, ps) {
  hap1 <- matrix(rbinom(n * length(ps), 1, ps), ncol = n)
  hap2 <- matrix(rbinom(n * length(ps), 1, ps), ncol = n)
  diploid_genotypes <- cbind(hap1,hap2)
  return(diploid_genotypes)
}
generate_genotype <- function(k) {
  if (k == 1) {
    # For k = 1, just generate a single genotype randomly
    return(list(ref.test =matrix( sample( c(0,1), 2, replace = TRUE ), nrow = 1, ncol = 2 ),ps=0.5))
  } else {
    ps <- gen.neut.sfs(k, 100)
    ref.test <- sim.ref.varfreq.diploid(1, ps)
    return(list(ref.test=ref.test,ps=ps))
  }
}

#simulates a reference panel of n genotypes with alleles at frequencies
#in the vector ps, all in linkage equilibrium
generate_genotype_panel <- function(k,ps) {
  if (k == 1) {
    # For k = 1, just generate a single genotype randomly
    return(list(ref.test =matrix( sample( c(0,1), 2, replace = TRUE ), nrow = 1, ncol = 2 ),ps=0.5))
  } else {
    ref.test <- sim.ref.varfreq.diploid(1, ps)
    return(list(ref.test=ref.test,ps=ps))
  }
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

############################### simulate the genotype and reads for the target individuals with the default parameters being k = 100.
# Generate target genotype, its allele frequency file and the snp_matrix file
result <- generate_genotype(k=100)
target_genotype_100_loci <- result$ref.test
loci_frequency <- result$ps
snp_matrix_100_loci <- simulate_legend_file(k=100)

# Sample indices uniformly at random to form the target genotype's at k=50 loci
sampled_indices <- sample(1:nrow(snp_matrix_100_loci), 50, replace=FALSE)
target_genotype_50_loci <- target_genotype_100_loci[sampled_indices, ]
#indices we extracted -- now ordered and unlisted
indices <- sort(unlist(sampled_indices,use.names=FALSE))
snp_matrix_sampled_50_loci <- snp_matrix_100_loci[indices,]
# generate the reads for the target genotype when k = 100
reads_50 <- simulate_reads_for_target(target_genotype_50_loci,eps=0.02,number_of_reads=8,k=50)
loci_frequency_50 <- loci_frequency[indices]


############################### simulate different scenarios
kvals <- c(5,50)
mvals <- 1:4
nsim <- 100

################## Starting with the cases where k is set to 50 (varying n, varying number of reads)
current_directory <- getwd()

##################  Vary the number of loci according to kvals while n is set to 100, p=0.5, eps=0.02, number of reads = 2
loci_to_sample <- c(5)
# Initialize lists to store the results
target_genotypes <- list()
snp_matrices <- list()
loci_freq_matrices <- list()

# we already have 50 loci defined
target_genotypes[["target_genotype_50_loci"]] <- target_genotype_50_loci
snp_matrices[["snp_matrix_sampled_50_loci"]] <- snp_matrix_sampled_50_loci
loci_freq_matrices[["loci_frequency_50_loci"]] <- loci_frequency_50

# Loop through each number of loci to sample
for (num_loci in loci_to_sample) {
  # Sample indices to satisfy the num_loci
  sampled_indices <- sample(1:nrow(target_genotype_100_loci), num_loci, replace=FALSE)
  
  # Combine sampled hets and homs into target genotype
  target_genotypes[[paste0("target_genotype_", num_loci, "_loci")]] <- target_genotype_100_loci[sampled_indices, ]
  
  indices <- sort(unlist(sampled_indices,use.names=FALSE))
  # Extract corresponding snp matrix
  snp_matrices[[paste0("snp_matrix_sampled_", num_loci, "_loci")]] <- snp_matrix_100_loci[indices,]
  
  # Extract corresponding allele frequencies for selected subsetted loci
  loci_freq_matrices[[paste0("loci_frequency_", num_loci, "_loci")]] <- loci_frequency[indices]
}

#simulate reads for k=5 (we already have reads for k=50)
# Define the list of k values
k_values <- c(5)

# Initialize an empty list to store the reads
reads_list <- list()

#add reads generated for 50 loci earlier 
reads_list[[as.character(50)]] <- reads_50

# Loop through each k value and simulate reads
for (k in k_values) {
  target_name <- paste0("target_genotype_", k, "_loci")
  reads_list[[as.character(k)]] <- simulate_reads_for_target(target_genotypes[[target_name]], eps = 0.02, number_of_reads = 8, k = k)
}

directory_name <- "copy_number"
directory <- file.path(current_directory, directory_name)
dir.create(directory)
setwd(directory)

for (i in 1:2) {
  k <- kvals[i]
  target_name <- paste0("target_genotype_", k, "_loci")
  snp_matrix_name <- paste0("snp_matrix_sampled_", k, "_loci")
  locifreq_name <- paste0("loci_frequency_",k,"_loci")
  getsims(n = 5000, k = k, eps = 0.02, nsim = nsim, number_of_reads = 2,
          snp_matrix = snp_matrices[[snp_matrix_name]], order_number = i, target_genotype_k = target_genotypes[[target_name]],
          reads = reads_list[[as.character(k)]],mvals,loci_frequency = loci_freq_matrices[[locifreq_name]])
}


############################### Run IBDGem on inputs from above scenarios
system("perfect_LD.sh")

############################### Analyze the output of IBDGem

################## Functions to get the likelihood ratios under both LD and non-LD modes and plot results
generate_main_directory <- function(list_name, kval, base_directory, i,x,nval) {
  main_directories <- file.path(base_directory, paste0("n",nval,"_k", kval, "_r2_order", i,"_m",x))
  return(main_directories)
}

run_summary <- function(base_directory,xvals,kval,n){
  if (kval == 5){
    i <- 1
  }
  else{
    i <- 2
  }
  list_of_dfs <- c()
  
  for (x in xvals){
    main_directory <- generate_main_directory(list_name = xvals,kval,base_directory,i,x,n)
    # List all subdirectories within the main directory
    subdirectories <- list.dirs(main_directory, recursive = FALSE)
    
    # Function to find the "output_LD" sub-subdirectory within a subdirectory
    
    find_output_LD <- function(directory) {
      output_LD_path <- file.path(directory,paste0(basename(directory),"_redundant"), "output_LD")
      if (file.exists(output_LD_path) && file.info(output_LD_path)$isdir) {
        return(output_LD_path)
      } else {
        return(NULL)
      }
    }
    
    # Function to find the "output_non_LD" sub-sub-directory within a sub-directory 
    
    find_output_non_LD_redundant <- function(directory) {
      output_non_LD_path <- file.path(directory, paste0(basename(directory),"_redundant"),"output_non_LD")
      if (file.exists(output_non_LD_path) && file.info(output_non_LD_path)$isdir) {
        return(output_non_LD_path)
      } else {
        return(NULL)
      }
    }
    
    find_output_non_LD_non_redundant <- function(directory) {
      output_non_LD_path <- file.path(directory,paste0(basename(directory),"_no_redundant"), "output_non_LD")
      if (file.exists(output_non_LD_path) && file.info(output_non_LD_path)$isdir) {
        return(output_non_LD_path)
      } else {
        return(NULL)
      }
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
    
    summary_data_LD_redundant <- list()
    summary_data_non_LD_redundant <- list()
    summary_data_non_LD_no_redundant <- list()
    
    # Loop through each subdirectory
    for (subdir in subdirectories) {
      output_LD_redundant <- find_output_LD(subdir)
      output_non_LD_redundant <- find_output_non_LD_redundant(subdir)
      output_non_LD_no_redundant <- find_output_non_LD_non_redundant(subdir)
      if (!is.null(output_LD_redundant)) {
        summary_data_LD_redundant[[subdir]] <- read_summary_file(output_LD_redundant)
      }
      if (!is.null(output_non_LD_redundant)) {
        summary_data_non_LD_redundant[[subdir]] <- read_summary_file(output_non_LD_redundant)
      }
      if (!is.null(output_non_LD_no_redundant)) {
        summary_data_non_LD_no_redundant[[subdir]] <- read_summary_file(output_non_LD_no_redundant)
      }
    }
    
    result_df <- data.frame(log10LR_LD_redundant = numeric(nsim),log10LR_non_LD_redundant = numeric(nsim), log10LR_non_LD_non_redundant = numeric(nsim))
    
    
    for (s in 1:length(summary_data_LD_redundant)){
      redundant_LD_LIBD0 = as.numeric(summary_data_LD_redundant[[s]][1,][4])
      redundant_LD_LIBD2 = as.numeric(summary_data_LD_redundant[[s]][1,][6])
      redundant_nonLD_LIBD0 = as.numeric(summary_data_non_LD_redundant[[s]][1,][4])
      redundant_nonLD_LIBD2 = as.numeric(summary_data_non_LD_redundant[[s]][1,][6])
      non_redundant_nonLD_LIBD0 = as.numeric(summary_data_non_LD_no_redundant[[s]][1,][4])
      non_redundant_nonLD_LIBD2 = as.numeric(summary_data_non_LD_no_redundant[[s]][1,][6])
      result_df[s,1] = log10(redundant_LD_LIBD2/redundant_LD_LIBD0)
      result_df[s,2] = log10(redundant_nonLD_LIBD2/redundant_nonLD_LIBD0)
      result_df[s,3] = log10(non_redundant_nonLD_LIBD2/non_redundant_nonLD_LIBD0)
    }
    df_name <- paste0("LR_", x)
    assign(df_name, result_df)
    list_of_dfs[[df_name]] <- result_df 
    print('finished here')
  }
  return(list_of_dfs)
}

compress_and_modify <- function(df) {
  df_compressed <- rbind(df, df, df)
  df_compressed$log10LR_LD_redundant[101:200] <- df_compressed$log10LR_non_LD_redundant[101:200]
  df_compressed$log10LR_LD_redundant[201:nrow(df_compressed)] <- df_compressed$log10LR_non_LD_non_redundant[201:nrow(df_compressed)]
  df_compressed <- df_compressed[, -c(2,3), drop = FALSE]
  return(df_compressed)
}

plot_results <- function(xvals, df, xlab, ylab, main,xlim=NULL) {
  plot(xvals, type="n",ylim = range(c(0,df)), bty = "n", pch = 19, xlim = xlim,
       xlab = xlab, ylab = ylab,las=1,xaxt="n")
  axis(1, at = seq(from = floor(min(xvals)), to = ceiling(max(xvals)), by = 1))
  mtext(main, adj = 0, cex = 1.2)
  for(i in 1:length(xvals)){
    points(rep(xvals[i], nsim), df[1:100,i], col = "#7570b3", pch = 19)
  }
  for(i in 1:length(xvals)){
    points(rep(xvals[i], nsim), df[201:nrow(df),i], col = "black", pch = 19)
  }
  for(i in 1:length(xvals)){
    points(rep(xvals[i], nsim), df[101:200,i], col = "#d95f02", pch = 19)
  }
  lines(xvals, log10(colMeans(10^df[1:100,])), col = "#7570b3", lty = 2  )
  lines(xvals, colMeans(df[1:100,]), col = "#7570b3", lty = 1  )
  
}

################## When we vary the number of loci

pdf("perfect_LD.pdf", width = 8, height = 6)
par(mfrow = c(2,2), mar = c(4.1, 4.1, 1.1, 1.1), mgp = c(2.2, 0.8, 0))

### plot for k=5
base_directory <- "/copy_number/n5000_k5_r2_order1"
df_list <- run_summary(base_directory,mvals,5,n=5000)
compressed_dataframes <- lapply(df_list, compress_and_modify)
combined_df_perfectLD_k5 <- do.call(cbind, compressed_dataframes)

plot_results(mvals, combined_df_perfectLD_k5,
             "length of perfect-LD blocks", "log likelihood ratio (base 10)", "D",xlim= range(c(1,mvals)))

### plot for k=50
base_directory <- "/copy_number/n5000_k50_r2_order2"
df_list <- run_summary(base_directory,mvals,50,n=5000)
compressed_dataframes <- lapply(df_list, compress_and_modify)
combined_df_perfectLD_k50 <- do.call(cbind, compressed_dataframes)

plot_results(mvals, combined_df_perfectLD_k50,
             "length of perfect-LD blocks", "log likelihood ratio (base 10)", "C",xlim= range(c(1,mvals)))

legend("topleft", pch = c(19,19,26, 26, 19), col = c("black", "#7570b3", "#7570b3", "#7570b3", "#d95f02"), lty = c(0,0,1,2,0), legend = c("Standard LR", "IBDGem LR (simulated)", "mean log(IBDGem LR)", "mean IBDGem LR", "LR assuming LE"), bty = "n",cex=0.9)

dev.off()





