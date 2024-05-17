#5/14/24
#Code to simulate and analyze haploid individuals
#under the IBDGem model.


#########################################################
#Table 1 and example

set.seed(8675309)

n <- 5
k <- 10
p <- 1/2
eps <- 0.02

#make a matrix holding a reference panel, individuals in rows
sim.ref <- function(n,k,p){
  matrix(rbinom(n*k, 1, p), nrow = n)
}

PD.I <- (1-eps)^k #probability that the haploid of interest generates reads that match
l <- 0:k
PD.Ls <- eps^l * (1-eps)^(k-l) # probability of generating reads that match given L (0:k) mismatches
P.Ls <- dbinom(l, k, p) #probability of l differences with POI
PD.U <- sum(PD.Ls * P.Ls) #probability of data from unrelated ind
PDU.12 <- (eps*p + (1-eps)*(1-p))^k

#probability of observing the data from a randomly
#drawn individual from the reference panel
PDR <- function(ref, eps){
  Ls <- rowSums(ref)
  PD.Ri <- eps^Ls * (1-eps)^(ncol(ref)-Ls)
  mean(PD.Ri)
}

PDR.log10 <- function(ref, eps){
  Ls <- rowSums(ref)
  PD.Ri.log10 <- Ls*log10(eps) + (ncol(ref) - Ls)*log10(1-eps) 
  log10(mean(10^PD.Ri.log10))
}


ref <- sim.ref(n, k, p) #simulate a reference panel
PD.R <- PDR(ref, eps) #compute probability of data given each panel member

PD.I
PD.U
PD.R
PD.I/PD.U #forensic LR
PD.I/PD.R #IBDGem LR

##################################################
#Figure 1

#Takes a reference panel size n, a number of loci k,
#sequencing error rate eps, frequency of the allele not 
#possessed by the individual of interest p, and number of simulations
#to run nsim.
#Computes the probability of observing 1 read at each
#site matching the individual of interest under three conditions:
#first, conditional on the reads coming from the individual of interest (PD.I)
#second, conditional on the reads coming a random unrelated individual
#from the population (PD.U),
#and third, conditional on the reads coming from a random individual
#from the reference panel (PDRs) for each of nsim simulated reference panels.
#Returns the likelihood ratio for unrelated individual from the population,
#and the nsim likelihood ratios for individual from the ref panel.
getsims <- function(n, k, eps, p = 0.5, nsim = 100){
  l <- 0:k
  PD.I <- (1-eps)^k
  PD.Ls <- eps^l * (1-eps)^(k-l) # probability of generating reads that match given L (0:k) mismatches
  P.Ls <- dbinom(l, k, p) #probability of l differences with POI
  PD.U <- sum(PD.Ls * P.Ls) #probability of data from unrelated ind
  PDRs <- numeric(nsim)
  for(i in 1:nsim){
    ref <- sim.ref(n, k, p)
    PDRs[i] <- PDR(ref, eps) 
  }
  return(c( PD.I/PD.U, PD.I/PDRs ))
}

#A version of getsims that works on a log10 scale for numerical stability
getsims.log10 <- function(n, k, eps, p = 0.5, nsim = 100){
  l <- 0:k
  PD.I.log10 <- k * log10(1-eps)
  PD.Ls.log10 <- l * log10(eps) + (k-l)*log10(1-eps)  # probability of generating reads that match given L (0:k) mismatches
  P.Ls.log10 <- dbinom(l, k, p, log = TRUE)/log(10) #probability of l differences with POI
  PD.U.log10 <- log10(sum(10^(PD.Ls.log10 + P.Ls.log10))) #probability of data from unrelated ind
  PDRs.log10 <- numeric(nsim)
  for(i in 1:nsim){
    ref <- sim.ref(n, k, p)
    PDRs.log10[i] <- PDR.log10(ref, eps) 
  }
  return(c( PD.I.log10 - PD.U.log10, PD.I.log10 - PDRs.log10 ))
}

# Function to simulate data
run_simulations <- function(vals, param = "n",epsilon,kval) {
  lr.f <- numeric(length(vals))
  lr.n <- matrix(ncol = length(vals), nrow = n_sim)
  
  for(i in 1:length(vals)){
    if (param == "n") {
      n <- vals[i]
      sim <- getsims.log10(n, k = kval, eps = epsilon, p = 0.5, nsim = n_sim)
    } else if (param == "p") {
      p <- vals[i]
      sim <- getsims.log10(n = 100, k = kval, eps = epsilon, p = p, nsim = n_sim)
    } else if (param == "k") {
      k <- vals[i]
      sim <- getsims.log10(n = 100, k = k, eps = epsilon, p = 0.5, nsim = n_sim)
    } else if (param == "eps") {
      eps <- vals[i]
      sim <- getsims.log10(n = 100, k = kval, eps = eps, p = 0.5, nsim = n_sim)
    } else {
      stop("Invalid parameter specified")
    }
    
    lr.f[i] <- sim[1]
    lr.n[,i] <- sim[2:length(sim)]
  }
  
  return(list(lr.f = lr.f, lr.n = lr.n))
}

# Function to plot results
plot_results <- function(xvals, lr.f, lr.n, xlab, ylab, main,xlim=NULL) {
  plot(xvals, lr.f, ylim = range(c(0, lr.n)), bty = "n", pch = 19, 
       xlab = xlab, ylab = ylab, las = 1,xlim=xlim)
  mtext(main, adj = 0, cex = 1.2)
  for (i in 1:length(xvals)) {
    points(rep(xvals[i], nsim), lr.n[, i], col = "#7570b3", pch = 19)
  }
  lines(xvals, log10(colMeans(10^lr.n)), col = "#7570b3", lty = 2)
  lines(xvals, colMeans(lr.n), col = "#7570b3", lty = 1)
  points(xvals, lr.f, pch = 19, col = "black")
}

# Values of n, p, k, and epsilon to plot in figures
nvals <- c(10, 50, 100, 200, 400, 600, 800, 1000)
pvals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
kvals <- c(1, 25, 50, 75, 100)
epsvals <- c(1e-4, 1e-3, 1e-2, .02, .05, 0.1, 0.2)
n_sim <- 100

pdf("Fig1.pdf", width = 8, height = 6)
par(mfrow = c(2, 2), mar = c(4.1, 4.1, 1.1, 1.1), mgp = c(2.2, 0.8, 0))

# Plot reference panel size results
result_nvals <- run_simulations(nvals, param = "n",eps=0.02,kval=50)
plot_results(nvals, result_nvals$lr.f, 
             result_nvals$lr.n, 
             "reference panel size", "log likelihood ratio (base 10)", "A")

# Plot allele frequency results
result_pvals <- run_simulations(pvals, param = "p",eps=0.02,kval=50)
plot_results(1 - pvals, result_pvals$lr.f, 
             result_pvals$lr.n, 
             "frequency of individual-of-interest allele", 
             "log likelihood ratio (base 10)", "B",xlim = c(0,1))

legend("topright", pch = c(19,19,26,26), col = c("black", "#7570b3", "#7570b3", "#7570b3"), lty = c(0,0,1,2), legend = c("Standard LR", "IBDGem LR (simulated)", "mean log(IBDGem LR)", "mean IBDGem LR"), bty = "n")

# Plot number of loci results
result_kvals <- run_simulations(kvals, param = "k",eps=0.02,kval=50)
plot_results(kvals, result_kvals$lr.f, 
             result_kvals$lr.n, 
             "number of loci", "log likelihood ratio (base 10)", "C")

# Plot error rate results
result_epsvals <- run_simulations(epsvals, param = "eps",eps=0.02,kval=50)
plot_results(log10(epsvals), result_epsvals$lr.f, 
             result_epsvals$lr.n, 
             "log sequencing error rate (base 10)", 
             "log likelihood ratio (base 10)", "D")

dev.off()


##################################################
#Figure 2 --- similar to figure 1, but with smaller numbers of loci

set.seed(8675309)

nvals <- c(1:16)*25
pvals <- c(0.1, 0.3, 0.5, 0.7, 0.9)
kvals <- c(1:8)
epsvals <- c(1e-4, 1e-3, 1e-2, .02, .05, 0.1, 0.2)

pdf("Fig2.pdf", width = 8, height = 6)
par(mfrow = c(2, 2), mar = c(4.1, 4.1, 1.1, 1.1), mgp = c(2.2, 0.8, 0))

# Plot reference panel size results
result_nvals <- run_simulations(nvals, param = "n",eps=0.001,kval=6)
plot_results(nvals, result_nvals$lr.f, 
             result_nvals$lr.n, 
             "reference panel size", "log likelihood ratio (base 10)", "A",xlim=c(0,400))

# Plot allele frequency results
result_pvals <- run_simulations(pvals, param = "p",eps=0.001,kval=6)
plot_results(1 - pvals, result_pvals$lr.f, 
             result_pvals$lr.n, 
             "frequency of individual-of-interest allele", 
             "log likelihood ratio (base 10)", "B",xlim = c(0,1))

legend("topright", pch = c(19,19,26,26), col = c("black", "#7570b3", "#7570b3", "#7570b3"), lty = c(0,0,1,2), legend = c("Standard LR", "IBDGem LR (simulated)", "mean log(IBDGem LR)", "mean IBDGem LR"), bty = "n")

# Plot number of loci results
result_kvals <- run_simulations(kvals, param = "k",eps=0.001,kval=6)
plot_results(kvals, result_kvals$lr.f, 
             result_kvals$lr.n, 
             "number of loci", "log likelihood ratio (base 10)", "C")

# Plot error rate results
result_epsvals <- run_simulations(epsvals, param = "eps",eps=0.001,kval=6)
plot_results(log10(epsvals), result_epsvals$lr.f, 
             result_epsvals$lr.n, 
             "log sequencing error rate (base 10)", 
             "log likelihood ratio (base 10)", "D")

dev.off()

#Figure 3---haploid LD

#A version of getsims for the 'perfect LD blocks of m loci' model
#in addition to likelihood ratios based on a reference panel and based 
#on assuming unlinked loci, computes the correct LR accounting for LD.
getsims.log10.mcopies <- function(n, k, eps, m = 1, p = 0.5, nsim = 100){
  l <- 0:k
  PD.I.log10 <- k * m * log10(1-eps)
  PD.Ls.log10 <- m * l * log10(eps) + m*(k-l)*log10(1-eps)  # probability of generating reads that match given L (0:k) mismatches
  P.Ls.log10 <- dbinom(l, k, p, log = TRUE)/log(10) #probability of l differences with POI
  PD.U.log10 <- log10(sum(10^(PD.Ls.log10 + P.Ls.log10))) #probability of data from unrelated ind
  
  P.Ls.log10.unlinked <- dbinom(m*l, m*k, p, log = TRUE)/log(10) 
  PD.U.log10.unlinked <- log10(sum(10^(PD.Ls.log10 + P.Ls.log10.unlinked)))
  PDRs.log10 <- numeric(nsim)
  for(i in 1:nsim){
    ref <- sim.ref(n, k, p)
    ref.rep <- matrix(rep(ref, m), nrow = nrow(ref))
    PDRs.log10[i] <- PDR.log10(ref, eps) 
  }
  return(c( PD.I.log10 - PD.U.log10, PD.I.log10 - PD.U.log10.unlinked, PD.I.log10 - PDRs.log10 ))
}

#simulate a reference panel based on the "caterpillar" complete-LD model
sim.ref.caterpillar <- function(n, k){
  hapnum <- sample(0:k, n , replace=T)
  makehap <- function(hn){
    c(rep(1, hn), rep(0,k-hn))
  }
  t(sapply(hapnum, makehap))
}

#A version of getsims() for the caterpillar complete-LD model. 
#In addition to likelihood ratios  based on a reference panel and based 
#on assuming unlinked loci, also computes the correct 
#population-level likelihood ratio based on haplotype frequencies.
getsims.log10.caterpillar <- function(n, k, eps, nsim = 100){
  l <- 0:k
  PD.I.log10 <- k * log10(1-eps)
  PD.Ls.log10 <-  l * log10(eps) + (k-l)*log10(1-eps)  # probability of generating reads that match given L (0:k) mismatches
  P.Ls.log10 <- log10(rep(1/(k+1), k+1)) #probability of l differences with POI
  PD.U.log10 <- log10(sum(10^(PD.Ls.log10 + P.Ls.log10))) #probability of data from unrelated ind
  i <- 1:k
  PD.U.log10.unlinked <- sum(log10(i*(1-eps)/(k+1) + eps*(k+1-i)/(k+1) ))
  PDRs.log10 <- numeric(nsim)
  for(i in 1:nsim){
    ref <- sim.ref.caterpillar(n, k)
    PDRs.log10[i] <- PDR.log10(ref, eps) 
  }
  return(c( PD.I.log10 - PD.U.log10, PD.I.log10 - PD.U.log10.unlinked, PD.I.log10 - PDRs.log10 ))
}


set.seed(8675309)

#function to simulate the data
run_simulations_2 <- function(vals,kval=NULL,n,param) {
  if (param == "m"){
    lr.f <- numeric(length(vals))
    lr.unlinked <- numeric(length(vals))
    lr.n <- matrix(ncol = length(vals), nrow = nsim)
    for(i in 1:length(vals)){
      m <- vals[i]
      sim <- getsims.log10.mcopies(n = n, k = kval, eps = 0.02, p = 0.5, m = vals[i], nsim = nsim)
      lr.f[i] <- sim[1]
      lr.unlinked[i] <- sim[2]
      lr.n[,i] <- sim[3:length(sim)]
    }
    
    return(list(lr.f = lr.f, lr.unlinked = lr.unlinked,lr.n = lr.n))    
  }
  else{
    lr.f <- numeric(length(vals))
    lr.unlinked <- numeric(length(vals))
    lr.n <- matrix(ncol = length(vals), nrow = nsim)
    for(i in 1:length(vals)){
      k <- vals[i]
      sim <- getsims.log10.caterpillar(n = n, k = k, eps = 0.02, nsim = nsim)
      lr.f[i] <- sim[1]
      lr.unlinked[i] <- sim[2]
      lr.n[,i] <- sim[3:length(sim)]
    }
    return(list(lr.f = lr.f, lr.unlinked = lr.unlinked,lr.n = lr.n))    
  }
}

# Function to plot results
plot_results_2 <- function(xvals, lr.f, lr.n, lr.unlinked,xlab, ylab, main,xlim=NULL) {
  plot(xvals, lr.f, ylim = range(c(0,lr.n, lr.unlinked)), bty = "n", pch = 19, 
       xlab = xlab, ylab = ylab, las = 1,xlim=xlim,xaxt="n")
  mtext(main, adj = 0, cex = 1.2)
  for (i in 1:length(xvals)) {
    points(rep(xvals[i], nsim), lr.n[, i], col = "#7570b3", pch = 19)
  }
  lines(xvals, log10(colMeans(10^lr.n)), col = "#7570b3", lty = 2)
  lines(xvals, colMeans(lr.n), col = "#7570b3", lty = 1)
  points(xvals, lr.unlinked, pch = 19, col = "#d95f02")
  points(xvals, lr.f, pch = 19, col = "black")
}


pdf("Fig3.pdf", width = 8, height = 6)
par(mfrow = c(2,2), mar = c(4.1, 4.1, 1.1, 1.1), mgp = c(2.2, 0.8, 0))

mvals <- 1:4
results_mvals_50 <- run_simulations_2(mvals,n=100,k=50,param="m")

plot_results_2(mvals, results_mvals_50$lr.f, 
               results_mvals_50$lr.n, 
               results_mvals_50$lr.unlinked,
               "length of perfect-LD blocks", "log likelihood ratio (base 10)", "A")
axis(1, at = c(1,2,3,4))


#Repeat, but with lower k
results_mvals_5 <- run_simulations_2(mvals,n=100,k=5,param="m")
plot_results_2(mvals, results_mvals_5$lr.f, 
               results_mvals_5$lr.n, 
               results_mvals_5$lr.unlinked,
               "length of perfect-LD blocks", "log likelihood ratio (base 10)", "B")
axis(1, at = c(1,2,3,4))


#Panel with caterpillar/graded sim
kvals <- c((1:10)*10)

#caterpillar, changing k
results_kvals <- run_simulations_2(kvals,n=20,kval=NULL,param="k")
plot_results_2(kvals, results_kvals$lr.f, 
               results_kvals$lr.n, 
               results_kvals$lr.unlinked,
               "number of linked loci", "log likelihood ratio (base 10)", "C",xlim=c(0,max(kvals)))

legend("topleft", pch = c(19,19,26, 26, 19), col = c("black", "#7570b3", "#7570b3", "#7570b3", "#d95f02"), lty = c(0,0,1,2,0), legend = c("Standard LR", "IBDGem LR (simulated)", "mean log(IBDGem LR)", "mean IBDGem LR", "LR assuming LE"), bty = "n")
axis(1, at = seq(0, 100, by = 20))

#caterpillar, changing k, higher n
results_kvals_higher_n <- run_simulations_2(kvals,n=100,kval=NULL,param="k")
plot_results_2(kvals, results_kvals_higher_n$lr.f, 
               results_kvals_higher_n$lr.n, 
               results_kvals_higher_n$lr.unlinked,
               "number of linked loci", "log likelihood ratio (base 10)", "D",xlim=c(0,max(kvals)))

axis(1, at = seq(0, 100, by = 20))

dev.off()



