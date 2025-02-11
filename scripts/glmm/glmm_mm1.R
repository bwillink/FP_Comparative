# Load packages
x <-
  c(
    "coda",
    "ape",
    "geiger",
    "phangorn",
    "phytools",
    "MCMCglmm", 
    "MASS"
  )

lapply(x, function(y) {
  # load package
  try(require(y, character.only = T), silent = T)
})

# Read in trees
all_tres<-read.tree(file ="../data/processed/HiSSE/Posterior_trees.trees", keep.multi = T)

# sample 600 posterior trees
tres <- random.trees<-sample(all_tres,size=600)

# adjust rounding so R recognizes the trees as ultrametric
for (i in 1:length(tres)) {
  temp <- nnls.tree(cophenetic(tres[[i]]), tres[[i]], rooted = TRUE)
  tres[[i]] <- temp
}

# Read in lit data
dat<-read.csv("../data/processed/glmm/Eco_lit_dat.csv", sep = ",", header = T, na.strings = c("", "NA"), colClasses = rep("factor", 4))
# taxa to row names
rownames(dat)<-dat$Taxon

### general data description

## data prep
# make ambiguous habitats NA
dat$Hab[dat$Hab == "C/O"] <- NA
dat$Hab <- droplevels(dat$Hab)

# drop species missing either ecological variable 
dat <- dat[complete.cases(dat[,2:3]) == T,] 

# check taxon labels match
matching <- name.check(tres[[1]], dat)

# prune trees to match data
PrunTres<-list()
for (i in 1:length(tres)){
  tempTre<-drop.tip(tres[[i]], matching$tree_not_data)
  PrunTres[[i]]<-tempTre
}

### Run the model

# we need a matrix of phylogenetic distances to account for shared evolutionary history
inv.phylo <- inverseA(PrunTres[[1]], nodes = "TIPS", scale = TRUE)

# residual covariance matrix 1/J*(I+J), I and J are identity and unit matrices respectively
IJ <- (1/3) * (diag(2) + matrix(1, 2, 2))

# kronecker prior for fixed effect close to being flat for the 2 way:
#(0 vs.1 , 1 vs. 2 and 0 vs. 2) marginal probabilities within each fixed effect
pr1 <- list(
  B = list(mu = rep(0, 8), V = diag(8) * (1.7 + pi ^ 2 / 3)),
  R = list(V = IJ, fix = 1),
  G = list(G1 = list(
    V = IJ,
    nu = 1000,
    alpha.mu = rep(0, 2),
    alpha.V = diag(2)
  ))
)

#Run the model without intercept to estimate the probability of each state in each combination of latitude and habitat
#starting value
mm1.start <-
  MCMCglmm(
    FemState ~ at.level(Hab, 'O'):at.level(Lat, 'Trp'):trait +
      at.level(Hab, 'O'):at.level(Lat, 'Tmp'):trait +
      at.level(Hab, 'C'):at.level(Lat, 'Trp'):trait +
      at.level(Hab, 'C'):at.level(Lat, 'Tmp'):trait - 1,
    random = ~ us(trait):Taxon,
    rcov = ~ us(trait):units,
    ginverse = list(Taxon = inv.phylo$Ainv),
    family = "categorical",
    data = dat,
    prior = pr1,
    pl = TRUE,
    slice = TRUE,
    nitt = 1000,
    burnin = 0,
    thin = 1,
    verbose = T
  )

# starting point after pre burn in
mm1.multiphyMCMC <- mm1.start

# run the model accounting for phylogenetic uncertainty 
for (i in 1:500) {
  IN.tree <- inverseA(PrunTres[[i]], nodes = "TIPS")
  start <- list(
    Liab = mm1.multiphyMCMC$Liab[1,],
    R = list(R1 = matrix(
      ncol = 2, nrow = 2, mm1.multiphyMCMC$VCV[1, 5:8]
    )),
    G = list(G1 = matrix(
      ncol = 2, nrow = 2, mm1.multiphyMCMC$VCV[1, c(1:4)]
    ))
  )
  
  mm1.multiphyMCMC <-
    MCMCglmm(
      FemState ~ at.level(Hab, 'O'):at.level(Lat, 'Trp'):trait +
        at.level(Hab, 'O'):at.level(Lat, 'Tmp'):trait +
        at.level(Hab, 'C'):at.level(Lat, 'Trp'):trait +
        at.level(Hab, 'C'):at.level(Lat, 'Tmp'):trait -
        1,
      random = ~ us(trait):Taxon,
      rcov = ~ us(trait):units,
      ginverse = list(Taxon = inv.phylo$Ainv),
      family = "categorical",
      data = dat,
      prior = pr1,
      pl = TRUE,
      slice = TRUE,
      nitt = 1000,
      thin = 1,
      burnin = 999,
      start = start,
      verbose = TRUE
    )
  
  if (i > 1) {
    mm1.start$VCV[i - 1, ] <- mm1.multiphyMCMC$VCV[1, ]
    mm1.start$Sol[i - 1, ] <- mm1.multiphyMCMC$Sol[1, ]
    mm1.start$Liab[i - 1, ] <- mm1.multiphyMCMC$Liab[1, ]
  }
  
  print(i)
}

write.table(mm1.start$VCV[1:500,], file = "../output/glmm/mm1_VCV.tsv", quote = F, row.names = "", sep = "\t")
write.table(mm1.start$Sol[1:500,], file = "../output/glmm/mm1_Sol.tsv", quote = F, row.names = "", sep = "\t")
write.table(mm1.start$Liab[1:500,], file = "../output/glmm/mm1_Liab.tsv", quote = F, row.names = "", sep = "\t")
