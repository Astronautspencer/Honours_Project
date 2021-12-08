args <- commandArgs(trailingOnly = TRUE) 
mn <- as.numeric(args[1])
library(spatstat)
library(INLA)
library(RandomFields)
library(rgeos)
### --------------------------------------------------------------------------->> Data simulations
### the data is simulated directly using Log Gaussian Cox Process.
### 5 sets of random point patterns, 10 replications.
source("data.sim.R")
source("multi.sim.R")
source("mesh.sim.R")
set.seed(180808)
pdf(paste0("demo2-ppp-", mn, ".pdf"))
ppp <- c(5, 4, 1, 10, 1, 3,
         5, 4, 1, 2, 1, 3,
         5, 2, 1, 10, 1, 3,
         5, 2, 1, 2, 1, 3)
pppmat <- matrix(ppp, ncol = 6, byrow = TRUE)
nsim <- multi.sim(pppmat)
dev.off()

### --------------------------------------------------------------------------->> Mesh Construction 
### Mesh size influences the computational time need to fit the model,
### define the study region
domain <- 3 * cbind(c(0,1,1,0,0), c(0,0,1,1,0))
mval <- c(c(0.8, 1), 0.005, c(0.05, -0.1),
          c(0.8, 1), 0.1, c(0.05, -0.1),
          c(0.8, 1), 0.1, c(0.3, 0.8),
          c(0.8, 1), 0.1, c(0.5, 1),
          c(0.5, 1), 0.005, c(0.05, -0.1),
          c(0.5, 1), 0.1, c(0.05, -0.1),
          c(0.5, 1), 0.1, c(0.3, 0.8),
          c(0.5, 1), 0.1, c(0.5, 1),
          c(0.6, 0.8), 0.005, c(0.05, -0.1),
          c(0.6, 0.8), 0.1, c(0.05, -0.1),
          c(0.6, 0.8), 0.1, c(0.3, 0.8),
          c(0.6, 0.8), 0.1, c(0.5, 1),
          c(0.4, 0.6), 0.005, c(0.05, -0.1),
          c(0.4, 0.6), 0.1, c(0.05, -0.1),
          c(0.4, 0.6), 0.1, c(0.3, 0.8),
          c(0.4, 0.6), 0.1, c(0.5, 0.8),
          c(0.3, 0.6), 0.005, c(0.05, -0.1),
          c(0.3, 0.6), 0.1, c(0.05, -0.1),
          c(0.3, 0.6), 0.1, c(0.3, 0.8),
          c(0.3, 0.6), 0.1, c(0.5, 1),
          c(0.3, 0.4), 0.00, c(0.05, -0.1),
          c(0.3, 0.4), 0.1, c(0.05, -0.1),
          c(0.3, 0.4), 0.1, c(0.3, 0.8),
          c(0.3, 0.4), 0.1, c(0.5, 1),
          c(0.2, 0.4), 0.005, c(0.05, -0.1),
          c(0.2, 0.4), 0.1, c(0.05, -0.1),
          c(0.2, 0.4), 0.1, c(0.3, 0.8),
          c(0.2, 0.4), 0.1, c(0.5, 0.8),
          c(0.2, 0.3), 0.005, c(0.05, -0.1),
          c(0.2, 0.3), 0.1, c(0.05, -0.1),
          c(0.2, 0.3), 0.1, c(0.3, 0.8),
          c(0.2, 0.3), 0.1, c(0.5, 1),
          c(0.1, 0.2), 0.1, c(0.2, -0.1),
          c(0.1, 0.2), 0.1, c(0.2, 0.1),
          c(0.1, 0.2), 0.1, c(0.3, 0.8),
          c(0.1, 0.2), 0.1, c(0.5, 1))
mesh.mat <- matrix(mval, ncol = 5, byrow = TRUE)
nsim.mesh <- mesh.sim(domain, 
                      c(mesh.mat[mn,1],mesh.mat[mn,2]), mesh.mat[mn,3], c(mesh.mat[mn,4],mesh.mat[mn,5]))

### --------------------------------------------------------------------------->> Dual Mesh and Weights
tmp <- list()
for (l in 1:nrow(pppmat)) {
  tmp[[l]] <- vector("list", length = unique(pppmat[,1]))
}

### get, nv, the number of observations
### get, n, the number of nsimulated points 
nsim.nv <- c()
nsim.n <- tmp
nsim.dmesh <- c()
nsim.weight <- c()

for (x in 1:nrow(pppmat)){
  for (y in 1:unique(pppmat[,1])) {
    ## get the vertices for each mesh 
    nsim.nv <- nsim.mesh$mesh$n
    nsim.n[[x]][[y]] <- rep(nrow(nsim[[x]]$locations[[y]]), length(nsim.nv))
    source("book.mesh.dual.R")
    ## create the dual mesh polygons
    nsim.dmesh <- book.mesh.dual(nsim.mesh$mesh)
    ## convert domain polygon into a Spatial Polygons
    bdomainSP <- SpatialPolygons(list(Polygons(list(Polygon(domain)), '0')))
    ## compute intersection between each polygon
    nsim.weight <- sapply(1:length(nsim.dmesh),
                          function(i) {
                            if (gIntersects(nsim.dmesh[i, ], bdomainSP))
                              return(gArea(gIntersection(nsim.dmesh[i, ], bdomainSP)))
                            else return(0)
                          })
  }
}

### --------------------------------------------------------------------------->> Projection Matrices
### --------------------------------------------------------------------------->> Setup the SPDE model
nsim.y.pp <- tmp
nsim.e.pp <- tmp
nsim.dmat <- tmp
nsim.lmat <- tmp
nsim.A <- tmp
nsim.spde <- tmp
nsim.stk <- tmp
nsim.res <- tmp
nsim.time <- tmp
transf <- tmp
for (p in 1:nrow(pppmat)){
  for (q in 1:unique(pppmat[,1])){
    ## define a vector of ones of the observations and zeros for the mesh nodes
    nsim.y.pp[[p]][[q]] <- rep(0:1, c(nsim.nv, nsim.n[[p]][[q]]))
    ## define the exposure vector 
    nsim.e.pp[[p]][[q]] <- c(nsim.weight, rep(0, nsim.n[[p]][[q]]))
    ## projection matrix
    nsim.dmat[[p]][[q]] <- Diagonal(nsim.nv, rep(1, nsim.nv))
    nsim.lmat[[p]][[q]] <- inla.spde.make.A(mesh = nsim.mesh$mesh, loc = nsim[[p]]$locations[[q]])
    nsim.A[[p]][[q]] <- rbind(nsim.dmat[[p]][[q]], nsim.lmat[[p]][[q]])
    ## set the SPDE
    nsim.spde[[p]][[q]] <- inla.spde2.pcmatern(nsim.mesh$mesh,
                                               prior.range = c(0.05, 0.01),
                                               prior.sigma = c(1, 0.01))
    ## set up the data stack
    nsim.stk[[p]][[q]] <- inla.stack(tag = "nsim",
                                     data = list(y = nsim.y.pp[[p]][[q]], e = nsim.e.pp[[p]][[q]]),
                                     A = list(1, nsim.A[[p]][[q]]),
                                     effects = list(list(b0 = rep(1, nsim.nv + nsim.n[[p]][[q]])),
                                                    list(i = 1:nsim.nv)))
    ## Fitting the model
    ## Computational time
    nsim.time[[p]][[q]] <-
      system.time(nsim.res[[p]][[q]] <- inla(y ~ 0 + b0 + f(i, model = nsim.spde[[p]][[q]]),
                                             family = 'poisson',
                                             data = inla.stack.data(nsim.stk[[p]][[q]]),
                                             control.predictor = list(A = inla.stack.A(nsim.stk[[p]][[q]])),
                                             control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE),
                                             E = inla.stack.data(nsim.stk[[p]][[q]])$e))
    transf[[p]][[q]] <- inla.spde.result(inla = nsim.res[[p]][[q]], name = 'i', 
                                         spde = nsim.spde[[p]][[q]], do.transf = TRUE) ## to user scale
  }
}

### --------------------------------------------------------------------------->> Final output CSV file

list_fixed <- tmp
marg_fixed <- tmp
marg_kappa <- tmp
marg_variance <- tmp
marg_range <- tmp
list_dic <- tmp
list_waic <- tmp
for (f in 1:nrow(pppmat)){
  for (g in 1:unique(pppmat[,1])){
    list_fixed[[f]][[g]] <- nsim.res[[f]][[g]]$summary.fixed
    marg_fixed[[f]][[g]] <- nsim.res[[f]][[g]]$marginals.fixed[[1]]
    marg_kappa[[f]][[g]] <- transf[[f]][[g]]$marginals.kappa[[1]]
    marg_variance[[f]][[g]] <- transf[[f]][[g]]$marginals.variance.nominal[[1]]
    marg_range[[f]][[g]] <- transf[[f]][[g]]$marginals.range.nominal[[1]]
    list_dic[[f]][[g]] <- nsim.res[[f]][[g]]$dic$dic
    list_waic[[f]][[g]] <- nsim.res[[f]][[g]]$waic$waic
  }
}

### --------------------------------------------------------------------------->> CSV summary.fixed
fixed <- data.frame(matrix(unlist(list_fixed), ncol = 7, byrow = TRUE))
colnames(fixed) <- names(nsim.res[[1]][[1]]$summary.fixed)
fixed$True <- rep(pppmat[,2], each=unique(pppmat[,1]))
fixed$Set <- rep(paste0("Set", 1:nrow(pppmat)), each=unique(pppmat[,1]))
fixed$Pattern <- paste0("Pattern", 1:(nrow(pppmat)*unique(pppmat[,1])))
fixed$Mesh <- rep(paste0("Mesh", mn), (nrow(pppmat)*unique(pppmat[,1])))
write.csv(fixed, paste0('demo-sum-fixed-', mn, '.csv'))


### --------------------------------------------------------------------------->> CSV marginals.fixed
mf <- lapply(rapply(marg_fixed, enquote, how="unlist"), eval)
margfixed <- data.frame(do.call('rbind', mf))
margfixed$Set <- rep(paste0("Set", 1:nrow(pppmat)), each=nrow(margfixed)/nrow(pppmat))
margfixed$Pattern <- rep(paste0("Pattern", 1:(nrow(pppmat)*unique(pppmat[,1]))), 
                         each=nrow(margfixed)/(nrow(pppmat)*unique(pppmat[,1])))
margfixed$Mesh <- rep(paste0("Mesh", mn), (nrow(pppmat)*unique(pppmat[,1])))
write.csv(margfixed, paste0('demo-marg-fixed-', mn, '.csv'))


### --------------------------------------------------------------------------->> CSV marginals.log.kappa
mk <- lapply(rapply(marg_kappa, enquote, how="unlist"), eval)
margkappa <- data.frame(do.call('rbind', mk))
margkappa$Set <- rep(paste0("Set", 1:nrow(pppmat)), each=nrow(margkappa)/nrow(pppmat))
margkappa$Pattern <- rep(paste0("Pattern", 1:(nrow(pppmat)*unique(pppmat[,1]))), 
                         each=nrow(margkappa)/(nrow(pppmat)*unique(pppmat[,1])))
margkappa$Mesh <- rep(paste0("Mesh", mn), (nrow(pppmat)*unique(pppmat[,1])))
write.csv(margkappa, paste0('demo-marg-kappa-', mn, '.csv'))


### --------------------------------------------------------------------------->> CSV marginals.log.variance.nominal
mv <- lapply(rapply(marg_variance, enquote, how="unlist"), eval)
margvar <- data.frame(do.call('rbind', mv))
margvar$Set <- rep(paste0("Set", 1:nrow(pppmat)), each=nrow(margvar)/nrow(pppmat))
margvar$Pattern <- rep(paste0("Pattern", 1:(nrow(pppmat)*unique(pppmat[,1]))), 
                       each=nrow(margvar)/(nrow(pppmat)*unique(pppmat[,1])))
margvar$Mesh <- rep(paste0("Mesh", mn), (nrow(pppmat)*unique(pppmat[,1])))
write.csv(margvar, paste0('demo-marg-var-', mn, '.csv'))


### --------------------------------------------------------------------------->> CSV marginals.log.range.nominal
mr <- lapply(rapply(marg_range, enquote, how="unlist"), eval)
margrange <- data.frame(do.call('rbind', mr))
margrange$Set <- rep(paste0("Set", 1:nrow(pppmat)), each=nrow(margrange)/nrow(pppmat))
margrange$Pattern <- rep(paste0("Pattern", 1:(nrow(pppmat)*unique(pppmat[,1]))), 
                         each=nrow(margrange)/(nrow(pppmat)*unique(pppmat[,1])))
margrange$Mesh <- rep(paste0("Mesh", mn), (nrow(pppmat)*unique(pppmat[,1])))
write.csv(margrange, paste0('demo-marg-range-', mn, '.csv'))


### --------------------------------------------------------------------------->> CSV DIC
dic <- data.frame(matrix(unlist(list_dic), ncol = 1, byrow = TRUE))
colnames(dic) <- "DIC"
dic$Set <- rep(paste0("Set", 1:nrow(pppmat)), each=unique(pppmat[,1]))
dic$Pattern <- paste0("Pattern", 1:(nrow(pppmat)*unique(pppmat[,1])))
dic$Mesh <- rep(paste0("Mesh", mn), (nrow(pppmat)*unique(pppmat[,1])))
write.csv(dic, paste0('demo-dic-', mn, '.csv'))


### --------------------------------------------------------------------------->> CSV WAIC
waic <- data.frame(matrix(unlist(list_waic), ncol = 1, byrow = TRUE))
colnames(waic) <- "WAIC"
waic$Set <- rep(paste0("Set", 1:nrow(pppmat)), each=unique(pppmat[,1]))
waic$Pattern <- paste0("Pattern", 1:(nrow(pppmat)*unique(pppmat[,1])))
waic$Mesh <- rep(paste0("Mesh", mn), (nrow(pppmat)*unique(pppmat[,1])))
write.csv(waic, paste0('demo-waic-', mn, '.csv'))


### --------------------------------------------------------------------------->> CSV computational time
time.df <- data.frame(matrix(unlist(nsim.time), ncol = 5, byrow = TRUE))
colnames(time.df) <- names(nsim.time[[1]][[1]])
time.df$Set <- rep(paste0("Set", 1:nrow(pppmat)), each=unique(pppmat[,1]))
time.df$Pattern <- paste0("Pattern", 1:(nrow(pppmat)*unique(pppmat[,1])))
time.df$Mesh <- rep(paste0("Mesh", mn), (nrow(pppmat)*unique(pppmat[,1])))
write.csv(time.df, paste0('demo-time-', mn, '.csv'))
