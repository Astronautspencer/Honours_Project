
data.sim <- function(nsim, beta0, sigma2x, kappa, nu, nw, seed = 1234, plot = TRUE) {
  argg <- c(as.list(environment()), list())
  set.seed(seed) 
  sim <- rLGCP('matern', nsim = nsim, beta0, var = sigma2x, 
               scale = 1 / kappa, nu = nu, 
               win = owin(c(0,nw), c(0,nw)))
  if (nsim == 1) {
    sim <- list(sim)
    }
  locs <- vector('list', length = nsim)
  for (i in 1:nsim) { 
    locs[[i]] <- cbind(sim[[i]]$x, sim[[i]]$y)
    if(plot) {
      plot(attr(sim[[i]], 'Lambda'), 
           main = paste("The simulated spatial field \n and point pattern", i))
      points(cbind(sim[[i]]$x, sim[[i]]$y)[ ,2:1], pch = 16)
    }
  }
  return(list(locations = locs, summary = sim, arguments = argg))
}
