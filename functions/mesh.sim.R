mesh.sim <- function(domain, edge.length, dist.point, extension) {
  argg <- c(as.list(environment()), list())
  timing <- system.time(mesh <- inla.mesh.2d(loc.domain = domain, max.edge = edge.length, 
                       cutoff = dist.point, offset = extension))
  return(list(mesh=mesh, time=timing, arguments = argg))
}
