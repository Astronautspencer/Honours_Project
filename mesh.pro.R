mesh.pro <- function(domain, pmat){
  mesh_list <- vector(mode = "list", length = nrow(pmat))
  for (row in 1:nrow(pmat)) {
    paras <- pmat[row, ]
    mesh_list[[row]] <- mesh.sim(domain, c(paras[1], paras[2]), paras[3], c(paras[4], paras[5]))
  }
  return(mesh_list)
}