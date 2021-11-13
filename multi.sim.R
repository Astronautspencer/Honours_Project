
multi.sim <- function (allval, seed = 123, simplify = FALSE) {
  sim_list <- vector(mode = "list", length = nrow(allval))
  for (row in 1:nrow(allval)) {
    mypar <- allval[row, ]
    sim_list[[row]] <- data.sim(mypar[1], mypar[2], mypar[3], mypar[4], mypar[5], mypar[6])
  }
  return(sim_list)
}