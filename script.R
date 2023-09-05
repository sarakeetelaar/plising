
estimate_all_models <- function(data, p, N, graph_type, n_cores, rep) {
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  param_grid <- expand.grid(
    rep = 1:rep,
    i = 1:length(N),
    j = 1:length(p),
    k = 1:length(graph_type)
  )
  
  complete_parameters = generate_all_parameters_complete(p, data)
  random_parameters = generate_all_parameters_random(p, data)
  smallworld_parameters = generate_all_parameters_smallworld(p, data)
  est <- foreach(idx = seq_len(nrow(param_grid)), .errorhandling = "pass",
                 .packages=c("dplyr", "plising", "psychonetrics")) %dopar% {
    row_params <- param_grid[idx,]
    r <- row_params[["rep"]]
    i <- row_params[["i"]]
    j <- row_params[["j"]]
    k <- row_params[["k"]]

    if (graph_type[k] == "random") {
      params = random_parameters[[j]]
    }
    else if (graph_type[k] == "smallworld") {
      params = smallworld_parameters[[j]]
    }
    else {
      params = complete_parameters[[j]]
    }
    single_simulation(params = params, N = N[i], seed = idx)
  }

  stopCluster(cl)
  
  return(est)
}
library(doParallel)
library(foreach)
library(dplyr)
library(psychonetrics)
load("data/Wenchuan.rda")
wenchuan <- data_preparation(Wenchuan)
data <- wenchuan

# complete graph
p = 5:15
N = c(50, 100, 500, 1000)
n_cores = 10
rep = 100

full_results_complete <- estimate_all_models(data, p, N, "complete", n_cores, rep)
full_results_random <- estimate_all_models(data, p, N, "random", n_cores, rep)
full_results_smallworld <- estimate_all_models(data, p, N, "smallworld", n_cores, rep)

saveRDS(full_results_complete, "full_results_complete.RDS")
saveRDS(full_results_random, "full_results_random.RDS")
saveRDS(full_results_smallworld, "full_results_smallworld.RDS")


