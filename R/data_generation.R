data_generation = function(N=1000, graph, mu, sigma, no.reps=1e3) {
  data = list()
  for (rep in no.reps) {
    data[[rep]] = IsingSampler::IsingSampler(N, graph, mu, sigma, method="direct")
  }
  return(data)
}
