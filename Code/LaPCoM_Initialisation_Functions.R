########################################################################################################################
# INITIALISATION FUNCTIONS
########################################################################################################################

# ----------------------------------------------------------------------------------------------------------------------
# FUNCTIONS TO CALCULATE THE DISTANCE BETWEEN TWO NETWORKS (FROM NATURE PAPER - REF?)
# ----------------------------------------------------------------------------------------------------------------------
entropy = function(a) {
  a = a[which(a > 0)]
  return(-sum(a * log(a)))
} # end entropy function

node_distance = function(gra, net_type) {
  nd_n = length(igraph::V(gra))
  
  if (nd_n == 1) {
    ret = 1
  } # end nd_n if statement
  if (nd_n > 1) {
    a = Matrix::Matrix(0, nrow = nd_n, ncol = nd_n, sparse = T)
    if (net_type == "count") {
      nd_m = igraph::distances(gra, algorithm = c("dijkstra"), weights = NA)
    } else if (net_type == "binary") {
      nd_m = igraph::distances(gra)
    } else {
      cat("Uh oh. There was a problem in the node_distance function.")
    } # end net_Type check
    
    nd_m[which(nd_m == "Inf")] = nd_n
    quem = setdiff(intersect(nd_m, nd_m), 0)
    
    for (j in (1:length(quem))) {
      l = which(nd_m == quem[j]) / nd_n
      lines = floor(l) + 1
      posm1 = which(l == floor(l))
      
      if (length(posm1) > 0) {
        lines[posm1] = lines[posm1] - 1
      } # end posm1 if statement
      
      a[1:nd_n, quem[j]] = hist(lines, plot = F, breaks = (0:nd_n))$counts
    } # end j for loop
    ret = (a / (nd_n - 1))
  } # end nd_n if statement
  return(ret)
} # end node_distance function

nnd = function(gra, net_type) {
  nnd_N = length(igraph::V(gra))
  nd = node_distance(gra, net_type)
  pdfm = colMeans(as.matrix(nd))
  nnd_norm = log(max(c(2, length(which(pdfm[1:(nnd_N - 1)] > 0)) + 1)))
  return(c(pdfm, max(c(0, entropy(pdfm) - entropy(as.matrix(nd)) / nnd_N)) / nnd_norm))
} # end nnd function

alpha_fun = function(gra, el) {
  alph_N = nrow(el)
  alph_r = sort(igraph::alpha_centrality(gra, exo = igraph::degree(gra) / (alph_N - 1), 
                                         alpha = 1 / alph_N)) / ((alph_N ^ 2))
  return(c(alph_r, max(c(0, 1 - sum(alph_r)))))
} # end alpha_fun function

nature_dist = function(el1, el2, gra1, gra2, w1, w2, w3, net_type) {
  first = 0
  second = 0
  third = 0
  nd_N = length(igraph::V(gra1))
  nd_M = length(igraph::V(gra2))
  PM = matrix(0, ncol = max(c(nd_N, nd_M)))
  
  if (w1 + w2 > 0) {
    pg = nnd(gra1, net_type)
    PM[1:(nd_N - 1)] = pg[1:(nd_N - 1)]
    PM[length(PM)] = pg[nd_N]
    
    ph = nnd(gra2, net_type)
    PM[1:(nd_M - 1)] = PM[1:(nd_M - 1)] + ph[1:(nd_M - 1)]
    PM[length(PM)] = PM[length(PM)] + ph[nd_M]
    
    PM = PM / 2
    
    first = sqrt(max(c((entropy(PM) - (entropy(pg[1:nd_N]) + entropy(ph[1:nd_M])) / 2) / log(2), 0)))
    second = abs(sqrt(pg[nd_N + 1]) - sqrt(ph[nd_M + 1]))
  } # end w1 + w2 if statement
  if (w3 > 0) {
    pg = alpha_fun(gra1, el1)
    ph = alpha_fun(gra2, el2)
    nd_m = max(c(length(pg), length(ph)))
    Pg = matrix(0, ncol = nd_m)
    Ph = matrix(0, ncol = nd_m)
    Pg[(nd_m - length(pg) + 1):nd_m] = pg
    Pg[(nd_m - length(ph) + 1):nd_m] = ph
    
    third = third + sqrt((entropy((Pg + Ph) / 2) - (entropy(pg) + entropy(ph)) / 2) / log(2)) / 2
    
    if (net_type == "count") {
      gra1 = igraph::complementer(gra1)
      el1 = cbind(igraph::as_edgelist(gra1), igraph::E(gra1)$weight)
      gra2 = igraph::complementer(gra2)
      el2 = cbind(igraph::as_edgelist(gra2), igraph::E(gra2)$weight)
    } else if (net_type == "binary") {
      gra1 = igraph::complementer(gra1)
      el1 = igraph::as_edgelist(gra1)
      gra2 = igraph::complementer(gra2)
      el2 = igraph::as_edgelist(gra2)
    } else {
      cat("Uh oh. There was a problem in the nature_dist function.")
    } # end net_type check
    
    pg = alpha_fun(gra1, el1)
    ph = alpha_fun(gra2, el2)
    nd_m = max(c(length(pg), length(ph)))
    Pg = matrix(0, ncol = nd_m)
    Ph = matrix(0, ncol = nd_m)
    Pg[(nd_m - length(pg) + 1):nd_m] = pg
    Ph[(nd_m - length(ph) + 1):nd_m] = ph
    
    third = third + sqrt((entropy((Pg + Ph) / 2) - (entropy(pg) + entropy(ph)) / 2) / log(2)) / 2
  } # end w3 if statement
  return(w1 * first + w2 * second + w3 * third)
} # end nature_dist function