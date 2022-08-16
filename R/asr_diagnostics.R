#' RERconverge
#'
#'
#' @docType package
#' @author
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib RERconverge
#' @name RERconverge


require(ape)
require(phytools)
require("castor")

#' @keywords internal
matrixToVector <- function(A) {
  v = A[col(A) != row(A)]
  return(v)
}

# if the rate models have the same number of free parameters, returns the one with more zeroes first
# mode can be "vector" or "index"
# vector returns (simpler model, more complex model) as vectors
# index returns the index of the simpler model in c(A, B) as either 1 or 2
#' @keywords internal
getSimplerRateModel <- function(A, B, mode = "vector") {
  if(mode == "vector") {
    a <- matrixToVector(A)
    b <- matrixToVector(B)
    if(max(a) > max(b)) {
      return(list(b, a))
    } else if (max(a) < max(b)) {
      return(list(a,b))
    } else {
      if(sum(a == 0) > sum(b == 0)) {
        return(list(a,b))
      } else if(sum(a == 0) < sum(b == 0)){
        return(list(b,a))
      }
      else {
        # warning("The two provided rate models are of equal simplicity")
        return(list(a,b))
      }
    }
  } else if(mode == "index") {
    # return 1 if A is simpler, return 2 if B is simpler
    a <- matrixToVector(A)
    b <- matrixToVector(B)
    if(max(a) > max(b)) {
      return(2)
    } else if (max(a) < max(b)) {
      return(1)
    } else {
      if(sum(a == 0) > sum(b == 0)) {
        return(1)
      } else if(sum(a == 0) < sum(b == 0)){
        return(2)
      }
      else {
        warning("The two provided rate models are of equal simplicity")
        return(1)
      }
    }
  }
  else {
    stop("invalid mode")
  }
}

# returns true if the simpler model is a special case of the more complex model
# returns false otherwise
#' @keywords internal
areNested <- function(A,B) {
  # get the simpler rate model & convert both to vectors with diagnonals removed
  V = getSimplerRateModel(A,B)
  # s - simpler model, c - more complex model
  s = V[[1]]
  c = V[[2]]

  # if s has only one rate
  if(length(unique(s)) <= 1) {
    if(sum(c == 0) == 0) {
      # if the complex model has no zeroes, return TRUE
      return(TRUE)
    }
    else {
      # otherwise, return FALSE (the simpler model cannot be all 0)
      return(FALSE)
    }
  }
  # if number of free params are equal, but one has more 0s than the other, they are not nested
  else if(max(s) == max(c) && sum(s == 0) > sum(c == 0)) {
    return(FALSE)
  }

  matchUp = vector(mode = "list", length = length(unique(s)))
  names(matchUp) = as.character(unique(s))

  for(rate in unique(s)) {
    matchUp[[as.character(rate)]] = c[s == rate]
  }

  # check for intersections between rows
  for(i in 1:(length(matchUp)-1)) {
    for(j in (i+1):length(matchUp)) {
      if(length(intersect(matchUp[[i]], matchUp[[j]])) > 0) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

#' Gets the matrix form of a rate model from its abbreviation
#' @param abbr The abbreviation for the rate model. Can be one of "ER", "ARD", or "SYM"
#' @param Nstates The number of phenotype states to include in the rate model
#' @return a rate model as an Nstates x Nstates matrix
#' @export
# abbr can be "ARD", "ER", or "SYM"
getMatrixFromAbbr <- function(abbr, Nstates) {
  if(abbr == "ER") {
    # generate an ER matrix
    rate_model = matrix(nrow = Nstates, ncol = Nstates)
    for(col in 1:Nstates) {
      for(row in 1:Nstates) {
        if(row == col) {
          rate_model[row,col] = 0
        } else{
          rate_model[row,col] = 1
        }
      }
    }

  } else if(abbr == "ARD"){
    # make a matrix describing the ARD model
    rate_model = matrix(nrow = Nstates, ncol = Nstates)
    cnt = 1
    for(row in 1:Nstates) {
      for(col in 1:Nstates) {
        if(row == col) {
          rate_model[col,row] = 0
        } else {
          rate_model[col, row] = cnt
          cnt = cnt + 1
        }
      }
    }
  } else if(abbr == "SYM") {
    rate_model = matrix(nrow = Nstates, ncol = Nstates)
    cnt = 1
    for(col in 1:Nstates) {
      for(row in col:Nstates) {
        if(row == col) {
          rate_model[col,row] = 0
        } else {
          rate_model[col, row] = cnt
          rate_model[row, col] = cnt
          cnt = cnt + 1
        }
      }
    }
  }
  return(rate_model)
}

# returns the degrees of freedom for a chi-squared distribution of likelihood ratios between two rate models
# assumes the rate models are nested
#' @keywords internal
getDegreesFreedom <- function(A,B) {
  return(abs(max(A) - max(B)))
}

# pass in the loglikelihoods of the limited and general cases
# limited and general should already be log likelihoods
#' compute the likelihood ratio from the log likeilhoods of two transition matrices
#' @param limited The log likelihood of the simpler model
#' @param general The log likelihood of the more complex model
#' @return the likelihood ratio between two rate models
#' @export
getLikelihoodRatio <- function(limited, general) {
  return(-2 * (limited - general))
}

# gets the null distribution for the LR between two rate models
# as described in Pagel 1994
# Q should be the MLE set of rate parameters fit on the data under the null rate model
#' Returns the p value of a likelihood ratio based on the simulated null distribution
#' @param LR the likelihood ratio between two rate models
#' @param tree the phylogenetic tree to obtain the null distribution on
#' @param Q the transition matrix fit on the null rate model
#' @param null_rm the null rate model (the simpler rate model)
#' @param alt_rm the alternative rate model (the more complex rate model)
#' @param nsims the number of simulations to use generate the null distribution
#' @param plot a boolean expressing whether or not to plot a histogram of likeilhood ratios representing the null distribution
#' @return the p value of a given likeihood ratio based on the null distribution that is generated
#' @export
getLRnullDistribution <- function(LR, tree, Q, null_rm, alt_rm, nsims = 100,
                                  plot = FALSE, ...) {
  ratios = c()
  for(i in 1:nsims) {
    sim = simulate_mk_model(tree, Q)
    null_fit = fit_mk(tree, Nstates = nrow(Q), tip_states = sim$tip_states, rate_model = null_rm,
                      root_prior = get_stationary_distribution(Q),...)
    alt_fit = fit_mk(tree, Nstates = nrow(Q), tip_states = sim$tip_states, rate_model = alt_rm,
                     root_prior = get_stationary_distribution(Q),...)
    lr = getLikelihoodRatio(null_fit$loglikelihood, alt_fit$loglikelihood)
    ratios = c(ratios, lr)
  }
  # probability that a larger loglikelihood will occur by chance
  p = sum(ratios > LR)/nsims
  if(plot) {
    hist(ratios, xlab = "Loglikelihood ratios", main = paste0("Loglikelihood ratio distribution\nLR = ", round(LR)))
    abline(v = LR, col = "red")
  }
  return(p)
}


#' Performs pairwise likelihood ratio comparisons on a list of rate models
#' @param rate_models a named list of rate models in matrix form
#' @param treesObj the trees object returned by readTrees
#' @param phenvals the named phenotype vector with names matching the tip labels on the tree
#' @param nsims the number of simulations to use in getLRnullDistribuutions when the rate models are not nested
#' @param nested_only a boolean indicating whether comparisons should only be made between nested models. If TRUE, this decreases runtime because simulations to get the null distribution are not performed.
#' @param return_type can be "pvals" or "ratios" indicating whether to return p values obtained from the likelihood ratios or the likelihood ratios themselves
#' @param ... additional parameters for the castor function, fit_mk used to fit the transition matrix
#' @return a square matrix of p-values or likelihood ratios for pairwise comparisons between rate models
#' @export
compareRateModels <- function(rate_models, treesObj, phenvals, nsims = 100,
                              nested_only = FALSE, return_type = "pvals",...) {
  # prune tree, order phenvals, and map states to ints
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  phenvals = phenvals[tree$tip.label]
  intlabels = map_to_state_space(phenvals)

  N = length(rate_models)
  if(N <= 1) {
    stop("Not enough rate models to compare")
  } else if(is.null(names(rate_models))) {
    stop("rate_models must have a names attribute")
  }
  # make a matrix to store results
  res = matrix(nrow = N, ncol = N, dimnames = list(names(rate_models), names(rate_models)))

  # make pairwise comparisons
  for(i in 1:(N-1)) {
    for(j in (i+1):N) {
      A = rate_models[[i]]
      B = rate_models[[j]]
      # fit under rate model A and get loglikelihood
      fitA = fit_mk(tree, Nstates = intlabels$Nstates,
                    tip_states = intlabels$mapped_states,
                    rate_model = A,...)
      # fit under rate model B and get loglikelihood
      fitB = fit_mk(tree, Nstates = intlabels$Nstates,
                    tip_states = intlabels$mapped_states,
                    rate_model = B,...)
      # get simpler model to put in the numerator
      n = getSimplerRateModel(A, B, mode = "index")
      # calculate likelihood ratio
      if(n == 1) {
        LR = getLikelihoodRatio(fitA$loglikelihood, fitB$loglikelihood)
      } else { # n == 2
        LR = getLikelihoodRatio(fitB$loglikelihood, fitA$loglikelihood)
      }
      if(areNested(A,B)) {
        if(return_type == "pvals") {
          # get degrees of freedom
          df = getDegreesFreedom(A,B)
          # calculate p-value from chi squared distribution
          p = pchisq(LR, df = df, lower.tail = FALSE)
          # enter in res matrix - more complex matrix corresponds to the column
          if(n == 1) { # A is simpler
            res[i,j] = round(p,5)
          } else { #n == 2, B is simpler
            res[j, i] = round(p,5)
          }
        } else if(return_type == "ratios") {
          if(n == 1) { # A is simpler
            res[i,j] = round(LR,5)
          } else { #n == 2, B is simpler
            res[j, i] = round(LR,5)
          }
        }
      } else { # models are not nested
        if(!nested_only) {
          if(return_type == "pvals") {
            # get p-value from Monte Carlo Methods
            if(n == 1) { # A is simpler
              p = getLRnullDistribution(LR, tree, Q = fitA$transition_matrix,
                                        null_rm = A, alt_rm = B, nsims = nsims)
              res[i,j] = round(p,5)
            } else { # n == 2, B is simpler
              p = getLRnullDistribution(LR, tree, Q = fitB$transition_matrix,
                                        null_rm = B, alt_rm = A, nsims = nsims)
              res[j,i] = round(p,5)
            }
          } else if(return_type == "ratios") {
            if(n == 1) { # A is simpler
              res[i,j] = round(LR,5)
            } else { # n == 2, B is simpler
              res[j,i] = round(LR,5)
            }
          }
        }
      }
    }
  }
  return(res)
}

# get relative entropies at each node
# entropies are in node order (order of anc liks returned by asr_mk_model)
# P is the "target" and M is the "estimate"
#' @keywords internal
getRelativeEntropies <- function(P, M) {
  entropies = c()
  for(i in 1:nrow(P)) {
    p = P[i,]
    m = M[i,]
    p[p == 0] = 1e-323
    m[m == 0] = 1e-323
    e = sum(p * log2(p/m))
    entropies = c(entropies, e)
  }
  return(entropies)
}

#' returns the states at each node corresponding to the state with the max likelihood
#' @param ancliks a table of ancestral likelihoods with rows in node order and columns corresponding to phenotype states
#' @param confidence_threshold the default is NULL, but if provided it will only obtain states from nodes whose max likelihood is greater than or equal to the confidence threshold
#' @return a vector of states for the internal nodes in order of the internal nodes
#' @export
getStatesAtNodes <- function(ancliks, confidence_threshold = NULL) {
  if(is.null(confidence_threshold)) {
    states = rep(0, length(ancliks[,1]))
    for (i in 1:length(states)) {
      states[i] = which.max(ancliks[i, ])
    }
    return(states)
  } else {
    states = rep(0, length(ancliks[,1]))
    for (i in 1:length(states)) {
      if(max(ancliks[i,]) >= confidence_threshold) {
        states[i] = which.max(ancliks[i, ])
      } else {
        states[i] = NA
      }
    }
    return(states)
  }
}


#' Visually compare the results of two rate models by plotting the differences between them on the phylogenetic tree
#' Relative entropy is calculated with A as the "target" and B as the "estimate"
#' Switching the order of rate model A and rate model B DOES give different relative entropies (relative entropy is not symmetric)
#' @param A the first rate model
#' @param B teh second rate model
#' @param treesObj the trees object returned by readTrees
#' @param phenvals the named phenotype vector with names matching the tip labels in the trees object
#' @param mode indicates how to represent differences between ancestral reconstructions. When mode is "entropy" the relative entropies between ancestral likelihoods at each node are plotted along the branches of the trees. When mode is "match", nodes are colored red when state assignments are different and red when they are the same.
#' @param cex a graphical argument representing the size of node labels and tip labels
#' @return a list of relative entropies or state assignemtns and ancestral likelihoods under each rate model
#' @export
visCompareTwoRateModels <- function(A, B, treesObj, phenvals, mode = "entropy",
                                    cex = 0.5,...) {
  # prune tree, order phenvals, and map states to ints
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  phenvals = phenvals[tree$tip.label]
  intlabels = map_to_state_space(phenvals)

  if(!isSymmetric(A)){
    reroot = FALSE
  } else {
    reroot = TRUE
  }

  arA = asr_mk_model(tree, tip_states = intlabels$mapped_states,
                     Nstates = intlabels$Nstates,
                     rate_model = A, reroot = reroot,...)$ancestral_likelihoods

  if(!isSymmetric(B)){
    reroot = FALSE
  } else {
    reroot = TRUE
  }

  arB = asr_mk_model(tree, tip_states = intlabels$mapped_states,
                     Nstates = intlabels$Nstates,
                     rate_model = B, reroot = reroot,...)$ancestral_likelihoods

  if(mode == "entropy") {
    relative_entropies = getRelativeEntropies(arA,arB)
    X = c(rep(0, length(tree$tip.label)), relative_entropies)[tree$edge[,2]]
    plotBranchbyTrait(tree, X, mode = "edge", edge.width = 2, cex = cex)
    names(relative_entropies) = 1:tree$Nnode + length(tree$tip.label)
    return(list(relative_entropies = sapply(relative_entropies,function(x){round(x,5)}), anc_liks1 = arA, anc_liks2 = arB))
  } else if (mode == "match") {
    statesA = getStatesAtNodes(arA)
    statesB = getStatesAtNodes(arB)
    res = cbind(statesA,statesB)
    rownames(res) = 1:tree$Nnode + length(tree$tip.label)
    mismatched_nodes = which(statesA != statesB) + length(tree$tip.label)
    correct_nodes = which(statesA == statesB) + length(tree$tip.label)
    plot(tree, cex = cex)
    nodelabels(node = mismatched_nodes, col = "black", frame = "circle", cex = cex,
               bg = "red")
    nodelabels(node = correct_nodes, col = "black", frame = "circle", cex = cex,
               bg = "lightblue")
    return(list(states = res, anc_liks1 = arA, anc_liks2 = arB))
  }
}

# Helper function for getBoxPlot
# number of incorrect inferrence / total number internal nodes
# sim - the simulation returned by simulate_mk_model
# recon - the reconstruction returned by asr_mk_model
#' @keywords internal
getPercentMatch <- function(sim, recon, confidence_threshold = NULL) {
  if(is.null(confidence_threshold)) {
    Nnodes = length(sim$node_states)
    states = getStatesAtNodes(recon$ancestral_likelihoods)
    return(sum(states == sim$node_states) / Nnodes)
  } else {
    states = getStatesAtNodes(recon$ancestral_likelihoods,
                       confidence_threshold = confidence_threshold)
    Nnodes = sum(!is.na(states))
    return(sum(states == sim$node_states, na.rm = TRUE) / Nnodes)
  }
}

# helper function that returns a boxplot for boxPlotTest function
#' @keywords internal
getBoxPlot <- function(tree, Q, rate_models, nsims,
                       root_prior = get_stationary_distribution(Q),
                       confidence_threshold) {
  N = length(rate_models)

  # make a matrix to store percent matches to generate the boxplot
  res = matrix(nrow = nsims, ncol = N, dimnames = list(NULL, names(rate_models)))

  for(i in 1:nsims) {
    print(i)
    # simulate evolution - ensure all states are represented in tips of the sim
    sim = simulate_mk_model(tree, Q)
    # sometimes, it might never be able to simulate every state at the tips, stop trying after nsims times
    count = 1
    while(count <= nsims && length(unique(sim$tip_states)) < nrow(Q)) {
      sim = simulate_mk_model(tree, Q)
      count = count + 1
    }

    # if nsims simulations do not generate all states at the tips, skip this rate model
    if(length(unique(sim$tip_states)) < nrow(Q)) {
      return(NULL)
    }

    # reconstruct the simulated tree under the different rate models
    for(j in 1:N) {
      rm = rate_models[[j]]

      # determine value of reroot
      if(is.character(rm)) {
        if(rm == "ARD") {reroot = FALSE} else {reroot = TRUE}
      } else {
        if(!isSymmetric(rm)){
          reroot = FALSE
        } else {
          reroot = TRUE
        }
      }

      if(reroot) { # if symmetrical, use asr_mk_model because it is faster
        recon = asr_mk_model(tree, tip_states = sim$tip_states, Nstates = nrow(Q),
                             rate_model = rm, reroot = reroot,
                             root_prior = root_prior)
      }
      else { # otherwise, use getAncLiks becaue it works on asymmetrical models
        fit = fit_mk(tree, tip_states = sim$tip_states, Nstates = nrow(Q),
                     rate_model = rm, root_prior = root_prior)

        tm = fit$transition_matrix

        liks = getAncLiks(tree, tipvals = sim$tip_states, Q = tm,
                          root_prior = root_prior)

        recon = list(ancestral_likelihoods = liks)
      }

      # get match_correct and store in results table
      if(sum(is.nan(liks)) > 0) { # if NaN produced, the rate model was wrong for the tips & can't calculate match_correct
        match_correct = NA
      } else {
        match_correct = getPercentMatch(sim, recon, confidence_threshold = confidence_threshold)
      }
      res[i,j] = match_correct
    }
  }
  # use ggplot2 to make a boxplot
  data = stack(as.data.frame(res))
  colnames(data) = c("percent_correct", "rate_model")
  plot <- ggplot(data, aes(x = rate_model, y = percent_correct, color = rate_model)) + geom_boxplot()
  return(plot)
}

# rate_models should be a named list of rate models
# the transition matrices returned are the Q used to SIMUlATE the data
# they were fit under the rate models to the real observed data
#' Compares prediction accuracy under different rate models and represents the results as a series of boxplots
#' @param treesObj the trees object returned by readTrees
#' @param phenvals the named phenotype vector with names matching the tip labels in the trees object
#' @param rate_models a list of rate models to compare
#' @param nsims the number of simulations to perform for each model of evolution (each rate model)
#' @param confidence_threshold if not NULL, only compares prediction accuracy for nodes whose maximum likelihood is greater than the confidence threshold. Must be between 0 and 1.
#' @return a list of boxplots for each rate model
#' @export
boxPlotTest <- function(treesObj, phenvals, rate_models, nsims = 100,
                        confidence_threshold = NULL, ...) {

  N = length(rate_models)

  # if not named, name the rate models numerically
  if(is.null(names(rate_models))) {
    names(rate_models) = as.character(1:N)
  }

  # make a list to store the Qs fit with each rate model
  TMs = vector(mode = "list", length = N)
  names(TMs) = names(rate_models)
  plots = vector(mode = "list", length = N)

  # prune tree, order phenvals by tiplabels, map to state space
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  phenvals = phenvals[tree$tip.label]
  intlabels = map_to_state_space(phenvals)

  # loop through the rate models as evolution models
  for(i in 1:N) {
    # get the next rate model
    rm = rate_models[[i]]

    # fit a transition matrix to simulate from
    Q = fit_mk(tree, Nstates = intlabels$Nstates,
               tip_states = intlabels$mapped_states, rate_model = rm,...)$transition_matrix
    TMs[[i]] = Q

    plot = getBoxPlot(tree, Q, rate_models = rate_models, nsims = nsims,
                      confidence_threshold = confidence_threshold)
    if(is.null(plot)) {
      plot = ggplot() + theme(text = element_text(size = 10))  + labs(title = paste("Simulations from", names(rate_models)[i], "evolution model failed to generate extant species of every category."))
    } else {
      plot = plot + labs(title = paste("Evolution Model:", names(rate_models)[i]))
    }
    plots[[i]] = plot
  }
  return(list(plots = plots, transition_matrices = TMs))
}

# helper function
#' @keywords internal
standardize <- function(v) {
  i = 1
  x = 1
  w = rep(0, length(v))
  while(x <= max(v) && i <= length(v)) {
    if(v[i] != 0 && !(i > 1 && (v[i] %in% v[1:(i-1)]))) {
      ii = which(v == v[i])
      # print(ii)
      w[ii] = x
      # print(w)
      x = x + 1
      i = i + 1
    }
    else {
      i = i + 1
    }
  }
  return(w)
}

# helper function
#' @keywords internal
vectToMatrix <- function(v,N) {
  m = matrix(rep(0, N*N), N)
  cnt = 1
  for(row in 1:N) {
    for(col in 1:N) {
      if(row != col) {
        m[row,col] = v[cnt]
        cnt = cnt + 1
      }
    }
  }
  return(m)
}

# helper function
#' @keywords internal
getClosestVals <- function(Q) {
  q = as.vector(Q[Q>=0])
  q = unique(q)
  # if there is only one unique value in Q, return this value as the two closest values
  if(length(q) == 1) {
    return(c(q[1],q[1]))
  }
  inds = order(q)
  dd = q[inds]
  dists = c()
  for(i in 1:(length(dd)-1)) {
    d = dd[i+1] - dd[i]
    dists = c(dists,d)
  }
  i = which.min(dists)[1] # use the first (arbitrarily) in case there is more than one
  return(q[inds[c(i, i+1)]])
}

# helper function
# fixes rate parameters to be in consecutive order
#' @keywords internal
cleanUp <- function(M) {
  m = unique(M[row(M) != col(M)]) # get the unique off diagonal elements
  m = m[m != 0] # get rid of zeroes
  m = m[order(m)] # place in increasing order
  M = apply(M, c(1,2),function(x){ifelse(x != 0, which(m == x), 0)})
  return(M)
}

# based off of algorithm from this paper:
# https://watermark.silverchair.com/msr128.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAsAwggK8BgkqhkiG9w0BBwagggKtMIICqQIBADCCAqIGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMq-Fd1lW6vh9qbd81AgEQgIICc84OGoWGpk22V34BKaHWMa5nlKipbxXcdrXOVTIGe6CYhnOBpy6uslY0yIHV0sHIfjlfLmbf8flF_Fmx19FyReu7Z4ERQopL_xgl8IH3L7VzSZezXut-JSI1d3FD0EcTWii4njZVmvOfdgGVP5Hz2t3ZvyhsLlno9SGr-sPJMEnMX2-jdwFWkKYKjAlbdAeFc-GZliSU0zqx108sA2MdXkZaIP-LM-kD3lB68Xxsu8vjhf6WixUdWPWx8MoVlAvJbtn8jxjsWzkXBvhDkzxrGQZzg-vCUYSgs4kZvv4HRHgW3Yny19d9_bN6QCBpOuJ-nbTF2EsKhBxpWWAHJUe4GuBnP_X1FzAVwE8m0vFuu8AtCWKLyad2-nns7IENx_pzAyi9NqqO3SuqFMkZc1XAd4qN5B-194ZhxEnPI3JS-jPrYWXJcTpqjgFHn-hTAJhPVCeuyry0nGIsjyi_URq-AoswNWE5xfcvyIrh6zJH90LDSB5cJrwBzSUUEpauvl2Zhj2FiYe3JjrMiB3rgbYfrvjFn2yvXmynwpX77YcR9Ljaby4rwuyNvYFRiHFgkeuFhqupHKE9BgnmyzP-7Mx-Z08lUJLaLQGC4Nypt1gb6Zd_afY2BfP58OKhjbcNt7DeGgsX9XYPJziiF_km97rwwfMXhRjAw3D9uOsvTFQ5laTXE37Bpo0YEXMWAPrxSdI07tPxbxHRy5UYO9ysgd6uAf_YHCVk0WPjRPDaF8dzhOm6ZoDI6ZopA0Iwyk-MqIrHgwMyXQ3cBa6Twke4M0CCA1Yvqcejpa59D575maM0MhIErOzY9kAtfwdiQHk45rtJBruX5w
#' iteratively searches for simpler rate models that still provide a good fit to the observed data at the tips of the tree
#' @param treesObj the trees object returned by readTrees
#' @param phenvals the named phenotype vector with names matching the tip labels in the trees object
#' @param pthreshold the more complex model is considered significantly better when the p value is less than pthrehsold
#' @param lthreshold the more compelx model is considered significantly better when the likelihood ratio is greater than the lthreshold, used when the more complex and simpler models are not nested
#' @param max_iterations stops the algorithm if the number of iterations exceeds max_iterations even if it isn't done in order to avoid extremely long runtimes
#' @param ... additional parameters ot pass to fit_mk, the castor function used to fit a transition matrix given a rate model, phylogenetic tree, and phenotype data
#' @return returns a list of the rate models generated during the iterative search for simpler rate models and a table of statistics including the likelihood ratios computed on each iteration
#' @export
searchRateModels <- function(treesObj, phenvals, pthreshold = 0.05,
                             lthreshold = 4,
                             max_iterations = 2000, ...) {
  dims = length(unique(phenvals))
  # make a list of models
  M = list(getMatrixFromAbbr("ARD", dims))
  # make a matrix to store p values and LRs
  tab = matrix(nrow = 0, ncol = 6,
               dimnames = list(NULL, c("prevIndex", "newIndex", "p_value", "likelihood_ratio", "prevAIC", "prevLogLik")))
  # make a hash table for quick lookup whether a model has been used yet
  H = new.env(hash = TRUE)
  getKey <- function(m) {
    v = standardize(as.vector(m))
    return(paste(as.character(v),collapse = ""))
  }
  # add ER and ARD to the hash table
  H[[getKey(M[[1]])]] = M[[1]]
  i = 1

  # prune tree, order phenvals by tiplabels, map to state space
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  phenvals = phenvals[tree$tip.label]
  intlabels = map_to_state_space(phenvals)

  # continue until you've visited every rate model in M
  while(i <= length(M)) {
    if(!is.null(max_iterations) && i > max_iterations) {
      print(paste("Stopping after", max_iterations, "iterations."))
      break
    }
    m = M[[i]] # get the next rate model
    print(m)

    # if only one free parameter is left, do not simplify further
    if(max(m) > 1) {
      # step 1: fit the transition matrix on the data
      fit = fit_mk(tree, intlabels$Nstates, intlabels$mapped_states,
                   rate_model = m,...)
      Q = fit$transition_matrix

      # step 2: get the new models
      # currently, if more than on parameter is equally close to another or equally close to zero, change each one at a time & change all of them
      min_pos = which(Q == min(Q[Q>=0]), arr.ind = TRUE)
      if(length(min_pos) > 1) {
        gz = vector(mode = "list", length = nrow(min_pos) + 1)
        for(k in 1:length(gz)) {
          if(k < length(gz)) {
            g = m
            g[min_pos[k,1], min_pos[k,2]] = 0
            # fix parameters to be in consecutive order
            g = cleanUp(g)
            gz[[k]] = g
          } else if(k == length(gz)) {
            g = m
            for(r in 1:nrow(min_pos)) {
              g[min_pos[r,1], min_pos[r,2]] = 0
            }
            # fix parameters to be in consecutive order
            g = cleanUp(g)
            gz[[k]] = g
          }
        }
      } else {
        g = m
        g[min_pos[1,1], min_pos[1,2]] = 0
        # fix parameters to be in consecutive order
        g = cleanUp(g)
        gz = list(g)
      }

      # generate ge
      ee = getClosestVals(Q) #returns the two closest values in Q
      e1 = which(Q == ee[1], arr.ind = TRUE) # get position(s) of first value
      e2 = which(Q == ee[2], arr.ind = TRUE) # get position(s) of second value

      # loop through the value with fewer positions in order to add fewer rate models (improve speed)
      if(nrow(e1) <= nrow(e2)) {
        if(nrow(e1) > 1) {
          rp = m[e2[1,1], e2[1,2]] # get the rate parameter for the value with more positions (uses first row arbitrarily)
          ge = vector(mode = "list", length = nrow(e1) + 1)
          for(k in 1:length(ge)) {
            if(k < length(ge)) {
              g = m
              g[e1[k,1], e1[k,2]] = rp
              # fix parameters to be in consecutive order
              g = cleanUp(g)
              ge[[k]] = g
            } else if(k == length(ge)) {
              g = m
              for(r in 1:nrow(e1)) {
                g[e1[r,1], e1[r,2]] = rp
              }
              # fix parameters to be in consecutive order
              g = cleanUp(g)
              ge[[k]] = g
            }
          }
        } else {
          g = m
          rp = m[e2[1,1], e2[1,2]]
          g[e1[1,1], e1[1,2]] = rp
          # fix parameters to be in consecutive order
          g = cleanUp(g)
          ge = list(g)
        }
      }
      else {
        if(nrow(e2) > 1) {
          rp = m[e1[1,1], e1[1,2]] # get the rate parameter for the value with more positions (uses first row arbitrarily)
          ge = vector(mode = "list", length = nrow(e2) + 1)
          for(k in 1:length(ge)) {
            if(k < length(ge)) {
              g = m
              g[e2[k,1], e2[k,2]] = rp
              # fix parameters to be in consecutive order
              g = cleanUp(g)
              ge[[k]] = g
            } else if(k == length(ge)) {
              g = m
              for(r in 1:nrow(e2)) {
                g[e2[r,1], e2[r,2]] = rp
              }
              # fix parameters to be in consecutive order
              g = cleanUp(g)
              ge[[k]] = g
            }
          }
        } else {
          g = m
          rp = m[e1[1,1], e1[1,2]]
          g[e2[1,1], e2[1,2]] = rp
          # fix parameters to be in consecutive order
          g = cleanUp(g)
          ge = list(g)
        }
      }

      # check whether each model in ge/gz have already been generated (are in the hash table)
      # if not add them, check whether m is NOT significanlty better, if not add them to M
      for(k in 1:length(gz)) {
        kz = getKey(gz[[k]])
        hz = H[[kz]]

        if(is.null(hz)) {
          H[[kz]] = gz[[k]] # add to hash table

          fitz = fit_mk(tree, intlabels$Nstates, intlabels$mapped_states,
                        rate_model = gz[[k]],...)

          LRz = getLikelihoodRatio(fitz$loglikelihood, fit$loglikelihood)

          if(areNested(gz[[k]], m)) {
            df = getDegreesFreedom(gz[[k]],m)
            p = pchisq(LRz, df = df, lower.tail = FALSE)
            # if m is not significantly better than gz
            if(p > pthreshold) {
              gzIndex = length(M) + 1
              M[[gzIndex]] = gz[[k]]
              tab = rbind(tab, c(i, gzIndex, p, LRz, fit$AIC, fit$loglikelihood))
            }
          } else{
            # if m is not significantly better than ge
            if(LRz < lthreshold) {
              gzIndex = length(M) + 1
              M[[gzIndex]] = gz[[k]]
              tab = rbind(tab, c(i, gzIndex, NA, LRz, fit$AIC, fit$loglikelihood))
            }
          }
        }
      }

      for(k in 1:length(ge)) {
        ke = getKey(ge[[k]])
        he = H[[ke]]
        # if they are not already in there, add them
        if(is.null(he)) {
          H[[ke]] = ge[[k]] # add to hash table

          fite = fit_mk(tree, intlabels$Nstates, intlabels$mapped_states,
                        rate_model = ge[[k]],...)
          LRe = getLikelihoodRatio(fite$loglikelihood, fit$loglikelihood)

          if(areNested(ge[[k]], m)) {
            df = getDegreesFreedom(ge[[k]],m)
            p = pchisq(LRe, df = df, lower.tail = FALSE)
            # if m is not significantly better than ge
            if(p > pthreshold) {
              geIndex = length(M) + 1
              M[[geIndex]] = ge[[k]]
              tab = rbind(tab, c(i, geIndex, p, LRe, fit$AIC, fit$loglikelihood))
            }
          } else{
            # if m is not significantly better than ge
            if(LRe < lthreshold) {
              geIndex = length(M) + 1
              M[[geIndex]] = ge[[k]]
              tab = rbind(tab, c(i, geIndex, NA, LRe, fit$AIC, fit$loglikelihood))
            }
          }
        }
      }

      # step 3: increment i
      i = i + 1
    } else {
      # step 3:
      i = i + 1
    }
  }
  return(list(models = M, stats = tab))
}

#' Filters a list of rate models based on a boolean expression expressing filtering criteria
#' @param models a list of rate models
#' @param criteria a function that is passed to sapply specifying the desired criteria to filter by. It should take one argument, a rate model, and should be a boolean expression expressing criteria for the rate model.
#' @return returns a subset of the list of rate models that match the criteria specified in the criteria function
#' @export
filterByCriteria <- function(models, criteria) {
  ii = sapply(models, criteria)
  return(models[ii])
}

