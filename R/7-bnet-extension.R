#' Compute univariate conditional distribution
#' 
#' Computes the f(x|markovBlanket(x)) function
#' 
#' @param bn A fitted bnet object
#' @param node The name of the node
#' @param evidence A named list with values of the MB
#' @param normalize Flag. If true the output function is normalized
#' 
#' @return List with function, min and max values
#' 
full_conditioned.bnet <- function(bn,
                                  node,
                                  evidence,
                                  normalize = F){
  
  # Check Bn is...
  if (!is.bnet(bn))
    stop("First argument is not an bnet object")
  
  if (!is.fitted.bnet(bn))
    stop("Bnet is not fitted")
  
  # Check node in variables
  if (!(node %in% variables(bn)))
    stop("Node not in model variables")
  
  # Get node descendants and parents (NAMES)
  node.ch <- child.bnet(bn, node)
  node.pa <- bn$variables[[node]]$parents
  node.mb <- setdiff(markovblanket.bnet(bn, node),node)
  
  # Evidence is a list with the values of the mb. It can include additional values
  evidence <- as.list(evidence)
  if (!all(node.mb %in% names(evidence)))
    stop("Missing evidence")
  
  # Select only mb variables
  evidence <- evidence[node.mb]
  
  # Now create density function
  unormalized.f <- function(x){
    # Vapply for vector argumetns
    return(
      vapply(x,{
        # Create column with the variable
        evidence[[node]] <- x
        return(
          prod( vapply( bn$variables[c(node,node.ch)], FUN = function(z){
            # Evaluate each node (node+children)
            Rbmop::evaluate.bmop(MIN = 1E-10, x = as.numeric(evidence[c(z$name,z$parents)]), object = z$prob)
          }, FUN.VALUE = numeric(1) ))
        )
      }, FUN.VALUE = numeric(1))
    )
  }
  
  # Get min and max values
  node.min <- min(bn$variables[[node]]$prob$knots[[1]])
  node.max <- max(bn$variables[[node]]$prob$knots[[1]])

  # Compute normalization term
  if (normalize) {
    normalization.cst <- integrate( unormalized.f,
                                    lower = node.min,
                                    upper = node.max)$value
    f <- function(x) unormalized.f(x)/normalization.cst
  }
  else{
    f <- unormalized.f
  }
    
  return(list(
    node = node,
    evidence = evidence,
    f  = f,
    min = node.min,
    max = node.max
  ))
}


#' Generate samples from nodes (with evidence)
#'
#' @param bn A fitted bnet object
#' @param evidence A named list with fixed values
#' @param N The number of samples to generate
#' @param burn.in Number of sampler to ensure chain convergence
#' @param thinnin Number of draws to be discarded between samples
#'
#' @return A named list with the complete sample
#' @export
mcmc_sampler.bnet <- function(bn,
                              evidence,
                              vars = names(bn$variables),
                              N = 1,
                              burn.in = 100,
                              thinnin = 10) {
  # Check Bn is...
  if (!is.bnet(bn))
    stop("First argument is not an bnet object")
  
  if (!is.fitted.bnet(bn))
    stop("Bnet is not fitted")
  
  # Check evidence
  evidence <- as.list(evidence)
  if (is.null(names(evidence)))
    warning("Evidence without names - will be ignored")

  # Get model variables, to generate and in evidence
  variables.in.evidence <- vars[vars %in% names(evidence)]
  variables.not.evidence <- vars[ !(vars %in% names(evidence)) ]
  evidence <- evidence[variables.in.evidence]
  
  # Collapse net
  bn.collapsed <- collapse.bnet(bn, evidence)  
  
  # Preprocess to speed up the process. Looks for isolated undetermined values that dont need MCMC
  isolated.vars <- mcmc_sampler.bnet.isolated(bn.collapsed, evidence, variables.not.evidence)
  
  # Remove isolated from variables.not.ev
  variables.not.evidence <- setdiff(variables.not.evidence, isolated.vars)
  
  # Prealloc DF
  df <- as.data.frame( matrix(NA, ncol = length(vars), nrow = N) )
  colnames(df) <- c(variables.in.evidence, isolated.vars, variables.not.evidence)
  
  # Initialize unassigned variables
  init.values <- vapply(variables.not.evidence, function(node){
    node.min <- min(bn.collapsed$variables[[node]]$prob$knots[[1]])
    node.max <- max(bn.collapsed$variables[[node]]$prob$knots[[1]])
    return(runif(1, min = node.min, max = node.max))
  }, FUN.VALUE = numeric(1))
  
  # Compute isolated functions
  isolated.funcs <- lapply(isolated.vars, function(node){
    return(full_conditioned.bnet(bn.collapsed,node,evidence,normalize = F))
  })
  names(isolated.funcs) <- isolated.vars
  
  # Sample isolated values
  isolated.values <- vapply(isolated.funcs, function(f){
                            return(metropolisSampler_auxiliar( f$f, min = f$min, max = f$max ))
                            }, FUN.VALUE = numeric(1))
  
  # Current state
  current.state <- as.list(c(evidence, isolated.values, init.values))
  names(current.state) <- c(variables.in.evidence,isolated.vars, variables.not.evidence)
  
  # Burnin
  for (i in 1:burn.in) {
    # Select a node to be upodated at random
    node.to.update <- sample( variables.not.evidence, 1)
    # Compute full conditioned
    f <- full_conditioned.bnet(bn.collapsed, 
                               node.to.update, 
                               current.state,normalize = F)
    # Generate sample
    current.state[[node.to.update]] <- metropolisSampler_auxiliar( f$f, min = f$min, max = f$max )
  }
  
  # Now we start
  if (N > 1) {
    for (i in 0:((N - 1)*thinnin - 1) ) {
      # Add sample
      if ( i %% thinnin == 0) {
        df[ floor(i/thinnin) + 1  ,] <- current.state
        
        # Isolated - update current state
        current.state[isolated.vars] <- vapply(isolated.funcs, function(f){
          return(metropolisSampler_auxiliar( f$f, min = f$min, max = f$max ))
        }, FUN.VALUE = numeric(1))
      }
      
      # Select a node to be upodated at random
      node.to.update <- sample( variables.not.evidence, 1)
      # Compute full conditioned
      f <- full_conditioned.bnet(bn.collapsed,
                                 node.to.update,current.state,normalize = F)
      # Generate sample
      current.state[[node.to.update]] <- metropolisSampler_auxiliar( f$f, min = f$min, max = f$max )
    }
  }
  df[N,] <- current.state
  
  return(df)
}

##
# Isolated undetermined values
##
mcmc_sampler.bnet.isolated <- function(bn,
                                       evidence,
                                       vars) {
  # look for isolated
  return(vars[vapply(vars, function(node){
    node.mb <- setdiff(markovblanket.bnet(bn, node),node)
    return(!any( vars %in% node.mb) )
  }, FUN.VALUE = logical(1))])
}


##
# Metropolis - Hastings auxiliar sample for f
##
metropolisSampler_auxiliar <- function(f, burn.in = 25, min = -Inf, max = Inf, sd = 1){
  
  # Initial value
  x.old <- runif(1,min,max)
  f.old <- f(x.old)
  
  # Loop
  for (i in 1:burn.in) {
    # Compute new values
    x.new <- x.old + rnorm(1, sd = sd)
    f.new <- f(x.new)
    
    # Acceptance ratio
    acc.ratio <- (f.new / f.old) * (pnorm(x.old, mean = x.new)/pnorm(x.new, mean = x.old))
    if ((acc.ratio >= 1) || (runif(1) < acc.ratio) ) {
      x.old <- x.new
      f.old <- f.new
    }
  } # END FOR LOOP
  return(x.old)
}