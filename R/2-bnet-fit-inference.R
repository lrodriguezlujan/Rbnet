.bnetenv<-new.env()
lockBinding(sym = ".bnetenv",env = environment())
assign(".bnetPars",value = list(searchBmop=FALSE,
                                intMethod="adaptative"), envir=.bnetenv  )
lockEnvironment(env = .bnetenv,bindings = ".bnetPars")

#' Set fitting parameters
#'
#' Set appropriate parameters to be used by fitting
#' functions.
#'
#'  @param ... see details.
#'  @details This function set parameters used for estimation.
#'           \code{N}:  The number of knots in every dimensions.
#'           \code{order}: The order of the B-spline in every dimensions.
#'           \code{alpha}: The exponent to compute the number of knots.
#'          \code{knotsMethod}: \code{"uniform"} or \code{"quantiles"} knots.
#'
#'           @export
bnetPar<-function(...){

  l <- list(...)

  old <- get(".bnetPars",envir = .bnetenv)
  if (length(l) == 0){
    return(old)
  }
  unlockBinding(sym = ".bnetPars",env = .bnetenv)
  for (a in names(l)){
    if (a %in% names(old)){
      if (a=="intmethod"){
        if (!(l[[a]] %in% c("adaptative","sampling"))){
          l[[a]]<-"adaptative"
          warning("intMethod parameter can be only adaptative or sampling,
                  set to adaptative by default")
        }
      }
      old[[a]] <- l[[a]]
    }
    else{
      warning(paste(a,"is not a Rbnet parameter"))
    }
  }
  assign(".bnetPars",old,envir = .bnetenv)
  lockBinding(sym = ".bnetPars",env = .bnetenv)
}




# Evaluate an object over points (INTERNAL FUNCTION)
#
# This is a generic functions
# @param x points over perform evaluations
# @param object R object to be evaluated \code{evaluate.foo} must be available
# for object of class \code{foo}
# @return array of evaluations
evaluate<-function(x,object,...){
  UseMethod("evaluate",object)
}


# (INTERNAL FUNCTION)
evaluate.evidence<-function(x,object,MIN=10^(-10)){
  if (is.null(dim(x))){
    dim(x)<-c(1,length(x))
  }
  if (dim(x)[1]>1){ return(apply(x,MARGIN = 1,FUN =
                                   evaluate.evidence,object=object,MIN=MIN))}
  if (x[1]==object){ return(1)}
  else{return(max(MIN,0))}
}


# 'Upper and Lower values for an object (INTERNAL FUNCTION)
#
#' This is a generic function, used by package bnet to compute lower and upper
#' boundaries for the \code{$prob} values of a node.
#' @param object an R object, e.g. a bmop object
#' @return numeric value
#'
lower <- function(object){ UseMethod("lower",object)}

#' Upper and Lower values for an object (INTERNAL FUNCTION)
#'
#' This is a generic function, used by package bnet to compute lower and upper
#' boundaries for the \code{$prob} values of a node.
#' @param object an R object, e.g. a bmop object
#' @return numeric value
#'
upper <- function(object){ UseMethod("upper",object)}


# internal lower for evidence (INTERNAL FUNCTION)
lower.evidence <- function(object){
  return(object)
}

# internal upper for evidence (INTERNAL FUNCTION)
upper.evidence <- function(object){
  return(object)
}


# fast fit, internal function (INTERNAL FUNCTION)
bmop_fit.bnet <- function(object,data,nodes,unif.variables = NULL, Mins = NULL,
                        Maxs=NULL, mc.cores = getOption("mc.cores",2), ...){
  # Subnode list
  nodes <- nodes[nodes %in% variables.bnet(object)]
  # Parallel mop learning
  mop.vars <- parallel::mclapply(nodes, function(a, object, unif.variables, data, Mins, Maxs, ...){
  #mop.vars <- lapply(nodes, function(a, object, unif.variables, data, Mins, Maxs, ...){
    # Select variables in the mop and data
    var <- c(a,object$variables[[a]]$parents)
    if (!is.element(el = a,set = unif.variables)) {
      return(
        list(
          # MOP
          prob = Rbmop::bmop_fit(data = as.data.frame(data[,var]),conditional = T,
                                 Min = as.numeric(Mins[var]),Max = as.numeric(Maxs[var]),
                                 ...),
          # Marginal
          marginal = Rbmop::bmop_fit(data = data[,a])
        ))# END RETURN
    }
    else {
      prob <- Rbmop::new_bmop(knots = Rbmop::generate_knots(data = NULL, N = 1, Min = Mins[var],Max = Maxs[var]),order = 1,ctrpoints = 1)
      prob <- Rbmop::normalize.bmop(prob)
      return(list(prob = prob))
    }
  }, object, unif.variables, data, Mins, Maxs, ..., mc.cores = mc.cores) # END parallel MCLapply
  #}, object, unif.variables, data, Mins, Maxs, ...) # END parallel MCLapply

  # Copy to object
  for (i in 1:length(nodes)) {
    object$variables[[ nodes[i] ]]$prob <- mop.vars[[i]]$prob
    object$variables[[ nodes[i] ]]$marginal <- mop.vars[[i]]$marginal
  }
  if (any( nodes %in% unif.variables)) object$info$fit$unif.variables <- unif.variables
  object$info$fitted <- TRUE
  return(object)
}



#' Fit a bnet object
#'
#' @param object a bnet obejct
#' @param data data.frame with names that include the variables of \code{object}
#' @param search logical, use penalized loglik search or not
#' @param ... additional parameters to be passed to the bmop fitting functions
#' @return a bnet object with fitted densities and conditional densities
#' @export
#' @examples
#' \dontrun{require(bnlearn)
#' data(gaussian.test)
#' bn<-hc(gaussian.test,maxp=2)
#' plot(bn)
#' bnet<-as.bnet(bn)
#' plot(bnet)
#' is.fitted.bnet(bnet)
#' bnet<-fit.bnet(bnet,gaussian.test[1:100,]) # just 100 observation
#'  are used to make computation lighter
#' is.fitted.bnet(bnet)
#' plot(bnet$variables$A$prob)
#' plot(bnet$variables$B$prob)}
#' X<-rnorm(100)
#' Y<-rnorm(100)
#' U<-rnorm(100,mean=X+Y)
#' V<-rnorm(100,mean=X)
#' Z<-rnorm(100,mean=U)
#' data<-data.frame(X,Y,U,V,Z)
#' mat<-matrix(nrow=5, c(0,0,1,1,0,
#'                      0,0,1,0,0,
#'                      0,0,0,0,1,
#'                      0,0,0,0,0))
#' bnet<-new_bnet(names(data),mat)
#' bnet<-fit.bnet(bnet,data)
fit.bnet<-function(object,data,nodes=NULL,Mins=NULL,Maxs=NULL,...){
  if (is.null(nodes)){
    nodes<-variables.bnet(object)
  }

  if (is.null(Mins)){
    Mins<-lapply(data,min)
    names(Mins)<-variables(object)
  }
  if (is.null(Maxs)){
    Maxs<-lapply(data,max)
    names(Maxs)<-variables(object)
  }
  if (dim(data)[1]==0){
 return(object)
  }
    return(bmop_fit.bnet(object =
                                  object,data = data,nodes=nodes,
                         Mins=Mins,Maxs=Maxs,...))

}


#' put evidence on a bnet object
#'
#' @param object a bnet object
#' @param evidence numeric with names or list, the value of evidence
#' @param store.old logical
#' @param propagate logical, see details
#' @return a bnet object with evidence
#' @details If  \code{propagate} is set to \code{TRUE}, the evidence will be
#' naively propagate, that is propagate just to the direct child.
#' @export
put_evidence.bnet<-function(object,evidence,store.old=T,propagate=T){
  evidence<-as.list(evidence)
  for (a in names(evidence)){
    if (is.null(object$variables[[a]]$name)){
      object$variables[[a]]$name<-a

    }
    if (store.old){
        object$variables[[a]]$prob.old<-object$variables[[a]]$prob
      if (propagate){
      for (n in child.bnet(object = object,nodes = a)){
        object$variables[[n]]$prob.old<-object$variables[[n]]$prob
      }
      }
    }
    object$variables[[a]]$prob<-evidence[[a]]
    class(object$variables[[a]]$prob)<-"evidence"
  }

  if (propagate){
    for (a in names(evidence)){
    for (n in child.bnet(object = object,nodes = a)){
      par<-object$variables[[n]]$parents
      if (!class(object$variables[[n]]$prob)=="evidence"){

      object$variables[[n]]$prob<-Rbmop::put_evidence.bmop(
        evd.pos = (2:(length(par)+1))[par==a] ,
        evidence = evidence[[a]]  ,
        object = object$variables[[n]]$prob  )
      }
    }
  }
  }
  object$store.old<-store.old
  return(object)
}


#' remove all evidence from a bnet object
#'
#' @param object a bnet object
#' @return a bnet object, the same as
#' \code{object} if old probabilities were not stored when evidence was set,
#' otherwise the restored bnet object
#' @export
clear_evidence.bnet<-function(object){
  if (is.null(object$store.old)){
    warning("No old probability, sure evidence was set?...")
    return(object)
  }
  if (!object$store.old){
    warning("Old probability was lost putting evidence...
            Next time try saving it")
    return(object)
  }
  else{
    for (a in variables(object)){
      object$variables[[a]]$prob<-object$variables[[a]]$prob.old
      object$variables[[a]]$prob.old<-NULL
    }
  }
  object$store.old<-NULL
}


#' compute probability for a bnet
#'
#' @param x vector of evidence
#' @param object
#' @param log logical, if \code{TRUE} the logarithm of the probabilities
#' is returned
#' @return numeric or \code{NULL} if partial evidence is provided
#' @export
probability.bnet<-function(x,object,log=F){
  if (!object$info$fitted){
    warning("not fitted bnet object, use fit.bnet to fit the
            Bayesian network, NULL is returned")
    return(NULL)}
  if (is.null(dim(x))){ dim(x)<-c(1,length(x))}
  if (is.null(names(x))){ names(x)<-variables.bnet(object)}
  x<-as.data.frame(x)
  vars<-variables.bnet(object)
  if (all(vars %in% names(x))){
    return(apply(x,MARGIN = 1,FUN = function(d){
      if (!log) {
      return(prod(sapply(object$variables,FUN = function(vari){
        return(evaluate(object = vari$prob,x = d[c(vari$name,vari$parents)],
                        MIN = 0))
      })))

      }
      else {
        return(sum(log((sapply(object$variables,FUN = function(vari){
          return(evaluate(object = vari$prob,x = d[c(vari$name,vari$parents)]))
        })))))
      }
    }))
  }
else {
  to.integrate<-var[!(var %in%  names(x))]
  warning("partial evidence is not yet implemented")
  return(NULL)
}
}

#' predict posterior values for a bnet object
#'
#' @param object a bnet object
#' @param newdata new data, evidence
#' @param targets targets variables
#' @param predictors predictors variables
#' @param method string "posteriorMean" or "posteriorMode"
#' @param optmethod \link{optim}
#' @param intmethod string "adaptative" or "sampling"
#' @return the values predicted
#' @export
predict.bnet<-function(object,newdata,targets=NULL,predictors=NULL,
                       method="posteriorMean",optmethod="Nelder-Mead",
                       intMethod=bnetPar()$intMethod,...){
  if (!object$info$fitted){
    warning("not fitted bnet object, use fit.bnet to fit the Bayesian
            network, NULL is returned")
    return(NULL)}
  if (is.null(targets)){
    if (is.null(object$targets)){
      targets<-variables.bnet(object)[1]
    }
    else{
    targets<-object$targets
    }
  }
  if (is.null(predictors)){
    if (is.null(object$predictors)){
      predictors<-variables.bnet(object)
      predictors<-predictors[!(predictors %in% targets) ]
    }
    else{
    predictors<-object$predictors
    }
  }
  nam<-names(newdata)
  for (a in targets){
    if (!(a %in% nam)){
      newdata[,a]<-rep(0,dim(newdata)[1])
    }
  }
  mb<-markovblanket.bnet(object = object,nodes=targets)
  if (!(all(mb %in% nam))){
    warning("please provide a dataset with the needed predictor variables,
            that is the markov blanket of the targets")
    return(NULL)
  }
    if (dim(newdata)[1]==1){
      part<-sapply(X = object$variables[targets],FUN=function(n){
        return(lower(n$prob))
      } )
      fin<-sapply(X = object$variables[targets],FUN=function(n){
        return(upper(n$prob))
      } )
      if (method=="posteriorMode"){
        if (is.null(dim(newdata))){
          dim(newdata)=c(1,length(newdata))
        }
        if (length(targets)==1){ optmethod="Brent" }
        ctr=list()
        ctr$fnscale=-1
        return(optim(lower =part ,upper =fin ,
                     method=optmethod,control=ctr,par=part,fn = function(x){
          newdata[,targets]<-x
          return(
            sum(sapply(object$variables,FUN =
                         function(a){
                           log(evaluate(
                             x = as.numeric(newdata[,c(a$name,a$parents)]) ,
                             object = a$prob))}))
          )
        })$par)
      }
      if (method=="posteriorMean"){
        child<-child.bnet(object,targets)
        #cc<-prod(sapply(object$variables[!(variables.bnet(object)
        #%in% child) ],FUN =function(a){ evaluate(MIN = 10^(-10),x =
          #             as.numeric(newdata[,c(a$name,a$parents)]) ,
        #object = a$prob)} ))
        densityi<-function(x){
          newdata[,targets]<-x
          return(prod(sapply(object$variables[c(targets,child)],
                             FUN = function(a){
                               evaluate(MIN = 10^(-10),
                                        x = as.numeric(
                                          newdata[,c(a$name,a$parents)]) ,
                                        object = a$prob)})))
        }
        if (intMethod=="adaptative"){
          R<-rep(0,times = length(targets))
          for (i in 1:length(targets)){
            konst<-
              cubature::adaptIntegrate(f =densityi,lowerLimit = part,
                                       upperLimit = fin )$integral
            R[i]<-cubature::adaptIntegrate(f =function(x){x[i]*densityi(x)},
                                           lowerLimit = part,upperLimit = fin )$
              integral/konst
          }
        }
        if (intMethod=="sampling"){
          D<-length(targets)
          xstart<-(part+fin)/2
          sample<-Rbmop::sampler_MH(N = 500,d = D,densit = densityi,
                                    h = 3,M=500,xstart=xstart,max=fin,min=part)
          if (D>1){colnames(sample)<-targets}
          R<-colMeans(sample)
        }
        return(R)
      }
    }
    else {
      l <- 1:dim(newdata)[1]
      dim(l) <- dim(newdata)[1]
      return(t(apply(l, MARGIN = 1, FUN = function(i){ return(
        predict.bnet(method = method, object = object, targets = targets,
                     predictors = predictors, newdata = newdata[i,],
                     intmethod = intmethod))
        })))

    }

}

#' Generate samples from a bnet
#'
#' @rdname sample.bnet
#'
#' @param bn A fitted bnet object
#' @param evidence A named list with fixed values (NULL for no evidence)
#' @param N The number of samples to generate
#' @param burn.in Number of sampler to ensure chain convergence
#' (uses \code{mcmc_sampler.bnet})
#' @param thinnin Number of draws to be discarded between samples
#' (uses \code{mcmc_sampler.bnet})
#' @param mc.cores Max number of parallel workers
#'
#' @importFrom igraph graph_from_edgelist components add.vertices
#' @importFrom dplyr bind_cols
#'
#' @return A data frame with one sample per row (Columns named)
#' @export
sample.bnet <- function(bn, evidence = NULL, N = 1,
                        burn.in = 100,
                        thinnin = 10,
                        mc.cores = getOption("mc.cores", default = 4)){

  # Check Bn is...
  if (!is.bnet(bn))
    stop("First argument is not an bnet object")

  if (!is.fitted.bnet(bn))
    stop("Bnet is not fitted")

  # Check evidence
  if (!is.null(evidence)) {
    evidence <- as.list(evidence)
    if (is.null(names(evidence))) {
      warning("Unnamed evidence - will be ignored")
      evidence <- NULL
    }
  }

  # If no evidence is provided use forward sampling
  if (length(evidence) == 0)
    return(forwardSampler.bnet(bn, N = N))

  # Get model variables, to generate and in evidence
  vars <- variables(bn)
  variables.in.evidence <- vars[vars %in% names(evidence)]
  variables.not.evidence <- vars[ !(vars %in% names(evidence)) ]
  evidence <- evidence[variables.in.evidence]

  # Collapse net
  bn <- collapse.bnet(bn, evidence)

  # Look for remaining connected components (and detect isolated nodes)
  l <- as.list(bn)
  isolated.nodes <- l$nodes[!(l$nodes %in% unique(as.character(l$arcs)))]
  g <- igraph::graph_from_edgelist(l$arcs, directed = T)
  # Add isolated nodes
  g <- igraph::add.vertices(g, length(isolated.nodes), name = isolated.nodes)
  # Find components
  g.components <- igraph::components( g , mode = "weak")
  # Generate samples for each component
  components.samples <- parallel::mclapply(1:g.components$no, function(comp, bn, evidence, comps){
    component.sampler(bn,
                      evidence,
                      names(comps$membership)[comps$membership == comp],
                      N = N,
                      burn.in = 100,
                      thinnin = 10)
  },bn, evidence, g.components, mc.cores = mc.cores)

  return(dplyr::bind_cols(components.samples))
}

#' Auxiliar function. per component sampler
#' @noRd
component.sampler <- function(bn, evidence, component, N, burn.in, thinnin){

  selector <- component %in% names(evidence)
  # Case A: All component members are in evidence
  if ( all(selector) ) {
    values <- as.numeric(evidence[component])
    df <- as.data.frame(matrix( rep(values, each = N), ncol = length(component)))
    colnames(df) <- component
    return( df )
  }
  # Case B: No evidence in the component
  else if (!any(selector)) {
    return(forwardSampler.bnet(bn, N = N, vars = component))
  }
  # Case C: Mixed evidence -> MCMC
  else{
    return(mcmcSampler.bnet(bn, evidence, N = N, vars = component, burn.in = burn.in, thinnin = thinnin))
  }
}


#' @rdname sample.bnet
#'
#' @export
forwardSampler.bnet <- function(bn, N = 1, vars = variables(bn) ){

  # Ancestral order
  bn <- orden.bnet(bn)

  # Variables to sample (in order)
  vars <- names(bn$variables)[names(bn$variables) %in% vars]

  # Prealloc DF and set current state
  df <- as.data.frame(matrix(NA, ncol = length(vars), nrow = N))
  colnames(df) <- vars
  current.state <- as.list(rep(NA,length(vars)))
  names(current.state) <- vars

  # Start sampling
  for (i in 1:N) {
    for (j in 1:length(vars)) {
      # Get var.evidence (there is no null values since parents are always in current state)
      parent.values <- current.state[bn$variables[[ vars[j] ]]$parents]

      if (length(parent.values) > 0) {
        # Indexes
        # +1: position 0 is reserved for the var. itself
        var.pos <- vapply(names(parent.values),
                          function(x){which(bn$variables[[ vars[j] ]]$parents == x) + 1},
                          FUN.VALUE = numeric(1) )

        mop <- Rbmop::put_evidence.bmop(bn$variables[[ vars[j] ]]$prob, as.numeric(parent.values), evd.pos = var.pos, normalize = F)
      }
      else
        mop <- bn$variables[[ vars[j] ]]$prob

      current.state[j] <- metropolisSampler_auxiliar( as.function(mop, MIN = 1E-10),
                                                      min = Rbmop::lower.bmop(mop),
                                                      max = Rbmop::upper.bmop(mop) )
    }
    df[i,] <- current.state
  }
  return(df)
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
mcmcSampler.bnet <- function(bn,
                              evidence,
                              vars = names(bn$variables),
                              N = 1,
                              burn.in = 100,
                              thinnin = 10,
                              collapse = F) {
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
  if (collapse)
    bn <- collapse.bnet(bn, evidence)

  # Prealloc DF
  df <- as.data.frame( matrix(NA, ncol = length(vars), nrow = N) )
  colnames(df) <- c(variables.in.evidence, variables.not.evidence)

  # Initialize unassigned variables
  init.values <- vapply(variables.not.evidence, function(node){
    node.min <- min(bn$variables[[node]]$prob$knots[[1]])
    node.max <- max(bn$variables[[node]]$prob$knots[[1]])
    return(runif(1, min = node.min, max = node.max))
  }, FUN.VALUE = numeric(1))

  # Current state
  current.state <- as.list(c(evidence, init.values))
  names(current.state) <- c(variables.in.evidence, variables.not.evidence)

  # Burnin
  for (i in 1:burn.in) {
    # Select a node to be upodated at random
    node.to.update <- sample( variables.not.evidence, 1)
    # Compute full conditioned
    f <- full_conditioned.bnet(bn,
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
      }

      # Select a node to be upodated at random
      node.to.update <- sample( variables.not.evidence, 1)
      # Compute full conditioned
      f <- full_conditioned.bnet(bn,
                                 node.to.update,current.state,normalize = F)
      # Generate sample
      current.state[[node.to.update]] <- metropolisSampler_auxiliar( f$f, min = f$min, max = f$max )
    }
  }
  df[N,] <- current.state
  return(df)
}


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

##
# Metropolis - Hastings auxiliar sample for f
##
#'@noRd
#'@importFrom truncnorm rtruncnorm dtruncnorm
metropolisSampler_auxiliar <- function(f, burn.in = 100, min = -Inf, max = Inf, sd = 2 ){

  # Initial value
  x.old <- runif(1,min,max)
  f.old <- f(x.old)

  # Loop
  for (i in 1:burn.in) {
    # Compute new values
    x.new <- truncnorm::rtruncnorm(1, a = min, b = max, mean = x.old, sd = sd)
    f.new <- f(x.new)

    # Acceptance ratio

    acc.ratio <- (f.new / f.old) * (truncnorm::dtruncnorm(x.old, a = min, b = max, mean = x.new, sd = sd)/
                                    truncnorm::dtruncnorm(x.new, a = min, b = max, mean = x.old, sd = sd))

    if (is.na(acc.ratio)) {
      acc.ratio <- 0.5 # FIXME
    }
    if ((acc.ratio >= 1) || (runif(1) < acc.ratio) ) {
      x.old <- x.new
      f.old <- f.new
    }
  } # END FOR LOOP
  return(x.old)
}
