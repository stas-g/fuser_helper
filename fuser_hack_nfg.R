L2fuse.cv <- function(X, Y, groups, G = NULL, lambda.vals = NULL, gamma.vals = NULL, cv.folds = 5, max.cores = 10, mc.flag = TRUE) {

  n <- length(Y)
  k <- length(unique(groups))

  if(is.null(G)) {
    G <- matrix(1, k, k)
    diag(G) <- 0
      }

  group.names <- sort(unique(groups))
  # numerical.groups <- 1 : k
  # names(numerical.groups) <- group.names
  # numerical.groups <- numerical.groups[groups]

  # Hyperparameters for the methods to loop over
  if(is.null(gamma.vals)) gamma.vals <- c(0, 1e-5, 9e-5, 1e-4, 0.001, 0.01, 0.1, 0.5, 0.75, 1, 2.5)
  if(is.null(lambda.vals)) lambda.vals <- c(0, 1e-5, 9e-5, 1e-4, 2e-4, 5e-4, 0.001, 0.01, 0.1, 0.5, 0.75, 1, 2.5)

  hyper.vals <- expand.grid(lambda.vals, gamma.vals)

  # Rough scaling of hyperparameters to stay within the same range
  #
  # hyper.scale.factor <- (4/5)*(n/k) # Training set size divided by number of groups
  # gamma.vals <- gamma.vals / hyper.scale.factor
  # lambda.vals <- lambda.vals / hyper.scale.factor
  # hyper.vals <- expand.grid(lambda.vals, gamma.vals)

  total.space <- length(gamma.vals)*length(lambda.vals)
  data.split <- split.data(length(Y), cv.folds, groups)

  objective <- array(NA, c(total.space, cv.folds, k)) #param combinations x folds x cell types

  for(cv.i in 1 : cv.folds) {
      train.indices <- data.split$train[[cv.i]]

      X.train <- X[train.indices, ]
      Y.train <- Y[train.indices]
      groups.train <- groups[train.indices]

      if(mc.flag) {
      beta.results <- mclapply(gamma.vals, FUN = function(gamma) {
        message(gamma)
          L2fuse.fit(X.train, Y.train, groups.train, G, lambda.vals, gamma)
      }, mc.cores = length(gamma.vals))
  } else {
        beta.results <- lapply(gamma.vals, FUN = function(gamma) {
            message(gamma)
            L2fuse.fit(X.train, Y.train, groups.train, G, lambda.vals, gamma)
      })
  }

      # Retrieve validation data and scale according to training data.
      valid.indices <- data.split$valid[[cv.i]]
      X.valid <- X[valid.indices,]
      Y.valid <- Y[valid.indices]
      groups.valid <- groups[valid.indices]

      m <- 0

      # Loop over hyperparameters and calculate RMSE for validation
      for(gamma.i in 1 : length(gamma.vals)) {
        cat('Gamma:', gamma.vals[gamma.i], ' | ')

        l <- length(lambda.vals)

        beta.tmp <- beta.results[[gamma.i]]
        Y.predict <- L2fuse.pred(beta.tmp, X.valid, groups.valid)
        Y.valid.m <- matrix(rep(Y.valid, l), ncol = l)

        for(k.i in 1 : k) {
              group.indices <- groups.valid == group.names[k.i]
              objective[l * m + (1 : l), cv.i, k.i] <- colMeans((Y.valid.m[group.indices, ] - Y.predict[group.indices, ])^2)
        }

        cat(round(objective[l * m + (1 : l), cv.i, k.i], digits=2), ' ')

        m <- m + 1
          }

        cat('\n')
  }

  #average over different groups
  objective.global <- apply(objective, 1 : 2, median)

  #finding optimal parameters
  dat <- data.frame(objective.global, lambda = hyper.vals[, 1], gamma = hyper.vals[, 2])
  dat$m <- rowMeans(dat[, 1 : 5])
  dat$sd <- apply(dat[, 1 : 5], 1, sd)

  opt.val <- dat[which.min(dat$m), ]
  # ind <- (dat$m < opt.val$m + opt.val$sd) & (dat$m > opt.val$m - opt.val$sd) & (dat$lambda >= opt.val$lambda) & (dat$gamma >= opt.val$gamma)
  # se1.val <- dat[min(which(ind)), ]

  #fitting full model with optimal parameters
  opt.beta <- L2fuse.fit(X, Y, groups, G, opt.val$lambda, opt.val$gamma)
  # se1.beta <- L2fuse.fit(X, Y, groups, G, se1.val$lambda, se1.val$gamma)

  return(list(opt.val = opt.val, errors = objective, opt.beta = opt.beta))
}


L2fuse.fit <- function(X, Y, id, G, lambda, gamma, extra = FALSE, intercept = TRUE, ...) {
    # Convert groups to numerical to be safe
    groups <- sort(unique(id))
    numerical.groups <- setNames(1 : length(groups), groups)
    numerical.id <- numerical.groups[id]

    block.data <- generateBlockDiagonalMatrices(X, Y, numerical.id, G, intercept = intercept)

    betas <- fusedL2DescentGLMNet(block.data$X, block.data$X.fused, block.data$Y,
      numerical.id, lambda, gamma, ...)

    if(intercept) {
      dimnames(betas) <- list(c(colnames(X), 'alpha'), groups, lambda)
      X <- cbind(X, rep(1, nrow(X)))
    } else {
      dimnames(betas) <- list(colnames(X), groups, lambda)
  }

    if(extra) {#if want extra stuff calculated
        #calculate fitted values
        fitted.vals <- matrix(NA, nrow = length(id), ncol = length(lambda), dimnames = list(rownames(X), lambda))
        for(l in 1 : length(lambda)) {
          for(p in groups) {
            ind <- id == p
            fitted.vals[ind, l] <- X[ind, , drop = FALSE] %*% betas[, p, l, drop = FALSE]
          }
        }

        #calculate residuals
        obs <- matrix(rep(Y, length(lambda)), ncol = length(lambda))
        resids <- obs - fitted.vals
        #calculate rsq
        rsq <- sapply(groups, FUN = function(p) {
          calc.mse(obs[id == p,], fitted.vals[id == p,], rsq = TRUE)
        })

        return(list(betas = betas, fitted.vals = fitted.vals, resids = resids, rsq = rsq))
      }

    return(betas)
  }


L2fuse.pred <- function(betas, X, id, intercept = TRUE) {
  if(intercept) X <- cbind(X, rep(1, nrow(X)))
  groups <- dimnames(betas)[[2]]

  if(length(dim(betas)) == 3) {
    preds <- matrix(NA, nrow = length(id), ncol = dim(betas)[3], dimnames = list(rownames(X), dimnames(betas)[[3]]))

    for(l in 1 : dim(betas)[3]) {
      for(p in groups) {
        ind <- id == p
        preds[ind, l] <- X[ind, , drop = FALSE] %*% betas[, p, l, drop = FALSE]
      }
  }
    } else {
      preds <- rep(NA, length(id))

        for(p in groups) {
          ind <- id == p
          preds[ind] <- X[ind, , drop = FALSE] %*% betas[, p, drop = FALSE]
        }
  }
  preds
}


fusedL2DescentGLMNet.nfg <- function (transformed.x, transformed.x.f, transformed.y, groups,
    lambda, gamma = 1, ...)
{
    transformed.x.f = transformed.x.f * sqrt(gamma * (dim(transformed.x)[1] +
        dim(transformed.x.f)[1]))
    transformed.x = rbind(transformed.x, transformed.x.f)
    group.names = sort(unique(groups))
    num.groups = length(group.names)
    glmnet.result = glmnet(transformed.x, transformed.y, standardize = FALSE, lambda = lambda,
        ...)
    beta.mat = array(NA, c(dim(transformed.x)[2]/num.groups,
        num.groups, length(lambda)))
    for (lambda.i in 1:length(lambda)) {
        coef.temp = coef.glmnet(glmnet.result, s = lambda[lambda.i])
        beta.mat[, , lambda.i] = coef.temp[2:length(coef.temp)]
    }
    return(beta.mat)
}

##
