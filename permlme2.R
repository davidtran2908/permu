mymodelparm <- function(model, coef., vcov., df, ...) 
  UseMethod("mymodelparm")

mymodelparm.default <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
{
  ### extract coefficients and their covariance matrix
  if(inherits(model, "glmmTMB")) beta <- try(coef.(model)[["cond"]]) else beta <- try(coef.(model))
  if (inherits(beta, "try-error"))
    stop("no ", sQuote("coef"), " method for ",
         sQuote("model"), " found!")
  
  if(inherits(model, "glmmTMB")) sigma <- try(Matrix::as.matrix(vcov.(model)[["cond"]])) else sigma <- try(vcov.(model))
  if (inherits(sigma, "try-error"))
    stop("no ", sQuote("vcov"), " method for ",
         sQuote("model"), " found!")       
  sigma <- as.matrix(sigma)
  
  if (any(length(beta) != dim(sigma))) 
    beta = na.omit(beta)
  # stop("dimensions of coefficients and covariance matrix don't match")
  
  ### determine degrees of freedom
  if (is.null(df)) {
    df <- 0
    ### check if a linear model was supplied
    if (inherits(model, "aov") || inherits(model, "lm") || inherits(model, "glm")) {
      class(model) <- "lm"
      df <- summary(model)$df[2]
    }
    if (inherits(model, "gls")) {
      dd <- model$dims
      df <- dd[["N"]] - dd[["p"]]
    }
    if (inherits(model, "glmerMod") || inherits(model, "glmmTMB")) {
      df <- summary(model)$AICtab["df.resid"]
    }		
    if (inherits(model, "parm")) df <- model$df
  } else {
    if (df < 0) stop(sQuote("df"), " is not positive")
  }
  
  ### try to identify non-estimable coefficients
  ### coef.aov removes NAs, thus touch coefficients 
  ### directly
  if(inherits(model, "glmmTMB")) ocoef <- coef.(model)[["cond"]] else ocoef <- coef.(model)
  if (inherits(model, "aov")) ocoef <- model$coefficients
  estimable <- rep(TRUE, length(ocoef))
  if (any(is.na(ocoef))) {
    estimable[is.na(ocoef)] <- FALSE
    beta <- ocoef[estimable]
    if (dim(sigma)[1]==length(estimable)) sigma <- sigma[estimable, estimable]
  }
  
  ### just in case...
  if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable))
    stop("could not extract coefficients and covariance matrix from ", 
         sQuote("model"))
  
  RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
  class(RET) <- "mymodelparm"
  RET
}

mymodelparm.aovlist <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
  stop("This function does not support objects of class ", sQuote("aovlist"))

mymodelparm.lme <- function(model, coef. = nlme::fixef, vcov. = vcov, df = NULL, ...)
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

mymodelparm.lmerMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

mymodelparm.glmerMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

mymodelparm.gls <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)	

mymodelparm.glmmTMB <- function(model, coef. = glmmTMB::fixef, vcov. = vcov, df = NULL, ...)
  mymodelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)	

##------------------------------------------------------------------------
## block-diagonal matrix addition, multiplication, and trace functions matrix_list
##------------------------------------------------------------------------

# turn matrix into a list of sub-matrices

sub_f <- function(x, fac, dim) {
  function(f) switch(dim,
                     row = x[fac==f, ,drop=FALSE],
                     col = x[ ,fac==f, drop=FALSE],
                     both = x[fac==f, fac==f, drop=FALSE])
}

matrix_list <- function(x, fac, dim) {
  if (is.vector(x)) {
    if (dim != "both") stop(paste0("Object must be a matrix in order to subset by ",dim,"."))
    x_list <- split(x, fac)
    lapply(x_list, function(x) diag(x, nrow = length(x)))
  } else {
    lapply(levels(fac), sub_f(x, fac, dim))
  }
}

# turn block-diagonal into regular matrix

unblock <- function(A, block = attr(A, "groups")) {
  
  if (is.null(block)) block <- factor(rep(names(A), times = sapply(A, function(x) dim(x)[1])))
  n <- length(block)
  mat <- matrix(0, n, n)
  for (i in levels(block)) {
    index <- i == block
    mat[index,index] <- A[[i]]
  }
  return(mat)
}


# sum of two conformable block-diagonal matrices

sum_blockblock <- function(A, B)
  mapply(function(a,b) a + b, a = A, b = B, SIMPLIFY = FALSE)


# generic matrix minus block-diagonal

matrix_minus_block <- function(A, B, block=attr(B, "groups")) {
  if (is.null(block)) block <- rep(names(B), times = sapply(B, function(x) dim(x)[1]))
  
  mat <- A
  for (i in unique(block)) {
    index <- i == block
    mat[index,index] <- mat[index, index] - B[[i]]
  }
  return(mat)
}


# block-diagonal minus generic matrix

block_minus_matrix <- function(A, B, block = attr(A, "groups")) {
  if (is.null(block))
    block <- rep(names(A), times = sapply(A, function(x) dim(x)[1]))
  
  mat <- -B
  for (i in unique(block)) {
    index <- i == block
    mat[index,index] <- mat[index, index] + A[[i]]
  }
  return(mat)
}

add_submatrices <- function(indices, small_mat, big_mat) {
  levs <- levels(indices)
  if (nlevels(indices) != length(small_mat)) stop("Levels of indices do not match entries of small_mat.")
  for (i in 1:length(levs)) {
    ind <- levs[i] == indices
    big_mat[ind,ind] <- big_mat[ind,ind] + small_mat[[i]]
  }
  big_mat
}

add_bdiag <- function(small_mats, big_mats, crosswalk) {
  small_indices <- lapply(split(crosswalk[[1]], crosswalk[[2]]), droplevels)
  big_indices <- unique(crosswalk)
  big_indices <- big_indices[[2]][order(big_indices[[1]])]
  small_mats_list <- split(small_mats, big_indices)
  Map(add_submatrices, indices = small_indices, small_mat = small_mats_list, big_mat = big_mats)
}

# sum of conformable diagonal matrix and block-diagonal matrix 

add_diag <- function(d, M) {
  diag(M) <- diag(M) + d
  M
}

add_diag_bdiag <- function(diag_mats, big_mats) {
  Map(add_diag, d = diag_mats, M = big_mats)
}

get_order <- function(A, B, block = attr(A, "groups")) { 
  if (is.null(names(A))) names(A) <- 1:length(A)
  A_names <- names(A)
  C_indx <- rep(A_names, times = sapply(A, function(x) dim(x)[1]))
  if (is.null(block)) block <- C_indx
  C <- B
  B_index <- rep(0, nrow(B))
  
  for (b in A_names) {
    ind <- block == b
    ind_C <- C_indx == b
    C[ind_C, ] <-  B[ind,]
    B_index[ind_C] <- which(ind)
  }
  
  return(list(C, B_index))
}

get_cor_grouping <- function(mod, levels = NULL) {
  if (!is.null(mod$groups)) {
    struct <- mod$modelStruct$corStruct
    if (is.null(struct)) struct <- mod
    mod_formula <- nlme::getGroupsFormula(struct)
    dat <- nlme::getData(mod)
    if (inherits(na.action(mod), "exclude")) dat <- dat[-as.integer(na.action(mod)),,drop=FALSE]
    grps <- stats::model.frame(mod_formula, data = dat)
    grps <- apply(grps, 1, paste, collapse = "/")
    if (is.null(levels)) levels <- unique(grps)
    grps <- factor(grps, levels = levels)
  } else if (!is.null(mod$modelStruct$corStruct)) {
    grps <- factor(rep("A",mod$dims$N))
  } else {
    grps <- factor(1:mod$dims$N)
  }
  
  grps
}

# Construct list of block-diagonal correlation matrices build_Sigma_mats

build_corr_mats <- function(mod) {
  
  if (is.null(mod$modelStruct$corStruct)) {
    return(NULL)
  } else {
    R_list <- nlme::corMatrix(mod$modelStruct$corStruct)
    grps <- get_cor_grouping(mod, levels = names(R_list))
    if (!is.list(R_list)) R_list <- list(A = R_list)
    attr(R_list, "groups") <- grps
    return(R_list)
  }
}

# Construct list of block-diagonal lowest-level var-cov matrices

get_sort_order <- function(mod) {
  groups <- mod$groups
  if (is.data.frame(groups)) {
    order(do.call(order, groups))
  } else if (!is.null(groups)) {
    order(order(groups))
  } else {
    1:mod$dims$N
  }
}

get_cor_grouping <- function(mod, levels = NULL) {
  if (!is.null(mod$groups)) {
    struct <- mod$modelStruct$corStruct
    if (is.null(struct)) struct <- mod
    mod_formula <- nlme::getGroupsFormula(struct)
    dat <- nlme::getData(mod)
    if (inherits(na.action(mod), "exclude")) dat <- dat[-as.integer(na.action(mod)),,drop=FALSE]
    grps <- stats::model.frame(mod_formula, data = dat)
    grps <- apply(grps, 1, paste, collapse = "/")
    if (is.null(levels)) levels <- unique(grps)
    grps <- factor(grps, levels = levels)
  } else if (!is.null(mod$modelStruct$corStruct)) {
    grps <- factor(rep("A",mod$dims$N))
  } else {
    grps <- factor(1:mod$dims$N)
  }
  
  grps
}

build_var_cor_mats <- function(mod, R_list = build_corr_mats(mod), sigma_scale = FALSE) {
  
  sigma <- if (sigma_scale) mod$sigma else 1
  
  if (is.null(R_list)) {
    
    # if there is no correlation structure,
    # then build block-diagonals with first available grouping variable
    
    if (is.null(mod$groups)) {
      
      # if there are no groups then make diagonal matrix-lists
      
      if (is.null(mod$modelStruct$varStruct)) {
        V_list <- as.list(rep(sigma^2, mod$dims$N))
      } else {
        sd_vec <- sigma / as.numeric(nlme::varWeights(mod$modelStruct$varStruct))
        V_list <- as.list(sd_vec^2)
      }
      grps <- factor(1:mod$dims$N)
      attr(V_list, "groups") <- grps
      names(V_list) <- levels(grps)
      
    } else {
      
      # if there are groups then make block-diagonal matrix-lists
      
      if (is.null(mod$modelStruct$varStruct)) {
        grps <- mod$groups[[1]]
        V_list <- tapply(rep(sigma^2, length(grps)),  grps, function(x) diag(x, nrow = length(x)))
      } else {
        sort_order <- get_sort_order(mod)
        sd_vec <- sigma / as.numeric(nlme::varWeights(mod$modelStruct$varStruct))[sort_order]
        V_list <- tapply(sd_vec^2, mod$groups[[1]], function(x) diag(x, nrow = length(x)))
      }
      attr(V_list, "groups") <- mod$groups[[1]]
    }
    
  } else {
    
    # if there is a correlation structure,
    # build block-diagonals according to its grouping structure
    
    if (is.null(mod$modelStruct$varStruct)) {
      V_list <- if (sigma_scale) lapply(R_list, function(x) x * mod$sigma^2) else R_list
    } else {
      sort_order <- get_sort_order(mod)
      sd_vec <- sigma / as.numeric(nlme::varWeights(mod$modelStruct$varStruct))[sort_order]
      sd_list <- split(sd_vec, attr(R_list, "groups"))
      V_list <- Map(function(R, s) tcrossprod(s) * R, R = R_list, s = sd_list)
    }
    
    attr(V_list, "groups") <- attr(R_list, "groups")
  }
  
  return(V_list)
}

# Create block-diagonal covariance structure from Z-design and Tau matrices

ZDZt <- function(D, Z_list) {
  lapply(Z_list, function(z) z %*% D %*% t(z))
}

# Construct list of block-diagonal matrices for each random effects grouping structure

build_RE_mats <- function(mod, sigma_scale = FALSE) {
  
  # Get random effects structure
  all_groups <- rev(mod$groups)
  
  if (length(all_groups) == 1) {
    
    D_mat <- as.matrix(mod$modelStruct$reStruct[[1]])
    if (sigma_scale) D_mat <- mod$sigma^2 * D_mat
    data <- nlme::getData(mod)
    Z_mat <- model.matrix(mod$modelStruct$reStruc, data[complete.cases(data), ])
    row.names(Z_mat) <- NULL
    Z_list <- matrix_list(Z_mat, all_groups[[1]], "row")
    ZDZ_list <- ZDZt(D_mat, Z_list)
    
    attr(ZDZ_list, "groups") <- all_groups[[1]]
    
  } else {
    if (sigma_scale) {
      D_list <- lapply(mod$modelStruct$reStruct, function(x) mod$sigma^2 * as.matrix(x))
    } else {
      D_list <- lapply(mod$modelStruct$reStruct, as.matrix)
    }
    data <- nlme::getData(mod)
    Z_mat <- model.matrix(mod$modelStruct$reStruc, data[complete.cases(data), ])
    Z_names <- sapply(strsplit(colnames(Z_mat), ".", fixed=TRUE), function(x) x[1])
    row.names(Z_mat) <- NULL
    Z_levels <- lapply(names(all_groups), function(x) Z_mat[,x==Z_names,drop=FALSE])
    Z_levels <- Map(matrix_list, x = Z_levels, fac = all_groups, dim = "row")
    ZDZ_lists <- Map(ZDZt, D = D_list, Z_list = Z_levels)
    # ZDZ_lists <- Map(function(x,fac) x[order(fac)], x = ZDZ_lists, fac = all_groups)
    
    for (i in 2:length(all_groups)) {
      ZDZ_lists[[i]] <- add_bdiag(small_mats = ZDZ_lists[[i-1]],
                                  big_mats = ZDZ_lists[[i]],
                                  crosswalk = all_groups[c(i-1,i)])
    }
    
    ZDZ_list <- ZDZ_lists[[i]]
    
    attr(ZDZ_list, "groups") <- all_groups[[i]]
    
  }
  
  ZDZ_list
  
}

build_Sigma_mats <- function(mod, invert = FALSE, sigma_scale = FALSE) UseMethod("build_Sigma_mats")

build_Sigma_mats.default <- function(mod, invert = FALSE, sigma_scale = FALSE) {
  mod_class <- paste(class(mod), collapse = "-")
  stop(paste0("Sigma matrices not available for models of class ", mod_class, "."))
}

build_Sigma_mats.gls <- function(mod, invert = FALSE, sigma_scale = FALSE) {
  
  # lowest-level covariance structure
  V_list <- build_var_cor_mats(mod, sigma_scale = sigma_scale)
  V_grps <- attr(V_list, "groups")
  
  if (invert) {
    V_list <- lapply(V_list, function(x) chol2inv(chol(x)))
    attr(V_list, "groups") <- V_grps
  }
  
  
  return(V_list)
}

build_Sigma_mats.lme <- function(mod, invert = FALSE, sigma_scale = FALSE) {
  
  if (inherits(mod, "nlme")) stop("not implemented for \"nlme\" objects")
  
  # lowest-level covariance structure
  V_list <- build_var_cor_mats(mod, sigma_scale = sigma_scale)
  
  # random effects covariance structure
  ZDZ_list <- build_RE_mats(mod, sigma_scale = sigma_scale)
  
  V_grps <- attr(V_list, "groups")
  
  # Check if lowest-level covariance structure is nested within RE structure
  ZDZ_grps <- attr(ZDZ_list, "groups")
  group_mapping <- tapply(ZDZ_grps, V_grps, function(x) length(unique(x)))
  nested <- all(group_mapping == 1L)
  
  if (nested) {
    Sigma_list <- add_bdiag(V_list, ZDZ_list, data.frame(V_grps, ZDZ_grps))
    Sigma_grps <- attr(ZDZ_list, "groups")
  } else {
    V_mat <- unblock(V_list, block = V_grps)
    ZDZ_mat <- unblock(ZDZ_list, block = ZDZ_grps)
    Sigma_list <- V_mat + ZDZ_mat
    Sigma_grps <- factor(rep("A", nrow(Sigma_list)))
  }
  
  if (invert) {
    Sigma_list <- lapply(Sigma_list, function(x) chol2inv(chol(x)))
  }
  
  attr(Sigma_list, "groups") <- Sigma_grps
  
  return(Sigma_list)
}

## thien phuc's suggestion 3
permmodels.tp3 <- function(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3),
                           test.statistic=c("Chisq", "F", "LR", "Wald"),  exact=FALSE, 
                           data=NULL, fo=NULL, prt=TRUE, ncore=3, seed) { 
  
  options(scipen=6)
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)	
  if (inherits(model, "glm") && !(test.statistic %in% c("LR", "Wald", "F"))) test.statistic <- "LR"
  if (any(inherits(model, "gls"), inherits(model, "lme"))) test.statistic <- "Chisq"
  if (type %in% c("I", "1")) test.statistic <- "F"	  
  if (class(model)[1]=="aovlist") stop("Plese use model 'lme' instead of 'aov'!")
  if (inherits(model, "glmerMod")) stop("This function is not applied to 'glmer' yet!")
  
  if (any(inherits(model, "lm"), inherits(model, "aov"), inherits(model, "glm"))) {
    if (!is.null(data)) mod_df <- data else {
      if (inherits(model, "glm")) mod_df <- as.data.frame(model$data) else mod_df <- as.data.frame(model.frame(model))
    }
    Terms <- terms(model)
    yname <- as.character(attr(Terms, "variables"))[[2]]    
    if (grepl("[,]", yname)) {
      yname <- unlist(strsplit(yname, "[,] "))[2]
      yname <- gsub("\\)", "", unlist(strsplit(yname, " - "))) # ynames
    }
    diff_yname <- setdiff(colnames(mod_df), yname)
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      mod_df[sample(1:nrow(mod_df)), yname, drop=FALSE]
    }, simplify = !(length(yname) > 1))
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x
        rfitmodel <- try(update(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      
    } else{	
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        if (!(type %in% c("I", "1"))) clusterEvalQ(cl, library(car))
        clusterExport(cl, c("mymodelparm", "mymodelparm.default", "model", "mod_df", "yname", "diff_yname", "type", "test.statistic"), envir = environment()) 
        permmod <- parLapplyLB(cl, permy, function(x) {
          if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x		
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      }else{  
        permmod <- mclapply(permy, function(x) {
          if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }    
  }
  
  if (any(inherits(model, "gls"), inherits(model, "lme"))) { 
    # fm$call$fixed
    # fm$call$random  
    model <- update(model, na.action=na.omit)
    mod_df <- as.data.frame(nlme::getData(model))
    #X <- model.matrix(model) ## thien phuc froze
    
    mod0 <- lme(y ~ 0 + Study + y0:Study, random=list(Study=~0+Group), weights=varIdent(form=~1|Study),
                na.action = na.omit, method = "REML", data=mod_df)
    X <- reStruct(list(y = ~ 0 + Study + y0:Study), data = mod_df) ## thien phuc
    X <- model.matrix(X, mod_df)
    
    #if (inherits(model, "gls")) beta <- coef(model) else beta <- nlme::fixef(model)
    if (inherits(model, "gls")) beta <- coef(mod0) else beta <- nlme::fixef(mod0) ## thien phuc
    
    # check for columns dropped from model
    col_names <- names(beta)
    #if (ncol(X) != length(col_names)) X <- X[,col_names,drop=FALSE]  ## thien phuc froze  
    
    #Terms <- terms(model)
    Terms <- terms(mod0) ## thien phuc
    yname <- as.character(attr(Terms, "variables"))[[2]]
    #y <- nlme::getResponse(model)
    y <- nlme::getResponse(mod0) ## thien phuc
    X_y <- cbind(X, y)
    #V_list <- build_Sigma_mats(model, sigma_scale =FALSE)
    V_list <- build_Sigma_mats(mod0, sigma_scale =FALSE) ## thien phuc
    ord_X_y <- get_order(A=V_list, B=X_y)
    X_y <- ord_X_y[[1]]
    mod_df <- mod_df[ord_X_y[[2]], ]
    xbeta <- as.vector(X_y[, col_names]%*%beta)
    #errors <- X_y[, "y"] - xbeta
    
    errors <- X_y[, "y"] - xbeta
    
    V_matrix <- Matrix::bdiag(V_list)
    Ut <- t(chol(V_matrix))
    wt <- solve(Ut)
    
    #Weighting the residuals.
    wterrors <- wt%*%errors    
    
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      as.data.frame(perm_y <- xbeta+as.vector(Ut%*%wterrors[rowindex]))
      #as.data.frame(perm_y <- xbeta-ztheta+as.vector(Ut%*%wterrors[rowindex])) ## long hao
      #as.data.frame(perm_y <- xbeta[rowindex]-ztheta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## long hao
    })
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        mod_df[, yname] <- x
        rfitmodel <- try(update(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      
    }else{
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        if (!(type %in% c("I", "1"))) clusterEvalQ(cl, library(car)) 
        clusterEvalQ(cl, library(nlme))
        clusterExport(cl, c("mymodelparm", "mymodelparm.lme", "mymodelparm.gls","mymodelparm.default", "model", "mod_df", "yname", "type", "test.statistic"), envir = environment()) 
        permmod <- parLapplyLB(cl, permy, function(x) {
          mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      }else{  
        permmod <- mclapply(permy, function(x) {
          mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }
  }
  
  if (any(inherits(model, "lmerMod"), inherits(model, "glmerMod"))) { 
    lmer_var_names <- all.vars(model@call$formula)
    lmer_var_names <- lmer_var_names[lmer_var_names!="pi"]
    mod_df <- as.data.frame(nlme::getData(model))
    mod_df <- na.omit(mod_df[, lmer_var_names])	
    model <- update(model, data=mod_df)	
    theta <- getME(model, "theta")
    fixef <- fixef(model) 
    Lambda <- getME(model, "Lambda")
    Lambdac <- Matrix::tcrossprod(Lambda)
    V <- getME(model, "Z")%*%Lambdac%*%getME(model, "Zt")+diag(dim(getME(model, "Z"))[1])
    Ut <- t(chol(V))
    wt <- solve(Ut)
    xbeta <- as.vector(getME(model, "X")%*%fixef)
    if (inherits(model, "glmerMod")) {
      errors <- slot(model, "resp")$family$linkfun(getME(model, "y")) - xbeta 
    }else errors <- getME(model, "y") - xbeta
    
    #Weighting the residuals.
    wterrors <- wt%*%errors
    
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      if (inherits(model, "glmerMod")) {
        perm_y <- slot(model, "resp")$family$linkinv(xbeta+as.vector(Ut%*%wterrors[rowindex])) 
        #perm_y <- slot(model, "resp")$family$linkinv(xbeta-ztheta+as.vector(Ut%*%wterrors[rowindex])) ## longhao
        #perm_y <- slot(model, "resp")$family$linkinv(xbeta[rowindex]-ztheta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## longhao
        if (slot(model, "resp")$family$family == "poisson") { 
          perm_y[perm_y<=0] <- 0
          perm_y <- round(perm_y, 0)
        }
        if (slot(model, "resp")$family$family == "binomial") {
          perm_y[perm_y<=0] <- 0
          perm_y[perm_y>=1] <- 1		
        }
        if (slot(model, "resp")$family$family == "gamma") {
          perm_y[perm_y<=0] <- 0		
        }
        as.data.frame(perm_y)
      }else{
        as.data.frame(perm_y <- xbeta+as.vector(Ut%*%wterrors[rowindex]))
        #as.data.frame(perm_y <- xbeta-ztheta+as.vector(Ut%*%wterrors[rowindex]))  ## long hao
        #as.data.frame(perm_y <- xbeta[rowindex]-ztheta[rowindex]+as.vector(Ut%*%wterrors[rowindex]))  ## long hao
      } 
    })
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        rfitmodel <- try(refit(model, x), TRUE)
        if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      
    }else{
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        clusterEvalQ(cl, library(lme4))
        clusterExport(cl, c("mymodelparm", "mymodelparm.lmerMod", "mymodelparm.default", "model", "type", "test.statistic"), envir = environment())       
        permmod <- parLapplyLB(cl, permy, function(x) {
          rfitmodel <- try(refit(model, x), TRUE)
          if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")){
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      }else{  
        permmod <- mclapply(permy, function(x) {
          rfitmodel <- try(refit(model, x), TRUE)
          if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }	
  }    
  permmod <- permmod[!(sapply(permmod,is.null))]
  nperm <- length(permmod)
  if (nperm==0) stop("permutation can't produce useful data!")
  if (!(inherits(model, "aov"))) {
    # tTable <- function(x) { # function to construct t-table for the model
    # if (class(model)[1]=="lmerMod"){
    # mp <- mymodelparm(x)
    # tTable <- cbind(mp$coef, sqrt(base::diag(mp$vcov)), mp$coef/sqrt(base::diag(mp$vcov)))
    # colnames(tTable) <- c("Estimate", "Std. Error", "t value")
    # return(round(tTable, 4))		
    # }else{
    # summ <- summary(x)
    # if (class(x)[1] %in% c("lme", "gls")) {
    # tTable <- summ$tTable
    # }else{
    # tTable <- coef(summ)
    # }
    # return(tTable)
    # }
    # }
    
    # model_tTable <- tTable(model)
    # Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
    
    # if (nrow(model_tTable)==1) {
    # Tper.p <- (sum(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(model_tTable[, Tvalue]),6))+!exact)/(nperm+!exact)
    # }else{
    # Tper.p <- (rowSums(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(model_tTable[, Tvalue]), 6))+!exact)/(nperm+!exact)
    # }
    
    model_tTable <- coef(summary(model))
    Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
    if (nrow(model_tTable)==1) {
      Tper.p <- (sum(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]),6))+!exact)/(nperm+!exact)
    }else{
      Tper.p <- (rowSums(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]), 6))+!exact)/(nperm+!exact)
    }
    
    if ("(Intercept)" %in% rownames(model_tTable)) Tper.p[1] <- NA
    COEFFICENTS <- round(cbind(model_tTable, "Perm_p_value"=Tper.p),4)
    if (prt) {
      cat("\nCoefficients of (fixed) effects:\n")
      print(COEFFICENTS)
      cat("\nNote: Perm_p_value of t test is obtained using", sQuote(nperm), "permutations.\n")
    }
  }else COEFFICENTS <- NULL
  
  if (type %in% c("I", "1")) {
    key_anova <- anova(model)
    Fvalue <- colnames(key_anova)[grep("F.value",colnames(key_anova))]
    if (class(model)[1] %in% c("glm")) Fvalue <- "Deviance"
    if (nrow(key_anova)==1) {
      Fper.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
    }else{
      Fper.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
    } 
    if ("(Intercept)" %in% rownames(key_anova)) Fper.p[1] <- NA	
    ANOVA <- round(cbind(key_anova, "Perm_p_value"=Fper.p),4)
  }else{
    key_anova <- car::Anova(model, type=type, test.statistic=test.statistic)
    if (class(model)[1] %in% c("lm", "aov")) stats_value <- "F value"  else stats_value <-  test.statistic 
    if (class(model)[1] %in% c("glm")) stats_value <- switch(test.statistic, LR = "LR Chisq", Wald = "Chisq", F = ifelse(type %in% c("2", "II"), "F value", "F values"))  
    if (nrow(key_anova)==1) {
      stats.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
    }else{
      stats.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
    }  
    if ("(Intercept)" %in% rownames(key_anova)) stats.p[1] <- NA
    ANOVA <- round(cbind(key_anova, "Perm_p_value"=stats.p),4) 
    attr(ANOVA, "heading") <- attr(key_anova, "heading")
  }
  
  ANOVA[is.na(ANOVA)] <-""
  if (prt) {
    cat("\nANOVA:\n")
    if (!is.null(attr(ANOVA, "heading"))) cat(paste(attr(ANOVA, "heading"), "\n", sep=""))
    print(ANOVA)
    cat(paste("\nNote: Perm_p_value of", test.statistic, "test is obtained using"), sQuote(nperm), "permutations.\n\n")
  }
  permlist <- vector("list", nperm)
  for (i in 1:nperm) {
    permlist[[i]] <- permmod[[i]][[1]]
  }
  return(invisible(list("permlist"=permlist, "COEFFICENTS"=COEFFICENTS, "ANOVA"=ANOVA)))
} 

## thien phuc's suggestion 3
permmodels.tp3 <- function(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3),
                           test.statistic=c("Chisq", "F", "LR", "Wald"),  exact=FALSE, 
                           data=NULL, fo=NULL, prt=TRUE, ncore=3, seed) { 
  
  options(scipen=6)
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)	
  if (inherits(model, "glm") && !(test.statistic %in% c("LR", "Wald", "F"))) test.statistic <- "LR"
  if (any(inherits(model, "gls"), inherits(model, "lme"))) test.statistic <- "Chisq"
  if (type %in% c("I", "1")) test.statistic <- "F"	  
  if (class(model)[1]=="aovlist") stop("Plese use model 'lme' instead of 'aov'!")
  if (inherits(model, "glmerMod")) stop("This function is not applied to 'glmer' yet!")
  
  if (any(inherits(model, "lm"), inherits(model, "aov"), inherits(model, "glm"))) {
    if (!is.null(data)) mod_df <- data else {
      if (inherits(model, "glm")) mod_df <- as.data.frame(model$data) else mod_df <- as.data.frame(model.frame(model))
    }
    Terms <- terms(model)
    yname <- as.character(attr(Terms, "variables"))[[2]]    
    if (grepl("[,]", yname)) {
      yname <- unlist(strsplit(yname, "[,] "))[2]
      yname <- gsub("\\)", "", unlist(strsplit(yname, " - "))) # ynames
    }
    diff_yname <- setdiff(colnames(mod_df), yname)
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      mod_df[sample(1:nrow(mod_df)), yname, drop=FALSE]
    }, simplify = !(length(yname) > 1))
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x
        rfitmodel <- try(update(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      
    } else{	
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        if (!(type %in% c("I", "1"))) clusterEvalQ(cl, library(car))
        clusterExport(cl, c("mymodelparm", "mymodelparm.default", "model", "mod_df", "yname", "diff_yname", "type", "test.statistic"), envir = environment()) 
        permmod <- parLapplyLB(cl, permy, function(x) {
          if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x		
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      }else{  
        permmod <- mclapply(permy, function(x) {
          if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }    
  }
  
  if (any(inherits(model, "gls"), inherits(model, "lme"))) { 
    # fm$call$fixed
    # fm$call$random  
    model <- update(model, na.action=na.omit)
    mod_df <- as.data.frame(nlme::getData(model))
    #X <- model.matrix(model) ## thien phuc froze
    
    mod0 <- lme(y ~ 0 + Study + y0:Study, random=list(Study=~0+Group), weights=varIdent(form=~1|Study),
                na.action = na.omit, method = "REML", data=mod_df)
    X <- reStruct(list(y = ~ 0 + Study + y0:Study), data = mod_df) ## thien phuc
    X <- model.matrix(X, mod_df)
    
    #if (inherits(model, "gls")) beta <- coef(model) else beta <- nlme::fixef(model)
    if (inherits(model, "gls")) beta <- coef(mod0) else beta <- nlme::fixef(mod0) ## thien phuc
    
    # check for columns dropped from model
    col_names <- names(beta)
    #if (ncol(X) != length(col_names)) X <- X[,col_names,drop=FALSE]  ## thien phuc froze  
    
    #Terms <- terms(model)
    Terms <- terms(mod0) ## thien phuc
    yname <- as.character(attr(Terms, "variables"))[[2]]
    #y <- nlme::getResponse(model)
    y <- nlme::getResponse(mod0) ## thien phuc
    X_y <- cbind(X, y)
    #V_list <- build_Sigma_mats(model, sigma_scale =FALSE)
    V_list <- build_Sigma_mats(mod0, sigma_scale =FALSE) ## thien phuc
    ord_X_y <- get_order(A=V_list, B=X_y)
    X_y <- ord_X_y[[1]]
    mod_df <- mod_df[ord_X_y[[2]], ]
    xbeta <- as.vector(X_y[, col_names]%*%beta)
    #errors <- X_y[, "y"] - xbeta
    
    errors <- X_y[, "y"] - xbeta
    
    V_matrix <- Matrix::bdiag(V_list)
    Ut <- t(chol(V_matrix))
    wt <- solve(Ut)
    
    #Weighting the residuals.
    wterrors <- wt%*%errors    
    
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      as.data.frame(perm_y <- xbeta+as.vector(Ut%*%wterrors[rowindex]))
      #as.data.frame(perm_y <- xbeta-ztheta+as.vector(Ut%*%wterrors[rowindex])) ## long hao
      #as.data.frame(perm_y <- xbeta[rowindex]-ztheta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## long hao
    })
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        mod_df[, yname] <- x
        rfitmodel <- try(update(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      
    }else{
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        if (!(type %in% c("I", "1"))) clusterEvalQ(cl, library(car)) 
        clusterEvalQ(cl, library(nlme))
        clusterExport(cl, c("mymodelparm", "mymodelparm.lme", "mymodelparm.gls","mymodelparm.default", "model", "mod_df", "yname", "type", "test.statistic"), envir = environment()) 
        permmod <- parLapplyLB(cl, permy, function(x) {
          mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      }else{  
        permmod <- mclapply(permy, function(x) {
          mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }
  }
  
  if (any(inherits(model, "lmerMod"), inherits(model, "glmerMod"))) { 
    lmer_var_names <- all.vars(model@call$formula)
    lmer_var_names <- lmer_var_names[lmer_var_names!="pi"]
    mod_df <- as.data.frame(nlme::getData(model))
    mod_df <- na.omit(mod_df[, lmer_var_names])	
    model <- update(model, data=mod_df)	
    theta <- getME(model, "theta")
    fixef <- fixef(model) 
    Lambda <- getME(model, "Lambda")
    Lambdac <- Matrix::tcrossprod(Lambda)
    V <- getME(model, "Z")%*%Lambdac%*%getME(model, "Zt")+diag(dim(getME(model, "Z"))[1])
    Ut <- t(chol(V))
    wt <- solve(Ut)
    xbeta <- as.vector(getME(model, "X")%*%fixef)
    if (inherits(model, "glmerMod")) {
      errors <- slot(model, "resp")$family$linkfun(getME(model, "y")) - xbeta 
    }else errors <- getME(model, "y") - xbeta
    
    #Weighting the residuals.
    wterrors <- wt%*%errors
    
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      if (inherits(model, "glmerMod")) {
        perm_y <- slot(model, "resp")$family$linkinv(xbeta+as.vector(Ut%*%wterrors[rowindex])) 
        #perm_y <- slot(model, "resp")$family$linkinv(xbeta-ztheta+as.vector(Ut%*%wterrors[rowindex])) ## longhao
        #perm_y <- slot(model, "resp")$family$linkinv(xbeta[rowindex]-ztheta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## longhao
        if (slot(model, "resp")$family$family == "poisson") { 
          perm_y[perm_y<=0] <- 0
          perm_y <- round(perm_y, 0)
        }
        if (slot(model, "resp")$family$family == "binomial") {
          perm_y[perm_y<=0] <- 0
          perm_y[perm_y>=1] <- 1		
        }
        if (slot(model, "resp")$family$family == "gamma") {
          perm_y[perm_y<=0] <- 0		
        }
        as.data.frame(perm_y)
      }else{
        as.data.frame(perm_y <- xbeta+as.vector(Ut%*%wterrors[rowindex]))
        #as.data.frame(perm_y <- xbeta-ztheta+as.vector(Ut%*%wterrors[rowindex]))  ## long hao
        #as.data.frame(perm_y <- xbeta[rowindex]-ztheta[rowindex]+as.vector(Ut%*%wterrors[rowindex]))  ## long hao
      } 
    })
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        rfitmodel <- try(refit(model, x), TRUE)
        if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      
    }else{
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        clusterEvalQ(cl, library(lme4))
        clusterExport(cl, c("mymodelparm", "mymodelparm.lmerMod", "mymodelparm.default", "model", "type", "test.statistic"), envir = environment())       
        permmod <- parLapplyLB(cl, permy, function(x) {
          rfitmodel <- try(refit(model, x), TRUE)
          if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")){
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      }else{  
        permmod <- mclapply(permy, function(x) {
          rfitmodel <- try(refit(model, x), TRUE)
          if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }	
  }    
  permmod <- permmod[!(sapply(permmod,is.null))]
  nperm <- length(permmod)
  if (nperm==0) stop("permutation can't produce useful data!")
  if (!(inherits(model, "aov"))) {
    # tTable <- function(x) { # function to construct t-table for the model
    # if (class(model)[1]=="lmerMod"){
    # mp <- mymodelparm(x)
    # tTable <- cbind(mp$coef, sqrt(base::diag(mp$vcov)), mp$coef/sqrt(base::diag(mp$vcov)))
    # colnames(tTable) <- c("Estimate", "Std. Error", "t value")
    # return(round(tTable, 4))		
    # }else{
    # summ <- summary(x)
    # if (class(x)[1] %in% c("lme", "gls")) {
    # tTable <- summ$tTable
    # }else{
    # tTable <- coef(summ)
    # }
    # return(tTable)
    # }
    # }
    
    # model_tTable <- tTable(model)
    # Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
    
    # if (nrow(model_tTable)==1) {
    # Tper.p <- (sum(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(model_tTable[, Tvalue]),6))+!exact)/(nperm+!exact)
    # }else{
    # Tper.p <- (rowSums(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(model_tTable[, Tvalue]), 6))+!exact)/(nperm+!exact)
    # }
    
    model_tTable <- coef(summary(model))
    Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
    if (nrow(model_tTable)==1) {
      Tper.p <- (sum(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]),6))+!exact)/(nperm+!exact)
    }else{
      Tper.p <- (rowSums(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]), 6))+!exact)/(nperm+!exact)
    }
    
    if ("(Intercept)" %in% rownames(model_tTable)) Tper.p[1] <- NA
    COEFFICENTS <- round(cbind(model_tTable, "Perm_p_value"=Tper.p),4)
    if (prt) {
      cat("\nCoefficients of (fixed) effects:\n")
      print(COEFFICENTS)
      cat("\nNote: Perm_p_value of t test is obtained using", sQuote(nperm), "permutations.\n")
    }
  }else COEFFICENTS <- NULL
  
  if (type %in% c("I", "1")) {
    key_anova <- anova(model)
    Fvalue <- colnames(key_anova)[grep("F.value",colnames(key_anova))]
    if (class(model)[1] %in% c("glm")) Fvalue <- "Deviance"
    if (nrow(key_anova)==1) {
      Fper.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
    }else{
      Fper.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
    } 
    if ("(Intercept)" %in% rownames(key_anova)) Fper.p[1] <- NA	
    ANOVA <- round(cbind(key_anova, "Perm_p_value"=Fper.p),4)
  }else{
    key_anova <- car::Anova(model, type=type, test.statistic=test.statistic)
    if (class(model)[1] %in% c("lm", "aov")) stats_value <- "F value"  else stats_value <-  test.statistic 
    if (class(model)[1] %in% c("glm")) stats_value <- switch(test.statistic, LR = "LR Chisq", Wald = "Chisq", F = ifelse(type %in% c("2", "II"), "F value", "F values"))  
    if (nrow(key_anova)==1) {
      stats.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
    }else{
      stats.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
    }  
    if ("(Intercept)" %in% rownames(key_anova)) stats.p[1] <- NA
    ANOVA <- round(cbind(key_anova, "Perm_p_value"=stats.p),4) 
    attr(ANOVA, "heading") <- attr(key_anova, "heading")
  }
  
  ANOVA[is.na(ANOVA)] <-""
  if (prt) {
    cat("\nANOVA:\n")
    if (!is.null(attr(ANOVA, "heading"))) cat(paste(attr(ANOVA, "heading"), "\n", sep=""))
    print(ANOVA)
    cat(paste("\nNote: Perm_p_value of", test.statistic, "test is obtained using"), sQuote(nperm), "permutations.\n\n")
  }
  permlist <- vector("list", nperm)
  for (i in 1:nperm) {
    permlist[[i]] <- permmod[[i]][[1]]
  }
  return(invisible(list("permlist"=permlist, "COEFFICENTS"=COEFFICENTS, "ANOVA"=ANOVA)))
} 

## permutation for blup
permblup <- function(model, nperm=4999, type=c("I", "II", "III", 1, 2, 3),
                     test.statistic=c("Chisq", "F", "LR", "Wald"),  exact=FALSE, 
                     data=NULL, fo=NULL, prt=TRUE, ncore=3, seed) { 
  
  options(scipen=6)
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)	
  if (inherits(model, "glm") && !(test.statistic %in% c("LR", "Wald", "F"))) test.statistic <- "LR"
  if (any(inherits(model, "gls"), inherits(model, "lme"))) test.statistic <- "Chisq"
  if (type %in% c("I", "1")) test.statistic <- "F"	  
  if (class(model)[1]=="aovlist") stop("Plese use model 'lme' instead of 'aov'!")
  if (inherits(model, "glmerMod")) stop("This function is not applied to 'glmer' yet!")
  
  if (any(inherits(model, "lm"), inherits(model, "aov"), inherits(model, "glm"))) {
    if (!is.null(data)) mod_df <- data else {
      if (inherits(model, "glm")) mod_df <- as.data.frame(model$data) else mod_df <- as.data.frame(model.frame(model))
    }
    Terms <- terms(model)
    yname <- as.character(attr(Terms, "variables"))[[2]]    
    if (grepl("[,]", yname)) {
      yname <- unlist(strsplit(yname, "[,] "))[2]
      yname <- gsub("\\)", "", unlist(strsplit(yname, " - "))) # ynames
    }
    diff_yname <- setdiff(colnames(mod_df), yname)
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      mod_df[sample(1:nrow(mod_df)), yname, drop=FALSE]
    }, simplify = !(length(yname) > 1))
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x
        rfitmodel <- try(update(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      
    } else{	
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        if (!(type %in% c("I", "1"))) clusterEvalQ(cl, library(car))
        clusterExport(cl, c("mymodelparm", "mymodelparm.default", "model", "mod_df", "yname", "diff_yname", "type", "test.statistic"), envir = environment()) 
        permmod <- parLapplyLB(cl, permy, function(x) {
          if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x		
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      }else{  
        permmod <- mclapply(permy, function(x) {
          if (length(yname) > 1) mod_df <- cbind(mod_df[, diff_yname], x) else mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }    
  }
  
  if (any(inherits(model, "gls"), inherits(model, "lme"))) { 
    # fm$call$fixed
    # fm$call$random  
    model <- update(model, na.action=na.omit)
    mod_df <- as.data.frame(nlme::getData(model))
    #X <- model.matrix(model) ## thien phuc froze
    
    X <- reStruct(list(y = ~ 0 + Study + y0:Study + Group), data = mod_df) ## thien phuc
    X <- model.matrix(X, mod_df)
    
    if (inherits(model, "gls")) beta <- coef(model) else beta <- nlme::fixef(model)
    # check for columns dropped from model
    col_names <- names(beta)
    #if (ncol(X) != length(col_names)) X <- X[,col_names,drop=FALSE]  ## thien phuc froze  
    
    Terms <- terms(model)
    yname <- as.character(attr(Terms, "variables"))[[2]]
    y <- nlme::getResponse(model)
    X_y <- cbind(X, y)
    V_list <- build_Sigma_mats(model, sigma_scale =FALSE)
    ord_X_y <- get_order(A=V_list, B=X_y)
    X_y <- ord_X_y[[1]]
    mod_df <- mod_df[ord_X_y[[2]], ]
    xbeta <- as.vector(X_y[, col_names]%*%beta)
    errors <- X_y[, "y"] - xbeta
    V_matrix <- Matrix::bdiag(V_list)
    Ut <- t(chol(V_matrix))
    wt <- solve(Ut)
    
    #Weighting the residuals.
    wterrors <- wt%*%errors    
    
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      #as.data.frame(perm_y <- xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex]))
      as.data.frame(perm_y <- xbeta+as.vector(Ut%*%wterrors[rowindex])) ## long hao
    })
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        mod_df[, yname] <- x
        rfitmodel <- try(update(model, data=mod_df), TRUE)
        if (class(rfitmodel)[1]==class(model)[1]) {
          # mp <- mymodelparm(rfitmodel)[1:2]		
          # if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          # coefT <- coef(summary(rfitmodel))
          
          #BLUP <- coef(rfitmodel)$Group ## thien phuc 
          blup.lmm <- function(model, PI=TRUE, pi.type="all", tau2=NULL, sigma2=NULL){   ## blup for random effects
            df <- model$data
            ni <- table(model$groups)
            n <- length(ni)
            
            y <- df$y
            
            xz <- reStruct(list(y = ~ 0 + Study + y0:Study + Group:Study), data = df) 
            xz <- model.matrix(xz, df)
            x <- xz[,-((3*n-3):(3*n))]
            z <- xz[,((3*n-3):(3*n))]
            
            if(is.null(tau2)==TRUE){
              tau2 <- VarCorr(model)
              tau2 <- as.numeric(tau2["Group", "Variance"])
              cluster <- getGroups(model) ## to which cluster each observation belongs
              sigma.y <- getVarCov(model, type = "marginal", individuals = levels(cluster)) ## response variance-covariance matrix
              sigma.y <- matrix(Matrix::bdiag(sigma.y), nrow=nrow(df))
              #sigma.y <- diag(diag(sigma.y))
            }else{
              if(is.null(sigma2)==TRUE){
                sigma <- model$sigma/varWeights(model$modelStruct$varStruct)
                sigma <- sigma[as.character(1:n)]
                sigma2 <- sigma^2
              }
              sigma.y <- tau2*z%*%t(z)+diag(rep(sigma2, ni))
            }    
            
            xgrp <- cbind(x, df$Group)
            
            sigma.beta.theta <- solve(t(xgrp)%*%solve(sigma.y)%*%xgrp)
            # sigma.beta.theta <- model$varFix
            # sigma.beta.theta <- rbind(cbind(sigma.beta.theta[-(n+1),-(n+1)], sigma.beta.theta[-(n+1),(n+1)]), sigma.beta.theta[(n+1),])
            mu.beta.theta <- sigma.beta.theta%*%t(xgrp)%*%solve(sigma.y)%*%y
            
            M <- cbind(-tau2*t(z)%*%solve(sigma.y)%*%x, (diag(n)-tau2*t(z)%*%solve(sigma.y)%*%z)%*%rep(1,n))
            BLUP <- M%*%mu.beta.theta+tau2*t(z)%*%solve(sigma.y)%*%y
            
            sigma.blup <- tau2*diag(1,n)-tau2^2*t(z)%*%solve(sigma.y)%*%z+M%*%sigma.beta.theta%*%t(M)
            if(PI==FALSE){
              return(split(data.frame(blup=BLUP,
                                      SD=sqrt(diag(sigma.blup))),seq(ni)) )
            }else{
              if(pi.type=="norm"){
                return(data.frame(blup=BLUP,
                                  low=BLUP-1.96*sqrt(diag(sigma.blup)),
                                  hi=BLUP+1.96*sqrt(diag(sigma.blup))))
              }
              if(pi.type=="kh"){
                return(data.frame(blup=BLUP,
                                  low=BLUP-qt(0.975,n-1)*sqrt(diag(sigma.blup)),
                                  hi=BLUP+qt(0.975,n-1)*sqrt(diag(sigma.blup))))
              }
              if(pi.type=="all"){
                return(data.frame(blup=BLUP,
                                  low.norm=BLUP-1.96*sqrt(diag(sigma.blup)),
                                  hi.norm=BLUP+1.96*sqrt(diag(sigma.blup)),
                                  low.kh=BLUP-qt(0.975,n-1)*sqrt(diag(sigma.blup)),
                                  hi.kh=BLUP+qt(0.975,n-1)*sqrt(diag(sigma.blup))))
              }
            }
          }
          BLUP.perm <- blup.lmm(rfitmodel, PI=FALSE, pi.type="all", tau2=NULL, sigma2=NULL)
          
          #return(list(mp, aT, coefT, BLUP.perm))
          return(BLUP.perm)          
        }
      })
      
    }else{
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        if (!(type %in% c("I", "1"))) clusterEvalQ(cl, library(car)) 
        clusterEvalQ(cl, library(nlme))
        clusterExport(cl, c("mymodelparm", "mymodelparm.lme", "mymodelparm.gls","mymodelparm.default", "model", "mod_df", "yname", "type", "test.statistic"), envir = environment()) 
        permmod <- parLapplyLB(cl, permy, function(x) {
          mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            # mp <- mymodelparm(rfitmodel)[1:2]		
            # if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            # coefT <- coef(summary(rfitmodel))
            
            #BLUP <- coef(rfitmodel)$Group ## thien phuc 
            blup.lmm <- function(model, PI=TRUE, pi.type="all", tau2=NULL, sigma2=NULL){   ## blup for random effects
              df <- model$data
              ni <- table(model$groups)
              n <- length(ni)
              
              y <- df$y
              
              xz <- reStruct(list(y = ~ 0 + Study + y0:Study + Group:Study), data = df) 
              xz <- model.matrix(xz, df)
              x <- xz[,-((3*n-3):(3*n))]
              z <- xz[,((3*n-3):(3*n))]
              
              if(is.null(tau2)==TRUE){
                tau2 <- VarCorr(model)
                tau2 <- as.numeric(tau2["Group", "Variance"])
                cluster <- getGroups(model) ## to which cluster each observation belongs
                sigma.y <- getVarCov(model, type = "marginal", individuals = levels(cluster)) ## response variance-covariance matrix
                sigma.y <- matrix(Matrix::bdiag(sigma.y), nrow=nrow(df))
                #sigma.y <- diag(diag(sigma.y))
              }else{
                if(is.null(sigma2)==TRUE){
                  sigma <- model$sigma/varWeights(model$modelStruct$varStruct)
                  sigma <- sigma[as.character(1:n)]
                  sigma2 <- sigma^2
                }
                sigma.y <- tau2*z%*%t(z)+diag(rep(sigma2, ni))
              }    
              
              xgrp <- cbind(x, df$Group)
              
              sigma.beta.theta <- solve(t(xgrp)%*%solve(sigma.y)%*%xgrp)
              # sigma.beta.theta <- model$varFix
              # sigma.beta.theta <- rbind(cbind(sigma.beta.theta[-(n+1),-(n+1)], sigma.beta.theta[-(n+1),(n+1)]), sigma.beta.theta[(n+1),])
              mu.beta.theta <- sigma.beta.theta%*%t(xgrp)%*%solve(sigma.y)%*%y
              
              M <- cbind(-tau2*t(z)%*%solve(sigma.y)%*%x, (diag(n)-tau2*t(z)%*%solve(sigma.y)%*%z)%*%rep(1,n))
              BLUP <- M%*%mu.beta.theta+tau2*t(z)%*%solve(sigma.y)%*%y
              
              sigma.blup <- tau2*diag(1,n)-tau2^2*t(z)%*%solve(sigma.y)%*%z+M%*%sigma.beta.theta%*%t(M)
              if(PI==FALSE){
                return(split(data.frame(blup=BLUP,
                                        SD=sqrt(diag(sigma.blup))),seq(ni)) )
              }else{
                if(pi.type=="norm"){
                  return(data.frame(blup=BLUP,
                                    low=BLUP-1.96*sqrt(diag(sigma.blup)),
                                    hi=BLUP+1.96*sqrt(diag(sigma.blup))))
                }
                if(pi.type=="kh"){
                  return(data.frame(blup=BLUP,
                                    low=BLUP-qt(0.975,n-1)*sqrt(diag(sigma.blup)),
                                    hi=BLUP+qt(0.975,n-1)*sqrt(diag(sigma.blup))))
                }
                if(pi.type=="all"){
                  return(data.frame(blup=BLUP,
                                    low.norm=BLUP-1.96*sqrt(diag(sigma.blup)),
                                    hi.norm=BLUP+1.96*sqrt(diag(sigma.blup)),
                                    low.kh=BLUP-qt(0.975,n-1)*sqrt(diag(sigma.blup)),
                                    hi.kh=BLUP+qt(0.975,n-1)*sqrt(diag(sigma.blup))))
                }
              }
            }
            BLUP.perm <- blup.lmm(rfitmodel, PI=FALSE, pi.type="all", tau2=NULL, sigma2=NULL)
            
            #return(list(mp, aT, coefT, BLUP.perm))     
            return(BLUP.perm)          
          }
        })
        stopCluster(cl)
      }else{  
        permmod <- mclapply(permy, function(x) {
          mod_df[, yname] <- x
          rfitmodel <- try(update(model, data=mod_df), TRUE)
          if (class(rfitmodel)[1]==class(model)[1]) {
            # mp <- mymodelparm(rfitmodel)[1:2]		
            # if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            # coefT <- coef(summary(rfitmodel))
            
            #BLUP <- coef(rfitmodel)$Group ## thien phuc 
            blup.lmm <- function(model, PI=TRUE, pi.type="all", tau2=NULL, sigma2=NULL){   ## blup for random effects
              df <- model$data
              ni <- table(model$groups)
              n <- length(ni)
              
              y <- df$y
              
              xz <- reStruct(list(y = ~ 0 + Study + y0:Study + Group:Study), data = df) 
              xz <- model.matrix(xz, df)
              x <- xz[,-((3*n-3):(3*n))]
              z <- xz[,((3*n-3):(3*n))]
              
              if(is.null(tau2)==TRUE){
                tau2 <- VarCorr(model)
                tau2 <- as.numeric(tau2["Group", "Variance"])
                cluster <- getGroups(model) ## to which cluster each observation belongs
                sigma.y <- getVarCov(model, type = "marginal", individuals = levels(cluster)) ## response variance-covariance matrix
                sigma.y <- matrix(Matrix::bdiag(sigma.y), nrow=nrow(df))
                #sigma.y <- diag(diag(sigma.y))
              }else{
                if(is.null(sigma2)==TRUE){
                  sigma <- model$sigma/varWeights(model$modelStruct$varStruct)
                  sigma <- sigma[as.character(1:n)]
                  sigma2 <- sigma^2
                }
                sigma.y <- tau2*z%*%t(z)+diag(rep(sigma2, ni))
              }    
              
              xgrp <- cbind(x, df$Group)
              
              sigma.beta.theta <- solve(t(xgrp)%*%solve(sigma.y)%*%xgrp)
              # sigma.beta.theta <- model$varFix
              # sigma.beta.theta <- rbind(cbind(sigma.beta.theta[-(n+1),-(n+1)], sigma.beta.theta[-(n+1),(n+1)]), sigma.beta.theta[(n+1),])
              mu.beta.theta <- sigma.beta.theta%*%t(xgrp)%*%solve(sigma.y)%*%y
              
              M <- cbind(-tau2*t(z)%*%solve(sigma.y)%*%x, (diag(n)-tau2*t(z)%*%solve(sigma.y)%*%z)%*%rep(1,n))
              BLUP <- M%*%mu.beta.theta+tau2*t(z)%*%solve(sigma.y)%*%y
              
              sigma.blup <- tau2*diag(1,n)-tau2^2*t(z)%*%solve(sigma.y)%*%z+M%*%sigma.beta.theta%*%t(M)
              if(PI==FALSE){
                return(split(data.frame(blup=BLUP,
                                        SD=sqrt(diag(sigma.blup))),seq(ni)) )
              }else{
                if(pi.type=="norm"){
                  return(data.frame(blup=BLUP,
                                    low=BLUP-1.96*sqrt(diag(sigma.blup)),
                                    hi=BLUP+1.96*sqrt(diag(sigma.blup))))
                }
                if(pi.type=="kh"){
                  return(data.frame(blup=BLUP,
                                    low=BLUP-qt(0.975,n-1)*sqrt(diag(sigma.blup)),
                                    hi=BLUP+qt(0.975,n-1)*sqrt(diag(sigma.blup))))
                }
                if(pi.type=="all"){
                  return(data.frame(blup=BLUP,
                                    low.norm=BLUP-1.96*sqrt(diag(sigma.blup)),
                                    hi.norm=BLUP+1.96*sqrt(diag(sigma.blup)),
                                    low.kh=BLUP-qt(0.975,n-1)*sqrt(diag(sigma.blup)),
                                    hi.kh=BLUP+qt(0.975,n-1)*sqrt(diag(sigma.blup))))
                }
              }
            }
            BLUP.perm <- blup.lmm(rfitmodel, PI=FALSE, pi.type="all", tau2=NULL, sigma2=NULL)
            
            #return(list(mp, aT, coefT, BLUP.perm))
            return(BLUP.perm)          
          }
        }, mc.cores=ncore)
      }
    }
  }
  
  if (any(inherits(model, "lmerMod"), inherits(model, "glmerMod"))) { 
    lmer_var_names <- all.vars(model@call$formula)
    lmer_var_names <- lmer_var_names[lmer_var_names!="pi"]
    mod_df <- as.data.frame(nlme::getData(model))
    mod_df <- na.omit(mod_df[, lmer_var_names])	
    model <- update(model, data=mod_df)	
    theta <- getME(model, "theta")
    fixef <- fixef(model) 
    Lambda <- getME(model, "Lambda")
    Lambdac <- Matrix::tcrossprod(Lambda)
    V <- getME(model, "Z")%*%Lambdac%*%getME(model, "Zt")+diag(dim(getME(model, "Z"))[1])
    Ut <- t(chol(V))
    wt <- solve(Ut)
    xbeta <- as.vector(getME(model, "X")%*%fixef)
    if (inherits(model, "glmerMod")) {
      errors <- slot(model, "resp")$family$linkfun(getME(model, "y")) - xbeta 
    }else errors <- getME(model, "y") - xbeta
    
    #Weighting the residuals.
    wterrors <- wt%*%errors
    
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      if (inherits(model, "glmerMod")) {
        #perm_y <- slot(model, "resp")$family$linkinv(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) 
        perm_y <- slot(model, "resp")$family$linkinv(xbeta+as.vector(Ut%*%wterrors[rowindex])) ## longhao
        if (slot(model, "resp")$family$family == "poisson") { 
          perm_y[perm_y<=0] <- 0
          perm_y <- round(perm_y, 0)
        }
        if (slot(model, "resp")$family$family == "binomial") {
          perm_y[perm_y<=0] <- 0
          perm_y[perm_y>=1] <- 1		
        }
        if (slot(model, "resp")$family$family == "gamma") {
          perm_y[perm_y<=0] <- 0		
        }
        as.data.frame(perm_y)
      }else{
        #as.data.frame(perm_y <- xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex]))
        as.data.frame(perm_y <- xbeta+as.vector(Ut%*%wterrors[rowindex]))  ## long hao
      } 
    })
    
    if (!is.null(fo)) {
      permmod <- lapply(permy, function(x) {
        rfitmodel <- try(refit(model, x), TRUE)
        if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
          mp <- mymodelparm(rfitmodel)[1:2]		
          if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
          coefT <- coef(summary(rfitmodel))
          return(list(mp, aT, coefT))
        }
      })
      
    }else{
      if (.Platform$OS.type=="windows") {
        cl <- makeCluster(ncore)	  
        clusterEvalQ(cl, library(lme4))
        clusterExport(cl, c("mymodelparm", "mymodelparm.lmerMod", "mymodelparm.default", "model", "type", "test.statistic"), envir = environment())       
        permmod <- parLapplyLB(cl, permy, function(x) {
          rfitmodel <- try(refit(model, x), TRUE)
          if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")){
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        })
        stopCluster(cl)
      }else{  
        permmod <- mclapply(permy, function(x) {
          rfitmodel <- try(refit(model, x), TRUE)
          if (class(rfitmodel)[1] %in% c("lmerMod", "glmerMod")) {
            mp <- mymodelparm(rfitmodel)[1:2]		
            if (type %in% c("I", "1")) aT <- anova(rfitmodel) else aT <- car::Anova(rfitmodel, type=type, test.statistic=test.statistic)
            coefT <- coef(summary(rfitmodel))
            return(list(mp, aT, coefT))
          }
        }, mc.cores=ncore)
      }
    }	
  }    
  permmod <- permmod[!(sapply(permmod,is.null))]
  # nperm <- length(permmod)
  # if (nperm==0) stop("permutation can't produce useful data!")
  # if (!(inherits(model, "aov"))) {
  #   # tTable <- function(x) { # function to construct t-table for the model
  #   # if (class(model)[1]=="lmerMod"){
  #   # mp <- mymodelparm(x)
  #   # tTable <- cbind(mp$coef, sqrt(base::diag(mp$vcov)), mp$coef/sqrt(base::diag(mp$vcov)))
  #   # colnames(tTable) <- c("Estimate", "Std. Error", "t value")
  #   # return(round(tTable, 4))		
  #   # }else{
  #   # summ <- summary(x)
  #   # if (class(x)[1] %in% c("lme", "gls")) {
  #   # tTable <- summ$tTable
  #   # }else{
  #   # tTable <- coef(summ)
  #   # }
  #   # return(tTable)
  #   # }
  #   # }
  #   
  #   # model_tTable <- tTable(model)
  #   # Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
  #   
  #   # if (nrow(model_tTable)==1) {
  #   # Tper.p <- (sum(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(model_tTable[, Tvalue]),6))+!exact)/(nperm+!exact)
  #   # }else{
  #   # Tper.p <- (rowSums(round(sapply(permmod, function(x) {mp <- x[[1]]; abs(mp$coef/sqrt(base::diag(mp$vcov)))}), 6) >= round(abs(model_tTable[, Tvalue]), 6))+!exact)/(nperm+!exact)
  #   # }
  #   
  #   model_tTable <- coef(summary(model))
  #   Tvalue <- colnames(model_tTable)[grep("value", colnames(model_tTable))][1] 
  #   if (nrow(model_tTable)==1) {
  #     Tper.p <- (sum(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]),6))+!exact)/(nperm+!exact)
  #   }else{
  #     Tper.p <- (rowSums(round(sapply(permmod, function(x) abs(x[[3]][, Tvalue])), 6) >= round(abs(model_tTable[, Tvalue]), 6))+!exact)/(nperm+!exact)
  #   }
  #   
  #   if ("(Intercept)" %in% rownames(model_tTable)) Tper.p[1] <- NA
  #   COEFFICENTS <- round(cbind(model_tTable, "Perm_p_value"=Tper.p),4)
  #   if (prt) {
  #     cat("\nCoefficients of (fixed) effects:\n")
  #     print(COEFFICENTS)
  #     cat("\nNote: Perm_p_value of t test is obtained using", sQuote(nperm), "permutations.\n")
  #   }
  # }else COEFFICENTS <- NULL
  # 
  # if (type %in% c("I", "1")) {
  #   key_anova <- anova(model)
  #   Fvalue <- colnames(key_anova)[grep("F.value",colnames(key_anova))]
  #   if (class(model)[1] %in% c("glm")) Fvalue <- "Deviance"
  #   if (nrow(key_anova)==1) {
  #     Fper.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
  #   }else{
  #     Fper.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, Fvalue]}), 6) >= round(key_anova[, Fvalue], 6))+!exact)/(nperm+!exact)
  #   } 
  #   if ("(Intercept)" %in% rownames(key_anova)) Fper.p[1] <- NA	
  #   ANOVA <- round(cbind(key_anova, "Perm_p_value"=Fper.p),4)
  # }else{
  #   key_anova <- car::Anova(model, type=type, test.statistic=test.statistic)
  #   if (class(model)[1] %in% c("lm", "aov")) stats_value <- "F value"  else stats_value <-  test.statistic 
  #   if (class(model)[1] %in% c("glm")) stats_value <- switch(test.statistic, LR = "LR Chisq", Wald = "Chisq", F = ifelse(type %in% c("2", "II"), "F value", "F values"))  
  #   if (nrow(key_anova)==1) {
  #     stats.p <- (sum(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
  #   }else{
  #     stats.p <- (rowSums(round(sapply(permmod, function(x) {aT <- x[[2]]; aT[, stats_value]}), 6) >= round(key_anova[, stats_value], 6))+!exact)/(nperm+!exact)
  #   }  
  #   if ("(Intercept)" %in% rownames(key_anova)) stats.p[1] <- NA
  #   ANOVA <- round(cbind(key_anova, "Perm_p_value"=stats.p),4) 
  #   attr(ANOVA, "heading") <- attr(key_anova, "heading")
  # }
  # 
  # ANOVA[is.na(ANOVA)] <-""
  # if (prt) {
  #   cat("\nANOVA:\n")
  #   if (!is.null(attr(ANOVA, "heading"))) cat(paste(attr(ANOVA, "heading"), "\n", sep=""))
  #   print(ANOVA)
  #   cat(paste("\nNote: Perm_p_value of", test.statistic, "test is obtained using"), sQuote(nperm), "permutations.\n\n")
  # }
  # permlist <- blup.lst <- vector("list", nperm)
  # for (i in 1:nperm) {
  #   #permlist[[i]] <- permmod[[i]][[1]]
  #   blup.lst[[i]] <- permmod[[i]][[4]]
  # }
  # blup.df <- as.data.frame(do.call(rbind, blup.lst))
  # pi.perm <- as.data.frame(t(apply(blup.df, 2, FUN=function(x) c(quantile(x, 0.025), quantile(x, 0.975)))))
  # colnames(pi.perm) <- c("low", "hi")
  #return(invisible(list("permlist"=permlist, "COEFFICENTS"=COEFFICENTS, "ANOVA"=ANOVA)))
  blup.lmm <- function(model, PI=TRUE, pi.type="all", tau2=NULL, sigma2=NULL){   ## blup for random effects
    df <- model$data
    ni <- table(model$groups)
    n <- length(ni)
    
    y <- df$y
    
    xz <- reStruct(list(y = ~ 0 + Study + y0:Study + Group:Study), data = df) 
    xz <- model.matrix(xz, df)
    x <- xz[,-((3*n-3):(3*n))]
    z <- xz[,((3*n-3):(3*n))]
    
    if(is.null(tau2)==TRUE){
      tau2 <- VarCorr(model)
      tau2 <- as.numeric(tau2["Group", "Variance"])
      cluster <- getGroups(model) ## to which cluster each observation belongs
      sigma.y <- getVarCov(model, type = "marginal", individuals = levels(cluster)) ## response variance-covariance matrix
      sigma.y <- matrix(Matrix::bdiag(sigma.y), nrow=nrow(df))
      #sigma.y <- diag(diag(sigma.y))
    }else{
      if(is.null(sigma2)==TRUE){
        sigma <- model$sigma/varWeights(model$modelStruct$varStruct)
        sigma <- sigma[as.character(1:n)]
        sigma2 <- sigma^2
      }
      sigma.y <- tau2*z%*%t(z)+diag(rep(sigma2, ni))
    }    
    
    xgrp <- cbind(x, df$Group)
    
    sigma.beta.theta <- solve(t(xgrp)%*%solve(sigma.y)%*%xgrp)
    # sigma.beta.theta <- model$varFix
    # sigma.beta.theta <- rbind(cbind(sigma.beta.theta[-(n+1),-(n+1)], sigma.beta.theta[-(n+1),(n+1)]), sigma.beta.theta[(n+1),])
    mu.beta.theta <- sigma.beta.theta%*%t(xgrp)%*%solve(sigma.y)%*%y
    
    M <- cbind(-tau2*t(z)%*%solve(sigma.y)%*%x, (diag(n)-tau2*t(z)%*%solve(sigma.y)%*%z)%*%rep(1,n))
    BLUP <- M%*%mu.beta.theta+tau2*t(z)%*%solve(sigma.y)%*%y
    
    sigma.blup <- tau2*diag(1,n)-tau2^2*t(z)%*%solve(sigma.y)%*%z+M%*%sigma.beta.theta%*%t(M)
    if(PI==FALSE){
      return(split(data.frame(blup=BLUP,
                              SD=sqrt(diag(sigma.blup))),seq(ni)) )
    }else{
      if(pi.type=="norm"){
        return(data.frame(blup=BLUP,
                          low=BLUP-1.96*sqrt(diag(sigma.blup)),
                          hi=BLUP+1.96*sqrt(diag(sigma.blup))))
      }
      if(pi.type=="kh"){
        return(data.frame(blup=BLUP,
                          low=BLUP-qt(0.975,n-1)*sqrt(diag(sigma.blup)),
                          hi=BLUP+qt(0.975,n-1)*sqrt(diag(sigma.blup))))
      }
      if(pi.type=="all"){
        return(data.frame(blup=BLUP,
                          low.norm=BLUP-1.96*sqrt(diag(sigma.blup)),
                          hi.norm=BLUP+1.96*sqrt(diag(sigma.blup)),
                          low.kh=BLUP-qt(0.975,n-1)*sqrt(diag(sigma.blup)),
                          hi.kh=BLUP+qt(0.975,n-1)*sqrt(diag(sigma.blup))))
      }
    }
  }
  BLUP <- blup.lmm(model, PI=FALSE, pi.type="all", tau2=NULL, sigma2=NULL)
  return(invisible(list("BLUP"=BLUP, "BLUP.perm.list"=permmod)))
} 

## permlmer thien phuc long hao
permlmer.tplh <- function(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed){
  
  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) stop("The model must be a lmer object!")
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) stop("Please check the response in your model!")
  
  c_deparse <- function (...) 
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")
  
  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1) 
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1) 
  
  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }  
  
  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(setequal(fixef0name, fixef1name))  ref <- TRUE else ref <- FALSE
  
  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  thetan[theta0name] <- theta0
  # Lambda1n <- getME(lmer1, "Lambda") 
  # Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  # Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  #V1n <- getME(lmer1, "Z")%*%Lambda1nc%*%getME(lmer1, "Zt")+diag(dim(getME(lmer1, "Z"))[1])
  Lambda1n <- getME(lmer0, "Lambda") ## thien phuc
  Lambda1n@x <- thetan[getME(lmer0, "Lind")]
  Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  #V1n <- getME(lmer0, "Z")%*%Lambda1nc%*%getME(lmer0, "Zt")+diag(dim(getME(lmer0, "Z"))[1])
  V1n <- unique(diag(Lambda1nc))*getME(lmer0, "Z")%*%getME(lmer0, "Zt")+getME(lmer0, "sigma")^2*diag(dim(getME(lmer0, "Z"))[1])
  Ut <- t(chol(V1n))
  wt <- solve(Ut)
  xbeta <- as.vector(getME(lmer0, "X")%*%fixef0)
  errors <- getME(lmer0, "y") - xbeta
  
  #Weighting the residuals.
  wterrors <- wt%*%errors
  
  # permute weighted resid, then unweighted it for 999 times
  # permResid <- matrix(0, length(wterrors), nperm)
  # for (i in 1:nperm) {
  # if(!missing(seed)) set.seed(seed+i)
  # permResid[, i] <- as.vector(Ut%*%sample(wterrors)) 
  # }
  # permy <- as.data.frame(xbeta+permResid)
  if (ref) {
    if(!missing(seed)) set.seed(seed)
    permy <- as.data.frame(xbeta+replicate(nperm, as.vector(Ut%*%sample(wterrors))))
  }else{   
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      #as.data.frame(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## thien phuc
      as.data.frame(xbeta+as.vector(Ut%*%wterrors[rowindex]))
    })
  }
  
  # Calculating the likelihood ratio test statistic for each permutation.    
  lrtest1 <- 2*(logLik(lmer1, REML=ref)-logLik(lmer0, REML=ref))
  lrtest1 <- ifelse(lrtest1 < 0, 0, lrtest1)
  
  # ttest1 <- summary(lmer1)$coefficients["Group", "t value"]
  # betahattest1 <- summary(lmer1)$coefficients["Group", "Estimate"]
  
  if (.Platform$OS.type=="windows") {
    cl <- makeCluster(ncore)	  
    clusterEvalQ(cl, library(lme4))
    clusterExport(cl, c("lmer0", "lmer1", "ref"), envir = environment()) 
    lrtest2 <- parLapplyLB(cl, permy, function(x) {
      LRT <- try(2*(logLik(refit(lmer1, x), REML=ref) - logLik(refit(lmer0, x), REML=ref)), TRUE)
      LRT <- ifelse(is.numeric(LRT), LRT, NA)
    })
    # ttest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
    #   TT <- try(summary(refit(lmer1, x))$coefficients["Group", "t value"], TRUE)
    #   TT <- ifelse(is.numeric(TT), TT, NA)
    # })
    # betahattest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
    #   BHT <- try(summary(refit(lmer1, x))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # })
    stopCluster(cl)
  }else{  
    lrtest2 <- mclapply(permy, function(x) {
      LRT <- try(2*(logLik(suppressMessages(refit(lmer1, x)), REML=ref) - logLik(suppressMessages(refit(lmer0, x)), REML=ref)), TRUE)
      LRT <- ifelse(is.numeric(LRT), LRT, NA)
    }, mc.cores=ncore)
    # ttest2 <- mclapply(permy, function(x) {
    #   TT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "t value"], TRUE)
    #   TT <- ifelse(is.numeric(TT), TT, NA)
    # }, mc.cores=ncore)
    # betahattest2 <- mclapply(permy, function(x) {
    #   BHT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # }, mc.cores=ncore)
    
  }
  
  #Calculating the p-values.  
  lrtest <- na.omit(unlist(lrtest2))
  lrtest <- ifelse(lrtest < 0, 0, lrtest)
  perm_p <- (sum(lrtest >= lrtest1) +1)/(length(lrtest) + 1)
  aod <- anova(lmer0, lmer1, refit=!ref)
  aod$'Perm-p' <- c(NA, perm_p)
  
  # ttest <- na.omit(unlist(ttest2))
  # ttest <- ifelse(ttest < 0, 0, ttest)
  # perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
  # aod$'Perm-p t test' <- c(NA, perm_p.t)
  # 
  # betahattest <- na.omit(unlist(betahattest2))
  # betahattest <- ifelse(betahattest < 0, 0, betahattest)
  # perm_p.betahat <- (sum(betahattest >= betahattest1) +1)/(length(betahattest) + 1)
  # aod$'Perm-p beta hat test' <- c(NA, perm_p.betahat)
  
  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}

permlmer.tplh2 <- function(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed, pvalue=FALSE, perc=FALSE){
  
  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) stop("The model must be a lmer object!")
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) stop("Please check the response in your model!")
  
  c_deparse <- function (...) 
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")
  
  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1) 
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1) 
  
  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }  
  
  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(setequal(fixef0name, fixef1name))  ref <- TRUE else ref <- FALSE
  
  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  thetan[theta0name] <- theta0
  # Lambda1n <- getME(lmer1, "Lambda") 
  # Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  # Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  #V1n <- getME(lmer1, "Z")%*%Lambda1nc%*%getME(lmer1, "Zt")+diag(dim(getME(lmer1, "Z"))[1])
  Lambda1n <- getME(lmer0, "Lambda") ## thien phuc
  Lambda1n@x <- thetan[getME(lmer0, "Lind")]
  #Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  Lambda1nc <- base::tcrossprod(Lambda1n)
  # V1n <- getME(lmer0, "Z")%*%Lambda1nc%*%getME(lmer0, "Zt")+diag(dim(getME(lmer0, "Z"))[1])
  V1n <- unique(diag(Lambda1nc))*getME(lmer0, "Z")%*%getME(lmer0, "Zt")+getME(lmer0, "sigma")^2*diag(dim(getME(lmer0, "Z"))[1])
  Ut <- t(chol(V1n))
  wt <- solve(Ut)
  xbeta <- as.vector(getME(lmer0, "X")%*%fixef0)
  errors <- getME(lmer0, "y") - xbeta
  
  #Weighting the residuals.
  wterrors <- wt%*%errors
  
  # permute weighted resid, then unweighted it for 999 times
  # permResid <- matrix(0, length(wterrors), nperm)
  # for (i in 1:nperm) {
  # if(!missing(seed)) set.seed(seed+i)
  # permResid[, i] <- as.vector(Ut%*%sample(wterrors)) 
  # }
  # permy <- as.data.frame(xbeta+permResid)
  if (ref) {
    if(!missing(seed)) set.seed(seed)
    permy <- as.data.frame(xbeta+replicate(nperm, as.vector(Ut%*%sample(wterrors))))
  }else{   
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      #as.data.frame(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## thien phuc
      as.data.frame(xbeta+as.vector(Ut%*%wterrors[rowindex]))
    })
  }
  
  # Calculating the likelihood ratio test statistic for each permutation.    
  ttest1 <- summary(lmer1)$coefficients["Group", "t value"]
  # betahattest1 <- summary(lmer1)$coefficients["Group", "Estimate"]
  
  if (.Platform$OS.type=="windows") {
    cl <- makeCluster(ncore)	  
    clusterEvalQ(cl, library(lme4))
    clusterExport(cl, c("lmer0", "lmer1", "ref"), envir = environment()) 
    ttest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
      TT <- try(summary(refit(lmer1, x))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    })
    # betahattest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
    #   BHT <- try(summary(refit(lmer1, x))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # })
    stopCluster(cl)
  }else{  
    ttest2 <- mclapply(permy, function(x) {
      TT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    }, mc.cores=ncore)
    # betahattest2 <- mclapply(permy, function(x) {
    #   BHT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # }, mc.cores=ncore)
    
  }
  
  #Calculating the p-values.  
  aod<-list(NULL)
  ttest <- na.omit(unlist(ttest2))
  # ttest <- ifelse(ttest < 0, 0, ttest) # thien phuc
  # perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
  # aod$'Perm-p' <- c(NA, perm_p.t)
  if(pvalue==TRUE){
    perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
    perm_p.t <- min(2*perm_p.t, 2-2*perm_p.t)
    aod$'Perm-p' <- c(NA, perm_p.t)
  }
  if(perc==TRUE){
    aod$'97.5-perc' <- quantile(ttest, 0.975)
    aod$'2.5-perc' <- quantile(ttest, 0.025)
  }
  
  # betahattest <- na.omit(unlist(betahattest2))
  # betahattest <- ifelse(betahattest < 0, 0, betahattest)
  # perm_p.betahat <- (sum(betahattest >= betahattest1) +1)/(length(betahattest) + 1)
  # aod$'Perm-p beta hat test' <- c(NA, perm_p.betahat)
  
  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}

permlmer.tplh3 <- function(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed){
  
  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) stop("The model must be a lmer object!")
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) stop("Please check the response in your model!")
  
  c_deparse <- function (...) 
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")
  
  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1) 
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1) 
  
  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }  
  
  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(setequal(fixef0name, fixef1name))  ref <- TRUE else ref <- FALSE
  
  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  thetan[theta0name] <- theta0
  Lambda1n <- getME(lmer1, "Lambda")
  Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  #V1n <- getME(lmer1, "Z")%*%Lambda1nc%*%getME(lmer1, "Zt")+diag(dim(getME(lmer1, "Z"))[1])
  V1n <- unique(diag(Lambda1nc))*getME(lmer1, "Z")%*%getME(lmer1, "Zt")+getME(lmer1, "sigma")^2*diag(dim(getME(lmer1, "Z"))[1])
  
  # Lambda1n <- getME(lmer0, "Lambda") ## thien phuc
  # Lambda1n@x <- thetan[getME(lmer0, "Lind")]
  # Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  # # V1n <- getME(lmer0, "Z")%*%Lambda1nc%*%getME(lmer0, "Zt")+diag(dim(getME(lmer0, "Z"))[1])
  # V1n <- unique(diag(Lambda1nc))*getME(lmer0, "Z")%*%getME(lmer0, "Zt")+getME(lmer0, "sigma")^2*diag(dim(getME(lmer0, "Z"))[1])
  
  Ut <- t(chol(V1n))
  wt <- solve(Ut)
  xbeta <- as.vector(getME(lmer0, "X")%*%fixef0)
  errors <- getME(lmer0, "y") - xbeta
  
  #Weighting the residuals.
  wterrors <- wt%*%errors
  
  # permute weighted resid, then unweighted it for 999 times
  # permResid <- matrix(0, length(wterrors), nperm)
  # for (i in 1:nperm) {
  # if(!missing(seed)) set.seed(seed+i)
  # permResid[, i] <- as.vector(Ut%*%sample(wterrors)) 
  # }
  # permy <- as.data.frame(xbeta+permResid)
  if (ref) {
    if(!missing(seed)) set.seed(seed)
    permy <- as.data.frame(xbeta+replicate(nperm, as.vector(Ut%*%sample(wterrors))))
  }else{   
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      #as.data.frame(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## thien phuc
      as.data.frame(xbeta+as.vector(Ut%*%wterrors[rowindex]))
    })
  }
  
  # Calculating the likelihood ratio test statistic for each permutation.    
  ttest1 <- summary(lmer1)$coefficients["Group", "t value"]
  # betahattest1 <- summary(lmer1)$coefficients["Group", "Estimate"]
  
  if (.Platform$OS.type=="windows") {
    cl <- makeCluster(ncore)	  
    clusterEvalQ(cl, library(lme4))
    clusterExport(cl, c("lmer0", "lmer1", "ref"), envir = environment()) 
    ttest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
      TT <- try(summary(refit(lmer1, x))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    })
    # betahattest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
    #   BHT <- try(summary(refit(lmer1, x))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # })
    stopCluster(cl)
  }else{  
    ttest2 <- mclapply(permy, function(x) {
      TT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    }, mc.cores=ncore)
    # betahattest2 <- mclapply(permy, function(x) {
    #   BHT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # }, mc.cores=ncore)
    
  }
  
  #Calculating the p-values.  
  aod<-list(NULL)
  ttest <- na.omit(unlist(ttest2))
  ttest <- ifelse(ttest < 0, 0, ttest)
  perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
  aod$'Perm-p' <- c(NA, perm_p.t)
  
  # betahattest <- na.omit(unlist(betahattest2))
  # betahattest <- ifelse(betahattest < 0, 0, betahattest)
  # perm_p.betahat <- (sum(betahattest >= betahattest1) +1)/(length(betahattest) + 1)
  # aod$'Perm-p beta hat test' <- c(NA, perm_p.betahat)
  
  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}


permlmer.tplh4 <- function(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed, pvalue=FALSE, perc=FALSE){
  
  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) stop("The model must be a lmer object!")
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) stop("Please check the response in your model!")
  
  c_deparse <- function (...) 
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")
  
  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1) 
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1) 
  
  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }  
  
  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(setequal(fixef0name, fixef1name))  ref <- TRUE else ref <- FALSE
  
  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  #thetan[theta0name] <- theta0 #thien phuc froze
  Lambda1n <- getME(lmer1, "Lambda")
  Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  #V1n <- getME(lmer1, "Z")%*%Lambda1nc%*%getME(lmer1, "Zt")+diag(dim(getME(lmer1, "Z"))[1])
  V1n <- unique(diag(Lambda1nc))*getME(lmer1, "Z")%*%getME(lmer1, "Zt")+getME(lmer1, "sigma")^2*diag(dim(getME(lmer1, "Z"))[1])
  
  # Lambda1n <- getME(lmer0, "Lambda") ## thien phuc
  # Lambda1n@x <- thetan[getME(lmer0, "Lind")]
  # Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  # # V1n <- getME(lmer0, "Z")%*%Lambda1nc%*%getME(lmer0, "Zt")+diag(dim(getME(lmer0, "Z"))[1])
  # V1n <- unique(diag(Lambda1nc))*getME(lmer0, "Z")%*%getME(lmer0, "Zt")+getME(lmer0, "sigma")^2*diag(dim(getME(lmer0, "Z"))[1])
  
  Ut <- t(chol(V1n))
  wt <- solve(Ut)
  #xbeta <- as.vector(getME(lmer0, "X")%*%fixef0)
  xbeta <- as.vector(getME(lmer0, "X")%*%fixef1[-(getME(lmer1, "l_i")+1)])
  #errors <- getME(lmer0, "y") - xbeta
  errors <- getME(lmer1, "y") - xbeta
  
  #Weighting the residuals.
  wterrors <- wt%*%errors
  
  # permute weighted resid, then unweighted it for 999 times
  # permResid <- matrix(0, length(wterrors), nperm)
  # for (i in 1:nperm) {
  # if(!missing(seed)) set.seed(seed+i)
  # permResid[, i] <- as.vector(Ut%*%sample(wterrors)) 
  # }
  # permy <- as.data.frame(xbeta+permResid)
  if (ref) {
    if(!missing(seed)) set.seed(seed)
    permy <- as.data.frame(xbeta+replicate(nperm, as.vector(Ut%*%sample(wterrors))))
  }else{   
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      #as.data.frame(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## thien phuc
      as.data.frame(xbeta+as.vector(Ut%*%wterrors[rowindex]))
    })
  }
  
  # Calculating the likelihood ratio test statistic for each permutation.    
  ttest1 <- summary(lmer1)$coefficients["Group", "t value"]
  # betahattest1 <- summary(lmer1)$coefficients["Group", "Estimate"]
  
  if (.Platform$OS.type=="windows") {
    cl <- makeCluster(ncore)	  
    clusterEvalQ(cl, library(lme4))
    clusterExport(cl, c("lmer0", "lmer1", "ref"), envir = environment()) 
    ttest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
      TT <- try(summary(refit(lmer1, x))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    })
    # betahattest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
    #   BHT <- try(summary(refit(lmer1, x))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # })
    stopCluster(cl)
  }else{  
    ttest2 <- mclapply(permy, function(x) {
      TT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    }, mc.cores=ncore)
    # betahattest2 <- mclapply(permy, function(x) {
    #   BHT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # }, mc.cores=ncore)
    
  }
  
  #Calculating the p-values.  
  aod<-list(NULL)
  ttest <- na.omit(unlist(ttest2))
  # ttest <- ifelse(ttest < 0, 0, ttest) # thien phuc
  # perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
  # aod$'Perm-p' <- c(NA, perm_p.t)
  if(pvalue==TRUE){
    perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
    perm_p.t <- min(2*perm_p.t, 2-2*perm_p.t)
    aod$'Perm-p' <- c(NA, perm_p.t)
  }
  if(perc==TRUE){
    aod$'97.5-perc' <- quantile(ttest, 0.975)
    aod$'2.5-perc' <- quantile(ttest, 0.025)
  }
  
  # betahattest <- na.omit(unlist(betahattest2))
  # betahattest <- ifelse(betahattest < 0, 0, betahattest)
  # perm_p.betahat <- (sum(betahattest >= betahattest1) +1)/(length(betahattest) + 1)
  # aod$'Perm-p beta hat test' <- c(NA, perm_p.betahat)
  
  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}

permlmer.tplh5 <- function(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed, pvalue=FALSE, perc=FALSE){
  
  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) stop("The model must be a lmer object!")
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) stop("Please check the response in your model!")
  
  c_deparse <- function (...) 
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")
  
  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1) 
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1) 
  
  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }  
  
  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(setequal(fixef0name, fixef1name))  ref <- TRUE else ref <- FALSE
  
  # thetan <- rep(0, length(theta1))
  # names(thetan) <- theta1name
  # thetan[theta0name] <- theta0
  # Lambda1n <- getME(lmer1, "Lambda") 
  # Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  # Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  #V1n <- getME(lmer1, "Z")%*%Lambda1nc%*%getME(lmer1, "Zt")+diag(dim(getME(lmer1, "Z"))[1])
  # Lambda1n <- getME(lmer0, "Lambda") ## thien phuc
  # Lambda1n@x <- thetan[getME(lmer0, "Lind")]
  # Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  # V1n <- getME(lmer0, "Z")%*%Lambda1nc%*%getME(lmer0, "Zt")+diag(dim(getME(lmer0, "Z"))[1])
  # V1n <- unique(diag(Lambda1nc))*getME(lmer0, "Z")%*%getME(lmer0, "Zt")+getME(lmer0, "sigma")^2*diag(dim(getME(lmer0, "Z"))[1])
  # Ut <- t(chol(V1n))
  # wt <- solve(Ut)
  xbeta <- as.vector(getME(lmer0, "X")%*%fixef0)
  errors <- getME(lmer0, "y") - xbeta
  
  #Weighting the residuals.
  #wterrors <- wt%*%errors
  
  # permute weighted resid, then unweighted it for 999 times
  # permResid <- matrix(0, length(wterrors), nperm)
  # for (i in 1:nperm) {
  # if(!missing(seed)) set.seed(seed+i)
  # permResid[, i] <- as.vector(Ut%*%sample(wterrors)) 
  # }
  # permy <- as.data.frame(xbeta+permResid)
  if (ref) {
    if(!missing(seed)) set.seed(seed)
    #permy <- as.data.frame(xbeta+replicate(nperm, as.vector(Ut%*%sample(wterrors)))) ## thien phuc froze
    permy <- as.data.frame(xbeta+replicate(nperm, diag(2*rbinom(length(errors), 1, 0.5)-1)%*%errors))
  }else{   
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      # rowindex <- sample(1:length(errors))
      # #as.data.frame(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## thien phuc froze
      # as.data.frame(xbeta+as.vector(Ut%*%wterrors[rowindex]))
      as.data.frame(xbeta+diag(2*rbinom(length(errors), 1, 0.5)-1)%*%errors)
    })
  }
  
  # Calculating the likelihood ratio test statistic for each permutation.    
  ttest1 <- summary(lmer1)$coefficients["Group", "t value"]
  # betahattest1 <- summary(lmer1)$coefficients["Group", "Estimate"]
  
  if (.Platform$OS.type=="windows") {
    cl <- makeCluster(ncore)	  
    clusterEvalQ(cl, library(lme4))
    clusterExport(cl, c("lmer0", "lmer1", "ref"), envir = environment()) 
    ttest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
      TT <- try(summary(refit(lmer1, x))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    })
    # betahattest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
    #   BHT <- try(summary(refit(lmer1, x))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # })
    stopCluster(cl)
  }else{  
    ttest2 <- mclapply(permy, function(x) {
      TT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    }, mc.cores=ncore)
    # betahattest2 <- mclapply(permy, function(x) {
    #   BHT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # }, mc.cores=ncore)
    
  }
  
  #Calculating the p-values.  
  aod<-list(NULL)
  ttest <- na.omit(unlist(ttest2))
  # ttest <- ifelse(ttest < 0, 0, ttest) # thien phuc
  # perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
  # aod$'Perm-p' <- c(NA, perm_p.t)
  if(pvalue==TRUE){
    perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
    perm_p.t <- min(2*perm_p.t, 2-2*perm_p.t)
    aod$'Perm-p' <- c(NA, perm_p.t)
  }
  if(perc==TRUE){
    aod$'97.5-perc' <- quantile(ttest, 0.975)
    aod$'2.5-perc' <- quantile(ttest, 0.025)
  }
  
  # betahattest <- na.omit(unlist(betahattest2))
  # betahattest <- ifelse(betahattest < 0, 0, betahattest)
  # perm_p.betahat <- (sum(betahattest >= betahattest1) +1)/(length(betahattest) + 1)
  # aod$'Perm-p beta hat test' <- c(NA, perm_p.betahat)
  
  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}

permlmer.tplh6 <- function(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed, pvalue=FALSE, perc=FALSE){
  
  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) stop("The model must be a lmer object!")
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) stop("Please check the response in your model!")
  
  c_deparse <- function (...) 
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")
  
  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1) 
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1) 
  
  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }  
  
  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(setequal(fixef0name, fixef1name))  ref <- TRUE else ref <- FALSE
  
  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  thetan[theta0name] <- theta0
  Lambda1n <- getME(lmer1, "Lambda") 
  Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  # Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  #V1n <- getME(lmer1, "Z")%*%Lambda1nc%*%getME(lmer1, "Zt")+diag(dim(getME(lmer1, "Z"))[1])
  Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  # V1n <- getME(lmer0, "Z")%*%Lambda1nc%*%getME(lmer0, "Zt")+diag(dim(getME(lmer0, "Z"))[1])
  V1n <- unique(diag(Lambda1nc))*getME(lmer1, "Z")%*%getME(lmer1, "Zt")+getME(lmer1, "sigma")^2*diag(dim(getME(lmer1, "Z"))[1])
  Ut <- t(chol(V1n))
  wt <- solve(Ut)
  xbeta <- as.vector(getME(lmer1, "X")%*%fixef1)
  errors <- getME(lmer1, "y") - xbeta
  
  #Weighting the residuals.
  wterrors <- wt%*%errors
  
  # permute weighted resid, then unweighted it for 999 times
  # permResid <- matrix(0, length(wterrors), nperm)
  # for (i in 1:nperm) {
  # if(!missing(seed)) set.seed(seed+i)
  # permResid[, i] <- as.vector(Ut%*%sample(wterrors)) 
  # }
  # permy <- as.data.frame(xbeta+permResid)
  if (ref) {
    if(!missing(seed)) set.seed(seed)
    permy <- as.data.frame(xbeta+replicate(nperm, as.vector(Ut%*%sample(wterrors))))
  }else{   
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      #as.data.frame(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## thien phuc
      as.data.frame(xbeta+as.vector(Ut%*%wterrors[rowindex]))
    })
  }
  
  # Calculating the likelihood ratio test statistic for each permutation.    
  ttest1 <- summary(lmer1)$coefficients["Group", "t value"]
  # betahattest1 <- summary(lmer1)$coefficients["Group", "Estimate"]
  
  if (.Platform$OS.type=="windows") {
    cl <- makeCluster(ncore)	  
    clusterEvalQ(cl, library(lme4))
    clusterExport(cl, c("lmer0", "lmer1", "ref"), envir = environment()) 
    ttest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
      TT <- try(summary(refit(lmer1, x))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    })
    # betahattest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
    #   BHT <- try(summary(refit(lmer1, x))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # })
    stopCluster(cl)
  }else{  
    ttest2 <- mclapply(permy, function(x) {
      TT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    }, mc.cores=ncore)
    # betahattest2 <- mclapply(permy, function(x) {
    #   BHT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # }, mc.cores=ncore)
    
  }
  
  #Calculating the p-values.  
  aod<-list(NULL)
  ttest <- na.omit(unlist(ttest2))
  # ttest <- ifelse(ttest < 0, 0, ttest) # thien phuc
  # perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
  # aod$'Perm-p' <- c(NA, perm_p.t)
  if(pvalue==TRUE){
    perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
    perm_p.t <- min(2*perm_p.t, 2-2*perm_p.t)
    aod$'Perm-p' <- c(NA, perm_p.t)
  }
  if(perc==TRUE){
    aod$'97.5-perc' <- quantile(ttest, 0.975)
    aod$'2.5-perc' <- quantile(ttest, 0.025)
  }
  
  # betahattest <- na.omit(unlist(betahattest2))
  # betahattest <- ifelse(betahattest < 0, 0, betahattest)
  # perm_p.betahat <- (sum(betahattest >= betahattest1) +1)/(length(betahattest) + 1)
  # aod$'Perm-p beta hat test' <- c(NA, perm_p.betahat)
  
  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}

permlmer.tplh7 <- function(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed, pvalue=FALSE, upper=FALSE, lower=FALSE, perc=FALSE, theta0, data){
  
  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) stop("The model must be a lmer object!")
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) stop("Please check the response in your model!")
  
  c_deparse <- function (...) 
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")
  
  df <- data ## thien phuc
  df$y<-data$y-theta0*data$Group 
  lmer0<-update(lmer0, data=df)
  lmer1<-update(lmer1, data=df)
  
  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1) 
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1) 
  
  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }  
  
  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(setequal(fixef0name, fixef1name))  ref <- TRUE else ref <- FALSE
  
  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  thetan[theta0name] <- theta0
  # Lambda1n <- getME(lmer1, "Lambda") 
  # Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  # Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  #V1n <- getME(lmer1, "Z")%*%Lambda1nc%*%getME(lmer1, "Zt")+diag(dim(getME(lmer1, "Z"))[1])
  Lambda1n <- getME(lmer0, "Lambda") ## thien phuc
  Lambda1n@x <- thetan[getME(lmer0, "Lind")]
  #Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  Lambda1nc <- base::tcrossprod(Lambda1n)
  # V1n <- getME(lmer0, "Z")%*%Lambda1nc%*%getME(lmer0, "Zt")+diag(dim(getME(lmer0, "Z"))[1])
  V1n <- unique(diag(Lambda1nc))*getME(lmer0, "Z")%*%getME(lmer0, "Zt")+getME(lmer0, "sigma")^2*diag(dim(getME(lmer0, "Z"))[1])
  Ut <- t(chol(V1n))
  wt <- solve(Ut)
  xbeta <- as.vector(getME(lmer0, "X")%*%fixef0)
  errors <- getME(lmer0, "y") - xbeta
  
  #Weighting the residuals.
  wterrors <- wt%*%errors
  
  # permute weighted resid, then unweighted it for 999 times
  # permResid <- matrix(0, length(wterrors), nperm)
  # for (i in 1:nperm) {
  # if(!missing(seed)) set.seed(seed+i)
  # permResid[, i] <- as.vector(Ut%*%sample(wterrors)) 
  # }
  # permy <- as.data.frame(xbeta+permResid)
  if (ref) {
    if(!missing(seed)) set.seed(seed)
    permy <- as.data.frame(xbeta+replicate(nperm, as.vector(Ut%*%sample(wterrors))))
  }else{   
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      #as.data.frame(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## thien phuc
      as.data.frame(xbeta+as.vector(Ut%*%wterrors[rowindex]))
    })
  }
  
  # Calculating the likelihood ratio test statistic for each permutation.    
  ttest1 <- summary(lmer1)$coefficients["Group", "t value"] 
  # betahattest1 <- summary(lmer1)$coefficients["Group", "Estimate"]
  
  if (.Platform$OS.type=="windows") {
    cl <- makeCluster(ncore)	  
    clusterEvalQ(cl, library(lme4))
    clusterExport(cl, c("lmer0", "lmer1", "ref"), envir = environment()) 
    ttest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
      TT <- try(summary(refit(lmer1, x))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    })
    # betahattest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
    #   BHT <- try(summary(refit(lmer1, x))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # })
    stopCluster(cl)
  }else{  
    ttest2 <- mclapply(permy, function(x) {
      TT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    }, mc.cores=ncore)
    # betahattest2 <- mclapply(permy, function(x) {
    #   BHT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # }, mc.cores=ncore)
    
  }
  
  #Calculating the p-values.  
  aod<-list(NULL)
  ttest <- na.omit(unlist(ttest2))
  # ttest <- ifelse(ttest < 0, 0, ttest) # thien phuc
  # perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
  # aod$'Perm-p' <- c(NA, perm_p.t)
  if(pvalue==TRUE){
    perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
    if(upper==TRUE){
      aod$'Perm-p' <- 1-perm_p.t
    }else{
      if(lower==TRUE){
        aod$'Perm-p' <- perm_p.t
      }else{
        perm_p.t <- min(2*perm_p.t, 2-2*perm_p.t)
        aod$'Perm-p' <- c(NA, perm_p.t)
      }
    }
    
  }
  if(perc==TRUE){
    aod$'97.5-perc' <- quantile(ttest, 0.975)
    aod$'2.5-perc' <- quantile(ttest, 0.025)
  }
  
  # betahattest <- na.omit(unlist(betahattest2))
  # betahattest <- ifelse(betahattest < 0, 0, betahattest)
  # perm_p.betahat <- (sum(betahattest >= betahattest1) +1)/(length(betahattest) + 1)
  # aod$'Perm-p beta hat test' <- c(NA, perm_p.betahat)
  
  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}


permlmer.ci <- function(lmer0, lmer1, nperm = 999, ncore=3, n.grid.l=5, n.grid.u=5, data){
  
  mod <- summary(lmer1) ## beta
  mod <- mod$coefficients
  mod[1:2] <- mod["Group", c("Estimate", "Std. Error")]
  
  c.l <- mod[1]-seq(1.96, 4, length.out=n.grid.l)*mod[2]
  c.u <- mod[1]+seq(1.96, 4, length.out=n.grid.u)*mod[2]
  
  # p.l <- unlist(lapply(c.l, FUN=function(x) permlmer.tplh7(lmer0=lmer0, lmer1=lmer1, nperm = nperm, ncore=ncore, pvalue=TRUE, lower=TRUE, data=data, theta0=x))$'Perm-p')
  # p.u <- unlist(lapply(c.u, FUN=function(x) permlmer.tplh7(lmer0=lmer0, lmer1=lmer1, nperm = nperm, ncore=ncore, pvalue=TRUE, upper=TRUE, data=data, theta0=x))$'Perm-p')
  # 
  # if(min(p.l)>0.025){
  #   low <- c.l[n.grid.l]
  #   print("cannot find the satisfactory lower bound")
  # }else{
  #   low <- c.l[which(c.l==max(c.l) & p.l<=0.025)]
  # }
  # if(min(p.u)>0.025){
  #   hi <- c.u[n.grid.u]
  #   print("cannot find the satisfactory upper bound")
  # }else{
  #   hi <- c.u[which(c.u==min(c.u) & p.u<=0.025)]
  # }
  
  low <- c.l[n.grid.l]
  hi <- c.u[n.grid.u]
  for (i in 1:(n.grid.l-1)) {
    p.l<-permlmer.tplh7(lmer0=lmer0, lmer1=lmer1, nperm = nperm, ncore=ncore, pvalue=TRUE, lower=TRUE, data=data, theta0=c.l[i])$'Perm-p'   
    if(p.l<=0.025){
      low<-c.l[i]
      break
    }
  }
  print(i)
  for (i in 1:(n.grid.u-1)) {
    p.u<-permlmer.tplh7(lmer0=lmer0, lmer1=lmer1, nperm = nperm, ncore=ncore, pvalue=TRUE, upper=TRUE, data=data, theta0=c.u[i])$'Perm-p'   
    if(p.u<=0.025){
      hi<-c.u[i]
      break
    }
  }
  print(i)
  return(c(low, hi))
  
  
}
permlmer.tplh7 <- function(lmer0, lmer1, nperm = 999, ncore=3, plot=FALSE, seed, pvalue=FALSE, upper=FALSE, lower=FALSE, perc=FALSE, theta0, data){
  
  if (any(!inherits(lmer0, "lmerMod"), !inherits(lmer1, "lmerMod"))) stop("The model must be a lmer object!")
  if (!setequal(getME(lmer0, "y"), getME(lmer1, "y"))) stop("Please check the response in your model!")
  
  c_deparse <- function (...) 
    paste(deparse(..., width.cutoff = 500), collapse = "")
  lmernames <- vapply(as.list(sys.call()[-1L]), c_deparse, "")
  
  df <- data ## thien phuc
  df$y<-data$y-theta0*data$Group 
  lmer0<-update(lmer0, data=df)
  lmer1<-update(lmer1, data=df)
  
  theta0 <- getME(lmer0, "theta")
  theta1 <- getME(lmer1, "theta")
  fixef0 <- fixef(lmer0)
  fixef1 <- fixef(lmer1) 
  theta0name <- names(theta0)
  theta1name <- names(theta1)
  fixef0name <- names(fixef0)
  fixef1name <- names(fixef1) 
  
  term0name <- attr(terms(lmer0),"term.labels")
  term1name <- attr(terms(lmer1),"term.labels")
  term0in1 <- rep(FALSE, length(term0name))
  names(term0in1) <- term0name
  for (i in term0name) {
    for (j in term1name){
      if (length(setdiff(unlist(strsplit(i, "\\:")), unlist(strsplit(j, "\\:"))))==0) {
        term0in1[i] <- TRUE
        break
      }
    }
  }  
  
  # if(!setequal(intersect(theta0name, theta1name), theta0name) || !setequal(intersect(fixef0name, fixef1name), fixef0name)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(!setequal(intersect(theta0name, theta1name), theta0name) || !all(term0in1)) stop(paste("The model", lmernames[1], "must be nested within the model", lmernames[2]))
  
  if(setequal(fixef0name, fixef1name))  ref <- TRUE else ref <- FALSE
  
  thetan <- rep(0, length(theta1))
  names(thetan) <- theta1name
  thetan[theta0name] <- theta0
  # Lambda1n <- getME(lmer1, "Lambda") 
  # Lambda1n@x <- thetan[getME(lmer1, "Lind")]
  # Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  #V1n <- getME(lmer1, "Z")%*%Lambda1nc%*%getME(lmer1, "Zt")+diag(dim(getME(lmer1, "Z"))[1])
  Lambda1n <- getME(lmer0, "Lambda") ## thien phuc
  Lambda1n@x <- thetan[getME(lmer0, "Lind")]
  #Lambda1nc <- Matrix::tcrossprod(Lambda1n)
  Lambda1nc <- base::tcrossprod(Lambda1n)
  # V1n <- getME(lmer0, "Z")%*%Lambda1nc%*%getME(lmer0, "Zt")+diag(dim(getME(lmer0, "Z"))[1])
  V1n <- unique(diag(Lambda1nc))*getME(lmer0, "Z")%*%getME(lmer0, "Zt")+getME(lmer0, "sigma")^2*diag(dim(getME(lmer0, "Z"))[1])
  Ut <- t(chol(V1n))
  wt <- solve(Ut)
  xbeta <- as.vector(getME(lmer0, "X")%*%fixef0)
  errors <- getME(lmer0, "y") - xbeta
  
  #Weighting the residuals.
  wterrors <- wt%*%errors
  
  # permute weighted resid, then unweighted it for 999 times
  # permResid <- matrix(0, length(wterrors), nperm)
  # for (i in 1:nperm) {
  # if(!missing(seed)) set.seed(seed+i)
  # permResid[, i] <- as.vector(Ut%*%sample(wterrors)) 
  # }
  # permy <- as.data.frame(xbeta+permResid)
  if (ref) {
    if(!missing(seed)) set.seed(seed)
    permy <- as.data.frame(xbeta+replicate(nperm, as.vector(Ut%*%sample(wterrors))))
  }else{   
    if(!missing(seed)) set.seed(seed)  
    permy <- replicate(nperm, {
      rowindex <- sample(1:length(errors))
      #as.data.frame(xbeta[rowindex]+as.vector(Ut%*%wterrors[rowindex])) ## thien phuc
      as.data.frame(xbeta+as.vector(Ut%*%wterrors[rowindex]))
    })
  }
  
  # Calculating the likelihood ratio test statistic for each permutation.    
  ttest1 <- summary(lmer1)$coefficients["Group", "t value"] 
  # betahattest1 <- summary(lmer1)$coefficients["Group", "Estimate"]
  
  if (.Platform$OS.type=="windows") {
    cl <- makeCluster(ncore)	  
    clusterEvalQ(cl, library(lme4))
    clusterExport(cl, c("lmer0", "lmer1", "ref"), envir = environment()) 
    ttest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
      TT <- try(summary(refit(lmer1, x))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    })
    # betahattest2 <- parLapplyLB(cl, permy, function(x) { ## thien phuc
    #   BHT <- try(summary(refit(lmer1, x))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # })
    stopCluster(cl)
  }else{  
    ttest2 <- mclapply(permy, function(x) {
      TT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "t value"], TRUE)
      TT <- ifelse(is.numeric(TT), TT, NA)
    }, mc.cores=ncore)
    # betahattest2 <- mclapply(permy, function(x) {
    #   BHT <- try(summary(suppressMessages(refit(lmer1, x)))$coefficients["Group", "Estimate"], TRUE)
    #   BHT <- ifelse(is.numeric(BHT), BHT, NA)
    # }, mc.cores=ncore)
    
  }
  
  #Calculating the p-values.  
  aod<-list(NULL)
  ttest <- na.omit(unlist(ttest2))
  # ttest <- ifelse(ttest < 0, 0, ttest) # thien phuc
  # perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
  # aod$'Perm-p' <- c(NA, perm_p.t)
  if(pvalue==TRUE){
    perm_p.t <- (sum(ttest >= ttest1) +1)/(length(ttest) + 1)
    if(upper==TRUE){
      aod$'Perm-p' <- 1-perm_p.t
    }else{
      if(lower==TRUE){
        aod$'Perm-p' <- perm_p.t
      }else{
        perm_p.t <- min(2*perm_p.t, 2-2*perm_p.t)
        aod$'Perm-p' <- c(NA, perm_p.t)
      }
    }
    
  }
  if(perc==TRUE){
    aod$'97.5-perc' <- quantile(ttest, 0.975)
    aod$'2.5-perc' <- quantile(ttest, 0.025)
  }
  
  # betahattest <- na.omit(unlist(betahattest2))
  # betahattest <- ifelse(betahattest < 0, 0, betahattest)
  # perm_p.betahat <- (sum(betahattest >= betahattest1) +1)/(length(betahattest) + 1)
  # aod$'Perm-p beta hat test' <- c(NA, perm_p.betahat)
  
  if (plot) {
    dev.new()
    plot (density(c(lrtest1, lrtest), kernel = "epanechnikov"), col="blue", lwd=2, xlab = "", main = "LR Test's density kernels")
    abline(v=lrtest1, col="red")
  }
  return(aod)
}


permlmer.ci <- function(lmer0, lmer1, nperm = 999, ncore=3, n.grid.l=5, n.grid.u=5, data){
  
  mod <- summary(lmer1) ## beta
  mod <- mod$coefficients
  mod[1:2] <- mod["Group", c("Estimate", "Std. Error")]
  
  c.l <- mod[1]-seq(1.96, 4, length.out=n.grid.l)*mod[2]
  c.u <- mod[1]+seq(1.96, 4, length.out=n.grid.u)*mod[2]
  
  # p.l <- unlist(lapply(c.l, FUN=function(x) permlmer.tplh7(lmer0=lmer0, lmer1=lmer1, nperm = nperm, ncore=ncore, pvalue=TRUE, lower=TRUE, data=data, theta0=x))$'Perm-p')
  # p.u <- unlist(lapply(c.u, FUN=function(x) permlmer.tplh7(lmer0=lmer0, lmer1=lmer1, nperm = nperm, ncore=ncore, pvalue=TRUE, upper=TRUE, data=data, theta0=x))$'Perm-p')
  # 
  # if(min(p.l)>0.025){
  #   low <- c.l[n.grid.l]
  #   print("cannot find the satisfactory lower bound")
  # }else{
  #   low <- c.l[which(c.l==max(c.l) & p.l<=0.025)]
  # }
  # if(min(p.u)>0.025){
  #   hi <- c.u[n.grid.u]
  #   print("cannot find the satisfactory upper bound")
  # }else{
  #   hi <- c.u[which(c.u==min(c.u) & p.u<=0.025)]
  # }
  
  low <- c.l[n.grid.l]
  hi <- c.u[n.grid.u]
  for (i in 1:(n.grid.l-1)) {
    p.l<-permlmer.tplh7(lmer0=lmer0, lmer1=lmer1, nperm = nperm, ncore=ncore, pvalue=TRUE, lower=TRUE, data=data, theta0=c.l[i])$'Perm-p'   
    if(p.l<=0.025){
      low<-c.l[i]
      break
    }
  }
  print(i)
  for (i in 1:(n.grid.u-1)) {
    p.u<-permlmer.tplh7(lmer0=lmer0, lmer1=lmer1, nperm = nperm, ncore=ncore, pvalue=TRUE, upper=TRUE, data=data, theta0=c.u[i])$'Perm-p'   
    if(p.u<=0.025){
      hi<-c.u[i]
      break
    }
  }
  print(i)
  return(c(low, hi))
  
  
}