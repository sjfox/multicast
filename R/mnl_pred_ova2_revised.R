library(MNLpred)

mnl_pred_ova2_revised <- function (model, data, x, xvals, z = NULL, z_value = NULL, 
                                   xvari, scenname, scenvalue, nsim = 1000, seed = "random", 
                                   probs = c(0.025, 0.975)) 
{
  # browser()
  output <- list()
  if (!missing(xvari)) {
    warning("The argument 'xvari' is deprecated; please use 'x' instead.\n\n", 
            call. = FALSE)
    x <- xvari
  }
  if (!missing(scenname)) {
    warning("The argument 'scenname' is deprecated; please use 'z' instead.\n\n", 
            call. = FALSE)
    z <- scenname
  }
  if (!missing(scenvalue)) {
    warning("The argument 'scenvalue' is deprecated; please use 'z_value' instead.\n\n", 
            call. = FALSE)
    z_value <- scenvalue
  }
  if (is.null(model) == TRUE) {
    stop("Please supply a model")
  }
  if (sum(grepl("multinom", model$call)) == 0) {
    stop("Please supply a multinom()-model")
  }
  if (is.null(data) == TRUE) {
    stop("Please supply a data set")
  }
  if (is.null(x) == TRUE | is.character(x) == FALSE) {
    stop("Please supply a character string of your x-variable of interest")
  }
  if (is.null(model$Hessian) == TRUE) {
    stop("There is no Hessian matrix. Please specify Hess = TRUE in your multinom() call.")
  }
  variables <- as.character(attr(model$terms, "variables"))[-1]
  if (!(x %in% variables) == TRUE) {
    stop("x-variable is not an independent variable in the model. There might be a typo.")
  }
  if (is.null(z) == FALSE & is.character(z) == FALSE) {
    stop("Please supply a character string of your scenario of interest")
  }
  if (is.null(z) == FALSE) {
    if (!(z %in% variables) == TRUE) {
      stop("The scenario variable is not an independent variable in the model. There might be a typo.")
    }
  }
  iv <- variables[2:length(variables)]
  output[["IV"]] <- iv
  if (length(iv) > 1) {
    if (sum(apply(data[, iv], 2, class) %in% c("numeric", 
                                               "integer")) < ncol(data[, iv])) {
      stop("Please supply data that consists of numeric values. The package can not handle factor or character variables, yet. For workarounds, please take a look at the github issues (https://github.com/ManuelNeumann/MNLpred/issues/1). The problem will hopefully be fixed with the 0.1.0 release.")
    }
  }
  else {
    if (class(eval(parse(text = paste0("data$", iv)))) %in% 
        c("numeric", "integer") == FALSE) {
      stop("Please supply data that consists of numeric values. The package can not handle factor or character variables, yet. For workarounds, please take a look at the github issues (https://github.com/ManuelNeumann/MNLpred/issues/1). The problem will hopefully be fixed with the 0.1.0 release.")
    }
  }
  dv <- variables[1]
  output[["DV"]] <- dv
  data_redux <- na.omit(data[, c(dv, iv)])
  obs <- nrow(data_redux)
  output[["Observations"]] <- obs
  coefmatrix <- coef(model)
  ncoef <- ncol(coefmatrix)
  mu <- as.vector(t(coef(model)))
  varcov <- MASS::ginv(model$Hessian)
  if (seed != "random") {
    set.seed(seed = seed)
  }
  S <- mvrnorm(nsim, mu, varcov)
  output[["S"]] <- S
  if (is.null(by) == TRUE) {
    by <- abs(min(eval(parse(text = paste0("data$", x))), 
                  na.rm = TRUE) - max(eval(parse(text = paste0("data$", 
                                                               x))), na.rm = TRUE))
  }
  variation <- xvals
  output[["Variation"]] <- variation
  nseq <- length(variation)
  if (nseq == 1) {
    stop("Please supply a dataset or a x-variable with variation")
  }
  output[["nVariation"]] <- nseq
  categories <- sort(unique(eval(parse(text = paste0("data$", 
                                                     dv)))))
  J <- length(categories)
  if (J < 3) {
    stop("Please supply a dataset with a dependent variable that has a sufficient number of outcomes (> 2)")
  }
  output[["ChoiceCategories"]] <- categories
  output[["nChoices"]] <- J
  ninteraction <- sum(grepl(":", model$coefnames))
  
  X <- matrix(NA, ncol = ncoef, nrow = obs)
  colnames(X) <- model$coefnames
  X[, 1] <- 1
  X[, 2:(length(iv) + 1)] <- as.matrix(data_redux[, iv])
  ovacases <- array(NA, c(dim(X), nseq))
  ovacases[, , ] <- X
  varidim <- which(colnames(X) == x)
  for (i in 1:nseq) {
    ovacases[, varidim, i] <- variation[i]
  }
  if (is.null(z) == FALSE) {
    scendim <- which(colnames(X) == z)
    for (i in 1:nseq) {
      ovacases[, scendim, i] <- z_value
    }
  }
  if (ninteraction != 0) {
    interactionterms <- which(grepl(":", model$coefnames) == 
                                TRUE)
    for (i in c(interactionterms)) {
      firstint <- gsub(":.*", "", model$coefnames[i])
      secondint <- gsub(".*:", "", model$coefnames[i])
      intdim1 <- which(colnames(X) == firstint)
      intdim2 <- which(colnames(X) == secondint)
      for (j in 1:nseq) {
        ovacases[, i, j] <- ovacases[, intdim1, j] * 
          ovacases[, intdim2, j]
      }
    }
  }
  ovaV <- array(NA, c(obs, nsim, nseq, J))
  pb_multiplication <- txtProgressBar(min = 0, max = nseq, 
                                      initial = 0)
  # browser()
  cat("Multiplying values with simulated estimates:\n")
  for (i in 1:nseq) {
    ovaV[, , i, 1] <- apply(matrix(0, nrow = nsim, ncol = ncol(X)), 
                            1, function(s) ovacases[, , i] %*% s)
    for (k in 2:J) {
      coefstart <- (k - 2) * ncoef + 1
      coefend <- (k - 1) * ncoef
      element <- parse(text = paste0("ovaV[,, i,", k, "] <- apply(S[, ", 
                                     coefstart, ":", coefend, "], 1, function(s) ovacases[,, i] %*% s)"))
      eval(element)
    }
    setTxtProgressBar(pb_multiplication, i)
  }
  Sexp <- rowSums(exp(ovaV), dims = 3L)
  P <- array(NA, c(nsim, J, nseq))
  pb_link <- txtProgressBar(min = 0, max = nseq, initial = 0)
  cat("\nApplying link function:\n")
  for (l in 1:nseq) {
    for (m in 1:J) {
      P[, m, l] <- colMeans(exp(ovaV[, , l, m])/Sexp[, 
                                                     , l])
      if (sum(is.na(P[, m, l])) != 0) {
        stop("Stop")
      }
    }
    setTxtProgressBar(pb_link, l)
  }
  output[["P"]] <- P
  if (is.null(z_value) == TRUE) {
    plotdat <- data.frame(iv = rep(variation, J), categories = rep(categories, 
                                                                   each = length(variation)), mean = NA, lower = NA, 
                          upper = NA)
  }
  else {
    plotdat <- data.frame(iv = rep(variation, J), categories = rep(categories, 
                                                                   each = length(variation)), scen = rep(z_value, each = length(variation)), 
                          mean = NA, lower = NA, upper = NA)
  }
  start <- 1
  for (i in 1:J) {
    end <- i * length(variation)
    plotdat[c(start:end), "mean"] <- colMeans(P[, i, ])
    plotdat[c(start:end), "lower"] <- apply(P[, i, ], 2, 
                                            quantile, probs = probs[1])
    plotdat[c(start:end), "upper"] <- apply(P[, i, ], 2, 
                                            quantile, probs = probs[2])
    start <- end + 1
  }
  if (is.null(z) == TRUE) {
    colnames(plotdat)[1:2] <- c(x, dv)
  }
  else {
    colnames(plotdat)[1:3] <- c(x, dv, z)
  }
  output[["plotdata"]] <- plotdat
  cat("\nDone!\n\n")
  return(output)
}
