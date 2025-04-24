library(tibble)
library(tidyr)
library(dplyr)

titerset_long <- titerset_b %>% 
  mutate(strain = "B_victoria") %>% 
  rows_append(titerset_h1n1 %>% 
                mutate( strain = "H1N1")) %>% 
  rows_append(titerset_h3n2 %>% 
                mutate( strain = "H3N2")) %>% 
  select(SubjectID = PID, Pre = baseline_titers, Post = outcome_titers, Strain = strain)
  
titerset_format <- FormatTiters(titerset_long)

FormatTiters <- function(titers, log2Transform = TRUE, fcMinZero = TRUE)
{
  if(log2Transform) {
    message("- Log transforming Pre and Post columns")
    trans <- log2
  } else {
    trans <- identity
  }
  if(fcMinZero) {
    message("- Setting any negative log fold changes to 0")
  }
  ## class_titers <- class(titers)
  ## if(length(class_titers) > 1 || class_titers[1] != "data.frame") {
  ##     message("- `titers` argument is not a data frame. Strange behavior has been observed with tibbles so coercing to data frame with `as.data.frame(titers)`")
  ##     titers <- as.data.frame(titers)
  ## }
  titer_list <- list()
  strains <- sort(unique(titers$Strain))
  for(i in seq_along(strains)) {
    result <- titers %>%
      dplyr::filter(Strain == strains[i]) %>%
      dplyr::mutate(Post = trans(Post), Pre = trans(Pre)) %>%
      dplyr::mutate(FC = Post - Pre) %>%
      dplyr::arrange(Pre, FC) %>%
      distinct()
    if(fcMinZero) {
      result <- result %>% dplyr::mutate(FC = pmax(FC, 0))
    }
    rownames(result) <- NULL    # remove rownames which may cause problems later
    titer_list[[strains[i]]] <- result
  }
  return(titer_list)
}



Calculate_maxRBA <- function(dat_list, subjectCol = "SubjectID",
                             method = c("exp", "lm"), yMinZero = FALSE,
                             scoreFun = max, discretize = c(0.2, 0.3),
                             normalize = FALSE, scaleResiduals = FALSE,
                             responseLabels = paste0(c("low", "moderate", "high"),
                                                     "Responder"), na_action = "na.fail",
                             ...) {
  method <- match.arg(method)
  if(length(unique(lapply(dat_list, dim))) != 1) {
    stop("Each data frame in `dat_list` must have the same dimensions")
  }
  if(method == "exp" && scaleResiduals) {
    warning("Scaling of residuals is not implemented for method == 'exp'.")
    scaleResiduals <- FALSE
  }
  ## Inverse Normal Transform
  .INT <- function(x, na.last = "keep",
                   ties.method = c("average", "first", "last",
                                   "random", "max", "min"), ...) {
    ties.method <- match.arg(ties.method)
    xranks <- rank(x, na.last = na.last, ties.method = ties.method)
    tempp <- (xranks - 0.5)/length(xranks)
    return(qnorm(tempp, ...))
  }
  residuals_list <- model_list <- list()
  ## Calculate residual for each strain
  for(i in seq_along(dat_list)) {
    dat <- dat_list[[i]]
    ## Check if arranged in order of Pre column
    ord <- order(dat$Pre)
    if(!all(ord == 1:nrow(dat))) {
      stop("`dat_list[[", i, "]]` is not ordered by 'Pre' column. Use only output from `FormatTiters`!")
    }
    if(method == "lm") {
      model <- lm(data = dat, formula = "FC ~ Pre", na.action = na_action, ...)
    } else if (method == "exp") {
      model <- nls(data = dat, formula = "FC ~ exp(a + b * Pre)",
                   start = list(a = 0, b = 0),
                   na.action = na_action, ...)
    }
    if(yMinZero && method == "lm") {
      residuals <- residuals(model)
      setToZero <- fitted(model) < 0
      residuals[setToZero] <- 0
    } else {
      residuals <- residuals(model)
    }
    names(residuals) <- dat[[subjectCol]]
    if(scaleResiduals) {
      cis <- stats::predict(model, na.action = na_action,
                            newdata = dat, se.fit = FALSE,
                            level = 0.95, interval = "confidence")
      intervals <- apply(cis, 1, function(row) { return(row["upr"] - row["lwr"]) })
      residuals <- residuals / intervals^2    
    }
    residuals_list[[i]] <- residuals[order(names(residuals))]
    model_list[[i]] <- model
  }
  ##:ess-bp-start::browser@nil:##
  if(length(residuals_list) > 1) {
    residual_mat <- Reduce(cbind, residuals_list)
  } else {
    residual_mat <- as.matrix(residuals_list[[1]])
  }
  colnames(residual_mat) <- names(dat_list)
  if(normalize) {
    residual_mat <- apply(residual_mat, 2, .INT)
  }
  maxRBA <- apply(residual_mat, 1, scoreFun, na.rm = TRUE)
  ## Calculated discretized metrics
  disList <- vector(mode = "list", length = length(discretize))
  names(disList) <- discretize
  for(dis in discretize) {
    tmp <- rep(NA, nrow(dat))
    names(tmp) <- names(maxRBA)
    lowR <- maxRBA <= quantile(maxRBA, dis, na.rm = TRUE)
    highR <- maxRBA >= quantile(maxRBA, 1 - dis, na.rm = TRUE)
    modR <- !(lowR | highR | is.na(maxRBA))
    tmp[lowR] <- 0
    tmp[modR] <- 1
    tmp[highR] <- 2
    disList[[as.character(dis)]] <- factor(tmp, labels = responseLabels, levels = 0:2)
  }
  disList <- setNames(disList, paste0("maxRBA_d",
                                      as.character(as.numeric(names(disList))*100)))
  names(model_list) <- names(dat_list)
  return(c(models = list(model_list), residualMatrix = list(residual_mat),
           maxRBA = list(maxRBA), disList))
}

titerset_format_2 <- Calculate_maxRBA(titerset_format)
View(titerset_format_2$maxRBA)

# Assuming titerset_format_2$maxRBA is a named vector
long_df <- enframe(titerset_format_2$maxRBA, name = "PID", value = "maxRBA") %>% 
  mutate(PID = as.double(PID))

# View result
print(long_df)

titerset_comb_corrected_seroneg %>% 
  left_join(long_df) %>% view()




corrected_FC_matrix <- titerset_format_2$residualMatrix
corrected_FC_df <- corrected_FC_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "PID") %>% 
  mutate(PID = as.double(PID))


titerset_comb_corrected_seroneg %>% 
  left_join(corrected_FC_df %>% select(PID, maxRBA_calcFC = B_victoria)) %>% 
  view()

