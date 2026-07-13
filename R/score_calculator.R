#' @import gamlss.dist

utils::globalVariables(c("mu", "sigma", "nu", "tau"))

#devtools::use_data(diamonds, overwrite = TRUE)
internal_env <- new.env()
df_name <- load("data/all_para_tables.RData", internal_env)

# initialize
sex_cats <- c("boys", "girls")
var_names <- c("bmi", "dbp", "glu", "hdl", "height", "homa", "insu", "MetS_shifted", "sbp", "trg", "waist", "crp")
par_cats <- c("cutoffp1", "mu", "sigma", "nu", "tau")



#' @title 2D Interpolation Function Constructor
#' @description Creates a closure that performs 2D interpolation over a grid using pracma::interp2.
#'
#' @param x Numeric vector. Grid values along the x-axis.
#' @param y Numeric vector. Grid values along the y-axis.
#' @param Z Matrix. Grid values corresponding to (x, y).
#' @param method Character. Interpolation method; default is "linear".
#'
#' @return A function that takes `x` and `y` (e.g. for age and height) vectors and returns interpolated values for `z`.
#' @keywords internal
approxfun2 <- function(x, y, Z, method = "linear") {
  force(x)
  force(y)
  force(Z)
  force(method)

  function(age, height) {
    #handle NAs
    none_na <- !(is.na(age) | is.na(height))
    rs <- rep(NA,length(none_na))
    if (sum(none_na)>0)
      rs[none_na] <- pracma::interp2(x, y, Z, age[none_na], height[none_na], method)
    return(rs)
  }
}




# vector for the distributions per variable
var_distribution <- list()

# fit splines for each parameter, variable and sex
approx_param_functions <- list()
for (sex in sex_cats) {
  for (vname in var_names) {
    # get table
    curr_tab <- get(paste("par", vname, sex, sep="_"), envir= internal_env)

    # remember distribution type
    if (is.null(var_distribution[[sex]])) var_distribution[[sex]] <- list()
    var_distribution[[sex]][[vname]] <- as.character(unique(curr_tab[["dist"]])[1])

    age_values    <- curr_tab[["age"]]
    height_values <- curr_tab[["height"]]
    curr_cats <- intersect(par_cats, colnames(curr_tab))
    for (param in curr_cats){
      parameters <- curr_tab[[param]]
      if (!is.null(height_values)) {
        parameter_matrix <- tapply(curr_tab[[param]], list(curr_tab$height, curr_tab$age), FUN = identity)
        approx_param_functions[[paste(vname,sex,param, sep="_")]] <- approxfun2(unique(age_values), unique(height_values), parameter_matrix)
        rm(parameter_matrix)
      } else {
        approx_param_functions[[paste(vname,sex,param, sep="_")]] <- approxfun(age_values, parameters, ties=mean)
      }
    }
  }
}



#' @title Create a Sex- and Age-Specific Interpolation Function
#'
#' @description Constructs a closure that returns interpolated values based on sex and age.
#' Internally selects from precomputed interpolation functions.
#'
#' @param vname Character. Variable name (e.g., "bmi", "waist").
#' @param param Character. Distribution parameter (e.g., "mu", "sigma").
#'
#' @return A function with signature `(sex, age)` returning interpolated values.
#' @keywords internal
mkfun_sex_age <- function(vname, param) {
  function(sex, age) {
    ifelse(is.na(sex) | is.na(age),
           NA,
          ifelse(sex=="m",
                 approx_param_functions[[paste(vname,"boys", param, sep="_")]](age),
                 approx_param_functions[[paste(vname,"girls",param, sep="_")]](age))
           )

  }
}

#' @title Create a Sex-, Age-, and Height-Specific Interpolation Function
#'
#' @description Constructs a closure for interpolating values that depend on sex, age, and height
#' Chooses appropriate internal spline function based on `vname` and `param`.
#'
#' @param vname Character. Variable name (e.g., "bmi", "waist").
#' @param param Character. Distribution parameter (e.g., "mu", "sigma").
#'
#' @return A function with signature `(sex, age, height)` returning interpolated values.
#' @keywords internal
mkfun_sex_age_height <- function(vname, param) {
  function(sex, age, height) {
    ifelse(is.na(sex) | is.na(age) | is.na(height),
           NA,
          ifelse(sex=="m",
                 approx_param_functions[[paste(vname,"boys", param, sep="_")]](age, height),
                 approx_param_functions[[paste(vname,"girls",param, sep="_")]](age, height))
    )
  }
}


for (vname in var_names) {
  all_cols <- lapply(sex_cats, function(x) colnames(get(paste("par", vname, x, sep="_"), envir= internal_env)))
  cols_for_all_sexes <- Reduce(intersect, all_cols)
  curr_cats <- intersect(par_cats, cols_for_all_sexes)
  for (param in curr_cats) {

    fname <- paste(vname, param, sep = "_")

    # Use a local scope to capture vname and param
    local({
      vname_ <- vname
      param_ <- param

      if (!"height" %in% cols_for_all_sexes) {
        approx_param_functions[[fname]] <<- mkfun_sex_age(vname_, param_)
      } else {
        approx_param_functions[[fname]] <<- mkfun_sex_age_height(vname_, param_)
      }
    })
  }
}

#approx_param_functions$dbp_boys_sigma(5, 120)
#approx_param_functions$dbp_mu(c("f","m","f","m"),c(5,5,5,5), c(120,120,121,121))
#approx_param_functions$bmi_sigma(c("f","m"),c(5,5))


#' @title Percentile calculator for the hurdle model
#'
#' @description
#'
#' @param y Numeric vector. Observed values of the variable to score.
#' @param p1 Numeric vector. First part of the hurdle model, representing the cutoff percentile from the logistic model.
#'
#' @return A numeric vector representing the calculated percentile based on the hurdle model
#' @keywords internal
p_hurdle <- function(y, p1, mu, sigma, nu, tau) {

  p <- p1 + (1 - p1) * pGB1(y,
                            mu = mu,
                            sigma = sigma,
                            nu = nu,
                            tau = tau)

  p[y <= 0] <- p1[y <= 0]

  p
}





#' @title Calculate Scores for Children and Young Adults
#'
#' @description Computes age-, sex- and (for `sbp` and `dbp`) height-specific percentiles or z-scores for anthropometric and metabolic variable
#' using IDEFICS study reference data and Biomarkers4Pediatrics collaboration (for `crp`).
#'
#' @param variable Character. The variable to assess. Must be one of the supported variables ("waist", "bmi", "hdl", "sbp", "dbp", "trg", "homa", "glu", "height","insu", or "crp").
#' @param sex Character vector. Same length as `age`, `height`, and `values`. Accepts "f" for female or "m" for male.
#' @param age Numeric vector. Ages of the children in years. Must be between 1 and 22 for `crp`, and between 2 and 11 otherwise.
#' @param height Numeric vector or NULL. Required for height-dependent models (`sbp` and `dbp`). Defaults to NULL.
#' @param values Numeric vector. Observed values of the variable to score.
#' @param return_values Character vector. Specifies which outputs to return. Options include "percentile", "z.score".
#'
#' @return A named list containing the requested scores. Each element is a numeric vector of the same length as `values`.
#'
#' @examples
#' get_scores(
#'   variable = "crp",
#'   sex = c("f","f"),
#'   age = c(3, 12),
#'   values = c(4, 0.5)
#'   )
#'
#' get_scores(
#'   variable="dbp",
#'   sex=c("f","m"),
#'   age=c(5,5),
#'   height=c(120,110),
#'   values=c(70,60)
#'   )
#'
#' @export
get_scores <- function(variable="waist", sex=c("f","m"), age=6:5, height=NULL, values=c(20,21), return_values=c("percentile","z.score"))  {
  if (length(variable)>1) warning("variable must be of length 1. Only the first value of variable is used.")

  sex_map <- c(f = "girls", m = "boys")
  dist <- var_distribution[[sex_map[sex][1] ]][[variable]]

  all_cols <- lapply(sex_cats, function(x) colnames(get(paste("par", variable, x, sep="_"), envir= internal_env)))
  cols_for_all_sexes <- Reduce(intersect, all_cols)
  curr_cats <- intersect(par_cats, cols_for_all_sexes)

  # assign parameters mu, sigma, nu, tau if used by current model
  if ("height" %in% cols_for_all_sexes)
    for (param in curr_cats) assign(param, approx_param_functions[[paste(variable,param,   sep="_")]](sex, age, height) )
  else
    for (param in curr_cats) assign(param, approx_param_functions[[paste(variable,param,   sep="_")]](sex, age) )

  if (dist == "BCCG") {
    #handle NAs
    none_na <- !(is.na(values) | is.na(mu) | is.na(sigma) | is.na(nu))

    percentiles <- rep(NA,length(none_na))
    percentiles[none_na] <- gamlss.dist::pBCCG(q = values[none_na], mu[none_na], sigma[none_na], nu[none_na])
  } else if (dist == "BCT") {

    #handle NAs
    none_na <- !(is.na(values) | is.na(mu) | is.na(sigma) | is.na(nu) | is.na(tau))

    percentiles <- rep(NA,length(none_na))
    percentiles[none_na] <- gamlss.dist::pBCT(q = values[none_na], mu[none_na], sigma[none_na], nu[none_na], tau[none_na])

  } else if (dist == "BCPE") {

    #handle NAs
    none_na <- !(is.na(values) | is.na(mu) | is.na(sigma) | is.na(nu) | is.na(tau))

    percentiles <- rep(NA,length(none_na))
    percentiles[none_na] <- gamlss.dist::pBCPE(q = values[none_na], mu[none_na], sigma[none_na], nu[none_na], tau[none_na])

  } else if (dist == "LO") {

    #handle NAs
    none_na <- !(is.na(values) | is.na(mu) | is.na(sigma))

    percentiles <- rep(NA,length(none_na))
    percentiles[none_na] <- gamlss.dist::pLO(q = values[none_na], mu[none_na], sigma[none_na])

  } else if (dist == "GB1") {
    values <- values/10 - 0.02

    #handle NAs
    none_na <- !(is.na(values) | is.na(mu) | is.na(sigma) | is.na(nu) | is.na(tau))

    percentiles <- rep(NA,length(none_na))
    percentiles[none_na] <- p_hurdle(y = values[none_na], p1 = cutoffp1[none_na],
                                     mu = mu[none_na], sigma = sigma[none_na],
                                     nu = nu[none_na], tau = tau[none_na])

  }

  # z.score requested?
  if ("z.score" %in% return_values) {
    zscores <- rep(NA,length(none_na))
    zscores[none_na] <- gamlss.dist::qNO(percentiles[none_na])
  }

  # build up return values as list
  rs <- list()
  if ("percentile" %in% return_values) rs$percentile <- percentiles
  if ("z.score" %in% return_values) rs$z.score <- zscores
  return(rs)
}




#' @title Calculate Metabolic Syndrome (MetS) Score
#'
#' @description Computes a composite MetS score from standardized z-scores of metabolic indicators.
#'
#' @param df A data frame containing the following columns: `waist_z.score`, `homa_z.score`, `sbp_z.score`, `dbp_z.score`, `trg_z.score`, and `hdl_z.score`.
#'           These are expected to be numeric vectors representing z-scores.
#'
#' @return A numeric vector representing the calculated MetS score for each row in the input data.
#'
#' @details The MetS score is computed as:\cr
#' `waist_z + homa_z + 0.5 × (sbp_z + dbp_z + trg_z - hdl_z)`
#'
#' @examples
#' df <- data.frame(
#'   waist_z.score = c(1.2, 0.5),
#'   homa_z.score = c(0.8, 0.6),
#'   sbp_z.score = c(0.7, 0.3),
#'   dbp_z.score = c(0.6, 0.4),
#'   trg_z.score = c(1.0, 0.9),
#'   hdl_z.score = c(-0.5, -0.2)
#' )
#' MetSScore(df)
#'
#' @export
MetSScore <- function(df) {
  if( is.null(df$waist_z.score) || is.null(df$homa_z.score) ||
      is.null(df$sbp_z.score) || is.null(df$dbp_z.score) || is.null(df$trg_z.score) || is.null(df$hdl_z.score) )
    stop("MatS score calculation requires the following columns not to be NULL: waist_z.score, homa_z.score, sbp_z.score, dbp_z.score, trg_z.score, hdl_z.score")
  MetS <- df$waist_z.score + df$homa_z.score +
    0.5*(df$sbp_z.score + df$dbp_z.score + df$trg_z.score - df$hdl_z.score)

  return(MetS)
}

#' @title Compute Action Levels
#'
#' @description Assigns monitoring or intervention levels based on variable percentiles using B4P thresholds for `crp` and standard IDEFICS thresholds
#' for all other variables.
#'
#' @param df A data frame containing percentile columns such as `waist_percentile`, `sbp_percentile`, `hdl_percentile`, etc.
#' @param lvl_name Character vector of level labels. Defaults to `c("none", "monit", "action")`.
#' @param perc_level Numeric vector of two percentiles used as cutoffs. Defaults to `c(0.9, 0.95)` for 90th and 95th percentile.
#' @param append Logical. If `TRUE`, appends action level columns to `df`. If `FALSE`, returns only the computed levels.
#' @param filter Character. Optional. Limits calculation to a specific domain: `"adiposity"`, `"blood_pressure"`, `"blood_lipids"`, `"blood_glu_insu"`, or `"overall"`.
#'
#' @return A list (or a modified data frame if `append = TRUE`) containing action level classifications for each domain.
#'
#' @details
#' Action levels are derived using `cut()` on percentile values. For example, a value > 95th percentile maps to `"action"`.
#' HDL is reversed (`1 - hdl_percentile`) since lower HDL values are riskier.
#'
#' @examples
#' df <- data.frame(waist_percentile = c(0.85, 0.96))
#' action_levels(df)
#'
# df <- data.frame(
#   sex = c("m", "m"),
#   hdl_percentile = c(0.1,0.5),
#   homa_percentile = c(0.4,0.9),
#   trg_percentile = c(0.6,0.5),
#   crp_percentile = c(0.95, 0.9),
#   waist_percentile = c(0.9,0.99),
#   sbp_percentile = c(0.8,0.01)
# )
#' action_levels(df)
#'
#' @export
action_levels <- function(df, lvl_name=c("none","monit","action"), perc_level=c(0.9, 0.95), append=FALSE, filter=NULL) {
  n <- nrow(df)

  rs <- list()

  # Helper function which maps a numeric percentile vector to ordered factor levels (e.g., "none", "monit", "action") using pre-defined cutoffs.
  perc_to_actlev <- function(perc) {
    if (length(perc) == 0) # catches NULL as well as 1-NULL
      ordered(rep(NA, n), levels=1:length(lvl_name), labels=lvl_name)
    else
      cut(perc, c(-Inf,perc_level,Inf), lvl_name, ordered_result = TRUE)
  }

  crp_perc_to_actlev <- function(perc, sex) {

    cutoffs <- ifelse(sex == "m", 0.935, 0.899) # hard-coded sex-specific cutoffs

    ordered(
      ifelse(perc > cutoffs, "Elevated", "Not elevated"),
      levels = c("Not elevated", "Elevated")
    )
  }

  if ("adiposity" %in% filter || "overall" %in% filter || (is.null(filter) && !is.null(df$waist_percentile)))
    rs$adiposity.action <- perc_to_actlev(df$waist_percentile)

  if ("crp" %in% filter || "overall" %in% filter || (is.null(filter) && !is.null(df$crp_percentile)))
    rs$crp.action <- crp_perc_to_actlev(df$crp_percentile, df$sex)

  if ("blood_pressure" %in% filter || "overall" %in% filter || (is.null(filter) && !(is.null(df$dbp_percentile) && is.null(df$sbp_percentile) ))) {
    dbp.action <- perc_to_actlev(df$dbp_percentile)
    sbp.action <- perc_to_actlev(df$sbp_percentile)
    rs$blood_pressure.action <- suppressWarnings(pmax(dbp.action,sbp.action, na.rm=T))
  }

  if ("blood_lipids" %in% filter || "overall" %in% filter || (is.null(filter) && !(is.null(df$trg_percentile) && is.null(df$hdl_percentile) ))) {
    trg.action <- perc_to_actlev(df$trg_percentile)
    hdl.action <- perc_to_actlev(1-df$hdl_percentile)
    rs$blood_lipids.action <- suppressWarnings(pmax(trg.action,hdl.action, na.rm=T))
  }

  if ("blood_glu_insu" %in% filter || "overall" %in% filter || (is.null(filter) && !(is.null(df$homa_percentile) && is.null(df$glu_percentile) ))) {
    homa.action <- perc_to_actlev(df$homa_percentile)
    glu.action <- perc_to_actlev(df$glu_percentile)
    rs$blood_glu_insu.action <- suppressWarnings(pmax(homa.action,glu.action, na.rm=T))
  }

  if ((is.null(filter) && sum( c(!is.null(df$waist_percentile), !is.null(df$dbp_percentile) || !is.null(df$sbp_percentile), !is.null(df$trg_percentile) || !is.null(df$hdl_percentile), !is.null(df$homa_percentile) || !is.null(df$glu_percentile) ) ) >= 3) || "overall" %in% filter) {
    rs$overall.action <-
      apply(cbind(rs$adiposity.action,rs$blood_pressure.action,rs$blood_lipids.action,rs$blood_glu_insu.action), # converts also factors to integer...
            1, # apply rowwise
            function(x) {
              # is level 1 / 2 / 3 reached or exceeded?
              compare_with_123 <- sapply(x, function(lev) lev >= 1:(length(perc_level)+1))
              #rownames(compare_with_123) <- lvl_name
              # which of the levels 1, 2, 3 are exceeded/reached at least 3 times
              levelcheck_123 <- apply(compare_with_123, 1, function(x) sum(x, na.rm=T)>=3)
              # return maximum level that is reached or exceeded at least three times
              levs <- (1:(length(perc_level)+1))[levelcheck_123]

              if (length(levs) == 0) NA_integer_ else max(levs)

            })
    rs$overall.action <- ordered(rs$overall.action, 1:(length(perc_level)+1), lvl_name)
  }

  if (!is.null(filter)) rs <- rs[names(rs) %in% paste0(filter,".action")]

  if (append) rs <- cbind(df,rs)

  return(as.data.frame(rs))
}





#' @title  Compute IDEFICS and B4P Scores for Multiple Variables
#'
#' @description Applies `get_scores()` to multiple anthropometric and metabolic variables in a data frame and optionally returns MetS and action levels.
#'
#' @param df A data frame with columns: `sex`, `age`, `height`, and observed values for any of the supported variables:
#'           `bmi`, `glu`, `hdl`, `height`, `homa`, `insu`, `trg`, `waist`, `sbp`, `dbp`, `crp`. Values for `height` are only required for `sbp` and `dbp`.
#' @param return_input Logical. If `TRUE`, includes original input columns in the result. Defaults to `FALSE`.
#' @param return_values Character vector. Specifies which scores to compute. Options include `"percentile"`, `"z.score"`, `"MetS"`, and `"action"`.
#'
#' @return A data frame with score columns named as `<variable>_percentile`, `<variable>_z.score`, etc.
#'         Includes additional columns like `MetS` or action levels if requested.
#'
#' @details
#' The function iteratively applies `get_scores()` to each recognized variable in the input and combines the results.
#' If `"MetS"` is requested, a Metabolic Syndrome score is computed and optionally transformed to percentiles/z-scores.
#'
#' @examples
#' df <- data.frame(
#'   sex = c("f", "m"),
#'   age = c(8, 15),
#'   height = c(120, 125),
#'   waist = c(55, 60),
#'   homa = c(1.2, 1.4),
#'   sbp = c(100, 105),
#'   dbp = c(65, 70),
#'   trg = c(0.9, 1.0),
#'   crp = c(6, 2),
#'   trg = c(0.9, 1.0),
#'   hdl = c(1.1, 1.0)
#'   )
#'
#' df <- data.frame(
#'   sex = c("f", "m"),
#'   age = c(15, 7),
#'   crp = c(6, 2)
#'   )
#'
#' ScoreCalc(df, return_values = c("percentile", "action"))
#'
#' @export
ScoreCalc <- function(df, return_input = F, return_values=c("percentile","z.score", "MetS", "action")) {

  # names of the variables to which get_scores will be applied
  vars <- c("bmi", "glu", "hdl", "height", "homa", "insu", "trg", "waist", "sbp", "dbp", "crp")

  ## define vector with parameters necessary to calculate the MetS-score
  #necparas <- c("waist", "homa", "sbp", "dbp", "trg", "hdl")

  # calculation of percentiles and z-scores

  # MetS score computation requires z-scores
  return_values_temp <- return_values
  if ("MetS" %in% return_values_temp && !"z.score" %in% return_values_temp)
    return_values_temp <- c(return_values_temp, "z.score")

  # Generate score columns with appropriate prefixes
  score_results <- lapply(vars, function(var_name) {
    if (var_name %in% names(df)) {
      scores <- as.data.frame(get_scores(var_name, df$sex, df$age, df$height, df[[var_name]], return_values_temp))
      # Prefix the column names
      colnames(scores) <- paste0(var_name, "_", colnames(scores))
      return(scores)
    } else {
      NULL
    }
  })

  # remove NULL
  score_results <- score_results[!sapply(score_results, is.null)]

  # bind columns
  score_results <- do.call(cbind, score_results)

  if ("MetS" %in% return_values) {
    tryCatch({
      MetS <- MetSScore(score_results)
      score_results <- cbind(score_results, MetS)

      # compute percentiles and/or z.scores for "+100"-shifted MetS score
      MetS_score_results <- as.data.frame(get_scores("MetS_shifted", df$sex, df$age, df$height, MetS+100, return_values))
      # Prefix the column names
      if (length(colnames(MetS_score_results)) > 0) colnames(MetS_score_results) <- paste0("MetS", "_", colnames(MetS_score_results))
      if (length(MetS_score_results)>0) score_results <- cbind(score_results, MetS_score_results)

      # remove z-scores if not requested...
      if (!"z.score" %in% return_values)
        score_results <- score_results[!endsWith(names(score_results), "_z.score")]
    },
    error = function(e) warning("MetS score could not be computed. Required variables: waist, homa, sbp, dbp, trg, hdl")
    )
  }

  if ("action" %in% return_values) {
    score_results$sex <- df$sex
    score_results <- action_levels(score_results, append = TRUE)
    score_results[["sex"]] <- NULL
  }

  # if input requested
  if (return_input) score_results <- cbind(df, score_results)

  return(score_results)
}
