#' @title Negative Hypergeometric Confidence Intervals
#' @description Core functions to compute confidence intervals for parameters of the Negative Hypergeometric distribution.
#' @keywords internal
NULL



#' @importFrom dplyr filter mutate arrange group_by ungroup slice pull bind_rows
#' @importFrom data.table as.data.table setorder rbindlist
#' @importFrom extraDistr pnhyper
#' @importFrom stats choose


###################################
#---------------------------------#
# Preliminary Functions           #
#---------------------------------#
###################################

ngh_pmf <- function(x, N, M, m) {
  # Conditions: check supports
  if (x > (N - M) || N < M || m > M) {
    return(0)
  }

  numerator <- choose(m + x - 1, m - 1) * choose(N - m - x, M - m)
  denominator <- choose(N, M)
  result <- numerator / denominator

  return(result)
}



ngh_cdf <- function(x, N, M, m, lower_tail = TRUE) {
  # m = # total successes (unknown) - (our notation: M)
  m_pmf = M

  # n = # total failures - (our notation: X = N - M)
  n_pmf = N - m_pmf

  # r = # fixed successes (our notation: m)
  r_pmf = m

  # x = # balls being drawn (our notation: n = m + x)
  x_pmf = r_pmf + x

  return(pnhyper(q = x_pmf, m = m_pmf, n = n_pmf, r = r_pmf, lower.tail = lower_tail))
}


sum_ngh_pmf <- function(N, M, m, min_x, max_x) {
  sum_pmf = 0
  for (x in min_x:max_x) {
    if (x >= 0 && x <= (N - M)) {
      sum_pmf = sum_pmf + ngh_pmf(x, N, M, m)
    }
  }
  return(sum_pmf)
}



###################################
#---------------------------------#
# M Unknown                       #
#---------------------------------#
###################################

###################################
#---------------------------------#
# Confidence Interval Function    #
#---------------------------------#
###################################

#' Compute Confidence Intervals for M (Unknown Successes)
#'
#' Wrapper function to compute confidence intervals for the number of successes (M) in the Negative Hypergeometric Distribution using different methods.
#'
#' @param N Total population size.
#' @param m Fixed number of successes to be observed.
#' @param conf_level Confidence level for the interval (default = 0.95).
#' @param method Method to compute the interval. Choices are "Analog to Clopper-Pearson", "CMC", "CG", or "Blaker".
#'
#' @return A data frame containing the confidence intervals for each observed number of failures observed before the mth success (x).
#'
#' @note
#' The function assumes \code{m \geq 1} and \code{m \leq N}.
#' Confidence intervals are computed for all appropriate values of the number of observed failures (\code{x}) based on the specified \code{N} and \code{m}.
#' The full set of valid \code{x} values is output without the need for additional parameter tuning.
#'
#' @examples
#' # Example 1: Using the Analog to Clopper-Pearson method
#' nhgCI_M(N = 30, m = 5, conf_level = 0.95, method = "Analog to Clopper-Pearson")
#'
#' # Example 2: Using the CMC method
#' nhgCI_M(N = 20, m = 4, conf_level = 0.90, method = "CMC")
#'
#' # Example 3: Using the CG method
#' nhgCI_M(N = 15, m = 3, conf_level = 0.99, method = "CG")
#'
#' # Example 4: Using the Blaker method
#' nhgCI_M(N = 25, m = 7, conf_level = 0.95, method = "Blaker")
#'
#' @export
nhgCI_M <- function(N, m, conf_level = 0.95, method = "Analog to Clopper-Pearson") {
  method <- match.arg(method, choices = c("Analog to Clopper-Pearson", "CMC", "CG", "Blaker"))

  if (m > N) {
    stop("m must be less than or equal to N.")
  }

  if (method == "Analog to Clopper-Pearson") {
    return(CI_cov_prob_vec(N, m, conf_level))
  } else if (method == "CMC") {
    return(cmc_ci_vec(N, m, conf_level))
  } else if (method == "CG") {
    return(minimal_cardinality_ci_vec(N, m, conf_level, procedure = "CG"))
  } else if (method == "Blaker") {
    return(blaker_ci_vec(N, m, conf_level))
  } else {
    stop("Invalid method choice. Must be one of 'Analog to Clopper-Pearson', 'CMC', 'CG', or 'Blaker'.")
  }
}



###################################
#---------------------------------#
# Analog to Clopper Pearson       #
#---------------------------------#
###################################

CI_cov_prob_vec <- function(N, m, conf_level = 0.95) {
  target_probability <- (1 - conf_level) / 2

  # Pre-allocate storage
  x_values      <- seq.int(0, N)
  lower_bounds  <- numeric(length(x_values))
  upper_bounds  <- numeric(length(x_values))

  # Loop once over all x
  for (i in seq_along(x_values)) {
    xi <- x_values[i]

    if (xi >= (N - m) && xi <= N) {
      # Special case
      lower_bounds[i] <- N - xi
      upper_bounds[i] <- N - xi
    } else {
      # Defaults
      lb <- m
      ub <- N

      # Find lower bound
      for (M_val in seq.int(m, N)) {
        area_left <- ngh_cdf(x = xi, N = N, M = M_val, m = m, lower_tail = TRUE)
        if (isTRUE(all.equal(area_left, target_probability)) ||
            (area_left > target_probability)) {
          lb <- M_val
          break
        }
      }

      # Find upper bound
      for (M_val in seq.int(N, m, by = -1)) {
        area_right <- ngh_cdf(x = xi - 1, N = N, M = M_val, m = m, lower_tail = FALSE)
        if (isTRUE(all.equal(area_right, target_probability)) ||
            (area_right > target_probability)) {
          ub <- M_val
          break
        }
      }

      lower_bounds[i] <- lb
      upper_bounds[i] <- ub
    }
  }

  # Build the final data frame
  results <- data.frame(
    x           = x_values,
    lower_bound = lower_bounds,
    upper_bound = upper_bounds
  )
  return(results)
}



###################################
#---------------------------------#
# Minimal Cardinality             #
#---------------------------------#
###################################

all_mc_ac_vec <- function(N, m, conf_level = 0.95) {
  # Initialize the final results
  results <- data.frame(
    M              = integer(),
    a              = integer(),
    b              = integer(),
    cardinality    = integer(),
    coverage_prob  = numeric()
  )

  # Start with these constraints
  min_a <- 0
  min_b <- 0

  # Loop over M from N down to 0
  for (M_val in seq(N, 0, by = -1)) {

    # We'll build temp_results for the current M_val
    temp_results <- data.frame(
      M              = integer(),
      a              = integer(),
      b              = integer(),
      cardinality    = integer(),
      coverage_prob  = numeric()
    )

    # Special End Case: When M < m
    if (M_val < m) {
      a_val <- N - M_val
      b_val <- N - M_val

      coverage_prob <- sum_ngh_pmf(N, M_val, m, a_val, b_val)
      cardinality   <- b_val - a_val + 1

      temp_results <- rbind(
        temp_results,
        data.frame(
          M              = M_val,
          a              = a_val,
          b              = b_val,
          cardinality    = cardinality,
          coverage_prob  = coverage_prob
        )
      )

    } else {
      # If M >= m, we build all (a,b) pairs with:
      #   a in [min_a, N - M_val]
      #   b in [min_b, N - M_val]
      #   b >= a
      # Then compute coverage in a vectorized fashion.

      ab_grid <- expand.grid(
        a = seq.int(min_a, N - M_val),
        b = seq.int(min_b, N - M_val)
      ) %>%
        filter(b >= a)

      # If there are no (a,b) pairs, skip
      if (nrow(ab_grid) > 0) {
        coverage_vec <- mapply(
          FUN = function(a_val, b_val) {
            sum_ngh_pmf(N, M_val, m, a_val, b_val)
          },
          ab_grid$a,
          ab_grid$b
        )

        temp_results <- data.frame(
          M              = M_val,
          a              = ab_grid$a,
          b              = ab_grid$b,
          cardinality    = ab_grid$b - ab_grid$a + 1,
          coverage_prob  = coverage_vec
        )
      }
    }

    # If M_val >= m, we apply the coverage/confidence filter
    if (M_val >= m) {
      temp_results <- temp_results %>%
        filter(coverage_prob >= conf_level & coverage_prob >= 0 & coverage_prob <= 1) %>%
        group_by(M) %>%
        slice_min(order_by = cardinality, with_ties = TRUE) %>%
        ungroup()

      # Update min_a, min_b if we found valid intervals
      if (nrow(temp_results) > 0) {
        min_a <- max(min_a, min(temp_results$a))
        min_b <- max(min_b, min(temp_results$b))
      }
    }

    # Append current iteration’s data to the final results
    results <- rbind(results, temp_results)
  }

  # Separate results into M >= m and M < m
  results_M_ge_m <- results %>% filter(M >= m)
  results_M_lt_m <- results %>% filter(M < m)

  # Combine and order
  filtered_results <- bind_rows(results_M_ge_m, results_M_lt_m) %>%
    arrange(desc(M))

  # Add x_set column
  filtered_results <- filtered_results %>%
    mutate(x_set = paste(a, b, sep = "-"))

  return(filtered_results)
}


minimal_cardinality_ci_vec <- function(N, m, conf_level = 0.95, procedure = "MST") {
  # Choose which minimal cardinality procedure
  if (procedure == "MST") {
    results <- mst_ac_vec(N, m, conf_level)
  } else if (procedure == "CG") {
    results <- cg_ac_vec(N, m, conf_level)
  } else if (procedure == "BK") {
    results <- bk_ac_vec(N, m, conf_level)
  } else {
    stop("Invalid procedure. Choose from 'MST', 'CG', or 'BK'.")
  }

  # We'll store rows in a list, then bind them once at the end
  x_values <- seq.int(0, N)
  row_list <- vector("list", length(x_values))

  for (i in seq_along(x_values)) {
    x_val <- x_values[i]

    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())

    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_ub <- first_occurrence$M
      ci_lb <- last_occurrence$M
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")

      row_list[[i]] <- data.frame(
        x     = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci    = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine all non-NULL rows
  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
  return(ci_results)
}



###################################
#---------------------------------#
# CG                              #
#---------------------------------#
###################################

cg_ac_vec <- function(N, m, conf_level = 0.95) {
  results <- all_mc_ac_vec(N, m, conf_level)

  row_list <- vector("list", length = N + 1)
  idx <- 1

  min_a <- 0
  min_b <- 0

  for (current_M in seq.int(N, 0, by = -1)) {
    subset_results <- results %>% filter(M == current_M)

    if (nrow(subset_results) == 1) {
      chosen_row <- subset_results
    } else {
      chosen_row <- subset_results %>%
        filter(a >= min_a, b >= min_b) %>%
        arrange(a, b) %>%
        slice(1)
    }

    if (nrow(chosen_row) > 0) {
      min_a <- max(min_a, chosen_row$a)
      min_b <- max(min_b, chosen_row$b)

      row_list[[idx]] <- chosen_row
      idx <- idx + 1
    }
  }

  final_results <- do.call(rbind, row_list[1:(idx - 1)])
  final_results <- final_results %>% arrange(desc(M))
  return(final_results)
}



###################################
#---------------------------------#
# Blaker                          #
#---------------------------------#
###################################

blaker_ac_vec <- function(N, m, conf_level = 0.95) {
  alpha <- 1 - conf_level

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1) FIRST PASS: Build 'results' with (M, x, min_tail_prob, pmf_x, acceptance_set)
  #    in a more vectorized way.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  results <- data.frame(
    M              = integer(),
    x              = integer(),
    min_tail_prob  = numeric(),
    pmf_x          = numeric(),   # we'll store PMF here to avoid repeated calls later
    acceptance_set = character()
  )

  for (M_val in seq.int(0, N)) {

    # If M_val < m, there's exactly one x = N - M_val, with no min_tail_prob needed
    if (M_val < m) {
      x_val <- N - M_val
      # Append a single row
      results <- rbind(
        results,
        data.frame(
          M              = M_val,
          x              = x_val,
          min_tail_prob  = NA,   # as in your original code
          pmf_x          = NA,   # won't be used if M_val < m
          acceptance_set = as.character(x_val)
        )
      )

    } else {
      # M_val >= m => we consider x from 0 to (N - M_val)
      x_seq <- seq.int(0, N - M_val)

      # Vectorized calls to ngh_cdf() for "left" and "right" tails:
      area_left_vec <- sapply(
        x_seq,
        function(xx) ngh_cdf(x = xx,    N = N, M = M_val, m = m, lower_tail = TRUE )
      )
      area_right_vec <- sapply(
        x_seq,
        function(xx) ngh_cdf(x = xx - 1, N = N, M = M_val, m = m, lower_tail = FALSE)
      )

      # Vector of min_tail_prob
      min_tail_prob_vec <- pmin(area_left_vec, area_right_vec)

      # Precompute PMF for each x as well
      pmf_vec <- sapply(
        x_seq,
        function(xx) ngh_pmf(x = xx, N = N, M = M_val, m = m)
      )

      # Build a temporary data frame for all x in [0, N - M_val]
      tmp_df <- data.frame(
        M              = M_val,
        x              = x_seq,
        min_tail_prob  = min_tail_prob_vec,
        pmf_x          = pmf_vec,
        acceptance_set = NA_character_
      )

      # Bind once per M_val
      results <- rbind(results, tmp_df)
    }
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2) SECOND PASS: Build final_results.
  #    We do exactly what your original code did, but we skip repeated PMF calls,
  #    because 'results' already has 'pmf_x'.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  final_results <- data.frame(
    M               = integer(),
    acceptance_set  = character(),
    a               = integer(),
    b               = integer(),
    cardinality     = integer(),
    coverage_prob   = numeric(),
    x_set           = character(),
    gap             = logical()
  )

  # We iterate from M=0 up to M=N, same as your code
  for (fixed_M in seq.int(0, N)) {

    acceptance_set <- c()
    x_set_str      <- ""

    if (fixed_M < m) {
      # The acceptance set is exactly { N - M }
      acceptance_set_num <- N - fixed_M

      a_val         <- acceptance_set_num
      b_val         <- acceptance_set_num
      cardinality   <- 1
      coverage_prob <- 0  # same as your code
      gap_val       <- FALSE
      acceptance_set_str <- as.character(acceptance_set_num)
      x_set_str     <- paste(a_val, b_val, sep = "-")

    } else {
      # M >= m
      # Subset 'results' to the rows for this M
      sub_df <- results %>% filter(M == fixed_M)

      # For each x in [0, N - M], run the original logic
      # "for (fixed_x in 0:(N - fixed_M)) { ... }"
      x_seq <- seq.int(0, N - fixed_M)

      for (fixed_x in x_seq) {
        # min_tail_prob for this x
        min_tail_prob_fixed_x <- sub_df$min_tail_prob[sub_df$x == fixed_x]

        # All x' with min_tail_prob <= min_tail_prob_fixed_x
        # (the original code calls 'results %>% filter(...)', we do it here in memory)
        possible_xs_df <- sub_df[sub_df$min_tail_prob <= min_tail_prob_fixed_x, ]

        # Sum pmf_x for these x'
        prob_sum <- sum(possible_xs_df$pmf_x, na.rm = TRUE)

        # If prob_sum > alpha, we add 'fixed_x' to acceptance_set
        if (prob_sum > alpha) {
          acceptance_set <- c(acceptance_set, fixed_x)
        }
      }

      acceptance_set <- unique(as.numeric(acceptance_set))  # just in case
      acceptance_set_str <- paste(acceptance_set, collapse = ",")

      # Calculate coverage_prob by summing precomputed pmf_x
      # Only for those x in acceptance_set
      coverage_prob <- 0
      if (length(acceptance_set) > 0) {
        # We'll match again on sub_df$x
        coverage_prob <- sum(
          sub_df$pmf_x[sub_df$x %in% acceptance_set],
          na.rm = TRUE
        )
      }

      # Build a, b, cardinality
      a_val       <- min(acceptance_set)
      b_val       <- max(acceptance_set)
      cardinality <- length(acceptance_set)

      # Check for a "gap"
      gap_val <- any(diff(sort(acceptance_set)) > 1)

      # Construct x_set
      if (!gap_val) {
        # No gap => single interval
        x_set_str <- paste(a_val, b_val, sep = "-")
      } else {
        # There is a gap => multiple intervals
        acceptance_set_sorted <- sort(acceptance_set)
        intervals <- c()
        start_int <- acceptance_set_sorted[1]

        for (i in seq.int(2, length(acceptance_set_sorted))) {
          if (acceptance_set_sorted[i] != acceptance_set_sorted[i - 1] + 1) {
            intervals <- c(intervals, paste(start_int, acceptance_set_sorted[i - 1], sep = "-"))
            start_int <- acceptance_set_sorted[i]
          }
        }
        # Add the last interval
        intervals <- c(intervals, paste(start_int, acceptance_set_sorted[length(acceptance_set_sorted)], sep = "-"))
        x_set_str <- paste(intervals, collapse = ", ")
      }
    }

    # Append to final_results
    final_results <- rbind(
      final_results,
      data.frame(
        M               = fixed_M,
        acceptance_set  = acceptance_set_str,
        a               = a_val,
        b               = b_val,
        cardinality     = cardinality,
        coverage_prob   = coverage_prob,
        x_set           = x_set_str,
        gap             = gap_val
      )
    )
  }

  # Finally, arrange in descending order of M
  final_results <- final_results %>%
    arrange(desc(M))

  return(final_results)
}


blaker_ci_vec <- function(N, m, conf_level = 0.95) {
  results <- blaker_ac_vec(N, m, conf_level)

  # If there's any gap, just return the message
  if (any(results$gap)) {
    return("Gaps present in acceptance sets")
  }

  x_values <- seq.int(0, N)
  row_list <- vector("list", length(x_values))

  for (i in seq_along(x_values)) {
    x_val <- x_values[i]

    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())

    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_ub <- first_occurrence$M
      ci_lb <- last_occurrence$M
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")

      row_list[[i]] <- data.frame(
        x     = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci    = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }

  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
  return(ci_results)
}




###################################
#---------------------------------#
# CMC                             #
#---------------------------------#
###################################

cmc_ac_vec <- function(N, m, conf_level = 0.95) {
  # Initialize the final results
  results <- data.frame(
    M              = integer(),
    a              = integer(),
    b              = integer(),
    cardinality    = integer(),
    coverage_prob  = numeric()
  )

  # Start with these constraints
  min_a <- 0
  min_b <- 0

  # Loop over M from N down to 0
  for (M_val in seq.int(N, 0, by = -1)) {
    # A temporary data frame for the current M_val
    temp_results <- data.frame(
      M              = integer(),
      a              = integer(),
      b              = integer(),
      cardinality    = integer(),
      coverage_prob  = numeric()
    )

    # Special End Case: M_val < m => we have exactly (a, b) = (N - M_val, N - M_val)
    if (M_val < m) {
      a_val <- N - M_val
      b_val <- a_val

      coverage_val  <- sum_ngh_pmf(N, M_val, m, a_val, b_val)
      cardinality   <- b_val - a_val + 1

      temp_results <- rbind(
        temp_results,
        data.frame(
          M              = M_val,
          a              = a_val,
          b              = b_val,
          cardinality    = cardinality,
          coverage_prob  = coverage_val
        )
      )

    } else {
      # M_val >= m => create a grid of all valid (a, b) pairs with:
      #   a in [min_a, N - M_val],
      #   b in [min_b, N - M_val],
      #   b >= a (to ensure non-decreasing).

      ab_grid <- expand.grid(
        a = seq.int(min_a, N - M_val),
        b = seq.int(min_b, N - M_val)
      ) %>%
        filter(b >= a)

      # If there's no valid (a, b), skip
      if (nrow(ab_grid) > 0) {
        # Vectorized coverage probability computation using mapply
        coverage_vec <- mapply(
          FUN = function(a_val, b_val) {
            sum_ngh_pmf(N, M_val, m, a_val, b_val)
          },
          ab_grid$a,
          ab_grid$b
        )

        # Build temp_results for all (a, b)
        temp_results <- data.frame(
          M              = M_val,
          a              = ab_grid$a,
          b              = ab_grid$b,
          cardinality    = ab_grid$b - ab_grid$a + 1,
          coverage_prob  = coverage_vec
        )
      }
    }

    # Filter out the sets with coverage probability >= conf_level
    # Among those, choose the acceptance curve with the highest a and the lowest b
    # (the same as "the one with the lowest coverage_prob that is still >= conf_level"
    # in your description).
    if (M_val >= m && nrow(temp_results) > 0) {
      temp_results <- temp_results %>%
        filter(coverage_prob >= conf_level & coverage_prob >= 0 & coverage_prob <= 1) %>%
        group_by(M) %>%
        filter(a == max(a)) %>%
        filter(b == min(b)) %>%
        ungroup()

      # Update min_a and min_b if we found any valid intervals
      if (nrow(temp_results) > 0) {
        min_a <- max(min_a, min(temp_results$a))
        min_b <- max(min_b, min(temp_results$b))
      }
    }

    # Append the current iteration’s data to the final results
    results <- rbind(results, temp_results)
  }

  # Arrange in descending order of M
  filtered_results <- results %>%
    arrange(desc(M))

  # Add a column "x_set" to show "a-b"
  filtered_results <- filtered_results %>%
    mutate(x_set = paste(a, b, sep = "-"))

  return(filtered_results)
}


cmc_ci_vec <- function(N, m, conf_level = 0.95) {
  results <- cmc_ac_vec(N, m, conf_level)

  x_values <- seq.int(0, N)
  row_list <- vector("list", length(x_values))

  for (i in seq_along(x_values)) {
    x_val <- x_values[i]

    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())

    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_ub <- first_occurrence$M
      ci_lb <- last_occurrence$M
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")

      row_list[[i]] <- data.frame(
        x     = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci    = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }

  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
  return(ci_results)
}




###################################
#---------------------------------#
# N Unknown                       #
#---------------------------------#
###################################

###################################
#---------------------------------#
# Confidence Interval Function    #
#---------------------------------#
###################################

#' Compute Confidence Intervals for N (Unknown Population Size)
#'
#' Wrapper function to compute confidence intervals for the population size (N) in the Negative Hypergeometric Distribution using different methods.
#'
#' @param M Number of successes in population.
#' @param m Fixed number of successes to be observed.
#' @param conf_level Confidence level for the interval (default = 0.95).
#' @param method Method to compute the interval. Choices are "Analog to Clopper-Pearson", "CMC", "CG", or "Blaker".
#' @param max_N Maximum value of N to consider (default = 250).
#'
#' @return A data frame containing the confidence intervals for each observed number of failures observed before the mth success (x).
#'
#' @note
#' The function assumes \code{m \geq 1} and \code{m \leq M}.
#' Note that the parameter \code{N} (population size) is unbounded above in the Negative Hypergeometric model.
#' As a result, the number of observed failures (\code{x}) is theoretically unbounded.
#' The \code{max_N} parameter sets a practical upper limit for \code{N} during computations.
#'
#' If the desired value of \code{x} does not appear in the output, try increasing \code{max_N}.
#' Be aware that increasing \code{max_N} may substantially increase computation time, especially for large \code{m} or when high accuracy is required.
#'
#' @examples
#' # Example 1: Using the Analog to Clopper-Pearson method
#' nhgCI_N(M = 10, m = 5, conf_level = 0.95, method = "Analog to Clopper-Pearson")
#'
#' # Example 2: Using the CMC method
#' nhgCI_N(M = 10, m = 5, conf_level = 0.95, method = "CMC")
#'
#' # Example 3: Using the CG method
#' nhgCI_N(M = 10, m = 5, conf_level = 0.95, method = "CG")
#'
#' # Example 4: Using the Blaker method
#' nhgCI_N(M = 10, m = 5, conf_level = 0.95, method = "Blaker")
#'
#' @export
nhgCI_N <- function(M, m, conf_level = 0.95, method = "Analog to Clopper-Pearson", max_N = 250) {
  method <- match.arg(method, choices = c("Analog to Clopper-Pearson", "CMC", "CG", "Blaker"))

  if (m > M) {
    stop("m must be less than or equal to M.")
  }

  if (method == "Analog to Clopper-Pearson") {
    return(CI_Analog_CP_N_Unknown_vec(M, m, conf_level, max_N))
  } else if (method == "CMC") {
    return(cmc_ci_N_unkown_vec(M, m, conf_level, max_N))
  } else if (method == "CG") {
    return(minimal_cardinality_ci_N_unkown_vec(M, m, conf_level, max_N, procedure = "CG"))
  } else if (method == "Blaker") {
    return(blaker_ci_N_unkown_vec(M, m, conf_level, max_N))
  } else {
    stop("Invalid method choice. Must be one of 'Analog to Clopper-Pearson', 'CMC', 'CG', or 'Blaker'.")
  }
}





###################################
#---------------------------------#
# Analog to Clopper Pearson       #
#---------------------------------#
###################################

CI_Analog_CP_N_Unknown_vec <- function(M, m, conf_level = 0.95, max_N = 1000) {
  target_probability <- (1 - conf_level) / 2
  max_x <- max_N - M

  x_values <- seq.int(0, max_x)
  lower_bounds <- rep(NA, length(x_values))
  upper_bounds <- rep(NA, length(x_values))

  previous_upper_bound <- 0

  for (i in seq_along(x_values)) {
    xi <- x_values[i]

    # Initialize lower_bound
    lb <- M + xi
    # Find lower bound
    for (N_val in seq.int(M + xi, M + max_x)) {
      area_right <- ngh_cdf(x = xi - 1, N = N_val, M = M, m = m, lower_tail = FALSE)

      if (isTRUE(all.equal(area_right, target_probability)) || (area_right > target_probability)) {
        lb <- N_val
        break
      }
    }

    # Initialize upper_bound
    ub <- lb
    # Find upper bound
    for (N_val in seq.int(lb, M + max_x)) {
      area_left <- ngh_cdf(x = xi, N = N_val, M = M, m = m, lower_tail = TRUE)

      if (isTRUE(all.equal(area_left, target_probability))) {
        ub <- N_val
        break
      } else if (area_left < target_probability) {
        ub <- N_val - 1
        break
      }
    }

    # Stop if the upper bound starts decreasing
    if (ub < previous_upper_bound) {
      # Mark the rest as NA and break
      break
    }
    previous_upper_bound <- ub

    lower_bounds[i] <- lb
    upper_bounds[i] <- ub
  }

  # Build final results
  results <- data.frame(
    x = x_values,
    lower_bound = lower_bounds,
    upper_bound = upper_bounds
  )

  # Filter out rows where upper_bound is NA
  results <- results[!is.na(results$upper_bound), ]
  return(results)
}



coverage_prob_ACP_N_unknown_vec <- function(M, N, m, conf_level = 0.95, max_N = 1000) {
  found_N_in_last_CI <- TRUE

  # Keep increasing max_N until N is no longer in the last CI
  while (found_N_in_last_CI) {
    ci_results <- CI_Analog_CP_N_Unknown_vec(M, m, conf_level, max_N)
    last_x_ci <- ci_results[nrow(ci_results), ]

    if (N >= last_x_ci$lower_bound && N <= last_x_ci$upper_bound) {
      max_N <- max_N + 100
    } else {
      found_N_in_last_CI <- FALSE
    }
  }

  # Then compute final coverage probability
  ci_results <- CI_Analog_CP_N_Unknown_vec(M, m, conf_level, max_N)

  covered_x <- ci_results %>%
    dplyr::filter(lower_bound <= N & upper_bound >= N) %>%
    dplyr::pull(x)

  if (length(covered_x) == 0) {
    return(data.frame(N = N, coverage_prob = NA, min_x = NA, max_x = NA))
  }

  min_x <- min(covered_x)
  max_x <- max(covered_x)

  total_prob <- sum(
    sapply(covered_x, function(x_val) ngh_pmf(x_val, N, M, m))
  )

  return(data.frame(N = N, coverage_prob = total_prob, min_x = min_x, max_x = max_x))
}




###################################
#---------------------------------#
# Minimal Cardinality             #
#---------------------------------#
###################################

minimal_cardinality_ci_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 1000, procedure = "MST") {
  # Choose the procedure.
  if (procedure == "MST") {
    results <- mst_ac_N_unknown_direct(M, m, conf_level, max_N)
  } else if (procedure == "CG") {
    results <- cg_ac_N_unknown_direct(M, m, conf_level, max_N)
  } else if (procedure == "BK") {
    results <- bk_ac_N_unknown_direct(M, m, conf_level, max_N)
  } else {
    stop("Invalid procedure. Choose from 'MST', 'CG', or 'BK'.")
  }

  dt_results <- as.data.table(results)
  setorder(dt_results, N)  # Ensure results are ordered by N

  # Determine max_x as the second-highest a (if available) or the highest.
  unique_a_values <- sort(unique(dt_results$a), decreasing = TRUE)
  if (length(unique_a_values) > 1) {
    max_x <- unique_a_values[2]
  } else {
    max_x <- unique_a_values[1]
  }

  ci_list <- vector("list", length = max_x + 1)
  idx <- 1
  for (x in 0:max_x) {
    # Find all intervals that contain x.
    subset_dt <- dt_results[a <= x & b >= x]
    if (nrow(subset_dt) == 0) next  # Skip if x is not covered.

    # The first occurrence gives the lowest N and the last gives the highest N.
    first_occurrence <- subset_dt[1]
    last_occurrence <- subset_dt[.N]

    ci_lb <- first_occurrence$N
    ci_ub <- if (m == 1) Inf else last_occurrence$N
    ci_str <- if (is.infinite(ci_ub)) paste0("[", ci_lb, ", ∞)") else paste0("[", ci_lb, ", ", ci_ub, "]")

    ci_list[[idx]] <- data.table(x = x, ci_lb = ci_lb, ci_ub = ci_ub, ci = ci_str)
    idx <- idx + 1
  }

  if (idx > 1) {
    ci_dt <- rbindlist(ci_list[1:(idx - 1)])
  } else {
    ci_dt <- data.table(x = integer(), ci_lb = numeric(), ci_ub = numeric(), ci = character())
  }

  return(as.data.frame(ci_dt))
}







###################################
#---------------------------------#
# CG                              #
#---------------------------------#
###################################

cg_ac_N_unknown_direct <- function(M, m, conf_level = 0.95, max_N = 1000) {
  # This function finds minimal‐cardinality acceptance curves following the CG procedure.
  # For each N from M to max_N, it searches for an (a, b) (with b = a + k - 1, where k is the candidate cardinality)
  # such that the coverage probability (computed by sum_ngh_pmf) meets the conf_level.
  # Among all candidates at the minimal viable cardinality, it chooses the one with the largest a, then largest b.
  # In addition, a and b are not allowed to decrease across N.

  results <- list()

  # Initialize state variables
  min_a <- 0
  min_b <- 0
  previous_cardinality <- NA  # for N=M, we start with candidate_card = 1

  for (current_N in M:max_N) {
    max_x <- current_N - M  # possible a, b range is 0:max_x
    candidate_card <- if (!is.na(previous_cardinality)) previous_cardinality else 1
    candidate_found <- FALSE
    best_candidate <- NULL

    # Loop over candidate cardinalities until a viable acceptance curve is found.
    while (candidate_card <= (max_x + 1) && !candidate_found) {
      # For candidate cardinality k, ensure a >= min_a and also a >= min_b - (k - 1) so that b >= min_b.
      a_min <- max(min_a, min_b - candidate_card + 1)
      a_max <- max_x - (candidate_card - 1)

      if (a_min > a_max) {
        candidate_card <- candidate_card + 1
        next
      }

      # For CG, choose the candidate with the largest a; if tied, choose the one with largest b.
      for (a_val in a_min:a_max) {
        b_val <- a_val + candidate_card - 1
        if (b_val < min_b) next

        cov_prob <- sum_ngh_pmf(current_N, M, m, a_val, b_val)
        if (cov_prob >= conf_level && cov_prob <= 1 && cov_prob >= 0) {
          if (is.null(best_candidate)) {
            best_candidate <- list(N = current_N, a = a_val, b = b_val,
                                   cardinality = candidate_card, coverage_prob = cov_prob)
          } else {
            # For CG, prefer larger a; if a is equal then prefer larger b.
            if (a_val > best_candidate$a ||
                (a_val == best_candidate$a && b_val > best_candidate$b)) {
              best_candidate <- list(N = current_N, a = a_val, b = b_val,
                                     cardinality = candidate_card, coverage_prob = cov_prob)
            }
          }
        }
      }

      if (!is.null(best_candidate)) {
        candidate_found <- TRUE
      } else {
        candidate_card <- candidate_card + 1
      }
    }  # end while over candidate_card

    if (!is.null(best_candidate)) {
      results[[length(results) + 1]] <- best_candidate
      min_a <- best_candidate$a
      min_b <- best_candidate$b
      previous_cardinality <- candidate_card
    }
  }  # end loop over N

  if (length(results) == 0) {
    return(data.frame(
      N = integer(),
      a = integer(),
      b = integer(),
      cardinality = integer(),
      coverage_prob = numeric(),
      x_set = character(),
      stringsAsFactors = FALSE
    ))
  }

  dt <- rbindlist(results, use.names = TRUE, fill = TRUE)
  dt[, x_set := paste(a, b, sep = "-")]
  setcolorder(dt, c("N", "a", "b", "cardinality", "coverage_prob", "x_set"))

  return(as.data.frame(dt))
}






###################################
#---------------------------------#
# Blaker                          #
#---------------------------------#
###################################

blaker_ac_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 250) {
  alpha <- 1 - conf_level

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 1) FIRST PASS: Build a data frame 'results' with:
  #    (N, x, min_tail_prob, pmf_x)
  #    Instead of looping repeatedly with rbind, we:
  #    - Create combinations for N in [M, max_N] and x in [0, N - M].
  #    - Compute min_tail_prob and pmf_x in a vectorized fashion.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  N_values <- seq.int(M, max_N)

  # Build all (N, x) pairs in one go
  combos <- do.call(
    rbind,
    lapply(N_values, function(N_val) {
      data.frame(
        N = N_val,
        x = seq.int(0, N_val - M)
      )
    })
  )

  # Compute min_tail_prob for each (N, x) pair
  # and also store pmf_x to avoid repeated calls later.
  # We'll do this in two steps for clarity.
  combos <- combos %>%
    rowwise() %>%
    mutate(
      area_left  = ngh_cdf(x, N, M, m, lower_tail = TRUE),
      area_right = ngh_cdf(x - 1, N, M, m, lower_tail = FALSE),
      min_tail_prob = min(area_left, area_right),
      pmf_x         = ngh_pmf(x, N, M, m)
    ) %>%
    ungroup()

  # We'll keep the "acceptance_set" column as NA here (to match your original structure)
  # though it's not used until the second pass.
  results <- combos %>%
    mutate(acceptance_set = NA_character_) %>%
    select(N, x, min_tail_prob, pmf_x, acceptance_set)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 2) SECOND PASS: Build 'final_results' by scanning each fixed N
  #    and determining acceptance sets. The logic is unchanged.
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  final_results <- data.frame(
    N              = integer(),
    acceptance_set = character(),
    a              = integer(),
    b              = integer(),
    cardinality    = integer(),
    coverage_prob  = numeric(),
    x_set          = character(),
    gap            = logical()
  )

  for (fixed_N in seq.int(M, max_N)) {
    # Acceptance set for this N
    acceptance_set <- integer(0)

    # Subset results for this N just once
    sub_df <- results %>% filter(N == fixed_N)

    # For each x in [0, fixed_N - M], check the condition
    x_seq <- seq.int(0, fixed_N - M)
    for (fixed_x in x_seq) {
      # min_tail_prob for x
      min_tp_fx <- sub_df$min_tail_prob[sub_df$x == fixed_x]

      # All x's with min_tail_prob <= min_tail_prob_fixed_x
      possible_xs <- sub_df$x[sub_df$min_tail_prob <= min_tp_fx]

      # Sum up the pmf of all these x's
      pmf_sum <- sum(sub_df$pmf_x[sub_df$x %in% possible_xs])

      # If sum of pmf is > alpha, include fixed_x in acceptance set
      if (pmf_sum > alpha) {
        acceptance_set <- c(acceptance_set, fixed_x)
      }
    }

    # Convert to numeric (unique), build acceptance_set string
    acceptance_set <- sort(unique(as.numeric(acceptance_set)))
    acceptance_set_str <- paste(acceptance_set, collapse = ",")

    # a, b, cardinality
    a_val         <- min(acceptance_set)
    b_val         <- max(acceptance_set)
    cardinality   <- length(acceptance_set)

    # Coverage probability: sum pmf_x for x in acceptance_set
    coverage_prob <- sum(sub_df$pmf_x[sub_df$x %in% acceptance_set])

    # Check if there is a gap
    gap_val <- any(diff(acceptance_set) > 1)

    # Build x_set (interval notation)
    x_set_str <- ""
    if (!gap_val) {
      # No gap => single interval
      x_set_str <- paste(a_val, b_val, sep = "-")
    } else {
      # Gap => multiple intervals
      intervals <- c()
      start_int <- acceptance_set[1]
      for (i in seq.int(2, length(acceptance_set))) {
        if (acceptance_set[i] != acceptance_set[i - 1] + 1) {
          intervals <- c(intervals, paste(start_int, acceptance_set[i - 1], sep = "-"))
          start_int <- acceptance_set[i]
        }
      }
      # Add last interval
      intervals <- c(intervals, paste(start_int, acceptance_set[length(acceptance_set)], sep = "-"))
      x_set_str <- paste(intervals, collapse = ", ")
    }

    # Add row to final_results
    final_results <- rbind(
      final_results,
      data.frame(
        N               = fixed_N,
        acceptance_set  = acceptance_set_str,
        a               = a_val,
        b               = b_val,
        cardinality     = cardinality,
        coverage_prob   = coverage_prob,
        x_set           = x_set_str,
        gap             = gap_val
      )
    )
  }

  final_results <- final_results %>%
    arrange(N)

  return(final_results)
}


blaker_ci_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 250) {
  results <- blaker_ac_N_unkown_vec(M, m, conf_level, max_N)

  # If there's any gap, just return
  if (any(results$gap)) {
    return("Gaps present in acceptance sets")
  }

  unique_a_values <- sort(unique(results$a), decreasing = TRUE)
  max_x <- if (length(unique_a_values) > 1) unique_a_values[2] else unique_a_values[1]

  x_values <- seq.int(0, max_x)
  row_list <- vector("list", length(x_values))

  for (i in seq_along(x_values)) {
    x_val <- x_values[i]

    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())

    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_lb <- first_occurrence$N
      ci_ub <- last_occurrence$N
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")

      row_list[[i]] <- data.frame(
        x = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }

  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
  return(ci_results)
}




###################################
#---------------------------------#
# CMC                             #
#---------------------------------#
###################################

cmc_ac_N_unknown_vec <- function(M, m, conf_level = 0.95, max_N = 250) {
  # Initialize the final results
  results <- data.frame(
    N              = integer(),
    a              = integer(),
    b              = integer(),
    cardinality    = integer(),
    coverage_prob  = numeric()
  )

  # Start with these constraints
  min_a <- 0
  min_b <- 0

  # Loop over N from M up to max_N
  for (N_val in seq.int(M, max_N)) {
    max_x <- N_val - M
    # Build all (a, b) pairs (with b >= a) in a single step
    ab_grid <- expand.grid(
      a = seq.int(min_a, min(min_a + 5, max_x)),
      # a = seq.int(min_a, min_a+1),
      b = seq.int(min_b, max_x)
    ) %>%
      filter(b >= a)

    # Compute coverage probability in a vectorized manner
    if (nrow(ab_grid) > 0) {
      coverage_vec <- mapply(
        FUN = function(a_val, b_val) {
          sum_ngh_pmf(N_val, M, m, a_val, b_val)
        },
        ab_grid$a,
        ab_grid$b
      )

      temp_results <- data.frame(
        N             = N_val,
        a             = ab_grid$a,
        b             = ab_grid$b,
        cardinality   = ab_grid$b - ab_grid$a + 1,
        coverage_prob = coverage_vec
      )
    } else {
      # If there's no valid (a,b) pair, skip
      temp_results <- data.frame(
        N             = integer(),
        a             = integer(),
        b             = integer(),
        cardinality   = integer(),
        coverage_prob = numeric()
      )
    }

    # Filter out sets with coverage_prob >= conf_level
    # Then pick the acceptance curve with the highest 'a' and the lowest 'b'
    temp_results <- temp_results %>%
      filter(coverage_prob >= conf_level & coverage_prob >= 0 & coverage_prob <= 1)

    # If no rows left, skip the group_by part
    if (nrow(temp_results) > 0) {
      temp_results <- temp_results %>%
        group_by(N) %>%
        filter(a == max(a)) %>%
        filter(b == min(b)) %>%
        ungroup()
    }

    # Update min_a, min_b if we found any valid intervals
    if (nrow(temp_results) > 0) {
      min_a <- max(min_a, min(temp_results$a))
      min_b <- max(min_b, min(temp_results$b))
    }

    # Append to main results
    results <- rbind(results, temp_results)
  }

  # Arrange by N in ascending order
  filtered_results <- results %>%
    arrange(N)

  # Add x_set column "a-b"
  filtered_results <- filtered_results %>%
    mutate(x_set = paste(a, b, sep = "-"))

  return(filtered_results)
}



cmc_ci_N_unkown_vec <- function(M, m, conf_level = 0.95, max_N = 250) {
  results <- cmc_ac_N_unknown_vec(M, m, conf_level, max_N)

  unique_a_values <- sort(unique(results$a), decreasing = TRUE)
  max_x <- if (length(unique_a_values) > 1) unique_a_values[2] else unique_a_values[1]

  x_values <- seq.int(0, max_x)
  row_list <- vector("list", length(x_values))

  for (i in seq_along(x_values)) {
    x_val <- x_values[i]

    first_occurrence <- results %>% filter(a <= x_val, x_val <= b) %>% slice(1)
    last_occurrence  <- results %>% filter(a <= x_val, x_val <= b) %>% slice(n())

    if (nrow(first_occurrence) > 0 && nrow(last_occurrence) > 0) {
      ci_lb <- first_occurrence$N
      ci_ub <- last_occurrence$N
      ci_str <- paste0("[", ci_lb, ", ", ci_ub, "]")

      row_list[[i]] <- data.frame(
        x = x_val,
        ci_lb = ci_lb,
        ci_ub = ci_ub,
        ci = ci_str,
        stringsAsFactors = FALSE
      )
    }
  }

  ci_results <- do.call(rbind, row_list[!sapply(row_list, is.null)])
  return(ci_results)
}


