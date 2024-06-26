#' @title Likelihood for SecSSE model
#' Logikelihood calculation for the SecSSE model given a set of parameters and
#' data, returning also the likelihoods along the branches
#'
#' @inheritParams default_params_doc
#'
#' @return A list containing: "output", observed states along evaluated time
#' points along all branches, used for plotting. "states" all ancestral states
#' on the nodes and "duration", indicating the time taken for the total
#' evaluation
#' @examples
#' set.seed(5)
#' phy <- ape::rphylo(n = 4, birth = 1, death = 0)
#' traits <- c(0, 1, 1, 0)
#' params <- secsse::id_paramPos(c(0, 1), 2)
#' params[[1]][] <- c(0.2, 0.2, 0.1, 0.1)
#' params[[2]][] <- 0.0
#' params[[3]][, ] <- 0.1
#' diag(params[[3]]) <- NA
#'
#' secsse_loglik_eval(parameter = params,
#'                    phy = phy,
#'                    traits = traits,
#'                    num_concealed_states = 2,
#'                    sampling_fraction = c(1, 1),
#'                    num_steps = 10)
#' @export
secsse_loglik_eval <- function(parameter,
                               phy,
                               traits,
                               num_concealed_states,
                               cond = "proper_cond",
                               root_state_weight = "proper_weights",
                               sampling_fraction,
                               setting_calculation = NULL,
                               loglik_penalty = 0,
                               is_complete_tree = FALSE,
                               num_threads = 1,
                               atol = 1e-8,
                               rtol = 1e-7,
                               method = "odeint::bulirsch_stoer",
                               num_steps = 100) {
  RcppParallel::setThreadOptions(numThreads = num_threads)
  lambdas <- parameter[[1]]
  mus <- parameter[[2]]
  parameter[[3]][is.na(parameter[[3]])] <- 0
  q_matrix <- parameter[[3]]

  check_input(traits,
              phy,
              sampling_fraction,
              root_state_weight,
              is_complete_tree)
  setting_calculation <- build_initStates_time(phy,
                                               traits,
                                               num_concealed_states,
                                               sampling_fraction,
                                               is_complete_tree,
                                               mus,
                                               traitStates = 
                                                 get_trait_states(parameter,
                                                                  num_concealed_states))
  eval_cpp(rhs = if (is.list(lambdas)) "ode_cla" else "ode_standard",
           ances = setting_calculation$ances,
           states = setting_calculation$states,
           forTime = setting_calculation$forTime,
           lambdas = lambdas,
           mus = mus,
           Q = q_matrix,
           method = method,
           atol = atol,
           rtol = rtol,
           is_complete_tree = is_complete_tree,
           num_steps = num_steps)
}

#' Plot the local probability along a tree
#' 
#' Plot the local probability along the tree, including the branches
#' 
#' @details This function will evaluate the log likelihood locally along
#' all branches and plot the result. When `num_steps` is left to `NULL`, all 
#' likelihood evaluations during integration are used for plotting. This may 
#' work for not too large trees, but may become very memory heavy for larger 
#' trees. Instead, the user can indicate a number of steps, which causes the 
#' probabilities to be evaluated at a distinct amount of steps along each branch
#' (and the probabilities to be properly integrated in between these steps). 
#' This provides an approximation, but generally results look very similar to 
#' using the full evaluation.
#' The function used for `prob_func` will be highly dependent on your system.
#' for instance, for a 3 observed, 2 hidden states model, the probability
#' of state A is `prob[1] + prob[2] + prob[3]`, normalized by the row sum.
#' `prob_func` will be applied to each row of the 'states' matrix (you can thus
#' test your function on the states matrix returned when
#' `'see_ancestral_states = TRUE'`). Please note that the first N columns of the
#' states matrix are the extinction rates, and the `(N+1):2N` columns belong to
#' the speciation rates, where `N = num_obs_states * num_concealed_states`.
#' A typical `prob_func` function will look like:
#' ```
#' my_prob_func <- function(x) {
#'   return(sum(x[5:8]) / sum(x))
#' }
#' ```
#' 
#' @inheritParams default_params_doc
#' 
#' @return ggplot2 object
#' @examples
#' set.seed(5)
#' phy <- ape::rphylo(n = 4, birth = 1, death = 0)
#' traits <- c(0, 1, 1, 0)
#' params <- secsse::id_paramPos(c(0, 1), 2)
#' params[[1]][] <- c(0.2, 0.2, 0.1, 0.1)
#' params[[2]][] <- 0.0
#' params[[3]][, ] <- 0.1
#' diag(params[[3]]) <- NA
#' #  Thus, we have for both, rates
#' # 0A, 1A, 0B and 1B. If we are interested in the posterior probability of
#' # trait 0,we have to provide a helper function that sums the probabilities of
#' # 0A and 0B, e.g.:
#' helper_function <- function(x) {
#'   return(sum(x[c(5, 7)]) / sum(x)) # normalized by total sum, just in case.
#' }
#'
#' out_plot <- plot_state_exact(parameters = params,
#'                              phy = phy,
#'                              traits = traits,
#'                              num_concealed_states = 2,
#'                              sampling_fraction = c(1, 1),
#'                              num_steps = 10,
#'                              prob_func = helper_function)
#' @export
plot_state_exact <- function(parameters,
                             phy,
                             traits,
                             num_concealed_states,
                             sampling_fraction,
                             cond = "proper_cond",
                             root_state_weight = "proper_weights",
                             is_complete_tree = FALSE,
                             method = "odeint::bulirsch_stoer",
                             atol = 1e-16,
                             rtol = 1e-16,
                             num_steps = 100,
                             prob_func = NULL,
                             verbose = FALSE) {
  if (is.null(prob_func)) {
    stop("need to set a probability function, check description to how")
  }

  eval_res <- secsse_loglik_eval(parameter = parameters,
                                 phy = phy,
                                 traits = traits,
                                 num_concealed_states =
                                  num_concealed_states,
                                 cond = cond,
                                 root_state_weight = root_state_weight,
                                 num_steps = num_steps,
                                 sampling_fraction = sampling_fraction,
                                 is_complete_tree = is_complete_tree,
                                 atol = atol,
                                 rtol = rtol,
                                 method = method)

  if (verbose) message("\nconverting collected likelihoods
                       to graph positions:\n")

  xs <- ape::node.depth.edgelength(phy)
  ys <- ape::node.height(phy)
  num_tips <- length(phy$tip.label)
  num_nodes <- (1 + num_tips):length(ys)

  nodes <- data.frame(x = xs, y = ys, n = c(1:num_tips, num_nodes))

  to_plot <- eval_res$output
 
  for_plot <- collect_branches(to_plot, nodes, prob_func, verbose)

  node_bars <- collect_node_bars(to_plot, nodes, prob_func, eval_res$states)

  if (verbose) message("\ngenerating ggplot object\n")

  focal_plot <- make_ggplot(for_plot, node_bars)
  return(focal_plot)
}


#' @importFrom rlang .data
#' @keywords internal
make_ggplot <- function(for_plot, node_bars) {
  ggplot_plot <- ggplot2::ggplot(for_plot) +
    ggplot2::geom_segment(ggplot2::aes(x = .data[["x0"]],
                                       y = .data[["y"]],
                                       xend = .data[["x1"]],
                                       yend = .data[["y"]],
                                       col = .data[["prob"]])) +
    ggplot2::geom_segment(data = node_bars,
                          ggplot2::aes(x = .data[["x"]],
                                       y = .data[["y0"]],
                                       yend = .data[["y1"]],
                                       xend = .data[["x"]],
                                       col = .data[["prob"]])
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_blank())

  return(ggplot_plot)
}

#' @keywords internal
collect_branches <- function(to_plot,
                             nodes,
                             prob_func,
                             verbose) {
  num_rows <- length(to_plot[, 1])

  for_plot <- matrix(nrow = num_rows, ncol = 6)
  for_plot_cnt <- 1
  if (verbose) pb <- utils::txtProgressBar(max = length(unique(to_plot[, 1])),
                                           style = 3)
  cnt <- 1
  for (parent in unique(to_plot[, 1])) {
    if (verbose) utils::setTxtProgressBar(pb, cnt)
    cnt <- cnt + 1

    to_plot2 <- subset(to_plot, to_plot[, 1] == parent)
    for (daughter in unique(to_plot2[, 2])) {
      indices <- which(to_plot2[, 2] == daughter)
      if (length(indices) > 0) {
        # we have a branch
        focal_branch <- to_plot2[indices, ]
        start_x <- nodes$x[which(nodes$n == parent)]
        end_x <- nodes$x[which(nodes$n == daughter)]
        y <- nodes$y[which(nodes$n == daughter)]

        bl <- end_x - start_x

        probs <- apply(focal_branch[, 4:length(focal_branch[1, ])],
                       1,
                       prob_func)
        for (s in 1:(length(focal_branch[, 1]) - 1)) {
          x0 <- start_x + bl - focal_branch[s, 3]
          x1 <- start_x + bl - focal_branch[s + 1, 3]
          ps <- probs[s]
          for_plot[for_plot_cnt, ] <- c(x0, x1, y, ps, parent, daughter)
          for_plot_cnt <- for_plot_cnt + 1
        }
      }
    }
  }
  colnames(for_plot) <- c("x0", "x1", "y", "prob", "p", "d")
  for_plot <- tibble::as_tibble(for_plot)
  return(for_plot)
}

#' @keywords internal
collect_node_bars <- function(to_plot,
                              nodes,
                              prob_func,
                              states) {
  node_bars <- matrix(nrow = length(unique(to_plot[, 1])), ncol = 4)
  node_bars_cnt <- 1
  for (parent in unique(to_plot[, 1])) {
    focal_data <- subset(to_plot, to_plot[, 1] == parent)
    daughters <- unique(focal_data[, 2])
    start_x <- nodes$x[which(nodes$n == parent)]
    y <- c()
    for (i in seq_along(daughters)) {
      y <- c(y, nodes$y[nodes$n == daughters[i]])
    }
    y <- sort(y)

    probs <- states[parent, ]
    rel_prob <- prob_func(probs)
    node_bars[node_bars_cnt, ] <- c(start_x, y, rel_prob)
    node_bars_cnt <- node_bars_cnt + 1
  }

  colnames(node_bars) <- c("x", "y0", "y1", "prob")
  node_bars <- tibble::as_tibble(node_bars)
  return(node_bars)
}
