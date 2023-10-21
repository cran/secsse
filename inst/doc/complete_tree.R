## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----sim_plot_tree------------------------------------------------------------
library(secsse)

spec_matrix <- c()
spec_matrix <- rbind(spec_matrix, c(0, 0, 0, 1))
spec_matrix <- rbind(spec_matrix, c(1, 1, 1, 1))
lambda_list <- secsse::create_lambda_list(state_names = c(0, 1),
                                          num_concealed_states = 2,
                                          transition_matrix = spec_matrix,
                                          model = "CR")

mu_vector <- secsse::create_mu_vector(state_names = c(0, 1),
                                      num_concealed_states = 2,
                                      model = "CR",
                                      lambda_list = lambda_list)

shift_matrix <- c()
shift_matrix <- rbind(shift_matrix, c(0, 1, 3))
shift_matrix <- rbind(shift_matrix, c(1, 0, 4))

q_matrix <- secsse::create_q_matrix(state_names = c(0, 1),
                                    num_concealed_states = 2,
                                    shift_matrix = shift_matrix,
                                    diff.conceal = FALSE)

# Set-up starting parameters
speciation_rate <- 0.8
extinction_rate <- 0.2
q_01 <- 0.1
q_10 <- 0.1
used_params <- c(speciation_rate, extinction_rate, q_01, q_10)

sim_lambda_list <- secsse::fill_in(lambda_list, used_params)
sim_mu_vector   <- secsse::fill_in(mu_vector, used_params)
sim_q_matrix    <- secsse::fill_in(q_matrix, used_params)

# Simulate and plot the tree

sim_tree_complete <- secsse::secsse_sim(lambdas = sim_lambda_list,
                                        mus = sim_mu_vector,
                                        qs = sim_q_matrix,
                                        crown_age = 5,
                                        num_concealed_states = 2,
                                        seed = 40,
                                        drop_extinct = FALSE)

if (requireNamespace("diversitree")) {
  traits_for_plot_complete <- data.frame(
    trait = as.numeric(sim_tree_complete$obs_traits),
    row.names = sim_tree_complete$phy$tip.label
  )
  diversitree::trait.plot(tree = sim_tree_complete$phy,
                          dat = traits_for_plot_complete,
                          cols = list("trait" = c("blue", "red")),
                          type = "p")
} else {
  plot(sim_tree_complete$phy)
}


## ----fitting_model_complete_tree----------------------------------------------
idparsopt <- 1:4 # our maximum rate parameter was 4 -> We are keeping
# concealed and examined traits the same for the MLE.
idparsfix <- c(0) # we want to keep all zeros at zero
initparsopt <- rep(0.1, 4)
initparsfix <- c(0.0) # all zeros remain at zero.
sampling_fraction <- c(1, 1)

idparslist <- list()
idparslist[[1]] <- lambda_list
idparslist[[2]] <- mu_vector
idparslist[[3]] <- q_matrix

complete_tree_ml_CR <- secsse_ml(phy = sim_tree_complete$phy,
                                 traits = sim_tree_complete$obs_traits,
                                 num_concealed_states = 2,
                                 idparslist = idparslist,
                                 idparsopt = idparsopt,
                                 initparsopt = initparsopt,
                                 idparsfix = idparsfix,
                                 parsfix = initparsfix,
                                 sampling_fraction = sampling_fraction,
                                 verbose = FALSE)

## ----complete_tree_res--------------------------------------------------------
CR_par_complete <- secsse::extract_par_vals(idparslist, complete_tree_ml_CR$MLpars)
complete_tree_ml_CR
CR_par_complete
spec_rates_complete <- CR_par_complete[1]
ext_rates_complete <- CR_par_complete[2]
Q_01_complete <- CR_par_complete[3]
Q_10_complete <- CR_par_complete[4]
spec_rates_complete
ext_rates_complete
Q_01_complete
Q_10_complete

## ----fitting_ml_reconstructed_tree--------------------------------------------

sim_tree_reconstructed <- secsse::secsse_sim(lambdas = sim_lambda_list,
                                             mus = sim_mu_vector,
                                             qs = sim_q_matrix,
                                             crown_age = 5,
                                             num_concealed_states = 2,
                                             seed = 40,
                                             drop_extinct = TRUE)

if (requireNamespace("diversitree")) {
  traits_for_plot_reconstructed <- data.frame(
    trait = as.numeric(sim_tree_reconstructed$obs_traits),
    row.names = sim_tree_reconstructed$phy$tip.label
  )
  diversitree::trait.plot(tree = sim_tree_reconstructed$phy,
                          dat = traits_for_plot_reconstructed,
                          cols = list("trait" = c("blue", "red")),
                          type = "p")
} else {
  plot(sim_tree_reconstructed$phy)
}

reconstructed_tree_ml <- secsse_ml(phy = sim_tree_reconstructed$phy,
                                   traits = sim_tree_reconstructed$obs_traits,
                                   num_concealed_states = 2,
                                   idparslist = idparslist,
                                   idparsopt = idparsopt,
                                   initparsopt = initparsopt,
                                   idparsfix = idparsfix,
                                   parsfix = initparsfix,
                                   sampling_fraction = sampling_fraction,
                                   verbose = FALSE,
                                   is_complete_tree = FALSE)


## ----reconstructed_tree_res_comparison----------------------------------------
reconstructed_tree_ml_CR <- reconstructed_tree_ml$ML
CR_par_reconstructed <- secsse::extract_par_vals(
  idparslist,
  reconstructed_tree_ml$MLpars
)
reconstructed_tree_ml
CR_par_reconstructed
spec_rates_reconstructed <- CR_par_reconstructed[1]
ext_rates_reconstructed <- CR_par_reconstructed[2]
Q_01_reconstructed <- CR_par_reconstructed[3]
Q_10_reconstructed <- CR_par_reconstructed[4]

knitr::kable(
  data.frame(
    Reconstructed_tree = c(
      spec_rates_reconstructed,
      ext_rates_reconstructed,
      Q_01_reconstructed,
      Q_10_reconstructed
    ),  
    Complete_tree = c(
      spec_rates_complete,
      ext_rates_complete,
      Q_01_complete,
      Q_10_complete
    ),
    Generating_parameters = c(
      speciation_rate,
      extinction_rate,
      q_01,
      q_10
    ),
    row.names = c(
      "Speciation rate",
      "Extinction rate",
      "Transition rate 01",
      "Transition rate 10"
    )
  )
)

