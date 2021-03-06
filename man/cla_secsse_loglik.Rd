% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cla_secsse_loglik.R
\name{cla_secsse_loglik}
\alias{cla_secsse_loglik}
\title{Likelihood for SecSSE model}
\usage{
cla_secsse_loglik(
  parameter,
  phy,
  traits,
  num_concealed_states,
  use_fortran = TRUE,
  methode = "ode45",
  cond = "proper_cond",
  root_state_weight = "proper_weights",
  sampling_fraction,
  run_parallel = FALSE,
  setting_calculation = NULL,
  setting_parallel = NULL,
  see_ancestral_states = FALSE,
  loglik_penalty = 0,
  is_complete_tree = FALSE,
  func = ifelse(is_complete_tree, "cla_secsse_runmod_ct_d", ifelse(use_fortran ==
    FALSE, cla_secsse_loglik_rhs, "cla_secsse_runmod")),
  reltol = 1e-12,
  abstol = 1e-16
)
}
\arguments{
\item{parameter}{list where the first is a table where lambdas across different modes of speciation are shown, the second mus and the third transition rates.}

\item{phy}{phylogenetic tree of class phylo, ultrametric, fully-resolved, rooted and with branch lengths.}

\item{traits}{vector with trait states, order of states must be the same as tree tips, for help, see vignette.}

\item{num_concealed_states}{number of concealed states, generally equivalent to number of examined states.}

\item{use_fortran}{Should the Fortran code for numerical integration be called? Default is TRUE.}

\item{methode}{Solver for the set of equations, default is "ode45".}

\item{cond}{condition on the existence of a node root: "maddison_cond","proper_cond"(default). For details, see vignette.}

\item{root_state_weight}{the method to weigh the states:"maddison_weights","proper_weights"(default) or "equal_weights". It can also be specified the root state:the vector c(1,0,0) indicates state 1 was the root state.}

\item{sampling_fraction}{vector that states the sampling proportion per trait state. It must have as many elements as trait states.}

\item{run_parallel}{should the routine to run in parallel be called?}

\item{setting_calculation}{argument used internally to speed up calculation. It should be leave blank (default : setting_calculation = NULL)}

\item{setting_parallel}{argument used internally to set a parallel calculation. It should be left blank (default : setting_parallel = NULL)}

\item{see_ancestral_states}{should the ancestral states be shown? Deafault FALSE}

\item{loglik_penalty}{the size of the penalty for all parameters; default is 0 (no penalty)}

\item{is_complete_tree}{whether or not a tree with all its extinct species is provided}

\item{func}{function to be used in solving the ODE system. Currently only for testing purposes.}

\item{reltol}{relative tolerance in integration}

\item{abstol}{absolute tolerance in integration}
}
\value{
The loglikelihood of the data given the parameters
}
\description{
Logikelihood calculation for the cla_SecSSE model given a set of parameters and data
}
\note{
To run in parallel it is needed to load the following libraries when windows: apTreeshape, doparallel and foreach. When unix, it requires: apTreeshape, doparallel, foreach and doMC
}
\examples{
rm(list=ls(all=TRUE))
library(secsse)
library(DDD)
library(deSolve)
#library(diversitree)
library(apTreeshape)
library(foreach)
set.seed(13)
phylotree <- ape::rcoal(12, tip.label = 1:12)
traits <- sample(c(0,1,2),ape::Ntip(phylotree),replace=TRUE)
num_concealed_states <- 3
sampling_fraction <- c(1,1,1)
methode <- "ode45"
phy <- phylotree
# the idparlist for a ETD model (dual state inheritance model of evolution)
# would be set like this:
idparlist <- cla_id_paramPos(traits,num_concealed_states)
lambd_and_modeSpe <- idparlist$lambdas
lambd_and_modeSpe[1,] <- c(1,1,1,2,2,2,3,3,3)
idparlist[[1]] <- lambd_and_modeSpe
idparlist[[2]][] <- 0
masterBlock <- matrix(4,ncol=3,nrow=3,byrow=TRUE) 
diag(masterBlock) <- NA
idparlist [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
# Now, internally, clasecsse sorts the lambda matrices, so they look like:
prepare_full_lambdas(traits,num_concealed_states,idparlist[[1]]) 
# which is a list with 9 matrices, corresponding to the 9 states (0A,1A,2A,0B,etc)
# if we want to calculate a single likelihood:
parameter <- idparlist
lambd_and_modeSpe <- parameter$lambdas
lambd_and_modeSpe[1,] <- c(0.2,0.2,0.2,0.4,0.4,0.4,0.01,0.01,0.01)
parameter[[1]] <- prepare_full_lambdas(traits,num_concealed_states,lambd_and_modeSpe) 
parameter[[2]] <- rep(0,9)
masterBlock <- matrix(0.07,ncol=3,nrow=3,byrow=TRUE) 
diag(masterBlock) <- NA
parameter [[3]] <- q_doubletrans(traits,masterBlock,diff.conceal = FALSE)
cla_secsse_loglik(parameter, phy, traits, num_concealed_states,
                 use_fortran = FALSE, methode = "ode45", cond = "maddison_cond",
                 root_state_weight = "maddison_weights", sampling_fraction,
                 run_parallel = FALSE, setting_calculation = NULL,
                 setting_parallel = NULL, see_ancestral_states = FALSE,
                 loglik_penalty = 0)
# LL = -37.8741
}
