# Metropolis Hastings Sampler using Normal Approximation Proposal

This repository contains the code for a Metropolis Hastings sampler using a normal approximation proposal. The code is written in the R programming language.

The Metropolis Hastings algorithm is a Markov Chain Monte Carlo (MCMC) method used to sample from a target distribution. The algorithm works by proposing a new state based on the current state and accepting or rejecting the proposed state based on the ratio of the target distribution at the proposed state to the target distribution at the current state.

The normal approximation proposal is a proposal distribution that approximates the target distribution using a normal distribution. The normal approximation proposal is used to generate new states in the Metropolis Hastings algorithm.

The sampler was developed as part of a course project for the "Bayesian Statistics" course at Utrecht University. Its functionality is demonstrated by applying it to integrate a Bayesian Data Integration approach. The full example can be reproduced using the provided scripts and the [reproducible evironment file](renv.lock). For more information, see the [project report](docs/main.pdf).

## Installation

To install the Metropolis Hastings sampler using a normal approximation proposal, you can download the `sampler.R` file from this repository and run it in R.

## Citation

The citation for this code and report can be found in the [CITATION.bib](CITATION.bib) file within this repository.
