## Copyright (C) 2021-2024 Robersy Sanchez <https://genomaths.com/>
## Author: Robersy Sanchez This file is part of the R package
## 'GenomAutomorphism'.  'GenomAutomorphism' is a free
## software: you can redistribute it and/or modify it under the
## terms of the GNU General Public License as published by the Free
## Software Foundation, either version 3 of the License, or (at
## your option) any later version.  This program is distributed in
## the hope that it will be useful, but WITHOUT ANY WARRANTY;
## without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for
## more details.  You should have received a copy of the GNU
## General Public License along with this program; if not, see
## <http://www.gnu.org/licenses/>.

#' @aliases automorphism_prob
#' @rdname automorphism_prob
#' @title Autmorphism Probability
#' @description
#' This function applies a Dirichlet-Multinomial Modelling (in Bayesian 
#' framework) to compute the posterior probability of each 
#' \emph{type of mutational event}. DNA bases are classified based on the 
#' physicochemical criteria used to ordering the set of codons: number of
#' hydrogen bonds (strong-weak, S-W), chemical type (purine-pyrimidine, Y-R),
#' and chemical groups (amino versus keto, M-K) (see reference 4). Preserved
#' codon positions are labeled with letter “H”.
#' 
#' As a result, mutational events are grouped by type of mutation covering all
#' the possible combinations of symbols: "Y", "R", "M", "W", "K", "S", and "H", 
#' for example:  "YSH", "RSH", "MSH", "WSH", "KSH", "SSH", "HSH", "YHH", "RHH", 
#' and so on. Insertion/deletion mutations are not considered.
#' 
#' ## Maximum Likelihood Estimation (MLE) for Dirichlet Parameters
#' 
#' Given the observed data \eqn{n = (n_1, ..., n_k)}, where \eqn{n_i} is the 
#' frequency for the mutation type \eqn{i}, we want to estimate the parameters
#' of the Dirichlet distribution, \eqn{\alpha = (\alpha_1, ..., \alpha_k)},
#' that maximize the marginal likelihood.
#' 
#' ## Marginal Likelihood
#' The marginal likelihood of the data under the Dirichlet-multinomial model is
#' given by
#' 
#' \deqn{P\left(n\middle|\alpha\right) = \frac{N!}{\prod_{i=1}^{k}n_i}
#' \frac{\Gamma\left(\sum_{i=1}^{k}\alpha_i\right)}{\Gamma
#' \left(N+\sum_{i=1}^{k}\alpha_i\right)}\prod_{i=1}^{k}\frac{\Gamma
#' \left(n_i+\alpha_i\right)}{\Gamma\left(\alpha_i\right)}}
#' 
#' where \eqn{N = \sum_{i=1}^{k}n_i}.
#' 
#' ## Optimization
#' To perform MLE, we maximize the log of this likelihood: 
#' \eqn{log(P\left(n\middle|\alpha\right))}. That is, We aim to maximize this 
#' log-likelihood with respect to \eqn{\alpha}. This is done numerically 
#' because there's no closed-form solution. Here, we use:
#' 
#' \deqn{Arg\max_{\alpha} {\{log\left(P\left(n\middle|\alpha\right)\right\}}}
#' 
#' with initial guess set as \eqn{\alpha_i = n_i + 1}.
#' 
#' ## Posterior Distribution
#' 
#' Let be \eqn{\theta} the vector of probabilities for each mutation type. It's 
#' the parameter we're estimating, which represents the probabilities of
#' observing each mutation type in the multinomial distribution. The conjugate
#' distribution in this context refers to the Dirichlet distribution, which is 
#' the prior distribution for \eqn{\theta}. When we have observed data, the
#' posterior distribution of \eqn{\theta} given this data is also a Dirichlet 
#' distribution due to the conjugate property.
#' 
#' The prior distribution, before observing any data, our belief about the 
#' distribution of \eqn{\theta} is given by:
#' 
#' \deqn{\theta\sim\ Dirichlet\left(\alpha_1,\alpha_k,\ldots,\alpha_k\right)}
#' 
#' where the \eqn{\alpha_i} are known as initial concentration parameters, 
#' which represents our initial belief or assumption about the distribution
#' of the mutation types. These parameters control how the probability mass is
#' distributed among the mutation types. We might set these initial parameters 
#' based on some prior knowledge or simply to provide a non-informative prior. 
#' 
#' Once we have the MLE for \eqn{\alpha}, the posterior distribution of 
#' \eqn{\theta} given the data updates these parameters to incorporate the
#' observed data:
#' 
#' \deqn{\theta\mid\ n\sim\ Dirichlet\left(\alpha_1^\ast+n1,\alpha_k^\ast+n_2,
#' \ldots,\alpha_k^\ast+n_k\right)} 
#' 
#' where \eqn{\alpha^\ast_i} are the MLE estimates, i.e., 
#' \eqn{\alpha_i^\ast+n_i} are our updated belief about the frequencies of
#' mutation types in the population sample. 
#' 
#' 
#' The expected values of \eqn{\theta_i}, given the data under the posterior 
#' Dirichlet distribution is:
#' 
#' \deqn{E\left[\theta_i\right]=\frac{\alpha_i^*+n_i}
#' {\sum_{j=1}^{k}\left(\alpha_j^*+n_j\right)}}
#' 
#' These expected values directly corresponds to the "posterior probability" of
#' observing mutation type \eqn{i}.
#' 
#' This approach provides a rigorous estimation of the Dirichlet parameters 
#' under the Dirichlet-multinomial model using MLE.
#' 
#' @param x An AutomorphismByCoef or anAutomorphismByCoefList-class object 
#' returned by function [automorphism_bycoef].
#' 
#' @param method,maxit,abstol Parameter values to pass into 
#' \code{\link[stats]{optim }}.
#' @returns A data frame with the posterior probabilities.
#' 
#' @examples
#' ## Load the data set
#' data("autby_coef", package = "GenomAutomorphism")
#' post_prob <- automorphism_prob(autby_coef[1:10])
#' head(post_prob,10)
#' @export
#' @seealso [automorphism_bycoef]

setGeneric(
    "automorphism_prob",
    function(x, ...) {
        standardGeneric("automorphism_prob")
    }
)


#' @aliases automorphism_prob
#' @rdname automorphism_prob
#' @param x An AutomorphismByCoefList-class object returned by function
#' [automorphism_bycoef].
#' @importFrom stats optim
#' @export
setMethod(
    "automorphism_prob", signature(x = "AutomorphismByCoef"),
    function(
        x, 
        method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN", "Brent"),
        maxit = 500,
        abstol = 10^-8) {
        
        x <- table(x$mut_type[ x$autm != 1 & x$autm != -1 ])
        x <- x[-1] ## remove indel mutations
        x <- sort(x, decreasing = TRUE)
        x <- data.frame(x)
        colnames(x) <- c("Codon", "Freq")
        
        ## Total number of observed mutations
        N <- sum(x$Freq)
        
        ## All possible codon combinations considering the given symbols, 
        ## excluding HHH
        symbols <- c("Y", "R", "M", "W", "K", "S", "H")
        all_combinations <- expand.grid(
            C1 = symbols,
            C2 = symbols,
            C3 = symbols
        )
        all_combinations <- all_combinations[!(all_combinations$C1 == "H" & 
                                        all_combinations$C2 == "H" & 
                                            all_combinations$C3 == "H"),]
        all_combinations$Codon <- paste(all_combinations$C1, 
                                        all_combinations$C2, 
                                        all_combinations$C3, sep = "")
        
        # Merge with data to include zeros for unobserved mutations
        full_data <- merge(all_combinations, x, by = "Codon", all.x = TRUE)
        
        ## fill NA with 0 for unobserved mutation types
        full_data$Freq[is.na(full_data$Freq)] <- 0  
        
        # Log-likelihood function for optimization
        log_likelihood <- function(alpha) {
            k <- length(alpha)
            sum_alpha <- sum(alpha)
            return(
                lgamma(sum_alpha) - lgamma(N + sum_alpha) + 
                    sum(lgamma(alpha + full_data$Freq) - lgamma(alpha))
            )
        }
        
        # Initial guess for alpha (using n_i + 1)
        initial_alpha <- full_data$Freq + 1
        
        # Optimization
        if (method == "L-BFGS-B") {
            optim_result <- optim(initial_alpha, 
                                log_likelihood, 
                                method = method, 
                                lower = rep(0.1, length(initial_alpha)),
                                control = list(
                                    fnscale = -1,
                                    maxit = maxit,
                                    pgtol = abstol))
        }
        else {
            optim_result <- optim(initial_alpha, 
                                log_likelihood, 
                                method = method, 
                                control = list(
                                    fnscale = -1,
                                    maxit = maxit,
                                    abstol = abstol))            
        }
        
        # Optimized alpha parameters
        alpha_mle <- optim_result$par
        
        # Compute posterior probabilities with the MLE alpha
        posterior_probs <- alpha_mle / sum(alpha_mle)
        
        # Combine with Codon data for clarity
        results <- data.frame(
            Codon = full_data$Codon,
            Posterior_Probability = posterior_probs
        )
        
        # Sort results by probability in descending order
        results <- results[order(-results$Posterior_Probability), ]
        rownames(results) <- NULL
        
        return(results)
        
    }
)



#' @aliases automorphism_prob
#' @rdname automorphism_prob
#' @param x An AutomorphismByCoefList-class object returned by function
#' [automorphism_bycoef].
#' @importFrom stats optim
#' @export
setMethod(
    "automorphism_prob", signature(x = "AutomorphismByCoefList"),
    function(
        x, 
        method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN", "Brent"),
        maxit = 500,
        abstol = 10^-8) {
        
        method <- match.arg(method)
        
        data <- unlist(x)
        
        return(automorphism_prob(
            x = data, 
            method = method,
            maxit = maxit,
            abstol = abstol))
    }
)

