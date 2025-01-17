#' probability_id_fathers
#'
#' \code{probability_id_fathers} saves a data frame recording the probabilities
#' of detecting 1 - max_fathers fathers given a specific number of potential
#' fathers, sample size, and paternal contribution mode.
#'
#' @param hatchlings_mu numeric value, the mean number of hatchlings produced in
#'    a clutch. Default value is 100.58.
#' @param hatchlings_sd numeric value, the standard deviation of the number of
#'    hatchlings produced in a clutch. Default value is 22.61.
#' @param max_fathers integer value, the maximum number of fathers that can
#'    fertilize the eggs from a single mother. Default value is 5.
#' @param n_sims integer value, the number of simulations to run. Default value
#'    is 10000.
#' @param sample_sizes vector of integer values, the sample size(s) to collect.
#'    Default value is c(32, 96).
#' @param paternal_contribution_modes vector of character values, the different
#'    paternal contribution modes to analyse. Default value is c('random',
#'    'exponential', 'dominant50', 'dominant70', 'dominant90', 'mixed_dominant').
#' @param min_clutch_size integer value, the number of hatchlings to be assigned
#'    to a clutch when 0 was returned from the normal distribution based on
#'    hatchlings_mu and hatchlings_sd. Default value is 10.
#'
#' @return saves object with probabilities of identifying the number of actual
#'    fathers given paternal contribution mode, sample size, and number of
#'    fathers
#'
#' @export
#'
#' @import dplyr
#' @import magrittr
#'
#' @examples
#' probability_id_fathers(hatchlings_mu = 100.58,
#'                        hatchlings_sd = 22.61,
#'                        max_fathers = 5,
#'                        n_sims = 1000,
#'                        sample_sizes = c(32, 96),
#'                        paternal_contribution_modes = c('random',
#'                                                        'exponential',
#'                                                        'dominant50',
#'                                                        'dominant70',
#'                                                        'dominant90',
#'                                                        'mixed_dominant'),
#'                        min_clutch_size = 10)

probability_id_fathers <- function(hatchlings_mu = 100.58,
                                   hatchlings_sd = 22.61,
                                   max_fathers = 5,
                                   n_sims = 10000,
                                   sample_sizes = c(32, 96),
                                   paternal_contribution_modes = c('random',
                                                                   'exponential',
                                                                   'dominant50',
                                                                   'dominant70',
                                                                   'dominant90',
                                                                   'mixed_dominant'),
                                   min_clutch_size = 10)

{

  ###### Error handling ########################################################

  # classes of variables
  if (!is.numeric(hatchlings_mu)) {stop('hatchlings_mu must be a numeric value.')}
  if (!is.numeric(hatchlings_sd)) {stop('hatchlings_sd must be a numeric value.')}
  if (!is.numeric(max_fathers)) {stop('max_fathers must be a numeric value.')}
  if (!is.numeric(n_sims)) {stop('n_sims must be a numeric value.')}
  if (sum(!is.numeric(sample_sizes)) > 0 & !is.factor(sample_size))
    {stop('sample_sizes must be numeric or factor values.')}
  if (!is.character(paternal_contribution_mode) &
      !is.factor(paternal_contribution_mode))
    {stop('paternal_contribution_mode must be a character or factor value.')}
  if (!is.numeric(min_clutch_size))
    {stop('min_clutch_size must be a numeric value.')}

  # acceptable values
  if (hatchlings_mu <= 0) {stop('hatchlings_mu must be greater than 0.')}
  if (hatchlings_sd <= 0) {stop('hatchlings_sd must be greater than 0.')}
  if (max_fathers < 2) {stop('max_fathers must be greater than 1.')}
  if (n_sims <= 0) {stop('n_sims must be greater than 0.')}
  if (sum(as.numeric(as.character(sample_sizes)) <= 0) > 0)
    {stop('sample_sizes must be greater than 0.')}
  if (sum(as.numeric(as.character(sample_sizes)) > 96) > 0)
    {stop('sample_sizes must be less than 97.')}
  if (sum(!(paternal_contribution_modes) %in% c('random', 'exponential',
                                                'dominant50', 'dominant70',
                                                'dominant90',
                                                'mixed_dominant')) > 0)
  {stop('paternal_contribution_mode(s) given are not recognized')}
  if (min_clutch_size <= 0) {stop('min_clutch_size must be greater than 0.')}

  # relational values
  if(min_clutch_size > (hatchlings_mu - 3*hatchlings_sd))
  {stop('min_clutch_size is too large given the distribution of hatchlings per
        clutch.')}

  ##############################################################################

  # dimensions
  nP <- length(paternal_contribution_modes)
  nF <- max_fathers
  nS <- length(sample_sizes)

  # pre-allocate data frame
  probabilities <- data.frame(Paternal_Contribution_Mode = NA,
                              Fathers_Actual = NA,
                              Sample_Size = NA,
                              Fathers_Observed = NA,
                              Probability = NA)

  # for each paternal contribution mode
  for (p in 1:nP) {

    # paternal contribution mode
    PCM <- paternal_contribution_modes[p]

    # for each number of Fathers_Actual:
    for (f in 1:max_fathers) {

      # if there's only one father, he's the full contributor
      if (f == 1) {

        contributions <- 1

      } else {

        # set contributions per father
        if (PCM == 'dominant90') {
          FC <- 0.90
          contributions <- c(FC, rep((1 - FC)/(f - 1), (f - 1)))

        } else if (PCM == 'dominant70') {
          FC <- 0.70
          contributions <- c(FC, rep((1 - FC)/(f - 1), (f - 1)))

        } else if (PCM == 'dominant50') {
          FC <- 0.50
          contributions <- c(FC, rep((1 - FC)/(f - 1), (f - 1)))

        } else if (PCM == 'exponential') {
          FC <- 0.5
          contributions <- 0.5^c(1:(f-1))
          contributions <- c(contributions, contributions[f-1])

        } else if (PCM == 'random') {
          contributions <- rep(1/f, f)

        } else if (PCM == 'mixed_dominant') {
          doms <- sample(c(0.50, 0.70, 0.90), size = n_sims, replace = TRUE)
          M1 <- matrix(doms, nrow = n_sims, ncol = 1)
          M2 <- matrix(rep((1 - doms) / (f - 1), f - 1), nrow = n_sims,
                       ncol = f - 1)
          probs <- cbind(M1, M2)

        }

      }

      # for each sample size
      for (s in 1:nS) {

        # pre-allocate vector of fathers observed (identified)
        fathers_observed <- rep(NA, n_sims)

        # initialize clutch sizes
        # pull numbers of hatchlings from normal distribution
        n_hatchlings <- stats::rnorm(n = n_sims,
                                     mean = hatchlings_mu,
                                     sd = hatchlings_sd)

        # set any clutches with 0 or fewer eggs to the minimum clutch size
        n_hatchlings[which(n_hatchlings <= 0)] <- min_clutch_size

        # for each simulation
        for (j in 1:n_sims) {

          # if mixed dominant, extract contributions (assuming f > 1)
          if (PCM == 'mixed_dominant' & f > 1) {

            contributions <- probs[j, ]

          }

          # make clutch with i fathers
          clutch <- sample(x = 1:f,
                           size = n_hatchlings[j],
                           replace = TRUE,
                           prob = contributions)

          # take sample of size j from clutch (or the whole clutch if < j)
          sample_size <- min(sample_sizes[s], length(clutch))
          samples <- sample(x = clutch,
                            size = sample_size,
                            replace = FALSE)

          # how many fathers were observed
          fathers_observed[j] <- length(unique(samples))

        }

        # # for troubleshooting, if needed
        # calculate index in data frame
        # index <- (i - 2)*(max(sample_sizes)) + j
        # print(index)

        # probabilities for each number of fathers
        props <- rep(NA, f)

        for (k in 1:f) {

          props[k] <- round(length(fathers_observed[fathers_observed == k]) /
                              length(fathers_observed), digits = 4)

        }

        # pre-allocate data frame
        DF <- data.frame(Paternal_Contribution_Mode = rep(PCM, times = f),
                         Fathers_Actual = rep(f, times = f),
                         Sample_Size = rep(sample_sizes[s], times = f),
                         Fathers_Observed = 1:f,
                         Probability = props)

        # attach to DF
        probabilities <- rbind(probabilities, DF)

      }

    }

  }

  # remove row of NAs at the top
  probabilities <- probabilities[-1, ]

  # save object
  save(probabilities,
       file = paste('output/probabilities_', n_sims, '.Rdata', sep = ''))

  # make prettier object
  probabilities_pretty <- probabilities %>%
    tidyr::pivot_wider(names_from = 'Fathers_Observed',
                       values_from = 'Probability') %>%
    dplyr::rename('1 Observed Father' = `1`,
                  '2 Observed Fathers' = `2`,
                  '3 Observed Fathers' = `3`,
                  '4 Observed Fathers' = `4`,
                  '5 Observed Fathers' = `5`) %>%
    dplyr::arrange(Paternal_Contribution_Mode, Sample_Size, Fathers_Actual)

  # save prettier object
  save(probabilities_pretty,
       file = paste('output/probabilities_pretty_', n_sims, '.Rdata', sep = ''))

  # make output
  output <- list(probabilities, probabilities_pretty)

  # return output
  return(output)

}
