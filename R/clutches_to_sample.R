#' clutches_to_sample
#'
#' \code{clutches_to_sample} Samples clutches across a whole breeding season to
#'    determine if all of the fathers that contributed were identified.
#'
#' @param n_sims integer value, the number of simulations to run. Default value
#'    is 10000.
#' @param pop_size integer value, the population size of all breeding adults.
#'    Default value is 100.
#' @param sample_size integer value, the sample size to collect. Default value
#'    is 32.
#' @param paternal_contribution_mode a character value defining the distribution
#'    of father contributions to fertilizing a single clutch. Potential values
#'    include random', 'exponential', dominant50', 'dominant70',
#'    'dominant90', 'mixed_dominant'). Default value is 'random'.
#' @param Fprob a numeric vector, the probabilities of eggs from 1 mother being
#'    fertilized by 1-max_fathers fathers.
#' @param Mprob a numeric vector, the probabilities of fathers fertilizing eggs
#'    from 1+ mothers
#' @param clutches_mu a  numeric value, the mean number of clutches a mother
#'    lays in one nesting season. Default value is 4.95.
#' @param clutches_sd a numeric value, the standard deviation of the number of
#'    clutches a mother lays in one nesting season. Default value is 2.09.
#' @param probs_id a data frame with columns
#'    "Paternal Contribution Mode",
#'    "Fathers_Actual" (number of contributing fathers),
#'    "Sample Size" (1 - 96),
#'    "Fathers_Observed" (number of contributing fathers observed), and
#'    "Probability" (the probability of identifying Fathers_Observed given
#'    Fathers_Actual).
#' @param scenario a character vector describing the distributions of Fprob and
#'    Mprob values. Options include 'uniform_F_no_M', 'uniform_F_uniform_M',
#'    'uniform_F_base_M', 'base_F_no_M', 'base_F_uniform_M', 'base_F_base_M'.
#' @param minimum_id a numeric value denoting the minimum acceptable
#'    proportion of fathers to ID to count as 'successful'. Default value is 1.
#'
#' @return returns a data frame with the proportion of simulations where all
#'    fathers are identified given the average and standard deviation of the
#'    numbers of hatchlings in a clutch, the maximum potential number of
#'    fathers, the number of simulations, the number of hatchlings to sample
#'    from clutches, and the paternal contribution mode.
#'
#' @export
#'
#' @import dplyr
#' @import magrittr
#' @import stats
#' @import lubridate
#'
#' @examples
#' probabilities_id_fathers <- probability_id_fathers(
#'                                         hatchlings_mu = 100.58,
#'                                         hatchlings_sd = 22.61,
#'                                         max_fathers = 5,
#'                                         n_sims = 1e+04,
#'                                         sample_sizes = c(32, 96),
#'                                         paternal_contribution_modes =
#'                                           'random',
#'                                         min_clutch_size = 10)[[1]]
#'
#' clutches_to_sample(n_sims = 100,
#'                    pop_size = 100,
#'                    sample_size = 32,
#'                    paternal_contribution_mode = 'random',
#'                    Fprob = c(0.463, 0.318, 0.157, 0.034, 0.028),
#'                    Mprob = c(1),
#'                    clutches_mu = 4.95,
#'                    clutches_sd = 2.09,
#'                    probs_id = probabilities_id_fathers,
#'                    scenario = 'uniform_F_uniform_M',
#'                    minimum_id = 1)

clutches_to_sample <- function(n_sims = 10000,
                               pop_size = 100,
                               sample_size = 32,
                               paternal_contribution_mode = 'random',
                               Fprob,
                               Mprob,
                               clutches_mu = 4.95,
                               clutches_sd = 2.09,
                               probs_id,
                               scenario,
                               minimum_id)

{

  ###### Error handling ########################################################

  # classes of variables
  if (!is.numeric(n_sims)) {stop('n_sims must be a numeric value.')}
  if (!is.numeric(pop_size)) {stop('pop_size must be a numeric value.')}
  if (!is.numeric(sample_size) & !is.factor(sample_size))
  {stop('sample_size must be a numeric or factor value.')}
  if (!is.character(paternal_contribution_mode) &
      !is.factor(paternal_contribution_mode))
  {stop('paternal_contribution_mode must be a character or factor value.')}
  if (!is.numeric(Fprob)) {stop('Fprob must be a numeric value.')}
  if (!is.numeric(Mprob)) {stop('Mprob must be a numeric value.')}
  if (!is.numeric(clutches_mu)) {stop('clutches_mu must be a numeric value.')}
  if (!is.numeric(clutches_sd)) {stop('clutches_sd must be a numeric value.')}
  if (!is.data.frame(probs_id)) {stop('probs_id must be a data frame.')}
  if (!is.character(probs_id[, 1]) & !is.factor(probs_id[, 1]))
  {stop('Paternal Contribution Mode in probs_id must be a character or
        factor value.')}
  if (!is.numeric(probs_id[, 2]))
  {stop('Fathers_Actual in probs_id must be a numeric value.')}
  if (!is.numeric(probs_id[, 3]) & !is.factor(probs_id[, 3]))
  {stop('Sample_Size in probs_id must be a numeric or factor value.')}
  if (!is.numeric(probs_id[, 4]))
  {stop('Fathers_Observed in probs_id must be a numeric value.')}
  if (!is.numeric(probs_id[, 5]))
  {stop('Probability in probs_id must be a numeric value.')}
  if (!is.character(scenario) & !is.factor(scenario))
  {stop('scenario must be a character or factor value.')}
  if (!is.numeric(minimum_id))
  {stop('minimum_id must be a numeric value.')}

  # acceptable values
  if (n_sims <= 0) {stop('n_sims must be greater than 0.')}
  if (pop_size <= 0) {stop('pop_size must be greater than 0.')}
  if (as.numeric(as.character(sample_size)) <= 0)
  {stop('sample_size must be greater than 0.')}
  if (!(paternal_contribution_mode) %in% c('random', 'exponential',
                                           'dominant50', 'dominant70',
                                           'dominant90', 'mixed_dominant'))
  {stop('paternal contribution mode(s) given not recognized.')}
  if (sum(Fprob < 0) > 0) {stop('Fprob values cannot be below 0.')}
  if (sum(Fprob > 1) > 0) {stop('Fprob values cannot be above 1.')}
  if (sum(Mprob < 0) > 0) {stop('Mprob values cannot be below 0.')}
  if (sum(Mprob > 1) > 0) {stop('Mprob values cannot be above 1.')}
  if (clutches_mu <= 0) {stop('clutches_mu must be greater than 0.')}
  if (clutches_sd <= 0) {stop('clutches_sd must be greater than 0.')}
  if (sum(!(unique(probs_id[, 1])) %in% c('random', 'exponential',
                                          'dominant50', 'dominant70',
                                          'dominant90',
                                          'mixed_dominant')) > 0)
  {stop('paternal contribution mode(s) given in probs_id not recognized.')}
  if (sum(probs_id[, 2] < 1) > 0)
  {stop('probs_id Fathers_Actual cannot be below 1.')}
  if (sum(as.numeric(as.character(probs_id[, 3])) < 0) > 0)
  {stop('probs_id Sample_Size cannot be below 0.')}
  if (sum(probs_id[, 4] < 0) > 0)
  {stop('probs_id Fathers_Observed cannot be below 0.')}
  if (sum(probs_id[, 5] < 0) > 0)
  {stop('probs_id Probability cannot be below 0.')}
  if (sum(probs_id[, 5] > 1) > 0)
  {stop('probs_id Probability cannot be above 1.')}
  if (!(scenario) %in% c('uniform_F_no_M',
                         'uniform_F_uniform_M',
                         'uniform_F_base_M',
                         'base_F_no_M',
                         'base_F_uniform_M',
                         'base_F_base_M'))
  {stop('scenario given not recognized.')}
  if (minimum_id < 0) {stop('minimum_id cannot be below 0.')}
  if (minimum_id > 1) {stop('minimum_id cannot be above 1.')}

  ##############################################################################

  # dimensions
  # max number of fathers that can fertilize eggs from a single mother
  maxFathers <- length(Fprob)
  # max number of mothers whose eggs a single father can fertilize
  maxMothers <- length(Mprob)

  # operational sex ratios
  OSRs <- seq(from = 0.05, to = 0.95, by = 0.05)
  nOSR <- length(OSRs)

  # proportion of clutches sampled
  propClutches <- seq(from = 0.05, to = 1, by = 0.05)
  nPC <- length(propClutches)

  # column names for probs_id
  Paternal_Contribution_Mode <- names(probs_id)[1]
  Fathers_Actual             <- names(probs_id)[2]
  Sample_Size                <- names(probs_id)[3]
  Fathers_Observed           <- names(probs_id)[4]
  Probability                <- names(probs_id)[5]

  # pre-allocate data frame for results
  DF2 <- data.frame(OSR = rep(OSRs, each = nPC),
                    PropClutches = rep(propClutches, times = nOSR),
                    Proportion = NA)

  # for each OSR population
  for (osr in 1:nOSR) {

    # make population of Fathers and Mothers
    nF <- as.integer(pop_size*OSRs[osr])
    nM <- as.integer(pop_size - nF)

    # for each proportion of clutches sampled
    for (pc in 1:nPC) {

      # initialize vector of whether or not all fathers were identified
      ID <- rep(NA, n_sims)

      # initialize number of clutches based on number of mothers and population
      # parameters
      nClutches <- matrix(round(stats::rnorm(n = nM*n_sims,
                                             mean = clutches_mu,
                                             sd = clutches_sd)),
                          nrow = nM,
                          ncol = n_sims)

      # make sure there aren't any negative or 0 clutches, replace with 1 clutch
      nClutches[nClutches < 1] <- 1

      # initialize number of fathers for each mother and each simulation
      nFathers <- matrix(sample(1:maxFathers,
                                size = nM*n_sims,
                                prob = Fprob,
                                replace = TRUE),
                         nrow = nM,
                         ncol = n_sims)

      # proportion of clutches sampled
      prop <- propClutches[pc]

      # for each simulation
      for (i in 1:n_sims) {

        # how many mothers eggs are fertilized by each father in the population
        # vector of 1s for populations with no polygyny
        nMothers <- sample(1:maxMothers,
                           size = nF,
                           prob = Mprob,
                           replace = TRUE)

        # make breeding pool of fathers
        BPf <- rep(1:nF, times = nMothers)

        # initialize clutches list
        clutches <- NA

        # for each mother
        for (m in 1:nM) {

          # if there are no fathers left, stop the loop for the simulation
          if (dplyr::n_distinct(stats::na.omit(BPf)) == 0) { break; break }

          # how many clutches for this mother
          nC_m <- nClutches[m, i]

          # how many fathers for this mother
          nF_m <- nFathers[m, i]

          # if there are not enough unique fathers left in the breeding pool for
          # this mother
          if (dplyr::n_distinct(BPf) < nF_m) {

            # change the number of fathers to how many unique fathers are left
            nF_m <- dplyr::n_distinct(BPf)

          }

          # who are the contributing fathers themselves,
          # sample from breeding pool without duplicates
          fathers_m <- sample(unique(BPf),
                              size = nF_m,
                              replace = FALSE)

          # updated breeding pool for fathers minus the ones that already bred
          BPf <- BPf[-match(fathers_m, BPf)]

          # if there's only 1 father
          if (nF_m == 1) {

            # append identified father to clutches nN_m times, since it will
            # automatically get identified
            clutches <- append(clutches, rep(list(fathers_m), times = nC_m))

          } else {

            # probability of identification of all possible fathers for this
            # mother, pulled from probs_id data frame given
            sub <- probs_id %>%
              dplyr::filter(Paternal_Contribution_Mode == paternal_contribution_mode,
                            Fathers_Actual == nF_m,
                            Sample_Size == sample_size,
                            Probability > 0)

            # if total fathers are automatically identified for any clutch
            if (nrow(sub) == 1) {

              clutches <- append(clutches, rep(list(fathers_m), times = nC_m))

            } else {

              # how many fathers were identified in each clutch for this mother?
              nF_id <- sample(sub$Fathers_Observed,
                              size = nC_m,
                              prob = sub$Proportion_Correct,
                              replace = TRUE)

              # if there are more than 1 clutch
              for (n in 1:nC_m) {

                # add the fathers identified from the clutches for this mother
                clutches <- append(clutches, list(sample(fathers_m,
                                                         size = nF_id[n],
                                                         replace = FALSE)))

              }

            }

          }

        }

        # remove NA from clutches
        clutches <- clutches[-1]

        # number of fathers actually represented across all clutches for this
        # population in this simulation
        num_fathers <- dplyr::n_distinct(unlist(clutches))

        # how many fathers represent the minimum required proportion?
        min_fathers <- minimum_id * as.integer(num_fathers)

        # number of clutches total
        num_clutches <- length(clutches)

        # how many clutches were sampled
        num_clutches_sampled <- round(num_clutches*prop)

        # if no clutches end up getting sampled, sample 1 clutch
        if (num_clutches_sampled < 1) { num_clutches_sampled <- 1 }

        # sample clutches
        indices <- sample(1:num_clutches,
                          size = num_clutches_sampled,
                          replace = FALSE)

        # WHICH fathers were identified, add to identified fathers vector
        identified_fathers <- unlist(clutches[indices])

        # were the minimum required proportion of fathers identified?
        ID[i] <- ifelse(
          (dplyr::n_distinct(identified_fathers) >= min_fathers), 1, 0
        )

      }

      # calculate index
      index <- (osr - 1)*nPC + pc

      # proportion of simulations where all fathers were identified
      all_fathers_ID <- mean(ID, na.rm = TRUE)

      # add ID to dataframe
      DF2$Proportion[index] <- all_fathers_ID

      # print progress while running
      update1 <- paste(lubridate::now(), ' - ', scenario, '- N', pop_size,
                       ' - sample size ', sample_size, ' - ',
                       paternal_contribution_mode, ' - ', n_sims,
                       ' sims - OSR ', OSRs[osr], ' - PC ', propClutches[pc],
                       ' - done!', sep = '')

      write(update1, file = 'progress.txt', append = TRUE)

    }

  }

  # return output
  return(DF2)

}
