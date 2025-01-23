#' hatchlings_to_sample
#'
#' \code{hatchlings_to_sample}  samples hatchlings from clutches to determine
#'    the confidence of identifying all of the fathers that contributed for
#'    populations that exhibit multiple paternity.
#'
#' @param hatchlings_mu numeric value, the mean number of hatchlings produced in
#'    a clutch. Default value is 100.58.
#' @param hatchlings_sd numeric value, the standard deviation of the number of
#'    hatchlings produced in a clutch. Default value is 22.61.
#' @param max_fathers integer value, the maximum number of fathers that mothers
#'    can mate with. Default value is 5.
#' @param n_sims integer value, the number of simulations to run. Default value
#'    is 10000.
#' @param sample_sizes vector of integer values, the sample size(s) to collect.
#'    Default value is c(32, 96).
#' @param paternal_contribution_mode a character value defining the distribution
#'    of paternal contributions to a single clutch. Potential values
#'    include random', 'exponential', dominant50', 'dominant70',
#'    'dominant90', 'mixed_dominant'). Default value is 'random'.
#'
#' @return creates and saves figures to plot confidence of identifying all
#'    fathers given different numbers of actual fathers and sample sizes.
#'    Creates and saves table of confidences for specific sample sizes.
#'
#' @export
#'
#' @import dplyr
#' @import magrittr
#'
#' @examples
#' hatchlings_to_sample(hatchlings_mu = 100.58,
#'                      hatchlings_sd = 22.61,
#'                      max_fathers = 5,
#'                      n_sims = 100,
#'                      sample_sizes = c(32, 96),
#'                      paternal_contribution_mode = 'random')
#'
hatchlings_to_sample <- function(hatchlings_mu = 100.58,
                                 hatchlings_sd = 22.61,
                                 max_fathers = 5,
                                 n_sims = 10000,
                                 sample_sizes = c(32, 96),
                                 paternal_contribution_mode = 'random')

{

  ###### Error handling ########################################################

  # classes of variables
  if (!is.numeric(hatchlings_mu))
    {stop('hatchlings_mu must be a numeric value.')}
  if (!is.numeric(hatchlings_sd))
    {stop('hatchlings_sd must be a numeric value.')}
  if (!is.numeric(max_fathers)) {stop('max_fathers must be a numeric value.')}
  if (!is.numeric(n_sims)) {stop('n_sims must be a numeric value.')}
  if (sum(!is.numeric(sample_sizes)) > 0 & !is.factor(sample_sizes))
    {stop('sample_sizes must be numeric or factor values.')}
  if (!is.character(paternal_contribution_mode) &
      !is.factor(paternal_contribution_mode))
  {stop('paternal_contribution_mode must be a character or factor value.')}

  # acceptable values
  if (hatchlings_mu <= 0) {stop('hatchlings_mu must be greater than 0.')}
  if (hatchlings_sd <= 0) {stop('hatchlings_sd must be greater than 0.')}
  if (max_fathers < 2) {stop('max_fathers must be greater than 1.')}
  if (n_sims <= 0) {stop('n_sims must be greater than 0.')}
  if (sum(as.numeric(as.character(sample_sizes)) <= 0) > 0)
    {stop('sample_sizes must be greater than 0.')}
  if (sum(as.numeric(as.character(sample_sizes)) > 96) > 0)
      {stop('sample_sizes must be less than 97.')}
  if (!(paternal_contribution_mode) %in% c('random', 'exponential',
                                           'dominant50', 'dominant70',
                                           'dominant90', 'mixed_dominant'))
  {stop('paternal_contribution_mode given is not recognized')}

  ##############################################################################


  # pre-allocate data frame
  proportion_correct <- data.frame(Paternal_Contribution_Mode =
                                     paternal_contribution_mode,
                                   Fathers = rep(2:max_fathers,
                                                 each = max(sample_sizes)),
                                   Sample_Size = rep(c(1:max(sample_sizes)),
                                                     times = (max_fathers - 1)),
                                   Proportion_Correct = NA)

  # for each number of fathers that contribute to a clutch:
  for (i in 2:max_fathers) {

    # set contributions per father based on paternal contribution mode
    if (paternal_contribution_mode == 'dominant90') {
      FC <- 0.90
      contributions <- c(FC, rep((1 - FC)/(i - 1), (i - 1)))
      title <- paste('Dominant (90%) paternal contribution mode', sep = '')

    } else if (paternal_contribution_mode == 'dominant70') {
      FC <- 0.70
      contributions <- c(FC, rep((1 - FC)/(i - 1), (i - 1)))
      title <- paste('Dominant (70%) paternal contribution mode', sep = '')

    } else if (paternal_contribution_mode == 'dominant50') {
      FC <- 0.50
      contributions <- c(FC, rep((1 - FC)/(i - 1), (i - 1)))
      title <- paste('Dominant (50%) paternal contribution mode', sep = '')

    } else if (paternal_contribution_mode == 'exponential') {
      FC <- 0.5
      contributions <- 0.5^c(1:(i-1))
      contributions <- c(contributions, contributions[i-1])
      title <- 'Exponential (1/2) paternal contribution mode'

    } else if (paternal_contribution_mode == 'random') {
      contributions <- rep(1/i, i)
      title <- 'Random paternal contribution mode'

    } else if (paternal_contribution_mode == 'mixed_dominant') {
      doms <- sample(c(0.50, 0.70, 0.90), size = n_sims, replace = TRUE)
      F1 <- matrix(doms, nrow = n_sims, ncol = 1)
      F2 <- matrix(rep((1 - doms) / (i - 1), i - 1),
                   nrow = n_sims, ncol = i - 1)
      probs <- cbind(F1, F2)
      title <- 'Mixed dominant paternal contribution mode'

    }

    # proportion_correct array
    prop_correct <- rep(NA, n_sims)

    # for each sample size
    for (j in 1:max(sample_sizes)) {

      # pre-allocate correct identifications of number of fathers
      correct <- rep(NA, n_sims)

      # initialize clutch sizes
      # pull numbers of hatchlings from normal distribution
      n_hatchlings <- stats::rnorm(n = n_sims,
                                   mean = hatchlings_mu,
                                   sd = hatchlings_sd)

      # set any clutches with 0 or fewer eggs to 10 eggs
      n_hatchlings[which(n_hatchlings <= 0)] <- 10

      for (k in 1:n_sims) {

        # if mixed dominant, extract contributions
        if (paternal_contribution_mode == 'mixed_dominant') {
          contributions <- probs[k, ]}

        # make clutch with i fathers
        clutch <- sample(x = 1:i,
                         size = n_hatchlings[k],
                         replace = TRUE,
                         prob = contributions)

        # take sample of size j from clutch (or the whole clutch if < j)
        sample_size <- min(j, length(clutch))
        samples <- sample(x = clutch,
                          size = sample_size,
                          replace = FALSE)

        # correct allocation of number of fathers?
        correct[k] <- length(unique(samples)) == i

        # print progress while running
        if ((n_sims/k) %% 10 == 0) {
          paste(Sys.time(), ' - ', i, ' max fathers', ' - sample size ',
                sample_sizes[j], ' - ', paternal_contribution_mode, ' - ',
                n_sims, ' sims - ', n_sims/i*100, '% done!', sep = '')

        }

      }

      # calculate index in data frame
      index <- (i - 2)*(max(sample_sizes)) + j
      # print(index) for troubleshooting

      # stick proportion in data frame
      proportion_correct[index, 4] <- mean(correct)

      # grab column means of probs for dominant mixed paternal contribution mode
      if (paternal_contribution_mode == 'mixed_dominant') {
        contributions2 <- colMeans(probs)
      } else { contributions2 <- contributions }

      # marginal contribution of last (least dominant) father
      proportion_correct[index, 5] <- contributions2[i]

    }

  }

  #### plot results

  # color-blind friendly color palette
  colors <- viridisLite::viridis(max_fathers)

  # plot results - proportion correct
  fig1 <- ggplot2::ggplot(proportion_correct,
                          ggplot2::aes(x = proportion_correct[, 3],
                                       y = proportion_correct[, 4],
                                       col = as.factor(proportion_correct[, 2]))) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = 2) +
    ggplot2::geom_path(lwd = 1) +
    ggplot2::labs(col = 'Number \n of Fathers') +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::ylab('Proportion Correct') +
    ggplot2::xlab('Hatchlings Sampled') +
    ggplot2::geom_vline(xintercept = c(sample_sizes), linetype = 3) +
    ggplot2::ggtitle(title)

  # What's our confidence if we sample our sample sizes of eggs?
  proportion_correct_samples <- proportion_correct %>%
    dplyr::filter(names(proportion_correct)[3] %in% sample_sizes)

  # prettier tibble
  proportion_correct_samples_pretty <- proportion_correct_samples %>%
    tidyr::pivot_wider(names_from = names(proportion_correct_samples)[3],
                       values_from = names(proportion_correct_samples)[4],
                       names_prefix = 'Sample Size ')

  # what will the code produce as output
  output <- list(fig1,
                 proportion_correct,
                 proportion_correct_samples,
                 proportion_correct_samples_pretty)

  return(output)

}
