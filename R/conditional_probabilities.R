#' conditional_probabilities
#'
#' @param probabilities a data frame that has columns for
#'
#'   character values of the different paternal contribution modes to analyse.
#'   Default values are c('random', 'exponential', 'dominant50', 'dominant70',
#'   'dominant90', 'mixed_dominant').
#'
#'   integer numbers of actual fathers for a clutch. Default values are from
#'   1 to 5.
#'
#'   integer numbers of hatchlings to sample from a clutch. Default values are
#'   32 and 96.
#'
#'   integer numbers of fathers that were observed.
#'
#'   the numeric probability of the number of fathers observed given the number
#'   of actual fathers. Values must be between 0 and 1, inclusive.
#'
#' @returns two data frames of conditional probabilities for the number of
#' actual fathers given the number of fathers identified. The first data frame
#' has columns for paternal contribution modes, sample sizes, the number of
#' fathers observed, the number of actual fathers, and a probability column that
#' is the new conditional probability. The second data frame is wider so that
#' the conditional probability values are displayed by rows of paternal
#' contribution mode and observed numbers of fathers with columns for actual
#' numbers of fathers.
#'
#' @export
#'
#' @import dplyr
#' @import magrittr
#'
#' @examples
#' output <- probability_id_fathers(hatchlings_mu = 100.58,
#'                                  hatchlings_sd = 22.61,
#'                                  max_fathers = 5,
#'                                  n_sims = 10000,
#'                                  sample_sizes = c(32, 96),
#'                                  paternal_contribution_modes =
#'                                     c('random', 'exponential',
#'                                       'dominant50', 'dominant70',
#'                                       'dominant90', 'mixed_dominant'),
#'                                  min_clutch_size = 10)
#'
#' probabilities <- output[[1]]
#'
#' conditional_probabilities(probabilities)

conditional_probabilities <- function(probabilities) {

  ###### Error handling ########################################################

  # classes of variables
  if (!is.character(probabilities[, 1]) & !is.factor(probabilities[, 1]))
    {stop('Paternal Contribution Mode in probabilities must be a character or
        factor value.')}
  if (!is.numeric(probabilities[, 2]))
    {stop('Fathers Actual in probabilities must be a numeric value.')}
  if (!is.numeric(probabilities[, 3]) & !is.factor(probabilities[, 3]))
    {stop('Sample_Size in probabilities must be a numeric value.')}
  if (!is.numeric(probabilities[, 4]))
    {stop('Fathers Observed in probabilities must be a numeric value.')}
  if (!is.numeric(probabilities[, 5]))
    {stop('Probability in probabilities must be a numeric value.')}

  # acceptable values
  if (sum(!(unique(probabilities[, 1])) %in% c('random', 'exponential',
                                               'dominant50', 'dominant70',
                                               'dominant90',
                                               'mixed_dominant')) > 0)
    {stop('paternal contribution mode(s) in probabilities not recognized.')}
  if (sum(probabilities[, 2] < 1) > 0)
    {stop('Fathers Actual in probabilities cannot be below 1.')}
  if (sum(probabilities[, 3] < 1) > 0)
    {stop('Sample_Size in probabilities cannot be below 1.')}
  if (sum(probabilities[, 4] < 1) > 0)
    {stop('Fathers Observed in probabilities cannot be below 1.')}
  if (sum(probabilities[, 4] > max(probabilities[, 2])) > 0)
    {stop('Fathers Observed in probabilities cannot be above the maximum value
        in Fathers Actual.')}
  if (sum(probabilities[, 5] < 0) > 0)
    {stop('Probability in probabilities cannot be below zero.')}
  if (sum(probabilities[, 5] > 1) > 0)
    {stop('Probability in probabilities cannot be above 1.')}
  if (names(probabilities)[1] != 'Paternal_Contribution_Mode')
    {stop('The first column of probabilities must be named
          "Paternal_Contribution_Mode".')}
  if (names(probabilities)[2] != 'Fathers_Actual')
    {stop('The second column of probabilities must be named "Fathers_Actual".')}
  if (names(probabilities)[3] != 'Sample_Size')
    {stop('The third column of probabilities must be named "Sample_Size".')}
  if (names(probabilities)[4] != 'Fathers_Observed')
    {stop('The fourth column of probabilities must be named
          "Fathers_Observed".')}
  if (names(probabilities)[5] != 'Probability')
    {stop('The fifth column of probabilities must be named "Probability".')}

  ##############################################################################

  # extract paternal contribution modes
  PCMs <- unique(probabilities[, 1])

  # extract max number of fathers
  max_fathers <- max(probabilities[, 2])

  # extract sample sizes
  sample_sizes <- unique(probabilities[, 3])

  # initialize dataframe
  conditional_probabilities <- data.frame()

  # for each paternal contribution mode
  for (p in 1:length(PCMs)) {

    # for each sample size
    for (s in 1:length(sample_sizes)) {

      # create subset from data
      # pull out paternal contribution mode and the sample size
      subset1 <- probabilities %>%
        dplyr::filter(probabilities[, 1] == PCMs[p] &
                        probabilities[, 3] == sample_sizes[s])

      # initialize dataframe
      DF <- data.frame(Paternal_Contribution_Mode = PCMs[p],
                       Sample_Size = sample_sizes[s],
                       Fathers_Observed = unlist(mapply(rep, 1:max_fathers,
                                                        max_fathers:1)),
                       Fathers_Actual = paste(unlist(mapply(seq, 1:max_fathers,
                                                            max_fathers)),
                                              ' Actual Father(s)',
                                              sep = ''),
                       Conditional_Probability = NA)

      # probabilities of 1:max_fathers actual fathers - assume uniform
      PA <- 1/max_fathers

      # probabilities of  1:max_fathers observed across all potential actual
      # numbers of fathers
      observed_probabilities <- subset1 %>%
        dplyr::group_by(names(subset1)[4]) %>%
        dplyr::summarize(total = sum(across(5)))

      # marginal probability for each number of actual fathers
      PBs <- observed_probabilities$total / sum(subset1[, 5])

      # index restart
      index <- 0

      # for i Fathers_Observed
      for (i in 1:max_fathers) {

        # for c actual fathers
        for (c in i:max_fathers) {

          # filter so that Fathers_Actual = c and Fathers_Observed = i
          subset2 <- subset1 %>%
            dplyr::filter(across(2) == c) %>%
            dplyr::filter(across(4) == i)

          # probability of i Fathers_Observed given c actual fathers
          PBA <- subset2[, 5]

          # calculate PAB (probability of actual fathers given fathers observed)
          PAB <- PA * PBA / PBs[i]

          # index
          index <- index + 1

          # # troubleshooting
          # print(index)

          # add the PAB to the data frame
          DF[index, 5] <- PAB

        }

      }

      # get rid of any NaN values if they exist in the Probability column
      if (sum(is.nan(DF[, 5])) > 1) {

        # replace values that are not numbers with NA
        DF[which(is.nan(DF[, 5])), ][, 5] <- 0

      }

      # round conditional probability to 3 digits
      DF[, 5] <- round(DF[, 5], 3)

      # add
      conditional_probabilities <- rbind(conditional_probabilities, DF)


    }

  }

  # make it a less obnoxiously long table
  conditional_probabilities_pretty <- conditional_probabilities %>%
    tidyr::pivot_wider(names_from = names(conditional_probabilities)[4],
                       values_from = names(conditional_probabilities)[5]) %>%
    dplyr::arrange(across(1), across(2), across(3))

  # put output objects together in a list
  output <- list(conditional_probabilities, conditional_probabilities_pretty)

  return(output)

}
