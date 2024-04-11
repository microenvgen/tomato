# - null hypothesis is that there is no difference between the means of the
#   "random" and "PCG" conditions (μ_random = μ_PCG)
# - alternative hypothesis is that there is a difference (μ_random ≠ μ_PCG)

## It's normal that there's more coverage if we understand that
## different taxa == more functional variability. That means the reason behind
## this might be that they are different taxa, not because they are actually
## functional groups

# APROX VALUES; var.equal~=TRUE, pero quizá más pequeña en pcg...
# mean_neg <- 0.812408888888889
# std_neg  <- 0.0312901019796617
# mean_pcg <- 0.83700142045454
# std_pcg  <- 0.02199178228350425
# 
# mean_100_neg <- 0.815604261363636
# std_100_neg  <- 0.0285377623997377
# mean_100_pcg <- 0.848864589235128
# std_100_pcg  <- 0.0224532334281521

# Input
setwd("~/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation/")
random <- read.csv("results_random_combis_12.csv", skip = 1)
pcg    <- read.csv("results_PCG_combis.csv", skip = 1)
# setwd("~/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation_all20//")
# random <- read.csv("results_random_combis_20.csv", skip = 1)
# pcg    <- read.csv("results_PCG_combis_.csv", skip = 1)

# # non-pcg vs pcgs
# setwd("~/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation_NON-COHERENT_ONLY/")
# random <- read.csv("~/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation_NON-COHERENT_ONLY/results_non-PCG_combis_12.csv", skip = 1)
# pcg    <- read.csv("~/AAA/2023-05-05_Metagenoma_minimo/results_PCG_evaluation_COHERENT_ONLY/results_PCG_combis.csv", skip = 1)


# Perform the t-test (2-sample and Welch)
t.test(random[3], pcg[3], var.equal = T)
# p-value < 2.2e-16
t.test(random[3], pcg[3], var.equal = FALSE)
# p-value = 1.246e-10


# Perform the t-test again (2-sample and Welch)
# random100 <- read.csv("results_random_combis100.csv", skip = 2)
# pcg100    <- read.csv("results_PCG_combis100.csv", skip = 2)
t.test(random100[3], pcg100[3], var.equal = T)
# p-value < 2.2e-16
t.test(random100[3], pcg100[3], var.equal = FALSE)
# p-value < 2.2e-16


# Extra
t.test(pcg[3], pcg100[3])
# p-value = 4.631e-12


# Table
results <- matrix(ncol = 4, nrow = 0,
                  dimnames = list(c(),
                                  c("Combinations",
                                    "Mean # of KEGG Reactions",
                                    "Mean coverage % of minimal metagenome",
                                    "Standard deviation of % coverage")))

datasets <- list("random" = random,
                 "pcg" = pcg,
                 "random100" = random100,
                 "pcg100" = pcg100)

for (n in 1:length(datasets)) {
  results <- rbind(results,
                   c(names(datasets)[n],
                     mean(datasets[[n]][[2]]),
                     mean(datasets[[n]][[3]]),
                     sqrt(var(datasets[[n]][[3]]))
                     )
  )
}

print(results)

