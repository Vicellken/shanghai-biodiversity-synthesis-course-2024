library(tidyverse)
library(mobsim)
library(mobr)
library(cowplot)

# set a RNG seed 
set.seed(42)

# Set all communities to have 5000 individuals, 
# a lognormal SAD, and individuals randomly distributed in space
Jpool <- 2000

# We want to treatment to half the number of species
Spool_control <- 200
# treatment removes 50% of the species
Spool_treatment <- 0.5 * Spool_control

# simulate a meta-analysis with twenty studies
n_study <- 20

meta_sim <- tibble(
  Jpool = Jpool,
  S_control = Spool_control,
  S_treatment = Spool_treatment) %>% 
  # create the n_studies
  uncount(n_study, .remove = FALSE) %>% 
  # here, we are interested in examining variation in sample grain.
  # Draw some random quadrat sizes from a uniform distribution
  mutate(sample_grain = runif(n_study,
                              min = 0.01,
                              max = 0.1)) %>% 
  # create index to identify each study
  mutate(study = as.character(1:n_study)) %>% 
  # prepare to generate control and treatment communities for each study
  group_by(study) %>% 
  nest(data = c(Jpool, S_control, S_treatment,
                sample_grain)) %>% 
  # simulate poisson distribution of individuals for the controls and 
  # treatments within each study
  mutate(control_comm = map(data, 
                            ~sim_poisson_community(
                              s_pool = .x$S_control,
                              n_sim = .x$Jpool,
                              sad_type = 'lnorm',
                              sad_coef = list('meanlog' = log(.x$S_control/.x$Jpool),
                                              'sdlog' = 1))),
         treatment_comm = map(data, 
                              ~sim_poisson_community(
                                s_pool = .x$S_treatment,
                                n_sim = .x$Jpool,
                                sad_type = 'lnorm',
                                sad_coef = list('meanlog' = log(.x$S_treatment/.x$Jpool),
                                                'sdlog' = 1)
                              )
         )
  ) %>% 
  # to check our treatment worked calculate abundance and richness for the 
  # whole community (which is the scale at which we applied the treatment)
  mutate(control_comm_J = map(control_comm, ~nrow(.x$census)),
         control_comm_S = map(control_comm, ~n_distinct(.x$census$species)),
         treatment_comm_J = map(treatment_comm, ~nrow(.x$census)),
         treatment_comm_S = map(treatment_comm, ~n_distinct(.x$census$species))) %>% 
  # now, we want to get some samples from the controls and treatments
  # we'll just keep the site x species matrix
  mutate(control_samps = map2(control_comm, data, 
                              ~ sample_quadrats(comm = .x,
                                                # we'll take 5 samples
                                                n_quadrats = 5,
                                                # with the sample_grain for this study
                                                quadrat_area = .y$sample_grain,
                                                method = 'grid',
                                                plot = FALSE)$spec_dat),
         treatment_samps = map2(treatment_comm, data, 
                                ~ sample_quadrats(comm = .x,
                                                  # we'll take 5 samples
                                                  n_quadrats = 5,
                                                  # with the sample_grain for this study
                                                  quadrat_area = .y$sample_grain,
                                                  method = 'grid',
                                                  plot = FALSE)$spec_dat)) %>% 
  # let's calculate metrics for abundance, richness, and the ENS conversion
  # of PIE (for inferences about evenness)
  # note, because we're simulating what we'd typically find when doing a meta-analysis,
  # we're only going to do the calculations at a single scale
  mutate(control_J = map(control_samps, ~ rowSums(.x)),
         control_S = map(control_samps, ~ vegan::specnumber(.x, MARGIN = 1)),
         control_SPIE = map(control_samps, ~ mobr::calc_comm_div(.x, index = 'S_PIE', 
                                                                 extrapolate = FALSE,
                                                                 scales = 'alpha')$value),
         treatment_J = map(treatment_samps, ~ rowSums(.x)),
         treatment_S = map(treatment_samps, ~ vegan::specnumber(.x, MARGIN = 1)),
         treatment_SPIE = map(treatment_samps, ~ mobr::calc_comm_div(.x, index = 'S_PIE', 
                                                                     extrapolate = FALSE,
                                                                     scales = 'alpha')$value)) %>% 
  ungroup()

# calculate log-ratio effect sizes
effect_sizes <- meta_sim %>% 
  unnest(c(data, control_J, control_S, control_SPIE,
           treatment_J, treatment_S, treatment_SPIE)) %>% 
  mutate(J_LRR = log(treatment_J/control_J),
         S_LRR = log(treatment_S/control_S),
         SPIE_LRR = log(treatment_SPIE/control_SPIE))

# plot come effect sizes, let's start with total numbers of individuals
effect_sizes %>% 
  # calculate the mean and sd
  mutate(J_LRR_mean = mean(J_LRR),
         J_LRR_sd = sd(J_LRR)) %>% 
  ggplot() +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = J_LRR_mean - J_LRR_sd,
                ymax = J_LRR_mean + J_LRR_sd),
            fill = 'light grey') +
  geom_hline(aes(yintercept = J_LRR_mean)) +
  # and plot the known effect size on total abundance
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(aes(x = study, y = J_LRR)) +
  labs(x = 'Study',
       y = 'Effect size (LRR)',
       subtitle = 'Treatment effect on total abundance') +
  theme_bw()


# now species richness, note our known effect size
S_LRR_known <- log(0.5)

effect_sizes %>% 
  # calculate the mean and sd
  mutate(S_LRR_mean = mean(S_LRR),
         S_LRR_sd = sd(S_LRR)) %>% 
  ggplot() +
  geom_rect(aes(xmin = -Inf, xmax = Inf,
                ymin = S_LRR_mean - S_LRR_sd,
                ymax = S_LRR_mean + S_LRR_sd),
            fill = 'light grey') +
  geom_hline(aes(yintercept = S_LRR_mean)) +
  # let's plot the known effect size as a dashed line
  geom_hline(yintercept = S_LRR_known, linetype = 2) + 
  geom_point(aes(x = study, y = S_LRR)) +
  labs(x = 'Study',
       y = 'Effect size (LRR)',
       subtitle = 'Treatment effect on species richness') +
  theme_bw()


effect_sizes %>% 
  # calculate the mean and sd of the sample ESs
  mutate(S_LRR_mean = mean(S_LRR),
         S_LRR_sd = sd(S_LRR)) %>% 
  ggplot() +
  geom_hline(aes(yintercept = S_LRR_mean)) +
  # plot the known effect size as a dashed line
  geom_hline(yintercept = S_LRR_known, linetype = 2) + 
  # plot the sample effect sizes
  geom_point(aes(x = sample_grain, y = S_LRR)) +
  # to see if we've got any evidence for scale-dependence
  # we'll visualise a linear model with the effect sizes as a function of
  # grain size
  stat_smooth(aes(x = sample_grain, y = S_LRR),
              method = 'lm') +
  labs(x = 'Grain size',
       y = 'Effect size (LRR)',
       subtitle = 'Estimated treatment effect on species richness') +
  theme_bw()
