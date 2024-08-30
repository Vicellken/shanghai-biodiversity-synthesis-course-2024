library(tidyverse)
library(mobsim)
library(mobr)
library(cowplot)

#--------set parameters for community abundance, richness and evenness-----
community_params <- expand_grid(
  # number of species
  Spool = c(10, 40, 80), 
  # total number of individuals
  J = c(50, 500),
  # evenness of species relative (log) abundance
  # smaller values result in more even distributions
  # of species relative abundance
  sdlog = c(0.8,2)) %>% 
  mutate(
    # species mean (log) abundance
    meanlog = log(J/Spool)
  ) %>% 
  # create identifier
  mutate(region = 1:n()) 


#--------simulate communities-------------------
communities <- community_params %>% 
  group_by(region) %>%
  nest(data = c(Spool, J, sdlog, meanlog)) %>%
  # create some replicates (multiple samples with the same parameters)
  uncount(20, .remove = FALSE) %>% 
  ungroup() %>% 
  mutate(sim = 1:n()) %>% 
  # generate sample with individuals randomly distributed in space
  # i.e., spatial sampling is poisson process
  mutate(poisson_comm = map(data, 
                            ~sim_poisson_community(
                              s_pool = .x$Spool,
                              n_sim = .x$J,
                              sad_type = 'lnorm',
                              sad_coef = list('meanlog' = .x$meanlog,
                                              'sdlog'= .x$sdlog)))) %>% 
  # prepare stem map for visualisation 
  mutate(poisson_stem_map = map(poisson_comm, 
                                ~tibble(x = .x$census$x,
                                        y = .x$census$y,
                                        species = .x$census$species))) %>% 
  # convert sample to SAD for further calculations
  mutate(poisson_sad = map(poisson_comm, ~community_to_sad(.x)),
         # to 2d
         poisson_sad = map(poisson_sad, ~tibble(species = names(.x),
                                                N = as.numeric(.x)))) %>% 
  # calculate IBR
  mutate(poisson_ibr = map(poisson_sad, ~rarefaction(x = .x$N,
                                                     method = 'IBR')),
         # wrangle for plotting
         poisson_ibr = map(poisson_ibr, ~tibble(individuals = 1:length(.x),
                                                expected_richness = as.numeric(.x)))) 

# examine abundance, richness and evenness for the different communities
metrics <- communities %>% 
  select(sim, region, data, poisson_sad) %>% 
  unnest(c(poisson_sad)) %>% 
  group_by(sim, region) %>% 
  summarise(J = sum(N),
            S = calc_div(N, index = 'S'),
            PIE = calc_div(N, index = 'PIE'),
            S_PIE = calc_div(N, index = 'S_PIE')) %>% 
  ungroup() %>% 
  # put the known parameters back in for visualisation
  left_join(community_params)

# plot ibr
communities %>% 
  unnest(c(data, poisson_ibr)) %>% 
  select(sim, region, individuals, expected_richness) %>% 
  mutate(spatial_dist = 'random') %>% 
  left_join(community_params, by = 'region') %>% 
  mutate(region_label = paste0('S = ',Spool,
                               ', J = ', J,
                               ', sd = ', sdlog)) %>% 
  ggplot() +
  facet_wrap(~region_label, scales = 'free') + 
  geom_line(aes(x = individuals, y = expected_richness, 
                group = sim)) +
  labs(y = 'Expected number of species',
       x = 'Number of individuals') +
  theme(legend.position = 'none')

