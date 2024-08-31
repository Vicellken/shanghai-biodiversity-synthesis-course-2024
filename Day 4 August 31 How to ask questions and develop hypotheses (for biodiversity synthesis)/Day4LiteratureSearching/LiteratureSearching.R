remotes::install_github("elizagrames/litsearchr")
remotes::install_github("nealhaddaway/citationchaser")
install.packages("ggraph")
install.packages("igraph")

library(litsearchr)
library(citationchaser)
library(ggraph)
library(igraph)
library(tidyverse)

# Litsearchr ####
import <-
  import_results(file = "LitsearchrUrbanFish.ciw", verbose = TRUE)

#extract the keywords from literatures
keywords <- extract_terms(keywords = import[, "keywords"],
                method = "tagged",
                min_n = 1)
#extract the title terms from literatures
title_terms <-
  extract_terms(
    text = import[, "title"],
    method = "fakerake",
    min_freq = 2,
    min_n = 1
  )

terms <- unique(c(keywords, title_terms))

#get abstract terms with title and keywrods together
docs <- paste(terms, import[, "abstract"])

dfm <- create_dfm(elements = docs, features = terms)

g <- create_network(dfm, min_studies = 2)
#see the network betweens terms
ggraph(g, layout = "stress") +
  coord_fixed() +
  expand_limits(x = c(-3, 3)) +
  geom_edge_link(aes(alpha = weight)) +
  geom_node_point(shape = "circle filled", fill = "white") +
  geom_node_text(aes(label = name), hjust = "outward", check_overlap = TRUE) +
  guides(edge_alpha = "none")

strengths <- strength(g)

term_strengths <- data.frame(term = names(strengths),
           strength = strengths,
           row.names = NULL) %>%
  mutate(rank = rank(strength, ties.method = "min")) %>%
  arrange(strength) 

term_strengths

cutoff_fig <-
  ggplot(term_strengths, aes(x = rank, y = strength, label = term)) +
  geom_line() +
  geom_point() +
  geom_text(
    data = filter(term_strengths, rank > 5),
    hjust = "right",
    nudge_y = 20,
    check_overlap = TRUE
  )

cutoff_fig

cutoff_cum <- find_cutoff(g, method = "cumulative", percent = 0.8)

cutoff_cum

cutoff_fig +
  geom_hline(yintercept = cutoff_cum, linetype = "dashed")

get_keywords(reduce_graph(g, cutoff_cum))

cutoff_change <- find_cutoff(g, method = "changepoint", knot_num = 3)

cutoff_change

cutoff_fig +
  geom_hline(yintercept = cutoff_change, linetype = "dashed")

g_redux <- reduce_graph(g, cutoff_change[2])

selected_terms <- get_keywords(g_redux)

selected_terms

write.csv(selected_terms, file = "selected_terms.csv")


# Citationchaser ####

doi_list <- import[,"doi"] %>% 
  as.data.frame() 

doi_list <- paste(doi_list[, 1], collapse = ",")

write.table(doi_list, file = "doi_list.txt",sep=" ",row.names = FALSE)



