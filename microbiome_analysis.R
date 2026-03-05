library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
set.seed(3208940)

# Read tropical forest soil TSV file with OTU counts, delete sites with less than 558 samples and delete taxa without any occurences
tropical <- read_tsv("tropical.tsv", skip = 1) %>%
  rename(OTU_ID = "#OTU ID") %>% 
  pivot_longer(-OTU_ID, names_to = "site", values_to = "values") %>%
  group_by(site) %>% 
  mutate(total = sum(values)) %>% 
  filter(total > 558) %>%
  group_by(OTU_ID) %>% 
  mutate(total = sum(values)) %>% 
  filter(total != 0) %>% 
  ungroup %>% 
  select(-total)

# For setting the 558 boundary
# sort(unique(tropical$total))
#ggplot(tropical, aes(x=total)) + 
#  geom_histogram()

temperate <- read_tsv("temperate.tsv", skip = 1) %>%
  rename(OTU_ID = "#OTU ID") %>% 
  pivot_longer(-OTU_ID, names_to = "site", values_to = "values") %>%
  group_by(site) %>% 
  mutate(total = sum(values)) %>% 
  filter(total > 558) %>%
  group_by(OTU_ID) %>% 
  mutate(total = sum(values)) %>% 
  filter(total != 0) %>% 
  ungroup %>% 
  select(-total)

# For setting the 558 boundary
#sort(unique(temperate$total))

# custom functions to calculate alpha diversity based on Riffomonas Project's Pat Schloss: https://www.youtube.com/watch?v=wq1SXGQYgCs
richness <- function(x){
  sum(x>0)
}
shannon_ind <- function(x){
  - sum(x[x>0]/sum(x) * log (x[x>0]/sum(x)))
}
simp_ind <- function(x){
  sum((x * (x - 1)) / (sum(x) * (sum(x)-1) ))
}

tropical <- tropical %>% 
  group_by(site) %>% 
  summarise(sobs = richness(values),
            veg_sobs = specnumber(values),
            myshannon = shannon_ind(values),
            veg_shannon = diversity(values, index = "shannon"),
            simpson = simp_ind(values),
            invsimpson = 1/simpson,
            n = sum(values)) %>% 
  mutate(biome = "Tropical")

temperate <- temperate %>% 
  group_by(site) %>% 
  summarise(sobs = richness(values),
            veg_sobs = specnumber(values),
            myshannon = shannon_ind(values),
            veg_shannon = diversity(values, index = "shannon"),
            simpson = simp_ind(values),
            invsimpson = 1/simpson,
            n = sum(values)) %>% 
  mutate(biome = "Temperate")

forests <- bind_rows(temperate, tropical)

ggplot(forests, aes(x = biome, y = myshannon, fill = biome))+
  geom_violin(alpha = 0.7) +
  geom_point()

ggplot(forests, aes(x = biome, y = invsimpson, fill = biome))+
  geom_violin() +
  geom_point()
