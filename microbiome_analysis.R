library(tidyverse)
library(phyloseq)
library(vegan)
library(ggplot2)
library(DESeq2)
library(ggtree)
library(car)
library(microbiome)
library(ggpubr)

# Read taxonomy table of all OTUs present in soil forest samples, then transform to matrix for phyloseq
taxonomy <- read_tsv("taxonomy.tsv") %>%
  dplyr::rename(OTU_ID = "#OTU_ID") %>% 
  column_to_rownames("OTU_ID") %>% 
  as.matrix()

# Read phylogenetic tree of all OTUs present in Earth Microbiome Project (EMP) samples
tree <- read.tree("silva_123.97_otus.tre")

# Read tropical and temperate forest soil TSV files with OTU counts, delete sites with less than 558 samples and delete taxa without any occurence
tropical <- read_tsv("tropical.tsv", skip = 1) %>%
  dplyr::rename(OTU_ID = "#OTU ID") %>%
  pivot_longer(-OTU_ID, names_to = "site", values_to = "values") %>%
  group_by(site) %>% 
  mutate(total = sum(values)) %>% 
  filter(total > 558) %>%
  group_by(OTU_ID) %>% 
  mutate(total = sum(values)) %>% 
  filter(total != 0) %>% 
  ungroup %>% 
  select(-total) %>% 
  pivot_wider(names_from = "site", values_from = "values")

## Optional code to help setting the 558 boundary
# sort(unique(tropical$total))
# ggplot(tropical, aes(x=total)) + 
#  geom_histogram()

temperate <- read_tsv("temperate.tsv", skip = 1) %>%
  dplyr::rename(OTU_ID = "#OTU ID") %>% 
  pivot_longer(-OTU_ID, names_to = "site", values_to = "values") %>%
  group_by(site) %>% 
  mutate(total = sum(values)) %>% 
  filter(total > 558) %>%
  group_by(OTU_ID) %>% 
  mutate(total = sum(values)) %>% 
  filter(total != 0) %>% 
  ungroup %>% 
  select(-total) %>% 
  pivot_wider(names_from = "site", values_from = "values")

## For setting the 558 boundary
#sort(unique(temperate$total))

# Read metadata tsv file and filter for samples present in tropical and temperate matrices
tropical_samples <- as.data.frame(read_tsv("tropical_metadata.tsv")) %>% 
  dplyr::rename(Sample_ID = "#SampleID") %>%
  subset(Sample_ID %in% colnames(tropical), select = c("Sample_ID", "Description", "title")) %>% 
  mutate(biome = "Tropical")
rownames(tropical_samples) <- tropical_samples[,1]
tropical_samples[,1] <- NULL

temperate_samples <- as.data.frame(read_tsv("temperate_metadata.tsv"))%>% 
  dplyr::rename(Sample_ID = "#SampleID") %>%
  subset(Sample_ID %in% colnames(temperate), select = c("Sample_ID", "Description", "title")) %>% 
  mutate(biome = "Temperate")
rownames(temperate_samples) <- temperate_samples[,1]
temperate_samples[,1] <- NULL

# Combine datasets as phyloseq object
otu_combined <- full_join(tropical, temperate, by = "OTU_ID") %>% 
  column_to_rownames("OTU_ID") %>% 
  as.matrix()
otu_combined[is.na(otu_combined)] <- 0

combined_samples <- bind_rows(tropical_samples, temperate_samples)

ps <- phyloseq(
  otu_table(otu_combined, taxa_are_rows = TRUE),
  tax_table(taxonomy),
  sample_data(combined_samples),
  phy_tree(tree)
)

## To optionally plot a phylogenetic tree or heatmap
#ps_filtered = filter_taxa(ps, function(x) mean(x) > 0.5, TRUE)
#ps_filtered <- prune_taxa(taxa_sums(tropical_ps) > 10000, tropical_ps)
#ggtree(phy_tree(tropical_ps_filtered), layout="circular") + geom_tiplab2()
#plot_heatmap(ps, sample.label = "Sample_ID", species.label= "Phylum")

# Rarefying the data
sort(sample_sums(ps))
rare_depth = 10727
psrar <- rarefy_even_depth(ps,
                           sample.size = rare_depth,
                           rngseed = 42,
                           replace = FALSE,
                           trimOTUs = TRUE
                           )

## Uncomment to optionally plot rarefaction curves
# rarecurve(t(as(otu_table(physeq_filtered), "matrix")), step = 500, cex = 0.5, label = FALSE)
# rarecurve(t(as(otu_table(psrar), "matrix")), step = 500, cex = 0.5, label = FALSE)

# Calculating Bray-Curtis dissimilarities to visualize rarefaction. Based on Riffomonas Project & Pat Schloss: https://www.youtube.com/watch?v=ht3AX5uZTTQ
otu_matrix = t(as(otu_table(ps), "matrix"))
dist <- vegdist(otu_matrix, method="bray")
rare_dist <- avgdist(otu_matrix, dmethod="bray", sample = rare_depth)

dist_tibble <- dist %>% 
  as.matrix() %>% 
  as_tibble(rownames="sample") %>% 
  pivot_longer(-sample) %>% 
  filter(name < sample)

rare_dist_tibble <- rare_dist %>% 
  as.matrix() %>% 
  as_tibble(rownames="sample") %>% 
  pivot_longer(-sample) %>% 
  filter(name < sample)

inner_join(dist_tibble, rare_dist_tibble, by = c("sample", "name")) %>%
  select(sample, name, no_rare = value.x, rare = value.y) %>%
  mutate(biomes = case_when(
    sample %in% rownames(tropical_samples) & name %in% rownames(tropical_samples) ~ "Tropical-Tropical",
    sample %in% rownames(temperate_samples) & name %in% rownames(temperate_samples) ~ "Temperate-Temperate",
    TRUE ~ "Tropical-Temperate"
  )) %>%
  ggplot(aes(x = no_rare, y = rare, color = biomes)) +
  geom_point(alpha = 0.5) + 
  xlab("Bray-Curtis dissimilarities without rarefaction") +
  ylab("Bray-Curtis dissimilarities with rarefaction") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")

# Apply custom functions to calculate alpha diversity based on Riffomonas Project & Pat Schloss: https://www.youtube.com/watch?v=wq1SXGQYgCs
richness <- function(x){
  sum(x>0)
}

shannon_ind <- function(x){
  - sum(x[x>0]/sum(x) * log (x[x>0]/sum(x)))
}

simp_ind <- function(x){
  sum((x * (x - 1)) / (sum(x) * (sum(x)-1) ))
}

otu_df <- as.data.frame(otu_table(psrar)) %>%
  rownames_to_column(var = "OTU_ID") %>%
  pivot_longer(-OTU_ID, names_to = "site", values_to = "values") %>%
  group_by(site) %>% 
  summarise(sobs = richness(values),
            veg_sobs = specnumber(values),
            myshannon = shannon_ind(values),
            veg_shannon = diversity(values, index = "shannon"),
            simpson = simp_ind(values),
            invsimpson = 1/simpson,
            chao = estimateR(values)[2],
            n = sum(values)) %>% 
  mutate(biome = case_when(
    site %in% rownames(tropical_samples) ~ "Tropical",
    site %in% rownames(temperate_samples) ~ "Temperate"
  ))

# Testing normality, homoscedasticity and running unpaired Wilcoxon Test (https://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/)
mod_shannon <- lm(myshannon ~ biome, data = otu_df)
mod_sobs    <- lm(sobs ~ biome, data = otu_df)
mod_chao    <- lm(chao ~ biome, data = otu_df)

shapiro_shannon <- shapiro.test(residuals(mod_shannon))
shapiro_sobs    <- shapiro.test(residuals(mod_sobs))
shapiro_chao    <- shapiro.test(residuals(mod_chao))

print(shapiro_shannon)  # violates normality
print(shapiro_sobs)     # violates normality
print(shapiro_chao)     # violates normality

leveneTest(myshannon ~ biome, data = otu_df)  # violates homoscedasticity
leveneTest(sobs ~ biome, data = otu_df)     # violates homoscedasticity
leveneTest(chao ~ biome, data = otu_df)         # violates homoscedasticity

wilcox.test(sobs ~ biome, data = otu_df)         # p-value = 0.0006068
wilcox.test(myshannon ~ biome, data = otu_df)    # p-value = 1.318e-05
wilcox.test(invsimpson ~ biome, data = otu_df)   # p-value = 8.495e-05
wilcox.test(chao ~ biome, data = otu_df)         # p-value = 0.03048

# Plot Alpha diversity metrics side-by-side
alpha_plot <- function(df, y_var, y_label) {
  ggplot(df, aes(x = biome, y = .data[[y_var]], fill = biome)) +
    geom_violin(alpha = 0.5, show.legend = FALSE, adjust = 0.4, linewidth = 1) +
    geom_boxplot(show.legend = FALSE, outlier.shape = NA, alpha = 0.9, 
                 width = 0.15, fill = "white", linewidth = 1) +
    stat_summary(fun = median, show.legend = FALSE,
                 geom = "crossbar", alpha = 0.7, width = 1, linewidth = 0.75) +
    stat_compare_means(method = "wilcox.test", label.x = 1.3) +
    labs(x = NULL, y = y_label) +
    theme_classic()
}

F_sobs <- alpha_plot(otu_df, "sobs", "Σ OTU observed")
F_chao <- alpha_plot(otu_df, "chao", "Chao1 Index")
F_shan <- alpha_plot(otu_df, "myshannon", "Shannon Index")
F_simp <- alpha_plot(otu_df, "invsimpson", "Inverse Simpson Index")

(F_sobs | F_chao) / (F_shan | F_simp)
ggsave("AlphaDiversity.png", width = 10, height = 8)

# Show taxonomic composition by plotting relative abundance by biome
tax.table <- tax_table(ps)
ps_glom <- tax_glom(ps,taxrank = "Phylum", 
                       NArm = FALSE)

top_tib <- ps_glom %>%
  psmelt() %>%
  group_by(Phylum) %>% 
  summarise(total = sum(Abundance)) %>% 
  arrange(desc(total))

top15 <- top_tib$Phylum[1:15]

top15_plot <- ps_glom %>%
  psmelt() %>%
  mutate(Phylum = ifelse(Phylum %in% top15, as.character(Phylum), "Other")) %>%
  mutate(Kingdom = ifelse(Phylum %in% top15, as.character(Kingdom), "Archaea/Bacteria")) %>%
  mutate(Kingdom_Phylum = paste0(Kingdom, ": ", Phylum)) %>%
  mutate(Kingdom_Phylum = fct_reorder(Kingdom_Phylum, Abundance, .fun = sum, .desc = TRUE)) %>%
  ggplot(aes(x = biome, y = Abundance, fill = Kingdom_Phylum)) +
  labs(fill = "Kingdom: Phylum") +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_viridis_d(option = "turbo") +
  labs(x = "Biome", y = "Relative Abundance") + 
  theme_minimal()

top15_plot
ggsave("top15_plot.png", width = 10, height = 8)

# Beta diversity analysis based on:
# https://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
# https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/beta-diversity-metrics.html
# https://r-from-a-learners-perspective.readthedocs.io/en/latest/part5/

# Beta diversity: Unweighted UniFrac
ordu.unwt.uni <- ordinate(psrar , "PCoA", "unifrac", weighted=F)

eigenvalues_unwt <- ordu.unwt.uni$values$Relative_eig[1:10] * 100
barplot(eigenvalues_unwt, ylab = "% Variance Explained", names.arg = paste0("PC", 1:10))

unwt.unifrac <- plot_ordination(psrar, 
                              ordu.unwt.uni, color="biome") +
  ggtitle("Unweighted UniFrac") +
  geom_point(size = 2) +
  theme_minimal() +
  scale_color_viridis_d(option = "viridis")

unwt.unifrac + stat_ellipse()
ggsave("unwt_unifrac.png", width = 10, height = 8)

# Beta diversity: Weighted UniFrac
ordu.wt.uni <- ordinate(psrar , "PCoA", "unifrac", weighted=T)

eigenvalues <- ordu.wt.uni$values$Relative_eig[1:10] * 100
barplot(eigenvalues, ylab = "% Variance Explained", names.arg = paste0("PC", 1:10))

wt.unifrac <- plot_ordination(psrar, 
                              ordu.wt.uni, color="biome") +
  ggtitle("Weighted UniFrac") +
  geom_point(size = 2) +
  theme_minimal() +
  scale_color_viridis_d(option = "viridis")

wt.unifrac + stat_ellipse()
ggsave("wt_unifrac.png", width = 10, height = 8)

# Beta diversity: NMDS based on Bray-Curtis
# Ordination

ps_nmds <- ordinate(
  physeq = psrar, 
  method = "NMDS", 
  distance = "bray"
)

nmds_plot <- plot_ordination(
  physeq = psrar,
  ordination = ps_nmds,
  color = "biome",
  title = "NMDS of Temperate vs Tropical Microbial Communities") +
  scale_color_viridis_d(option = "viridis") +
  geom_point(aes(color = biome), size = 2)

nmds_plot
ggsave("nmds_plot.png", width = 10, height = 8)

# Running PERMANOVA via adonis and testing homogeneity via betadisper
metadf <- data.frame(sample_data(psrar))

unifrac.dist <- UniFrac(psrar, 
                        weighted = TRUE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)

permanova <- adonis2(unifrac.dist ~ biome, data = metadf)

permanova #was significant (p < 0.001)

psrar.disper <- betadisper(unifrac.dist, metadf$biome)
permutest(psrar.disper, pairwise = TRUE) #was not significant (p > 0.8)

# Differential Abundance using DESeq2
ds2 <- phyloseq_to_deseq2(ps, ~ biome)
dds <- DESeq(ds2, test = "Wald", fitType = "local", sfType = "poscounts")
