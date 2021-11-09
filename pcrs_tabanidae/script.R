# chamando as libs necessárias
pacman::p_load(tidyverse, magrittr, lubridate, iNEXT, vegan, permute, dplyr, networkD3)

df <- read_csv("../data/data_pcrs.csv")

df <- df %>%
  pivot_longer(!`Espécie`, names_to = "Localidade", values_to = "Abundância")

distinct(df, Localidade)

# top 5 specie

top_specie <- df %>% arrange(-Abundância) %>% slice(1:5)  


# Sankey


nodes <- data.frame(name = c(as.character(df$Espécie), as.character(df$Localidade)) %>% unique())


df$IDsource <- match(df$Espécie, nodes$name) - 1
df$IDtarget <- match(df$Localidade, nodes$name) - 1


ColourScal <- 'd3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF",
"#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'


sankeyNetwork(
  Links = df, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",
  Value = "Abundância", NodeID = "name",
  sinksRight = FALSE, colourScale = ColourScal, nodeWidth = 40, fontSize = 13, nodePadding = 20
)


# NMDS

rich <- df %>%
  group_by(Localidade) %>%
  filter(Abundância > 0) %>%
  summarise(rich = n_distinct(Espécie))

data_nmds <- df %>%
  group_by(Espécie, Localidade) %>%
  pivot_wider(names_from = Espécie, values_from = Abundância) %>%
  mutate(
    across(everything(), ~ replace_na(.x, 0))
  ) %>%
  left_join(rich, by = "Localidade")


run_nmds <- data_nmds

run_nmds$Localidade <- NULL
run_nmds$rich <- NULL

dist_bray <- vegdist(run_nmds, method = "bray", binary = TRUE)

dist_bray

nmds <- metaMDS(dist_bray)

scores(nmds) %>%
  as_tibble() %>%
  cbind(data_nmds) %>%
  as_tibble() %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(size = rich)) +
  stat_ellipse(geom = "polygon", aes(group = rich), alpha = 0.3) +
  annotate("text", x = -2, y = 0.95, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  theme_bw()

# Permanova

adonis(dist_bray~data_nmds$Localidade, permutations = 1000)

# Rich x Abund
rich_abund <- data_nmds %>% 
  pivot_longer(!c(Localidade, rich), names_to = "Espécie", values_to = "abund") %>% 
  filter(Espécie %in% top_specie$Espécie)

ggplot(rich_abund, aes(x = rich, y = abund)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15) +
  facet_wrap(~Espécie)

# Diversidade

abund <- read_csv("../data/data_pcrs.csv") %>%
    column_to_rownames(var = "Espécie")

resultados_tabanidae <- iNEXT(abund,
                              q = c(0, 1, 2),
                              datatype = "abundance",
                              endpoint = 400
)


# Resultado 0 = Species richness; 1 = Shannon diversity; 2 = Simpson diversity

resultados_tabanidae$AsyEst

ggiNEXT(resultados_tabanidae, type = 1, facet.var = 'site') + theme_light() + 
  facet_wrap(~site, ncol=3, scales = 'free') 
