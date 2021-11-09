# chamando as libs necessárias
pacman::p_load(tidyverse, magrittr, lubridate, iNEXT, vegan, RAM, permute,dplyr, nVennR)

df <- read_csv("../data/data_pcrs.csv")

df <- df %>%
  pivot_longer(!`Espécie`, names_to = "Localidade", values_to = "Abundância")

distinct(df, Localidade)

# top 5 specie

top_specie <- df %>% arrange(-Abundância) %>% slice(1:5)  

# NMDS

df$Abundância[df$Abundância == 0] <- NA

rich <- df  %>%  group_by(Localidade) %>% 
  filter(Abundância > 0) %>%
  summarise(rich = n_distinct(Espécie))

data_nmds <- df %>%
  group_by(Espécie, Localidade) %>%
  pivot_wider(names_from = Espécie, values_from = Abundância) %>%
  mutate(
    across(everything(), ~ replace_na(.x, 0))
  ) %>% left_join(rich, by = "Localidade")


run_nmds <- data_nmds

run_nmds$Localidade <- NULL
run_nmds$rich <- NULL

dist_bray <- vegdist(run_nmds, method = "bray", binary = TRUE)

dist_bray

nmds <- metaMDS(dist_bray)

scores(nmds)  %>%
  as_tibble() %>%
  cbind(data_nmds) %>%
  as_tibble()%>%
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

abund <- df %>%
  pivot_wider(names_from = Localidade, values_from = Abundância)  %>%
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
