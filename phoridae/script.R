# chamando as libs necessárias
pacman::p_load(tidyverse, magrittr, lubridate, iNEXT, vegan, RAM, permute)

df <- read_csv("../data/Phoridae.csv") %>% 
  pivot_longer(!c(points,area, rivers), names_to = "Espécie", values_to = "Abundância") 

df <- df %>% group_by(area, rivers, Espécie) %>%
  summarise(Abundância = sum(Abundância))

distinct(df, area)

# top 5 specie

top_specie <- df %>% arrange(-Abundância) %>% 
  filter(Abundância > 219)  

# NMDS

rich <- df  %>%  
  group_by(area) %>% 
  filter(Abundância > 0) %>%
  summarise(rich = n_distinct(Espécie))

data_nmds <- df %>%
  group_by(Espécie, area) %>%
  pivot_wider(names_from = Espécie, values_from = Abundância) %>%
  mutate(
    across(everything(), ~ replace_na(.x, 0))
  ) %>% left_join(rich, by = "area")


run_nmds <- data_nmds

run_nmds$area <- NULL
run_nmds$rivers <- NULL
run_nmds$rich <- NULL

dist_bray <- vegdist(run_nmds, method = "bray", binary = TRUE)

nmds <- metaMDS(dist_bray)

scores(nmds)  %>%
  as_tibble() %>%
  cbind(data_nmds) %>%
  as_tibble()%>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(size = rich, color = rivers)) +
  stat_ellipse(geom = "polygon", aes(group = rivers, color = rivers, fill = rivers), alpha = 0.3) +
  annotate("text", x = -2, y = 0.95, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  theme_bw()

# Permanova

adonis(dist_bray~data_nmds$area, permutations = 1000)

adonis(dist_bray~data_nmds$rivers, permutations = 1000)

# Rich x Abund

rich_abund <- data_nmds %>% 
  pivot_longer(!c(area, rich, rivers), names_to = "Espécie", values_to = "abund") %>% 
  filter(Espécie %in% top_specie$Espécie)

ggplot(rich_abund, aes(x = rich, y = abund)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15) +
  facet_wrap(~Espécie)


# Diversidade

abund <- df %>%
  group_by(Espécie, area) %>%
  summarise(Abundância = sum(Abundância)) %>%
  pivot_wider(names_from = area, values_from = Abundância)  %>%
  replace_na(list(`ARE` = 0, `COR` = 0, `DCO` = 0, `FOZ` =  0, `PEL` = 0,
                  `TUR` = 0)) %>%
  column_to_rownames(var = "Espécie")

resultados <- iNEXT(abund,
  q = c(0, 1, 2),
  datatype = "abundance",
  endpoint = 400
)

## Resultado

ggiNEXT(resultados, type = 1, facet.var = 'site') + theme_light() + 
  facet_wrap(~site, ncol=3, scales = 'free') 
