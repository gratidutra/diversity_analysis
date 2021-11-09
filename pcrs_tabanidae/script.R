# chamando as libs necessárias
pacman::p_load(tidyverse, magrittr, lubridate, iNEXT, vegan, RAM, permute,dplyr)

df <- read_csv("../data/data_pcrs.csv")

df <- df %>%
  pivot_longer(!`Espécie`, names_to = "Localidade", values_to = "Abundância")

distinct(df, Localidade)

# Diagrama de Venn

count <- df %>% 
  group_by(`Espécie`, Localidade) %>%
  summarise(total = n())

Pelotas <- count %>%
  filter(Localidade == 'Pelotas') %>%
  select(Localidade, `Espécie`)

Corrientes <- count %>%
  filter(Localidade == 'Corrientes') %>%
  select(Localidade, `Espécie`)

Turuçu <- count %>%
  filter(Localidade == 'Turuçu') %>%
  select(Localidade, `Espécie`)

RPPN <- count %>%
  filter(Localidade == 'RPPN') %>%
  select(Localidade, `Espécie`)

Lami <- count %>%
  filter(Localidade == 'Lami') %>%
  select(Localidade, `Espécie`)

Pacheca <- count %>%
  filter(Localidade == 'Pacheca') %>%
  select(Localidade, `Espécie`)

Taim <- count %>%
  filter(Localidade == 'Taim') %>%
  select(Localidade, `Espécie`)

Itapua <- count %>%
  filter(Localidade == 'Itapuã') %>%
  select(Localidade, `Espécie`)

Itapeva <- count %>%
  filter(Localidade == 'Itapeva') %>%
  select(Localidade, `Espécie`)

Lagoa <- count %>%
  filter(Localidade == 'Lagoa do Peixe') %>%
  select(Localidade, `Espécie`)



group.venn(list(Pelotas        = Pelotas$Espécie,
                Corrientes     = Corrientes$Espécie,
                Turuçu         = Turuçu$Espécie,
                RPPN           = RPPN$Espécie,
                Lami           = Lami$Espécie,
                Pacheca        = Pacheca$Espécie,
                Taim           = Taim$Espécie,
                Itapua       = `Itapuã`$Espécie,
                Itapeva        = Itapeva$Espécie,
                Lagoa = `Lagoa do Peixe`$Espécie), label=FALSE, 
           fill = c("lightpink", "lightblue", "green", "red", "yelow", "blue", "pink", "gray", "black", "lightgreen"),
           cat.pos = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
           cex = 1.8)


# top 5 specie

top_specie <- df %>% arrange(-Abundância) %>% slice(1:5)  


# NMDS

df$Abundância[df$Abundância == 0] <- NA

rich <- df  %>%  group_by(Localidade) %>%
  summarise(rich = n_distinct(Espécie))

data_nmds <- df %>%
  group_by(Espécie, Localidade) %>%
  summarise(abund = n()) %>%
  pivot_wider(names_from = Espécie, values_from = abund) %>%
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
  geom_point(aes(size = rich, color = Mês)) +
  stat_ellipse(geom = "polygon", aes(group = Mês, color = Mês, fill = Mês), alpha = 0.3) +
  annotate("text", x = -2, y = 0.95, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  theme_bw()

# Permanova

adonis(dist_bray~data_nmds$Localidade, permutations = 1000)

# Rich x Abund
rich_abund <- data_nmds %>% 
  pivot_longer(!c(Localidade, Mês, rich), names_to = "Espécie", values_to = "abund") %>% 
  filter(Espécie %in% top_specie$Espécie)

ggplot(rich_abund, aes(x = rich, y = abund)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15) +
  facet_wrap(~Espécie)



# Rich x Abund

data_rich <- df %>%
  group_by(Localidade) %>%
  summarise(rich = n_distinct(`Espécie`))




# Diversidade

abund <- df %>%
  group_by(Espécie, Localidade) %>%
  summarise(total = n()) %>%
  pivot_wider(names_from = Localidade, values_from = total)  %>%
  replace_na(list(`P1 - Argeu` = 0, `P2 - Hélio` = 0, `P3 - Paulo Cabeça` = 0,
                  `P4 - Paulo Vicentino` = 0, `P5 - Nelcivaldo` = 0)) %>%
  column_to_rownames(var = "Espécie")

resultados_tabanidae <- iNEXT(abund,
                              q = c(0, 1, 2),
                              datatype = "abundance",
                              endpoint = 400
)




## Resultado
resultados_tabanidae$AsyEst

ggiNEXT(resultados_tabanidae, type = 1, facet.var = 'order') + theme_light()
