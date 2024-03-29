# chamando as libs necessárias
pacman::p_load(tidyverse, magrittr, lubridate, iNEXT, vegan, RAM, permute)

df <- read_csv("data/data.csv")
df <-
  df %>%
  mutate(
    `Espécie` = case_when(
      `Espécie` == "Tabanus antarcticus" ~ "Tabanus antarticus",
      `Espécie` == "Tabaus antarcticus" ~ "Tabanus antarticus",
      `Espécie` == "Tabanus occidentalis var. dorsovittatus" ~ "Tabanus occidentalis",
      `Espécie` == "Tabanus occidentalis var. modestus" ~ "Tabanus occidentalis",
      `Espécie` == "Tabaus occidentalis var. modestus" ~ "Tabanus occidentalis",
      `Espécie` == "Tabanus occidentalis var.modestus" ~ "Tabanus occidentalis",
      `Espécie` == "Tabanus occidentalis var. ?dorsovittatus" ~ "Tabanus occidentalis",
      `Espécie` == "Tabanus rupripes" ~ "Tabanus rubripes",
      `Espécie` == "Tabanus ?fuscofasciatus" ~ "Tabanus fuscofasciatus",
      `Espécie` == "Tabanus occiedntalis" ~ "Tabanus occidentalis",
      TRUE ~ `Espécie`
    ),
    `Estação` = case_when(
      `Estação` == "Chuva" ~ "Rainy",
      TRUE ~ "Dry"
    ),
    Mês = case_when(
      Mês == "Agosto" ~ "August",
      Mês == "Fevereiro" ~ "February",
      Mês == "Junho" ~ "June",
      Mês == "Maio" ~ "May",
      Mês == "Novembro" ~ "November"
    )
  )

# Diagrama de Venn

count <- df %>%
  group_by(`Espécie`, Localidade) %>%
  summarise(total = n())

p1 <- count %>%
  filter(Localidade == "P1 - Argeu") %>%
  select(Localidade, `Espécie`)


p5 <- count %>%
  filter(Localidade == "P5-Neucivaldo") %>%
  select(Localidade, `Espécie`)


group.venn(list(Argeu = p1$Espécie, Neucivaldo = p5$Espécie),
  label = FALSE,
  fill = c("lightgray", "darkgray"),
  cat.pos = c(0, 0),
  cex = 2
)

# top 5 specie

top_specie <-
  df %>%
  group_by(Espécie) %>%
  summarise(total = n()) %>%
  arrange(-total) %>%
  slice(1:5)

# Coletas por mês

collect_by_month <-
  df %>%
  group_by(Mês, Espécie, `Estação`) %>%
  summarise(total = n()) %>%
  arrange(-total) %>%
  filter(Espécie %in% c('Tabanus antarticus','Stypommisa aripuana', 'Pityocera cervus'))

collect_by_month$Mês <-
  factor(collect_by_month$Mês,
    levels = c("November", "February", "May", "June", "August")
  )

ggplot(collect_by_month, aes(Mês, total, fill = `Estação`)) +
  geom_bar(binwidth = 1, stat = "identity") +
  theme_light() +
  scale_fill_manual(values = c("lightgray", "black")) +
  # scale_fill_gradient(low = 'lightgray', high = 'black') +
  facet_wrap(~Espécie) +
  ylab("Abundance") +
  xlab("Month") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_polar()

# Coletas por hora

collect_by_hour <-
  df %>%
  group_by(hour = hour(`Hora de captura`), Espécie) %>%
  summarise(total = n()) %>%
  arrange(-total) %>%
  filter(total > 10)

ggplot(collect_by_hour, aes(hour, total, fill = total)) +
  geom_bar(binwidth = 1, stat = "identity") +
  theme_light() +
  scale_fill_gradient(low = "lightgray", high = "black") +
  scale_x_continuous(breaks = seq(5, 17, by = 2)) +
  facet_wrap(~Espécie) +
  ylab("Abundance") +
  xlab("Hour") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_polar()

# Coletas por hora - região do cavalo

collect_by_horse <-
  df %>%
  group_by(hour = hour(`Hora de captura`), `Local de pouso`) %>%
  summarise(total = n()) %>%
  arrange(-total) %>%
  filter(total >= 8)

ggplot(collect_by_horse, aes(as.factor(hour), total, fill = total)) +
  geom_bar(binwidth = 1, stat = "identity") +
  theme_light() +
  scale_fill_gradient(low = "white", high = "black") +
  facet_wrap(~`Local de pouso`) +
  ylab("Abundance") +
  xlab("Hour") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# Coletas por estação

collect_by_station <-
  df %>%
  group_by(`Estação`, Espécie) %>%
  summarise(total = n()) %>%
  arrange(-total) %>%
  filter(total >= 12)

ggplot(collect_by_station, aes(`Estação`, total, fill = `Estação`)) +
  geom_bar(binwidth = 1, stat = "identity") +
  theme_light() +
  scale_fill_manual(values = c("lightgray", "black")) +
  facet_wrap(~Espécie) +
  ylab("Abundance") +
  xlab("Season") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# NMDS

rich <-
  df %>%
  group_by(Localidade, Mês) %>%
  summarise(rich = n_distinct(Espécie))

data_nmds <-
  df %>%
  group_by(Espécie, Localidade, Mês) %>%
  summarise(abund = n()) %>%
  pivot_wider(names_from = Espécie, values_from = abund) %>%
  replace(is.na(.), 0) %>%
  left_join(rich, by = c("Localidade", "Mês"))

run_nmds <-
  data_nmds

run_nmds$Localidade <- NULL
run_nmds$Mês <- NULL
run_nmds$rich <- NULL

dist_bray <- vegdist(run_nmds, method = "bray", binary = TRUE)

dist_bray

nmds <- metaMDS(dist_bray)

scores(nmds) %>%
  as_tibble() %>%
  cbind(data_nmds) %>%
  as_tibble() %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(size = rich, color = Mês)) +
  stat_ellipse(geom = "polygon", aes(group = Mês, color = Mês, fill = Mês), alpha = 0.3) +
  annotate("text", x = -2, y = 0.95, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  theme_bw()

# Permanova

adonis(dist.jac ~ data_nmds$Localidade, permutations = 1000)

# Rich x Abund
rich_abund <-
  data_nmds %>%
  pivot_longer(!c(Localidade, Mês, rich), names_to = "Espécie", values_to = "abund") %>%
  filter(Espécie %in% c(
    "Tabanus antarcticus", "Stypommisa aripuana",
    "Pityocera cervus", "Tabanus occidentalis",
    "Tabaus antarcticus"
  ))

ggplot(rich_abund, aes(x = rich, y = abund)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15) +
  facet_wrap(~Espécie)


# Temperatura e umidade

df_climatics <-
  df %>%
  group_by(Mês) %>%
  summarise(
    total = n(),
    temp = mean(`Temperatura (°C)`),
    um = mean(`Umidade (%)`)
  )

ggplot(df_climatics, aes(x = temp, y = total)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15)

ggplot(df_climatics, aes(x = um, y = total)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15)

# Rich x Abund

data_rich <-
  df %>%
  group_by(Localidade, Mês) %>%
  summarise(
    "Temp" = mean(`Temperatura (°C)`),
    "um" = mean(`Umidade (%)`),
    rich = n_distinct(`Espécie`)
  )

data_rich

# Rich x Temp

ggplot(data_rich, aes(x = Temp, y = rich)) +
  geom_point() +
  geom_smooth(method = "glm", family = "poisson", se = TRUE)

# Modelo 1

model1 <- glm(data = data_rich, rich ~ Mês * Temp * um, family = poisson)

summary(model1)

anova(model1)

# Diversidade

abund <-
  df %>%
  group_by(Espécie, Estação) %>%
  summarise(total = n()) %>%
  pivot_wider(names_from = Estação, values_from = total) %>%
  replace_na(list(`Dry` = 0, Rainy = 0)) %>%
  column_to_rownames(var = "Espécie")

resultados_tabanidae <-
  iNEXT(abund,
    q = c(0, 1, 2),
    datatype = "abundance",
    endpoint = 800
  )

## Resultado

ggiNEXT(resultados_tabanidae,
  facet.var = "site",
  color.var = "order"
) +
  theme_light()


# Partes cavalo -----------------------------------------------------------
df <-
  read_csv("partes_cavalo.csv")
df <-
  df %>%
  mutate(`Espécie` = case_when(
    `Espécie` == "Tabanus antarcticus" ~ "Tabanus antarticus",
    `Espécie` == "Tabaus antarcticus" ~ "Tabanus antarticus",
    `Espécie` == "Tabanus occidentalis var. dorsovittatus" ~ "Tabanus occidentalis",
    `Espécie` == "Tabanus occidentalis var. modestus" ~ "Tabanus occidentalis",
    `Espécie` == "Tabaus occidentalis var. modestus" ~ "Tabanus occidentalis",
    `Espécie` == "Tabanus occidentalis var.modestus" ~ "Tabanus occidentalis",
    `Espécie` == "Tabanus occidentalis var. ?dorsovittatus" ~ "Tabanus occidentalis",
    `Espécie` == "Tabanus rupripes" ~ "Tabanus rubripes",
    `Espécie` == "Tabanus ?fuscofasciatus" ~ "Tabanus fuscofasciatus",
    `Espécie` == "Tabanus occiedntalis" ~ "Tabanus occidentalis",
    TRUE ~ `Espécie`
  ))

unique(df$`Local de pouso_1`)

df <-
  df %>%
  mutate(local_do_pouso_simplificado = case_when(
    `Local de pouso_1` == "Face" |
      `Local de pouso_1` == "Pescoço" |
      `Local de pouso_1` == "Narina" ~ "Head",
    `Local de pouso_1` == "Barriga" |
      `Local de pouso_1` == "Aparelho genital" ~ "Ventral",
    `Local de pouso_1` == "Rabo" |
      `Local de pouso_1` == "Quartela traseira" ~ "Lumbar",
    `Local de pouso_1` == "Perna dianteira" |
      `Local de pouso_1` == "Perna traseira" ~ "Legs",
    TRUE ~ "Frontal"
  ))
unique(df$local_do_pouso_simplificado)

collect_by_horse <-
  df %>%
  group_by(local_do_pouso_simplificado) %>%
  summarise(total = n()) %>%
  arrange(-total) %>%
  filter(total >= 8)

ggplot(df, aes(as.factor(local_do_pouso_simplificado),
  fill = as.factor(local_do_pouso_simplificado)
)) +
  geom_bar() +
  theme_light() +
  scale_fill_grey(start = 0.25, end = 0.75) +
  ylab("Abundance") +
  xlab("Horse Region") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

