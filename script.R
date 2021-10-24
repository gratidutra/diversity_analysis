# chamando as libs necessárias
pacman::p_load(tidyverse, magrittr, lubridate, iNEXT, vegan, RAM)

df <- read_csv("data/data.csv")

# top 5 specie
top_specie <- df %>% group_by(Espécie) %>%
  summarise(total = n()) %>%
  arrange(-total) %>% slice(1:5)   

# Coletas por mês

collect_by_month <- df %>% group_by(Mês, Espécie) %>%
summarise(total = n()) %>%
arrange(-total) %>%  filter(total >= 12)

ggplot(collect_by_month, aes(Mês, total, fill = total)) +
  geom_bar(binwidth = 1, stat = "identity") +
  theme_light() +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  facet_wrap(~Espécie) +
  ylab("Abundance") +
  xlab("Month") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_polar()

# Coletas por hora

collect_by_hour <- df %>% group_by(hour = hour(`Hora de captura`), Espécie) %>%
summarise(total = n()) %>%
arrange(-total) %>% filter(total > 10)

ggplot(collect_by_hour, aes(hour, total, fill = total)) +
  geom_bar(binwidth = 1, stat = "identity") +
  theme_light() +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  scale_x_continuous(breaks = seq(5, 17, by = 2)) +
  facet_wrap(~Espécie) +
  ylab("Abundance") +
  xlab("Hour") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_polar()

# Coletas por hora - região do cavalo

collect_by_horse <- df %>% group_by(hour = hour(`Hora de captura`), `Local de pouso`) %>%
summarise(total = n()) %>%
arrange(-total) %>% filter(total >= 8)

ggplot(collect_by_horse, aes(as.factor(hour), total, fill = total)) +
  geom_bar(binwidth = 1, stat = "identity") +
  theme_light() +
  scale_fill_gradient(low = "lightblue", high = "blue") +
  facet_wrap(~`Local de pouso`) +
  ylab("Abundance") +
  xlab("Hour") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# NMDS

rich <- df  %>%  group_by(Localidade, Mês) %>%
summarise(rich = n_distinct(Espécie))

data_nmds <- df %>%
group_by(Espécie, Localidade, Mês) %>%
summarise(abund = n()) %>%
pivot_wider(names_from = Espécie, values_from = abund) %>%
mutate(
  across(everything(), ~ replace_na(.x, 0))
) %>% left_join(rich, by = c("Localidade", "Mês"))


run_nmds <- data_nmds

run_nmds$Localidade <- NULL
run_nmds$Mês <- NULL
run_nmds$rich <- NULL

nmds <- metaMDS(run_nmds)

scores(nmds)  %>%
as_tibble() %>%
cbind(data_nmds) %>%
as_tibble()%>%
ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(size = rich, color = Mês)) +
  stat_ellipse(geom = "polygon", aes(group = Mês, color = Mês, fill = Mês), alpha = 0.3) +
  annotate("text", x = -2, y = 0.95, label = paste0("stress: ", format(nmds$stress, digits = 4)), hjust = 0) +
  theme_bw()

# Rich x Abund
rich_abund <- data_nmds %>% 
  pivot_longer(!c(Localidade, Mês, rich), names_to = "Espécie", values_to = "abund") %>% 
  filter(Espécie %in% c("Tabanus antarcticus","Stypommisa aripuana",
                      "Pityocera cervus", "Tabanus occidentalis", 
                      "Tabaus antarcticus"))

ggplot(rich_abund, aes(x = rich, y = abund)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15) +
  facet_wrap(~Espécie)


# Temperatura e umidade

df_climatics <- df %>% group_by(Mês) %>%
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

ggplot(df_climatics, aes(Mês, temp)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_light()

# Diversidade

abund <- df %>%
group_by(Espécie, Localidade) %>%
summarise(total = n()) %>%
pivot_wider(names_from = Localidade, values_from = total)  %>%
replace_na(list(`P1 - Argeu` = 0, `P5-Neucivaldo` = 0)) %>%
column_to_rownames(var = "Espécie")

resultados_tabanidae <- iNEXT(abund,
  q = 0,
  datatype = "abundance",
  endpoint = 800
)

## Diagrama de Venn

count <- df %>% 
  group_by(`Espécie`, Localidade) %>%
  summarise(total = n())

p1 <- count %>%
  filter(Localidade == 'P1 - Argeu') %>%
  select(Localidade, `Espécie`)


p5 <- count %>%
  filter(Localidade == 'P5-Neucivaldo') %>%
  select(Localidade, `Espécie`)


group.venn(list(Argeu=p1$Espécie, Neucivaldo=p5$Espécie), label=FALSE, 
           fill = c("lightpink", "lightblue"),
           cat.pos = c(0, 0),
           cex = 1.8)



## Resultado
resultados_tabanidae$AsyEst

ggiNEXT(resultados_tabanidae, type = 1) + theme_light()
