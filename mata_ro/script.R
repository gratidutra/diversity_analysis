# chamando as libs necessárias
pacman::p_load(tidyverse, magrittr, lubridate, iNEXT, vegan, RAM, permute, dplyr)

df <- 
  read_csv("../data/data_mata.csv")

df <-
  df %>% mutate(
    `Espécie` =
      case_when(
        `Espécie` == "Tabanus antarcticus" ~ "Tabanus antarticus",
        `Espécie` == "Tabaus antarcticus" ~ "Tabanus antarticus",
        `Espécie` == "Tabanus occidentalis var. dorsovittatus" ~ "Tabanus occidentalis",
        `Espécie` == "Tabanus occidentalis var. modestus" ~ "Tabanus occidentalis",
        `Espécie` == "Tabanus occidentalis var.modestus" ~ "Tabanus occidentalis",
        `Espécie` == "Tabanus occidentalis var. ?dorsovittatus" ~ "Tabanus occidentalis",
        `Espécie` == "Tabanus rupripes" ~ "Tabanus rubripes",
        `Espécie` == "Tabanus ?fuscofasciatus" ~ "Tabanus fuscofasciatus",
        `Espécie` == "Tabanus occiedntalis" ~ "Tabanus occidentalis",
        `Espécie` == "Tabanus ?cicur" ~ "Tabanus circus",
        TRUE ~ `Espécie`
      ), 
    Localidade = case_when(
      Localidade == "P2- Hélio" ~ "P2 - Hélio",
      TRUE ~ Localidade
    ),
    Mês = case_when(
      Mês == "Janeiro" ~ "January",
      Mês == "Fevereiro" ~ "February",
      Mês == "Março" ~ "March",
      Mês == "Abril" ~ "April",
      Mês == "Maio" ~ "May",
      Mês == "Junho" ~ "June",
      Mês == "Julho" ~ "July",
      Mês == "Agosto" ~ "August",
      Mês == "Setembro" ~ "September",
      Mês == "Outubro" ~ "October",
      Mês == "Novembro" ~ "November",
      Mês == "Dezembro" ~ "December"
    )
  )

df <- 
  df %>% filter(`Espécie` != 'Família Syrphidae')

# Diagrama de Venn

count <- 
  df %>%
  group_by(`Espécie`, Localidade) %>%
  summarise(total = n())

p1 <- 
  count %>%
  filter(Localidade == "P1 - Argeu") %>%
  select(Localidade, `Espécie`)

p2 <- 
  count %>%
  filter(Localidade == "P2- Hélio") %>%
  select(Localidade, `Espécie`)

p3 <- 
  count %>%
  filter(Localidade == "P3 - Paulo Cabeça") %>%
  select(Localidade, `Espécie`)

p4 <-
  count %>%
  filter(Localidade == "P4 - Paulo Vicentino") %>%
  select(Localidade, `Espécie`)

p5 <- 
  count %>%
  filter(Localidade == "P5 - Nelcivaldo") %>%
  select(Localidade, `Espécie`)

group.venn(list(Argeu = p1$Espécie, `Hélio` = p2$Espécie, `Paulo Cabeça` = p3$Espécie, `Paulo Vicentino` = p4$Espécie, Neucivaldo = p5$Espécie),
  label = FALSE,
  fill = c("lightpink", "lightblue", "green", "red", "yelow"),
  cat.pos = c(1, 1, 1, 1, 1),
  cex = 1.8
)
# top 5 specie

top_specie <- df %>% group_by(Espécie) %>%
  summarise(total = n()) %>%
  arrange(-total) %>% slice(1:5)  

# Coletas por mês

collect_by_month <- 
  df %>% group_by(Mês, Espécie) %>%
  summarise(total = n()) %>%
  arrange(-total) %>%  filter(total >= 12)

collect_by_month$Mês <-
  factor(collect_by_month$Mês,
         levels = c("November", "December", "January", 
                    "February", "March", "April", "May",
                    "June", "July", "August", "September",
                    "October")
  )

ggplot(collect_by_month, aes(Mês, total, fill = total)) +
  geom_bar(binwidth = 1, stat = "identity") +
  theme_light() +
  scale_fill_gradient(low = "lightgray", high = "black") +
  facet_wrap(~Espécie) +
  ylab("Abundance") +
  xlab("Month") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_polar()

# Sexo

df %>% 
  group_by(Sexo, Localidade) %>%
  count() %>%
  ggplot(aes(x=Sexo, y=n, fill=Sexo)) +
  geom_col() +
  theme_light()+
  facet_wrap(~Localidade)


# NMDS

rich <- 
  df  %>%  
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

dist_bray <- 
  vegdist(run_nmds, method = "bray", binary = TRUE)

dist_bray

nmds <- 
  metaMDS(dist_bray)

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
  geom_smooth(method = "glm", alpha = .15) +
  facet_wrap(~Espécie)


# Temperatura, umidade. pressão, orvalho

df_climatics <- 
  df %>% 
  group_by(Mês) %>%
  summarise(
    total = n(),
    temp = mean(as.numeric(`Temperatura média (°C)`)),
    um = mean(as.numeric(`Umidade média (%)`)),
    or = mean(as.numeric(`Pto. Orvalho média (°C)`)),
    pre = mean(as.numeric(`Pressão média (hPa)`))
  )

df_climatics[4,3] <- 25.36


ggplot(df_climatics, aes(x = temp, y = total)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15)

ggplot(df_climatics, aes(x = um, y = total)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15)

ggplot(df_climatics, aes(x = or, y = total)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15)

ggplot(df_climatics, aes(x = pre, y = total)) +
  geom_point() +
  geom_smooth(method = "lm", alpha = .15)


# Rich x Abund

data_rich <- df %>%
  group_by(Localidade, Mês) %>%
  summarise(temp = mean(as.numeric(`Temperatura média (°C)`)),
            um = mean(as.numeric(`Umidade média (%)`)),
            or = mean(as.numeric(`Pto. Orvalho média (°C)`)),
            pre = mean(as.numeric(`Pressão média (hPa)`)),
            rich = n_distinct(`Espécie`))

data_rich[4,3] <- 25.36
data_rich

# Rich x Temp

ggplot(data_rich, aes(x = temp , y = rich)) +
  geom_point() +
  geom_smooth(method="glm", family="poisson", se=TRUE)


# Modelo 1

model1 <- glm(data = data_rich, rich ~ Mês*temp*um*or*pre, family = poisson)

summary(model1)

anova(model1)


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

# Resultado

resultados_tabanidae$AsyEst

ggiNEXT(resultados_tabanidae, type = 1, facet.var = 'order') + theme_light()
