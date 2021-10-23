df_model <- df |> group_by(Mês, Espécie) |> 
  summarise(total = n(),
            temp = mean(`Temperatura (°C)`),
            um = mean(`Umidade (%)`))

model1 <- glm(total ~ Mês*temp*um, data = df_model, family = 'poisson')

summary(model1)

anova(model1)



model2 <- glm.nb(total ~ Mês*temp*um, data = df_model)

summary(model2)

anova(model2)

exp(cbind(OR = coef(m2.1), confint(m2.1)))


ajuste = c('model1', 'model2')

aic    = c(AIC(model1), AIC(model2))

verossimilhança = c(logLik(model1),logLik(model2))

data.frame(ajuste, aic, verossimilhança)



hnp(model1, xlab = 'Percentil da N(0,1)', ylab = 'Resíduos', main = 'Gráfico Normal de Probabilidades')

hnp(model2, xlab = 'Percentil da N(0,1)', ylab = 'Resíduos', main = 'Gráfico Normal de Probabilidades')


m2.1 <- step(model2, direction = "both")

summary(m2.1)

anova(model2, m2.1)

exp(cbind(OR = coef(m2.1), confint(m2.1)))
