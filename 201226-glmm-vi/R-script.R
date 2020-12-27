library(tidyverse)
library(lme4)
library(sjPlot)

read_csv('zeb2_table.csv') -> data

glmer(ZEB2 ~ 0 + cell_group + offset(log(total_count)) + (1 | patient_assignment), data = data, family = poisson) -> m1

summary(m1)

data.frame(
  cell_group = fixef(m1) %>% names(),
  mean = fixef(m1),
  std = coef(summary(m1))[, 2]
) -> fixed_effects

plot_model(m1, transform = NULL)

ranef(m1) -> randoms
randoms$patient_assignment -> rand.interc
data.frame(
  patient_assignment = rownames(rand.interc),
  mean = randoms$patient_assignment[, 1]
) -> random_effects

randoms$patient_assignment[, 1]

(
  ggplot(random_effects, aes(x = mean, y = patient_assignment))
  + geom_point()
  + geom_errorbarh(aes(xmin = mean - re_std, xmax = mean + re_std))
)

VarCorr(m1)


fixed_effects %>% write_csv('lme4/fixed_effects.csv')
random_effects %>% write_csv('lme4/random_effects.csv')
VarCorr(m1) %>% as.data.frame %>% write_csv('lme4/varcorr.csv')


