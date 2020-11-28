setwd("C:/Users/vale/Google Drive/Blog/201111 - Cell count regression/Smillie et al Cell 2019/")

library(tidyverse)
library(ggrepel)
library(emmeans)

read_csv('cell_type_counts.csv') %>% mutate(batch = as.character(batch)) -> cell_type_counts

cell_type_counts

cell_type_counts %>% group_by(Health, batch, Location, Subject)

cell_type_counts %>% filter(Total > 0) -> df_a

df_a 

model0 <- glm(
  formula = cbind(Count, Other) ~ Cluster,
  family = binomial(link = 'logit'),
  data = df_a
)

emm0 <- emmeans(model0, specs = ~ Cluster)
emm0 %>%
  summary(infer = TRUE, type = 'response') %>%
  arrange(prob) -> cell_type_probs

cell_type_probs %>% head

cell_type_probs %>% pull(prob) %>% sum

df_a %>% pull(Count) %>% sum -> N
df_a %>% 
  group_by(Cluster) %>% 
  summarise(fraction = sum(Count) / N) %>% 
  arrange(fraction) %>% 
  as.data.frame -> cell_type_fractions

cell_type_fractions %>% left_join(cell_type_probs) -> cell_type_probs

(
  ggplot(aes(x = fraction, y = prob), data = cell_type_probs)
  + geom_abline(color = 'red')
  + geom_point()
  + geom_segment(aes(xend = fraction, y = asymp.LCL, yend = asymp.UCL))
  + scale_x_log10()
  + scale_y_log10()
  + theme_minimal()
  + labs(x = 'Fraction of all cells', y = 'Probability from GLM', subtitle = 'Cell type fraction vs probability')
)



df_a %>% filter(Health %in% c('Healthy', 'Inflamed')) -> df
formula = cbind(Count, Other) ~ Cluster * Health + Cluster * Location + Cluster * batch
model1 <- glm(formula = formula, family = 'binomial', data = df)

summary(model1)

emm1 <- emmeans(model1, specs = revpairwise ~ Health | Cluster)
emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> c_results

c_results %>% arrange(desc(odds.ratio))

(
  ggplot(aes(x = odds.ratio, y = -log10(p.value)), data = c_results)
  + geom_point()
  + geom_segment(aes(x = asymp.LCL, xend = asymp.UCL, yend = -log10(p.value)))
  + geom_text_repel(aes(label = Cluster), data = c_results %>% filter(p.value < 1e-50))
  + scale_x_log10()
  + theme_minimal()
  + labs(x = 'Inflamed / Healthy', title = 'Cell type proportion change')
)

(
  ggplot(aes(x = Health, y = Count / Total), data = df %>% filter(Cluster == 'Inflammatory Fibroblasts'))
  + geom_jitter(height = 0)
  + scale_y_log10()
  + theme_minimal()
  + labs(subtitle = 'Inflammatory Fibroblasts')
)

(
  ggplot(aes(x = Health, y = Count / Total), data = df %>% filter(Cluster == 'NKs'))
  + geom_jitter(height = 0)
  + scale_y_log10()
  + theme_minimal()
  + labs(subtitle = 'NKs')
)

(
  ggplot(aes(x = Health, y = Count / Total), data = df %>% filter(Cluster == 'MT-hi'))
  + geom_jitter(height = 0)
  + scale_y_log10()
  + theme_minimal()
  + labs(subtitle = 'MT-hi')
)

emm2 <- emmeans(model1, specs = ~ Cluster)
emm2 %>%
  summary(type = 'response') %>%
  select(Cluster, prob) -> mean_probs

c_results %>% left_join(mean_probs) -> m_results


(
  ggplot(aes(x = prob, y = odds.ratio, color = p.value < 0.001), data = m_results)
  + geom_point()
  + geom_text_repel(aes(label = Cluster), color = 'black', data = m_results %>% filter(abs(log(odds.ratio)) > log(2)))
  + scale_x_log10()
  + scale_y_log10()
  + theme_minimal()
  + labs(y = 'Inflamed / Healthy (odds ratio)', title = 'Cell type proportion change', x = 'Average abundance (probability)')
)

m_results %>% arrange(desc(odds.ratio))

















