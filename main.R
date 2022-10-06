## dave
## 2022-01-10

## explore peds data for madaline's project


library(tidyverse)
library(data.table)

#readxl::excel_sheets('Copy of Deidentified Sample Information for University of Vermont.xlsx')
indata <- readxl::read_xlsx('Copy of Deidentified Sample Information for University of Vermont.xlsx') %>%
  data.table::as.data.table()

## fill groups
indata[,grp:=Group[1], by=cumsum(!is.na(Group))]
LOD = indata[`Final titer`!=0, min(`Final titer`)] #20
indata[,fill_titer:=ifelse(`Final titer`==0, LOD/10, `Final titer`)]


ggplot(indata, aes(x=Age, y=fill_titer, group=1)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(floor(`FAM Ct`/5) ~floor(Age/5), scales = 'free_x') +
  scale_y_log10()

p_round <- function(x, precision) {
  lower <- floor(x/precision) * precision
  
  to_label <- function(x) {
    paste0(x, "-", x-1+precision)
  }
  
  factor(to_label(lower), levels = unique(to_label(sort(lower))))
}

data.frame(age=0:17 + 0.9,age_group = as.character(p_round(0:17 + 0.9, 4)))

breaks <- 4

apply_manual_breaks <- function(x) {
  x_string = case_when(0 < x & x < 1 ~ '<1',
            1<=x & x < 6 ~ '1-5',
            x>=6 & x<12 ~'6-11',
            x>=12 & x<18 ~'12-17',
            TRUE ~ 'Other')
  return(factor(x_string, levels = c('<1', '1-5', '6-11', '12-17', 'Other')))
}
data.frame(age=0:17 + 0.9,age_group = as.character(apply_manual_breaks(0:17 + 0.9)))


ggplot(indata, aes(x=`FAM Ct`, y=fill_titer, color=factor(p_round(Age, breaks)),fill=factor(p_round(Age, breaks)), shape=factor(p_round(Age, breaks)), group=p_round(Age, breaks))) +
  geom_point(size=5, alpha=0.8) +
  geom_hline(yintercept=LOD, linetype = 2) +
  scale_size_identity() +
  scale_shape_manual(values = 21:25) +
  scale_color_manual(values= rcartocolor::carto_pal(name='Prism')[c(1:5)]) +
  scale_fill_manual(values= rcartocolor::carto_pal(name='Prism')[c(1:5)]) +
  geom_smooth(method='lm', se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  labs(color = 'Age', fill = 'Age', shape = 'Age', y='Titer (FFU/ml)', x='RNA (Ct)')
ggsave('titer_vs_ct.pdf', height = 6, width=8)

ggplot(indata, aes(x=`FAM Ct`, y=fill_titer, color=apply_manual_breaks(Age),
                   fill=apply_manual_breaks(Age), shape=apply_manual_breaks(Age),
                   group=apply_manual_breaks(Age))) +
  geom_point(size=5, alpha=0.8) +
  geom_hline(yintercept=LOD, linetype = 2) +
  scale_size_identity() +
  scale_shape_manual(values = 21:25) +
  scale_color_manual(values= rcartocolor::carto_pal(name='Prism')[c(1:5)]) +
  scale_fill_manual(values= rcartocolor::carto_pal(name='Prism')[c(1:5)]) +
  geom_smooth(method='lm', se = FALSE) +
  scale_y_log10() +
  theme_bw() +
  labs(color = 'Age', fill = 'Age', shape = 'Age', y='Titer (FFU/ml)', x='RNA (Ct)')
ggsave('titer_vs_ct_manual_breaks_no_ci.pdf', height = 6, width=8)

ggplot(indata, aes(x=`FAM Ct`, y=fill_titer, color=apply_manual_breaks(Age),
                   fill=apply_manual_breaks(Age), shape=apply_manual_breaks(Age),
                   group=apply_manual_breaks(Age))) +
  geom_point(size=5, alpha=0.8) +
  geom_hline(yintercept=LOD, linetype = 2) +
  scale_size_identity() +
  scale_shape_manual(values = 21:25) +
  scale_color_manual(values= rcartocolor::carto_pal(name='Prism')[c(1:5)]) +
  scale_fill_manual(values= rcartocolor::carto_pal(name='Prism')[c(1:5)]) +
  geom_smooth(method='lm', se = TRUE, level=0.95) +
  scale_y_log10() +
  theme_bw() +
  labs(color = 'Age', fill = 'Age', shape = 'Age', y='Titer (FFU/ml)', x='RNA (Ct)')
ggsave('titer_vs_ct_manual_breaks_95_ci.pdf', height = 6, width=8)



### stats
no_age_model <- indata %>% select(ct = `FAM Ct`, titer = fill_titer, age = Age) %>% mutate( age_grp = p_round(age, breaks)) %>%
  lm(log(titer) ~ ct, .) 

## programmatic breaks
group_model <- indata %>% select(ct = `FAM Ct`, titer = fill_titer, age = Age) %>% mutate( age_grp = p_round(age, breaks)) %>%
  lm(log(titer) ~ ct + age_grp, .) 
group_model %>% summary()
multcomp::glht(group_model, linfct = multcomp::mcp(age_grp='Tukey')) %>% summary()
anova(no_age_model, group_model)

## manual breaks
group_model <- indata %>% select(ct = `FAM Ct`, titer = fill_titer, age = Age) %>% mutate( age_grp = apply_manual_breaks(age)) %>%
  lm(log(titer) ~ ct + age_grp, .) 
group_model %>% summary()
multcomp::glht(group_model, linfct = multcomp::mcp(age_grp='Tukey')) %>% summary()
anova(no_age_model, group_model)

#continuous
continuous_model <- indata %>% select(ct = `FAM Ct`, titer = fill_titer, age = Age) %>%
  lm(log(titer) ~ ct + age, .) 
continuous_model %>% summary()
anova(no_age_model, continuous_model)

## shapes and colors
display_carto_all()
ggplot(CJ(x=0:2, y=0:1 * 3) %>% mutate(z=x+y) %>% arrange(z), aes(x=x, y=y, shape=factor(z), color=factor(z), fill = factor(z))) +
  #scale_shape_identity() +
  scale_shape_manual(values = 21:26) +
  scale_color_manual(values = rcartocolor::carto_pal(name='Prism')[c(2:7)])
  geom_point(size=5) +
  scale_size_identity()
  