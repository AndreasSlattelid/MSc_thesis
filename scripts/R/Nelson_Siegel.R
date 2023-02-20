library(tidyverse)
library(latex2exp)
library(ggthemes)

t <- c(0.25, 0.75, 1.25, 1.60, 1.85)
P <- c(0.98, 0.85, 0.83, 0.77, 0.67)
df <- data.frame(t,P)

text <- c("T1" = parse(text = TeX('$T_{1}$')), 
          "T2" = parse(text = TeX('$T_{2}$')),
          "T3" = parse(text = TeX('$T_{3}$')), 
          "T4" = parse(text = TeX('$T_{4}$')), 
          "T5" = parse(text = TeX('$T_{5}$'))
)

df %>% 
  ggplot(aes(x = t, y = P)) + 
  geom_point(shape=4, size = 4) + 
  #geom_smooth(method=lm, se=FALSE, fullrange=TRUE) + 
  geom_smooth(method='lm', 
              formula= y~poly(x,4), se=FALSE, col="orange", lwd=1) +
  scale_x_continuous(breaks = t, 
                     labels = text) + 
  labs(title = "Estimated Term Structure", 
       x = "Time", 
       y = "Zero Coupon Bond") + 
  ylim(0,1) + 
  theme_minimal()

ggsave(filename = "Estimating_term_structure.png",
       plot = last_plot(),
       dpi = 300,
       height = 8,
       width = 12, units = "cm")
