library(tidyverse)
library(ggthemes)
library(latex2exp)



kappa <- 0.070
kappa_ESG <- c(kappa, 0.065, 0.060, 0.055, 0.060)
time <- c(0, 1, 2, 3, 4)


df <- data.frame(time, kappa_ESG)

text <- c("T0" = parse(text = TeX('$T_{0}$')),
          "T1" = parse(text = TeX('$T_{1}$')), 
          "T2" = parse(text = TeX('$T_{2}$')), 
          "T3" = parse(text = TeX('$T_{3}$')),
          "T4" = parse(text = TeX('$T_{4}$')), 
          "T5" = parse(text = TeX('$T_{5}$'))
)

df %>% 
  ggplot(aes(x = time, y = kappa_ESG)) +
  scale_x_continuous(breaks = c(0,1,2,3,4,5), 
                     labels = text) +
  geom_segment(aes(x=0,xend=5,y= kappa,yend=kappa), color = "orange") +
  geom_segment(aes(x=1,xend=2,y=kappa_ESG[2],yend=kappa_ESG[2]),color = "orange" ) + 
  geom_segment(aes(x=2,xend=3,y=kappa_ESG[3],yend=kappa_ESG[3]), color = "orange") + 
  geom_segment(aes(x=3,xend=4,y=kappa_ESG[4],yend=kappa_ESG[4]), color = "orange") + 
  geom_segment(aes(x=4,xend=5,y=kappa_ESG[5],yend=kappa_ESG[5]), color = "orange") + 
  labs(
    x = "Time", 
    y = "fixed rate process") + 
  ggtitle(TeX("Trajectory of $K^{ESG} = (K_{i}^{ESG}(\\omega))_{i=1}^{4}$")) +
  theme_pander() + 
  ylim(0.050, 0.075) + 
  annotate("text", x = 2.5, y = 0.071, label = TeX("$\\kappa$")) + 
  annotate("text", x = 1.6, y = 0.067, label = TeX("$K_{1}^{ESG}$")) + 
  annotate("text", x = 2.6, y = 0.062, label = TeX("$K_{2}^{ESG}$")) + 
  annotate("text", x = 3.6, y = 0.057, label = TeX("$K_{3}^{ESG}$")) + 
  annotate("text", x = 4.6, y = 0.062, label = TeX("$K_{4}^{ESG}$")) + 
  geom_point(aes(x=1,y=0.065),colour="black") + 
  geom_point(aes(x=2,y=0.065),colour="black", shape = 1)  + 
  geom_point(aes(x=2,y=0.060),colour="black") + 
  geom_point(aes(x=3,y=0.060),colour="black", shape = 1) + 
  geom_point(aes(x=3,y=0.055),colour="black") + 
  geom_point(aes(x=4,y=0.055),colour="black", shape = 1) + 
  geom_point(aes(x=4,y=0.060),colour="black") + 
  geom_point(aes(x=5,y=0.060),colour="black", shape = 1) 


#ggsave(filename = "SBM_ESG_path.png",
#       plot = last_plot(),
#       dpi = 300,
#       height = 8,
#       width = 12, units = "cm")

