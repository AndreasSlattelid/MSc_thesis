library(tidyverse)
library(ggthemes)
library(latex2exp)


OU_CPP <- function(z0, beta, sigma, lambda, mu, t_end){
  
  #dZ(t) = -betaZ(t)dt + sigma*dW(t) + dI(t)
  #I(t): CPP with intensity lambda, and jumps J_{k}~ Exp(mu)
  
  #t_end: for how long we want to simulate. 
  
  #stepsize 
  dt <- 1/360
  time <- seq(0, t_end, by = dt)
  
  n <- length(time)
  
  #draws number of jumps in each interval of length dt N(t)~ Pois(lambda*dt)
  N <- rpois(n, lambda*dt)
  
  
  #implementation of: dZ(t) = -beta*Z(t)dt + sigma*dW(t) + dI(t)
  z <- numeric(n)
  z[1] <- z0
  
  for (i in 2:n){
    dI <- sum(rexp(N[i], mu)) - sum(rexp(N[i-1], mu))
    z[i] <- z[i-1] - beta*z[i-1]*dt + sigma*rnorm(1)*dt + dI
  }
  
  x <- 100*exp(-z)
  
  #plot(time, x, type = "l", xlab = "Time in years", ylab = "X(t)")
  df <- data.frame(time, x)
  df
}


#set.seed(123)
df <- OU_CPP(z0 = -log(20/100), beta = -0.05, sigma = 0.00, lambda = 20, mu =150, t_end = 5)

relevant_times <- c(0.0, 1.0, 2.0, 3.0, 4.0)

#ESG-criterias for different simulations:
criteria1 <- c(17.6, 16.6, 15.6, 14.55) #simulation 1
criteria2 <- c(18, 17, 16, 15)          #simulation 2
criteria3 <- 2.5                        #simulation 3




#Use NA's as plot would not render properly with 0.
df_tmp <- df %>% 
  mutate(
    timepoints = ifelse(time %in% relevant_times, 
                          time, NA)
  ) %>% 
  mutate(
    timepoints_y = ifelse(time %in% relevant_times, 
                       x, NA)
  )

text <- c("T0" = parse(text = TeX('$T_{0}$')),
          "T1" = parse(text = TeX('$T_{1}$')), 
          "T2" = parse(text = TeX('$T_{2}$')), 
          "T3" = parse(text = TeX('$T_{3}$')),
          "T4" = parse(text = TeX('$T_{4}$'))
)

ESG_OU_path <- df_tmp %>% 
  ggplot(aes(x=time, y = x)) + 
  geom_line(color="orange") + 
  #geom_step(mapping = aes(x = timepoints, y = timepoints_y), color = "black") + 
  geom_point(mapping =  aes(x = timepoints, y = timepoints_y), 
             color = "blue") + 
  scale_x_continuous(breaks = relevant_times, 
                     labels = text) + 
  labs( 
    x = "Time", 
    y = "ESG-risk score") + 
  theme_pander()

ESG_OU_path

ggsave(filename = "ESG_OU_path.png",
       plot = last_plot(),
       dpi = 300,
       height = 8,
       width = 12, units = "cm")



ESG_plt_criteria1 <- df_tmp %>% 
  ggplot(aes(x=time, y = x)) + 
  geom_line(color="orange") + 
  #geom_step(mapping = aes(x = timepoints, y = timepoints_y), color = "black") + 
  geom_point(mapping =  aes(x = timepoints, y = timepoints_y), 
             color = "blue") + 
  geom_vline(xintercept = relevant_times[2:length(relevant_times)], 
             linetype="dashed") + 
  geom_segment(aes(x=0,xend=1,y=criteria1[1],yend=criteria1[1])) + 
  geom_segment(aes(x=1,xend=2,y=criteria1[2],yend=criteria1[2])) + 
  geom_segment(aes(x=2,xend=3,y=criteria1[3],yend=criteria1[3])) + 
  geom_segment(aes(x=3,xend=4,y=criteria1[4],yend=criteria1[4])) + 
  scale_x_continuous(breaks = relevant_times, 
                     labels = text) + 
  labs( 
       x = "Time", 
       y = "ESG-risk score") + 
  theme_pander()

ESG_plt_criteria1

ESG_plt_criteria2 <- df_tmp %>% 
  ggplot(aes(x=time, y = x)) + 
  geom_line(color="orange") + 
  #geom_step(mapping = aes(x = timepoints, y = timepoints_y), color = "black") + 
  geom_point(mapping =  aes(x = timepoints, y = timepoints_y), 
             color = "blue") + 
  geom_vline(xintercept = relevant_times[2:length(relevant_times)], 
             linetype="dashed") + 
  geom_segment(aes(x=0,xend=1,y=criteria2[1],yend=criteria2[1])) + 
  geom_segment(aes(x=1,xend=2,y=criteria2[2],yend=criteria2[2])) + 
  geom_segment(aes(x=2,xend=3,y=criteria2[3],yend=criteria2[3])) + 
  geom_segment(aes(x=3,xend=4,y=criteria2[4],yend=criteria2[4])) + 
  scale_x_continuous(breaks = relevant_times, 
                     labels = text) + 
  labs( 
    x = "Time", 
    y = "ESG-risk score") + 
  theme_pander()

ESG_plt_criteria2


ggsave(filename = "ESG_plt_criteria2.png",
       plot = last_plot(),
       dpi = 300,
       height = 8,
       width = 12, units = "cm")

ESG_plt_criteria3 <- df_tmp %>% 
  ggplot(aes(x=time, y = x)) + 
  geom_line(color="orange") + 
  #geom_step(mapping = aes(x = timepoints, y = timepoints_y), color = "black") + 
  geom_point(mapping =  aes(x = timepoints, y = timepoints_y), 
             color = "blue") + 
  geom_vline(xintercept = relevant_times[2:length(relevant_times)], 
             linetype="dashed") + 
  geom_segment(aes(x=0,xend=4,y = criteria3, yend = criteria3)) + 
  scale_x_continuous(breaks = relevant_times, 
                     labels = text) + 
  labs( 
    x = "Time", 
    y = "ESG-risk score") + 
  theme_pander()

ESG_plt_criteria3

ggsave(filename = "ESG_plt_criteria3.png",
       plot = last_plot(),
       dpi = 300,
       height = 8,
       width = 12, units = "cm")



