setwd("~/Desktop/Ongoing/PICCAgeing/Analysis&Data")
rm(list=ls())
load('PICCOA_2021.RData')
source('Models/seq/loglik.R')
library(dplyr)
library(ggplot2); library(ggpubr)
library(sjPlot)

# 0. Randomly select data to use as test stim
set.seed(1456)
condition = 'Decreasing'
randid = sample(d.dots$id[d.dots$condition==condition], 1)
testset <- d.dots[d.dots$id == randid, ]
red = as.numeric(gsub(",", "", substr(testset$trialstim, 2, 3)))
blue = testset$trialintensity + 154 
f =  f = 1 - 2*((red - min(red))/(max(red) - min(red))) 
#c = ifelse(testset$response==1, 1, -1)
pc0 = 1/(1 + exp(-f[1]))
c = rep(1, length(f)) # dummy responses; used to indicate the probability of choosing 1 

# 1.1. Specify median parameters (except params of interest, in this case lambdas)
betas = as.numeric(dots.fit[dots.fit$id == randid, 4:7])

# 1.2. Select range of params of interest to test
n = 10 # number of params to test (can get expensive quickly)
lf = seq(0, 1, length.out = n)
lc = seq(0, 1, length.out = n)
lambdas = expand.grid(lf, lc) # all combinations of lf and lc

# 2. Loop over lambdas
plots = list()
cpdf <- data.frame() # data frame for choice probabilities 
pcount = 1
for(l in 1:nrow(lambdas)) { 
  cat('lambda combo:', l, '/', nrow(lambdas), '\n')
  params = c(betas, lambdas[l,1], lambdas[l,2])
 # params = as.numeric(dots.fit[dots.fit$id == randid, 4:9]) # to test
  cp = loglik(f, c, params, return = 'cp')
 
  simresp <- c()
  for(t in 1:length(cp)) simresp[t] = sample(c(1, 0), 1, prob = c(cp[t], 1-cp[t]))
  
  # 2.1. Plot simulated responses 
  testset$simresp <- simresp
  thispicc <- 
    testset %>% 
    ungroup() %>% 
    mutate(timebin4 = factor(timebin4)) %>% 
    group_by(timebin4, trialintensity) %>% 
    summarise(psimBlue = mean(simresp)) %>% 
    filter(timebin4 == 1 | timebin4 == 4) %>% 
    ggplot(aes(x = trialintensity, y = psimBlue, color = timebin4)) + 
    geom_point() + 
    geom_line(stat="smooth",method = "glm", method.args = list(
      family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
    labs(x= '', y = '% Dots Judged as Blue', 
         title = paste0('lf = ', params[5],'\nlc = ', params[6]), 
         caption = paste0('data from subject: ', as.character(randid))) +
    scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), 
                        name=NULL,values=c("#0066CC", '#990000')) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                       limits = c(0, 1)) +
    scale_x_continuous(breaks = c(5, 96), 
                       labels = c('Very\nPurple', 'Very\nBlue')) + 
    theme_bw() + 
    theme(plot.caption=element_text(hjust = -0.2))
  
    plots[[pcount]] = thispicc
    pcount = pcount+1
    
    # 2.2. Save choice probabilities 
    thiscpdf <- data.frame(id = randid, cond=testset$condition, 
                       f = f, trial=1:length(f), cp = cp, 
                       lf = lambdas[l, 1], lc = lambdas[l,2])
    
    # 2.3. Do subject-level regression
    reg <- glm(simresp ~ trial0*colour0, family='binomial', data=testset)
    thiscpdf$b_trial0 = coef(reg)['trial0']
    thiscpdf$b_colour0 = coef(reg)['colour0']
    thiscpdf$b_trial0colour0 = coef(reg)['trial0:colour0']
    
    cpdf <- rbind(cpdf, thiscpdf)
  }

# 3. Plot choice probabilities by lambdas
fbins <- quantile(f, seq(0, 1, length.out = 10))
names(fbins) = paste0(round(seq(0, 1, length.out = 10), 2)*100, '%')
heatmap <- 
  cpdf %>% 
  mutate(fnames = sapply(f, function(x) names(fbins[fbins==x])), 
         fnames = factor(paste(fnames, ' blue'),
                         levels = paste(names(fbins), ' blue')), 
         timebin4 = ntile(trial, 4)) %>% 
  filter(f %in% fbins, timebin4 %in% c(1, 4)) %>%
  mutate(timebin4 = ifelse(timebin4==1, 'First 200 Trials', 'Last 200 Trials')) %>% 
  ggplot(aes(x = lc , y = lf, fill=cp)) +
  geom_tile() + 
  facet_grid(timebin4~fnames) + 
  scale_fill_gradient(low='purple', high = 'blue') + 
  labs(x =  expression('\u03BB'[c]), y =  expression('\u03BB'[F]), 
       fill = 'p(Choose Blue)', title=paste(condition, 'condition.')) + 
  theme_bw()

# 4. Plot regression coefficients per lambda value
reg2 <- lm(b_trial0colour0 ~ lc+lf, data=cpdf[cpdf$b_trial0colour0 > -1000,])
reg2plot <- plot_model(reg2, type='eff') 
for(p in 1:length(reg2plot)) { 
  xl = ifelse(p==1, expression('\u03BB'[c]), expression('\u03BB'[F]))
  reg2plot[[p]] <-
    reg2plot[[p]] + 
    labs(x=xl,  y = expression('\u03B2'['trial \u00D7 colour']), 
         title = 'Predicted values of\n\u03B2 trial \u00D7 colour')
    theme_bw() 
} 
reg2plot[[3]] <- 
  cpdf %>% 
  select(lc, lf, b_trial0colour0) %>% 
  filter(b_trial0colour0 > -75 & b_trial0colour0 < 75) %>% 
  melt(id.vars = 'b_trial0colour0') %>% 
  mutate(variable = ifelse(variable == 'lc', 
                           '\u03BBc', 
                           '\u03BBF')) %>% 
  ggplot(aes(x = value, y = b_trial0colour0, colour=variable)) + 
  stat_summary(fun = mean, geom='point', alpha=0.2) + 
  stat_summary(fun.data = mean_se, geom='errorbar', alpha=0.2, width = 0.1) +
  #geom_point(alpha =0.05) + 
  geom_smooth(method='lm') + 
  scale_color_brewer(palette = 'Dark2') + 
  labs(x = 'Parameter Value', y = expression('\u03B2'['trial \u00D7 colour']), 
       colour='', title='Actual values of \n\u03B2 trial \u00D7 colour', 
       caption = 'parameter space constrained OR \u2208 [-75, 75]') + 
  theme_bw() + 
  theme(plot.caption=element_text(hjust = 0))
  
# 5. plots
pdf(paste0('Simulate/lambdasims_',condition, '_',randid, '.pdf'), width = 6.66, height = 4)
for(p in 1:length(plots)) { 
  cat(p, '/', length(plots), '.\n')
  suppressWarnings(suppressMessages(print(plots[[p]])))
}
dev.off()

ggarrange(plotlist=reg2plot[3], 
  ggarrange(plotlist = reg2plot[1:2]), 
  nrow=2, legend = 'bottom') %>% 
  ggsave(filename = paste0('Simulate/predicted_actual_', condition, '_', randid, '.pdf'), 
         device = cairo_pdf, width = 6, height = 8)

heatmap %>% ggsave(filename=paste0('Simulate/heatmap_', condition, '_', randid, '.pdf'), 
                   device = cairo_pdf, width = 6, height = 4)
