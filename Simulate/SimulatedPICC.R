##### SETUP ####
# function to load packages or install if not already installed
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("dplyr", "ggplot2", "grid", 'gridExtra')
ipak(packages)
rm(packages, ipak)

#### SIMULATE STIMULI LISTS #### 
stimRandomizer <- function(colorfreq, numtrials){
  reds <- seq(from = 0, to = 100)
  greens <- rep(0, 100)
  blues <- seq(from = 255, to = 155)
  
  colors = list()
  
  for(i in seq(100)){
    colors[[length(colors)+1]] <- list(reds[i], greens[i], blues[i])
  }
  colorsignals <- colors[1:50] # blue spectrum
  colorsnoise <- colors[51: 100]

  colorsignals <- sample(size = as.integer(numtrials*colorfreq), x = colorsignals)
  colorsnoise <- sample(size = as.integer(numtrials*(1-colorfreq)), x = colorsnoise)
  colorstim <- sample(c(colorsignals, colorsnoise), size = numtrials)
  
  return(colorstim)
}

stablefreq <- rep(0.5, 16) # stable prevalence condition 
decreasefreq <- c(.5,.5,.5,.5,.4,.28,.14,.06,.06,.06,.06,.06,.06,.06,.06,.06) # changing prevalence condition 

SimSubjs <- 80 # 40 simulated subjects per condition 

SimData <- data.frame(id = numeric(), condition = numeric(), trial = numeric(), block = numeric(), 
                      red = numeric(), green = numeric(), blue = numeric(), 
                      colorfreq = numeric())
id = 1
for (sub in seq(SimSubjs)){
  print(paste0('Simulating data for subject: ', id))
  if(sub < 30){
    block = 1
    trial = 1
    for(b in stablefreq){ # subjects 1-30 = stable 
        thisblock <- stimRandomizer(b, 50)
          for(t in thisblock){
          SimData[nrow(SimData)+1, ] <- c(id, 0, trial, block, t[1], t[2], t[3], b)
          trial=trial+1
          }
        block = block+1
        } 
    }
    else{
      block = 1
      trial = 1
      for(b in decreasefreq){ # subjects 31-60 = decreasing
        thisblock <- stimRandomizer(b, 50)
        for(t in thisblock){
          SimData[nrow(SimData)+1, ] <- c(id, 1, trial, block, t[1], t[2], t[3], b)
          trial=trial+1
        }
        block = block+1
    }
    }
  id = id+1
}

rm(b, block, decreasefreq, id, SimSubjs, stablefreq, sub, t, thisblock, trial)

# write.csv(SimData, 'SimulatedData(noresp).csv', row.names = F)

#### ORIGINAL MODEL: SIMULATE RESPONSES ####
# Predict participants potential responses based on parameters in the model 
# H1 - Older adults are less sensitive to PICC due to greater effect of past response; more perseveration (higher Bc)
# H2 - Older adults are more sensitive to PICC due to greater effect of past stimulus; more outsourcing of control (higher BF)

SimChoices = function(f, beta, lambdaF, lambdaC){ # calculates probability of choice in model-based analysis below
  Fbar = c(0) # list of exponentially weighted past stimuli 
  Cbar = c(0) # list of exponentially weighted past responses 
  CP = c() # choice probabilities 
  resp = rep(1, length(f)) # dummy responses; used to indicate the probability of choosing 1 
  for(t in seq(length(f))){
    
    # predict choice probability from regression parameters
    dV = beta[1] + beta[2]*f[t] + beta[3]*Fbar[t] + beta[4]*Cbar[t] 
    P = 1/ (1 + exp(dV)) 
    CP[t] <- 1-P
    
    Fbar[t+1] <- lambdaF*Fbar[t] + f[t] # update list weighing it by decay parameter lambda
    Cbar[t+1] <- lambdaC*Cbar[t] + resp[t]
  }

  return(CP) # represents the probability of presssing 1
}

H0 <- c(-4, 15, -1.15, 0.5, 0.8, 0.5) # similar paramters to Wilson (2018) 
H1 <- c(-4, 15, -0.75, 2, 0.8, 0.5) # Bc > BF
H2 <- c(-4, 15, -1.75, 0.35, 0.8, 0.8) # Bc < BF  

SimData$f <- SimData$f <- 1 - 2*((SimData$red - min(SimData$red))/(max(SimData$red) - min(SimData$red))) # objective value of blueness; higher means more blue


for(id in unique(SimData$id)) {
  print(paste0('Simulating responses for subject ',id))
  
  
  b0 <- c(rnorm(1, mean=H0[1], sd=0.1), rnorm(1, mean=H0[2], sd=0.1), rnorm(1, mean=H0[3], sd=0.2),
          rnorm(1, mean=H0[4], sd=0.2), rnorm(1, mean=H0[5], sd=0.1), rnorm(1, mean=H0[6], sd=0.1))
  
  b1 <- c(rnorm(1, mean=H1[1], sd=0.1), rnorm(1, mean=H1[2], sd=0.1), rnorm(1, mean=H1[3], sd=0.2),
          rnorm(1, mean=H1[4], sd=0.2), rnorm(1, mean=H1[5], sd=0.1), rnorm(1, mean=H1[6], sd=0.1))
  
  b2 <- c(rnorm(1, mean=H2[1], sd=0.1), rnorm(1, mean=H2[2], sd=0.1), rnorm(1, mean=H2[3], sd=0.2),
          rnorm(1, mean=H2[4], sd=0.2), rnorm(1, mean=H2[5], sd=0.1), rnorm(1, mean=H2[6], sd=0.1))
  
  d.subject <- SimData[SimData$id == id, ]
  SimData[SimData$id == id,'simresp0'] <- ifelse(SimChoices(d.subject$f, b0[1:4], b0[5], b0[6]) > 0.5, 1, -1)
  SimData[SimData$id == id,'simrespH1']<- ifelse(SimChoices(d.subject$f, b1[1:4], b1[5], b1[6]) > 0.5, 1, -1)
  SimData[SimData$id == id,'simrespH2']<- ifelse(SimChoices(d.subject$f, b2[1:4], b2[5], b2[6]) > 0.5, 1, -1)
}
rm(id, d.subject)

#### ORIGINAL MODEL: VISUALISE SIMULATED RESPONSES ####
SimData$condition <- ifelse(SimData$condition == 0, 'Stable Prevalence Condition', 
                            'Decreasing Prevalence Condition')
SimData$simresp0 <- ifelse(SimData$simresp0 == -1, 0, 1)
SimData$simrespH1 <- ifelse(SimData$simrespH1 == -1, 0, 1)
SimData$simrespH2 <- ifelse(SimData$simrespH2 == -1, 0, 1)

SimData$trialintensity <- SimData$blue - 154
SimData$timebin4 <- as.factor(apply(as.matrix(SimData$trial), 2,
                                    cut, c(seq(0,800,200)), labels=FALSE)) # create timebins

desc.sim <- SimData %>% group_by(condition, trialintensity, timebin4) %>% 
  summarise(pBlue0 = length(simresp0[simresp0 == 1])/length(simresp0), 
            pBlueH1 = length(simrespH1[simrespH1 == 1])/length(simrespH1), 
            pBlueH2 = length(simrespH2[simrespH2 == 1])/length(simrespH2)) %>%
  rename(colour = trialintensity, trialbins = timebin4)
desc0.sim<- desc.sim[desc.sim$condition == 'Stable Prevalence Condition' & (desc.sim$trialbins==1 | desc.sim$trialbins ==4),]
desc1.sim <- desc.sim[desc.sim$condition == 'Decreasing Prevalence Condition' & (desc.sim$trialbins==1 | desc.sim$trialbins ==4),]

# PLOTS #
stable0 <- ggplot(desc0.sim, aes(x = colour, y = pBlue0, color = trialbins)) + geom_point() + ggtitle("H0: Similar to Wilson's Parameters\nStable Prevalence") + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  labs(y = '', x = '') + 
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 96), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.position = 'none', 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))

stableH1 <- ggplot(desc0.sim, aes(x = colour, y = pBlueH1, color = trialbins)) + geom_point() + ggtitle(expression(atop(paste("H1: Higher ", beta[c]), "Stable Prevalance"))) + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  labs(y = '', x = '') + 
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 96), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.position = 'none', 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))

stableH2 <- ggplot(desc0.sim, aes(x = colour, y = pBlueH2, color = trialbins)) + geom_point() +ggtitle(expression(atop(paste("H2: Higher ", beta[F]), "Stable Prevalance"))) + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  labs(y = '', x = '') + 
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 96), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.position = 'none', 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))

decrease0 <- ggplot(desc1.sim, aes(x = colour, y = pBlue0, color = trialbins)) + geom_point() + ggtitle("H0: Similar to Wilson's Parameters\nDecreasing Prevalence") + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  ylab('') + xlab('') +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 96), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.position = 'none', 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))

decreaseH1 <- ggplot(desc1.sim, aes(x = colour, y = pBlueH1, color = trialbins)) + geom_point() + ggtitle(expression(atop(paste("H1: Higher ", beta[c]), "Decreasing Prevalance"))) + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  ylab('') + xlab('') +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 96), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.position = 'none', 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))

decreaseH2 <- ggplot(desc1.sim, aes(x = colour, y = pBlueH2, color = trialbins)) + geom_point() + ggtitle(expression(atop(paste("H2: Higher ", beta[F]), "Decreasing Prevalance"))) + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  ylab('') + xlab('') +
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 96), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.justification = c(1, 0), legend.position = c(1, 0), 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))

grid.arrange(stableH1,stable0, stableH2, 
             decreaseH1, decrease0, decreaseH2, 
             nrow = 2, ncol = 3,
             #top = textGrob("X",x=0,hjust=0, gp=gpar(fontsize=30)),  
             left = textGrob("% Dots Judged as Blue",gp=gpar(fontsize=20), rot = 90))

rm(stable0, stableH1, stableH2, decrease0, decreaseH1, decreaseH2, b0, b1, b2, desc.sim, desc0.sim, desc1.sim)

##### MODIFIED MODEL: SIMULATE RESPONSES ####
# Simulate data based on a proposed modified model that incorporates RT 
# Uses the same parameters as observed in the young adults, but visualises the impact of high vs. low RT 
# High RTs = normal distribution with mean and sd RT from older adults, low RTs = same from young adults

SimChoicesNEW = function(f, beta, lambdaF, lambdaC, RT){ # calculates probability of choice in model-based analysis, by taking into account RT
  Fbar = c(0) # list of exponentially weighted past stimuli 
  Cbar = c(0) # list of exponentially weighted past responses 
  CP = c() # choice probabilities 
  resp = rep(1, length(f)) # dummy responses; used to indicate the probability of choosing 1 
  for(t in seq(length(f))){
    
    # predict choice probability from regression parameters
    dV = beta[1] + beta[2]*f[t] + beta[3]*Fbar[t] + beta[4]*Cbar[t] 
    P = 1/ (1 + exp(dV)) 
    CP[t] <- 1-P
    
    Fbar[t+1] <- (lambdaF^RT[t])*Fbar[t] + f[t] # update list weighing it by decay parameter lambda
    Cbar[t+1] <- (lambdaC^RT[t])*Cbar[t] + resp[t]
  }
  
  return(CP) # represents the probability of presssing 1
}

Params <- c(2.20, 11.7, -0.54, 0.33, 0.77, 0.41) # OA
# Params <- c(2.49, 11.70, -1.16, 0.57, 0.70, 0.53) # OA
# Params <- c(2.44, 12.13, -1.09, 0.66, 0.63, 0.44) # same parameters as YA; RTs to be simulated below 

# Low RT 
LowRTSim <- c()
for(id in unique(SimData$id)){
  d.subject <- SimData[SimData$id == id, ]
  
  print(paste0('simulating low RT responses for subject: ', id))
  
  RT <- rnorm(n=length(d.subject$trial), mean=0.287, sd = 0.458) # values taken from YA data 
  
  LowRTSim <- c(LowRTSim, ifelse(SimChoicesNEW(d.subject$f, Params[1:4], Params[5], Params[6], RT) > 0.5, 1, -1))
}
SimData$LowRTSim <- LowRTSim

# High RT 
HighRTSim <- c()
for(id in unique(SimData$id)){
  d.subject <- SimData[SimData$id == id, ]
  
  print(paste0('simulating high RT responses for subject: ', id))
  
  RT <- rnorm(n=length(d.subject$trial), mean=0.5703, sd = 1.140) # values taken from OA data 
  
  HighRTSim <- c(HighRTSim, ifelse(SimChoicesNEW(d.subject$f, Params[1:4], Params[5], Params[6], RT) > 0.5, 1, -1))
}
SimData$HighRTSim <- HighRTSim

#### MODIFED MODEL: VISUALISE SIMULATED RESPONSES ####
SimData$LowRTSim <- ifelse(SimData$LowRTSim == -1, 0, 1)
SimData$HighRTSim <- ifelse(SimData$HighRTSim == -1, 0, 1)

desc.sim2<- SimData %>% group_by(condition, trialintensity, timebin4) %>% 
  summarise(pBlueLow = length(LowRTSim[LowRTSim == 1])/length(LowRTSim), 
            pBlueHigh = length(HighRTSim[HighRTSim == 1])/length(HighRTSim)) %>%
  rename(colour = trialintensity, trialbins = timebin4)
desc0.sim2<- desc.sim2[desc.sim2$condition == 'Stable Prevalence Condition' & (desc.sim2$trialbins==1 | desc.sim2$trialbins ==4),]
desc1.sim2 <- desc.sim2[desc.sim2$condition == 'Decreasing Prevalence Condition' & (desc.sim2$trialbins==1 | desc.sim2$trialbins ==4),]

# PLOTS # 
LowRTStable <- ggplot(desc0.sim2, aes(x = colour, y = pBlueLow, color = trialbins)) + geom_point() + ggtitle("Simulated Responses when RT is Fast\nStable Prevalence") + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  labs(y = '', x = '') + 
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 100), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.position = 'none', 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))

LowRTDecrease <- ggplot(desc1.sim2, aes(x = colour, y = pBlueLow, color = trialbins)) + geom_point() + ggtitle("Simulated Responses when RT is Fast\nDecreasing Prevalence") + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  labs(y = '', x = '') + 
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 100), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.position = 'none', 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))

HighRTStable <- ggplot(desc0.sim2, aes(x = colour, y = pBlueHigh, color = trialbins)) + geom_point() + ggtitle("Simulated Responses when RT is Slow\nStable Prevalence") + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  labs(y = '', x = '') + 
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 100), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.position = 'none', 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))

HighRTDecrease <- ggplot(desc1.sim2, aes(x = colour, y = pBlueHigh, color = trialbins)) + geom_point() + ggtitle("Simulated Responses when RT is Slow\nDecreasing Prevalence") + 
  geom_line(stat="smooth",method = "glm", method.args = list(
    family="binomial"), se=FALSE,size=1.8,alpha=.7) + 
  labs(y = '', x = '') + 
  scale_colour_manual(labels=c("Initial 200 trials", "Final 200 trials"), name=NULL,values=c("#0066CC", '#990000')) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(0, 100), 
                     labels = c('Very Purple', 'Very Blue')) + 
  theme_bw() + theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), 
                     legend.position = 'none', 
                     axis.text.y = element_text(size = 10),
                     axis.text.x = element_text(size = 10, colour = 'black'), 
                     axis.title.y = element_text(size = 12)) + 
  ylim(c(0.00, 1.00))


grid.arrange(LowRTStable, LowRTDecrease, 
             HighRTStable, HighRTDecrease,
             nrow = 2, ncol = 2,
             #top = textGrob("X",x=0,hjust=0, gp=gpar(fontsize=30)),  
             left = textGrob("% Dots Judged as Blue",gp=gpar(fontsize=20), rot = 90))
