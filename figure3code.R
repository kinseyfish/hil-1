#partial deletions wormsizer
setwd("/Users/kinseyfisher/Documents/hil-1/Figure 3 material")
wormsizer <- read.csv("hil-1_wormsizer_backcrossed.csv", header = T)
wormsizer<-subset(wormsizer,pass=="TRUE")
wormsizer$Days_Starved <- as.character(wormsizer$Days_Starved)

wormsizer_plot<-ggplot(wormsizer,aes(x=Days_Starved,y=length))+
  #geom_boxplot(outlier.shape = NA)+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_dotplot(inherit.aes = TRUE,binaxis = "y",alpha=1,binwidth=2,stackdir = "center", dotsize = 1)+
  facet_grid(.~Strain) +
  theme_bw()+
  theme(legend.position="none")+
  stat_summary(
    fun = mean,
    geom = 'line',
    aes(group= Strain),
    position = position_dodge(width = 0.9),
    color = "black") +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=1, color="black")
wormsizer_plot

d1d8<-subset(wormsizer,Strain=="WT" | Strain=="hil-1(gk229)")
interaction<-lme(length~Strain*Days_Starved,random=~1|Replicate, data=d1d8) #the * notation is what we use for the interaction
summary(interaction)


#Full deletion wormsizer
wormsizer <- read.csv("hil1_fulldeletion_wormsizer_bx.csv", header = T)
wormsizer<-subset(wormsizer,pass=="TRUE")
wormsizer$Strain <- factor(wormsizer$Strain,levels = c("N2", "hil-1 (syb9165)"))
wormsizer$Day <- as.character(wormsizer$Day)
wormsizer_plot<-ggplot(wormsizer,aes(x=Day,y=length))+
  #geom_boxplot(outlier.shape = NA)+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_dotplot(inherit.aes = TRUE,binaxis = "y",alpha=1,binwidth=2,stackdir = "center", dotsize = 1)+
  facet_grid(.~Strain) +
  theme_bw()+
  theme(legend.position="none")+
  stat_summary(
    fun = mean,
    geom = 'line',
    aes(group= Strain),
    position = position_dodge(width = 0.9),
    color = "black") +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=1, color="black")
wormsizer_plot
d1d8<-subset(wormsizer,Strain=="N2" | Strain=="hil-1 (syb9165)")
interaction<-lme(length~Strain*Day,random=~1|Rep, data=d1d8) #the * notation is what we use for the interaction
summary(interaction)

#AID wormsizer
wormsizer <- read.csv("hil-1_AID_Wormsizer_KNAA.csv" , header = T)
wormsizer<-subset(wormsizer,pass=="TRUE")
wormsizer$Condition <- factor(wormsizer$Condition,levels = c("RY16", "RY16; TIR1", "RY16; TIR1 + aux", "RY17","RY17; TIR1", "RY17; TIR1 + aux"))
wormsizer$Day <- as.character(wormsizer$Day)
wormsizer$Replicate <- as.character(wormsizer$Replicate)
wormsizer_plot<-ggplot(wormsizer,aes(x=Day,y=length))+
  #geom_boxplot(outlier.shape = NA)+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_dotplot(inherit.aes = TRUE,binaxis = "y",alpha=1,binwidth=2,stackdir = "center", dotsize = 1)+
  facet_grid(.~Condition) +
  theme_bw()+
  theme(legend.position="none")+
  stat_summary(
    fun = mean,
    geom = 'line',
    aes(group= Condition),
    position = position_dodge(width = 0.9),
    color = "black") +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=1, color="black")
wormsizer_plot

d1d8<-subset(wormsizer,Condition=="RY16; TIR1" | Condition=="RY16; TIR1 + aux")
interaction<-lme(length~Condition*Day,random=~1|Replicate, data=d1d8) #the * notation is what we use for the interaction
summary(interaction)

d1d8<-subset(wormsizer,Condition=="RY17; TIR1" | Condition=="RY17; TIR1 + aux")
interaction<-lme(length~Condition*Day,random=~1|Replicate, data=d1d8) #the * notation is what we use for the interaction
summary(interaction)


#Survival partial deletions:
ss <- read.csv("hil-1_partialdeletions_survival_raw.csv")
ss$replicate <- as.character(ss$replicate)
ss$geno <- factor(ss$geno,levels = c("WT", "hil-1(gk229)", "hil-1(tm1442)"))

survival_plot<-ggplot()+
  geom_point(data=ss, aes(x = day, y = prop_alive, color = geno))+
  geom_smooth(data=ss, aes(x = day, y = prop_alive, color = geno), method = "glm", method.args = list(family = "quasibinomial"), level = 0)+
  scale_x_continuous(breaks = seq(1, 30, by = 2))+
  theme_classic(base_size = 15)+
  theme(aspect.ratio = 1)+
  ggtitle("survival rate ~ time")+labs(x="Days of L1 arrest",y="Proportion alive") 
survival_plot


result <- unique(subset(ss, select = "geno"))
for (i in 1:length(result$geno)) {
  model <- glm(data = ss, subset = (ss$geno == result$geno[i]),
               formula = prop_alive ~ day,
               family = "quasibinomial")
  half_life = - model$coefficients[1]/model$coefficients[2]
  goodness_of_fit <- 1 - model$deviance/model$null.deviance
  result$half_life[i] <- half_life
  result$goodness_of_fit[i] <- goodness_of_fit
}  
result

#round half_life to 2 digits
result$half_life <- round(result$half_life, 2)
result2 <- unique(subset(ss, select = c("geno", "replicate")))
result2$condition <- paste(result2$geno, result2$replicate, sep = "_")

#then generate model
for (i in 1:length(result2$condition)) {
  model <- glm(data = ss, subset = (ss$geno == result2$geno[i]
                                    & ss$replicate == result2$replicate[i]),
               formula = prop_alive ~ day,
               family = "quasibinomial")
  half_life = - model$coefficients[1]/model$coefficients[2]
  goodness_of_fit <- 1 - model$deviance/model$null.deviance
  result2$half_life[i] <- half_life
  result2$goodness_of_fit[i] <- goodness_of_fit
}  
result2

result2$half_life <- round(result2$half_life, 2)
result2$replicate <- as.factor(result2$replicate)
bartlett.test(half_life~geno, data = result2)
pairwise.t.test(result2$half_life,result2$geno, pool.sd = TRUE, p.adjust = "none")



#Survival full deletion:
ss <- read.csv("hil-1_partialdeletions_survival_raw.csv")
ss$replicate <- as.character(ss$replicate)
ss$strain <- factor(ss$strain,levels = c("N2", "hil-1 (syb9165)"))

survival_plot<-ggplot()+
  geom_point(data=ss, aes(x = day, y = prop_alive, color = strain))+
  geom_smooth(data=ss, aes(x = day, y = prop_alive, color = strain), method = "glm", method.args = list(family = "quasibinomial"), level = 0)+
  scale_x_continuous(breaks = seq(1, 30, by = 2))+
  theme_classic(base_size = 15)+
  theme(aspect.ratio = 1)+
  ggtitle("survival rate ~ time")+labs(x="Days of L1 arrest",y="Proportion alive") 
survival_plot


result <- unique(subset(ss, select = "strain"))
for (i in 1:length(result$strain)) {
  model <- glm(data = ss, subset = (ss$strain == result$strain[i]),
               formula = prop_alive ~ day,
               family = "quasibinomial")
  half_life = - model$coefficients[1]/model$coefficients[2]
  goodness_of_fit <- 1 - model$deviance/model$null.deviance
  result$half_life[i] <- half_life
  result$goodness_of_fit[i] <- goodness_of_fit
}  
result

#round half_life to 2 digits
result$half_life <- round(result$half_life, 2)
result2 <- unique(subset(ss, select = c("strain", "replicate")))
result2$condition <- paste(result2$strain, result2$replicate, sep = "_")

#then generate model
for (i in 1:length(result2$condition)) {
  model <- glm(data = ss, subset = (ss$strain == result2$strain[i]
                                         & ss$replicate == result2$replicate[i]),
               formula = prop_alive ~ day,
               family = "quasibinomial")
  half_life = - model$coefficients[1]/model$coefficients[2]
  goodness_of_fit <- 1 - model$deviance/model$null.deviance
  result2$half_life[i] <- half_life
  result2$goodness_of_fit[i] <- goodness_of_fit
}  
result2

result2$half_life <- round(result2$half_life, 2)
result2$replicate <- as.factor(result2$replicate)
bartlett.test(half_life~strain, data = result2)
pairwise.t.test(result2$half_life,result2$strain, pool.sd = TRUE, p.adjust = "none")



write_xlsx(ss, "hil-1_partialdeletions_survival_raw.xlsx")


#GFP degradation confirmation
Auxin <- read.csv("Auxin_Comparisons.csv", header = T)

ggplot()+
  geom_dotplot(Auxin, mapping = aes(x=Condition,y=Average.Intensity.Average),inherit.aes = TRUE,binaxis = "y",alpha=1,binwidth=5,stackdir = "center", dotsize = 3)+
  #geom_boxplot()+
  theme(axis.text.x = element_text(angle = 180)) +
  #scale_y_continuous(limits = c(600, 1100))+
  theme(text = element_text(size = 20))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  #facet_wrap(.~Strain)+
  labs(y = "Intensity")+
  stat_summary(Auxin, mapping = aes(x=Condition,y=Average.Intensity.Average), fun.y=mean, geom="point", shape=18,
               size=1, color="red")
Auxin_RY16 <- subset(Auxin, Condition == "RY16xRPL28_H20" |Condition == "RY16xRPL28_100uM_KNAA" )
anova <- aov(formula = Average.Intensity.Average ~ Condition, data = Auxin)
summary(anova)

tukey_results <- TukeyHSD(anova)
print(tukey_results)


