setwd("/Users/kinseyfisher/Documents/hil-1/Figure 2 material/")


#Figure 2c: Crosses to daf-2/daf-16 GFP expression

#RY16 is syb5748 and RY17 is syb5764
#CF1038 is mu86
#F_liq = fed, S_liq = starved
Crossesconditions <- read.csv("hil1_expression_24hrs_crosses.csv", header = T)
Crossesconditions$Condition <- factor(Crossesconditions$Condition, levels = c("RY16_F_liq","RY16xCF1038_F", "RY16xE1370_F", "RY16_S_liq", "RY16xCF1038_S", "RY16xE1370_S",
                                                                              "RY17_F_Liq","RY17xCF1038_F","RY17xE1370_F","RY17_S_Liq","RY17xCF1038_S", "RY17xE1370_S"))

ggplot()+
  geom_dotplot(Crossesconditions, mapping = aes(x=Condition,y=Average.Intensity.Average, color = Condition, fill = Condition),inherit.aes = TRUE,binaxis = "y",alpha=1,binwidth=10,stackdir = "center", dotsize = 1)+
  #geom_boxplot()+
  theme(axis.text.x = element_text(angle = 180)) +
  theme(text = element_text(size = 20))+
  theme_bw()+
  theme(aspect.ratio = 1) +
  labs(y = "Intensity")+
  stat_summary(Crossesconditions, mapping = aes(x=Condition,y=Average.Intensity.Average), fun.y=mean, geom="point", shape=18,
               size=1, color="black")
#Stats
Crossesconditions_RY16 <- subset(exp24hrs, Strain == "RY16"| Strain == "RY16xCF1038" | Strain == "RY16xE1370")
anova <- aov(formula = Average.Intensity.Average ~ Condition, data = Crossesconditions_RY16)
summary(anova)

tukey_results <- TukeyHSD(anova)
print(tukey_results)

Crossesconditions_RY17 <- subset(exp24hrs, Strain == "RY17"| Strain == "RY17xCF1038" | Strain == "RY17xE1370")
anova <- aov(formula = Average.Intensity.Average ~ Condition, data = Crossesconditions_RY17)
summary(anova)

tukey_results <- TukeyHSD(anova)
print(tukey_results)


#Figure 2d: timeseries expression (starved)

Starved <- read.csv("hil1_timeseries_expression.csv", header = T)

Results_Sum_S <- Starved %>%
  group_by(Condition, Timepoint) %>%
  summarise(
    Intensity = mean(Average.Intensity.Average),
    SD = sd(Average.Intensity.Average),
    N = n(),# Count of replicates
    SE = SD / sqrt(N),# Standard Error
    CI_lower = Intensity - 1.96 * SE,# Lower 95% CI
    CI_upper = Intensity + 1.96 * SE# Upper 95% CI
  )
TS_plot <- ggplot(data = Results_Sum_S, aes(x = Timepoint, y = Intensity, color = Condition, fill = Condition)) +
  geom_line(show.legend = FALSE) +# Line for mean intensity
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, color = NA) +# Shaded 95% CI
  scale_x_continuous(breaks = seq(1, 24, by = 1)) +
  theme_bw()

TS_plot


#Figure 2e: timeseries expression (starved then fed)
starvedthenfed <- read.csv("hil1_starvedthenfed.csv", header = T)


starvedthenfed$timepoint <- as.numeric(starvedthenfed$timepoint)

Results_Sum <- starvedthenfed %>%
  group_by(Strain, timepoint, Condition) %>%
  summarise(
    Intensity = mean(Average.Intensity),
    SD = sd(Average.Intensity),
    N = n(),
    SE = SD / sqrt(N),
    CI95 = 1.96 * SE    # 95% CI for normal approx
  )

TS_plot <- ggplot(Results_Sum, aes(x = timepoint, y = Intensity, color = Strain, fill = Strain, linetype = Condition)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = Intensity - CI95, ymax = Intensity + CI95),
              alpha = 0.2, color = NA, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(1, 24, by = 1)) +
  theme_bw()

TS_plot


#Figure 2f: Nutrient dilutions
dilutions <- read.csv("hil1_nutrientdilutions.csv", header = T)
dilutions <- na.omit(dilutions)
dilutions$Food.Dilution <- factor(dilutions$Food.Dilution,levels = c("25", "12.5","6.25", "3.13", "1.6", "EtOH", "Starved"))
ggplot(
  dilutions,
  aes(x = Food.Dilution, y = Average.Intensity)
) +
  geom_dotplot(inherit.aes = TRUE, binaxis = "y", alpha = 1, binwidth = 1.5, stackdir = "center", dotsize = 3) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 15,           # solid square
    size = 1,
    color = "red",
    position = position_dodge(width = 0.3)
  ) +
  facet_grid(.~Strain)+
  labs(
    x = "Food Dilution",
    y = "GFP Intensity"
  ) +
  theme_minimal()


dilutions_filtered <- dilutions %>% filter(Strain == "syb5748")
anova_result <- aov(Average.Intensity ~ Food.Dilution, data = dilutions_filtered)
summary(anova_result)
TukeyHSD(anova_result)

dilutions_filtered <- dilutions %>% filter(Strain == "syb5764")
anova_result <- aov(Average.Intensity ~ Food.Dilution, data = dilutions_filtered)
summary(anova_result)
TukeyHSD(anova_result)






