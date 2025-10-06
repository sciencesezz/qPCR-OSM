
getwd()
setwd("/mnt/data/home/sarahsczelecki/qPCR-OSM")
library(tidyverse)
library(ggplot2)
library(tidyr)

#copied from excel sheet that had done all the normalisation, so I just need to graph it.
data <- read.csv("/mnt/data/home/sarahsczelecki/input-files/cell-proliferation-for-R.csv")


head(data)
str(data)

data <- data %>%
  mutate(Condition = factor(Condition, 
                            levels = c("Control", "DNOSM", "COCSM", "COCUL", "OSF1", "OSF2")))
levels(data$Condition)

summary <- data %>%
  filter(Condition != "OSF1") %>%
  filter(Condition != "OSF2") %>%
  group_by(Condition) %>%
  summarise(
    mean_RelExp = mean(RelExp, na.rm = TRUE),
    sem_RelExp  = sd(RelExp, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

data <- data %>%
  filter(Condition != "OSF1") %>%
  filter(Condition != "OSF2") 

colours <- c(
  "Control" = "#A09E9F",  # blue
  "DNOSM"   = "#1E88E5",  # orange
  "COCSM"   = "#FFC107",  # green
  "COCUL"  = "#004D40")  # red
  
# Plot
bars <- ggplot() +
  # Mean bars
  geom_col(data = summary, aes(x = Condition, y = mean_RelExp, fill = Condition), alpha = 1) +
  # SEM error bars
  geom_errorbar(data = summary, 
                aes(x = Condition, ymin = mean_RelExp - sem_RelExp, ymax = mean_RelExp + sem_RelExp),
                width = 0.2, color = "black") +
  # Individual points
  geom_jitter(data = data, aes(x = Condition, y = RelExp),width = 0.15, size = 1, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +  # baseline = 1
 # facet_wrap(~ Target, scales = "free_y") +          # one panel per gene
  scale_fill_manual(values = colours) +
  theme_light(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  
  labs(y = "Fold Change", x = "Sample Group")
print(bars)

ggsave(filename = "cell-proliferation.png", plot = bars, width = 8, height = 8, dpi = 800)

#boxplot
box <- ggplot(data = data, aes(x = Condition, y = RelExp, fill = Condition)) +
  geom_boxplot(alpha = 1, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, color = "black") +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = colours) +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right") +
  labs(y = "Fold Change", x = "Sample Group")
print(box)
ggsave(filename = "cell-proliferation-box-leg.png", plot = box, width = 8, height = 8, dpi = 800)


##CCSM experiment
#copied from excel sheet that had done all the normalisation, so I just need to graph it.
data <- read.csv("/mnt/data/home/sarahsczelecki/input-files/cell-proliferation-CCSM-Exp-for-R.csv")

head(data)
str(data)

data <- data %>%
  mutate(Condition = factor(Condition, 
                            levels = c("Control", "DNOSM", "CCSM", "2xCCSM", "OSF")))
levels(data$Condition)

summary <- data %>%
  filter(Condition != "OSF") %>%
  group_by(Condition) %>%
  summarise(
    mean_RelExp = mean(RelExp, na.rm = TRUE),
    sem_RelExp  = sd(RelExp, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )

data <- data %>%
  filter(Condition != "OSF")

colours <- c(
  "Control" = "#A09E9F",  # blue
  "DNOSM"   = "#1E88E5",  # orange
  "CCSM"   = "#47CE27",  # green
  "2xCCSM"  = "#882E7F")  # red

# Plot
bars <- ggplot() +
  # Mean bars
  geom_col(data = summary, aes(x = Condition, y = mean_RelExp, fill = Condition), alpha = 1) +
  # SEM error bars
  geom_errorbar(data = summary, 
                aes(x = Condition, ymin = mean_RelExp - sem_RelExp, ymax = mean_RelExp + sem_RelExp),
                width = 0.2, color = "black") +
  # Individual points
  geom_jitter(data = data, aes(x = Condition, y = RelExp),width = 0.15, size = 1, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +  # baseline = 1
  # facet_wrap(~ Target, scales = "free_y") +          # one panel per gene
  scale_fill_manual(values = colours) +
  theme_light(base_size = 25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  
  labs(y = "Fold Change", x = "Sample Group")
print(bars)

ggsave(filename = "cell-proliferation-ccsm.png", plot = bars, width = 8, height = 8, dpi = 800)

#boxplot
box <- ggplot(data = data, aes(x = Condition, y = RelExp, fill = Condition)) +
  geom_boxplot(alpha = 1, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, color = "black") +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = colours) +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  labs(y = "Fold Change", x = "Sample Group")
print(box)
ggsave(filename = "cell-proliferation-ccsm-box.png", plot = box, width = 8, height = 8, dpi = 800)


#statistical analysis
#Check normality assumptions
#histogram - I think not as good for relatively small n
library(car)
#q-q
ggplot(data, aes(sample = RelExp)) +
  stat_qq() + stat_qq_line() +
  #facet_wrap(~Target, scales = "free") +
  theme_minimal()


# Output file path
output_file <- "cellproliferation-CCSM-stats.txt"

# Start redirecting output
sink(output_file)

#Shapiro Wilk Test for normality
shapiro.test(data$RelExp)

#loop for levene's and anova
#Check homogeneity of variances(homoscedacity)
#groups should have roughly equal variances

# Levene's test
print("Levene's Test:")
print(leveneTest(RelExp ~ Condition, data = data))
  
# ANOVA
anova_result <- aov(RelExp ~ Condition, data = data)
print("ANOVA Summary:")
print(summary(anova_result))
  
# Tukey post-hoc
print("Tukey HSD:")
print(TukeyHSD(anova_result))
#nonparametric 
print("Kruskal Wallis:")
print(kruskal.test(RelExp ~ Condition, data = data))
  
# Post-hoc pairwise Wilcoxon (Holm correction)
cat("\nPost-hoc Pairwise Wilcoxon (Holm):\n")
print(pairwise.wilcox.test(
data$RelExp, 
data$Condition, 
p.adjust.method = "holm"
  ))

# Stop redirecting output
sink()
  



