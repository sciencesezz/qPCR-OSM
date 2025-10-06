
getwd()
setwd("/mnt/data/home/sarahsczelecki/qPCR-OSM")
library(tidyverse)
library(ggplot2)
library(tidyr)

#copied from excel sheet that had done all the normalisation, so I just need to graph it.
data <- read.csv("/mnt/data/home/sarahsczelecki/input-files/inhibitors-for-R.csv")

head(data)
str(data)

exps <- unique(data$Experiment)

for (exp in exps) {
  
  #subset
  data_exp <- subset(data, Experiment == exp)
  
  #set levels
  s.levels <- unique(  data_exp$Condition)
  data_exp$Condition <- factor(  data_exp$Condition, levels = s.levels)
  
  #mean and SEM per experiment
  summary_exp <- data_exp %>%
    group_by(Condition) %>%
    summarise(
      mean_RelExp = mean(RelExp, na.rm = TRUE),
      sem_RelExp = sd(RelExp, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )

  #graph bars
  # Plot
  bars <- ggplot() +
    # Mean bars
    geom_col(data = summary_exp, aes(x = Condition, y = mean_RelExp, fill = Condition), alpha = 1) +
    # SEM error bars
    geom_errorbar(data = summary_exp, 
                  aes(x = Condition, ymin = mean_RelExp - sem_RelExp, ymax = mean_RelExp + sem_RelExp),
                  width = 0.2, color = "black") +
    # Individual points
    geom_jitter(data = data_exp, aes(x = Condition, y = RelExp),width = 0.15, size = 1, color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed") +  # baseline = 1
    # facet_wrap(~ Target, scales = "free_y") +          # one panel per gene
    scale_fill_brewer(palette = "PuRd") +
    theme_light(base_size = 25) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    labs(y = "Fold Change", x = "Sample Group")
  print(bars)
  
  ggsave(filename = paste0("cell-proliferation-", exp, ".png"), plot = bars, width = 8, height = 8, dpi = 800)
 
  #boxplot
  box <- ggplot(data = data_exp, aes(x = Condition, y = RelExp, fill = Condition)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, color = "black") +
    #geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_brewer(palette = "PuRd") +
    theme_bw(base_size = 25) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          #axis.ticks.x = element_blank(),
          legend.position = "none") +
    labs(y = "Fold Change", x = "Sample Group")
  print(box)
  ggsave(filename = paste0("cell-proliferation-box-", exp, ".png"), plot = box, width = 8, height = 8, dpi = 800)
   
}

data$logRelExp <- log(data$RelExp)

#statistical analysis
#Check normality assumptions
#histogram - I think not as good for relatively small n
library(car)
#q-q
ggplot(data, aes(sample = RelExp)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~Experiment, scales = "free") +
  theme_minimal()


# Output file path
output_file <- "inhibitors-stats-logtrans.txt"

# Start redirecting output
sink(output_file)

#Shapiro Wilk Test for normality
by(data$logRelExp, data$Experiment, shapiro.test)

#exps defined above
for (exp in exps) {
  
  cat("\n===== Experiment:", exp, "=====\n")
  
  # Subset the gene
 inhib_data <- subset(data, Experiment == exp)
  
  # Levene's test
  print("Levene's Test:")
  print(leveneTest(logRelExp ~ Condition, data = inhib_data))
  
  # ANOVA
  anova_result <- aov(logRelExp ~ Condition, data = inhib_data)
  print("ANOVA Summary:")
  print(summary(anova_result))
  
  # Tukey post-hoc
  print("Tukey HSD:")
  print(TukeyHSD(anova_result))
  #nonparametric 
  print("Kruskal Wallis:")
  print(kruskal.test(logRelExp ~ Condition, data = inhib_data))
  
  # Post-hoc pairwise Wilcoxon (Holm correction)
  cat("\nPost-hoc Pairwise Wilcoxon (Holm):\n")
  print(pairwise.wilcox.test(
    inhib_data$logRelExp, 
    inhib_data$Condition, 
    p.adjust.method = "BH"
  ))
  
}

# Stop redirecting output
sink()

