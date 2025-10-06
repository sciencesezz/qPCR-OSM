
getwd()
setwd("/mnt/data/home/sarahsczelecki/qPCR-OSM")
library(tidyverse)
library(ggplot2)
library(tidyr)

#copied from excel sheet that had done all the normalisation, so I just need to graph it.
data <- read.csv("/mnt/data/home/sarahsczelecki/input-files/Cell-prolieration-fig1a-for-R.csv")

head(data)
str(data)

data <- data %>%
  mutate(Condition = factor(Condition, 
                            levels = c("Control", "COCSM")))
levels(data$Condition)

data <- data %>%
  mutate(CellLine = factor(CellLine, 
                            levels = c("mOSE-T2", "mOSE-T2BR", "Ovcar-3")))
levels(data$CellLine)

summary <- data %>%
  group_by(Condition, CellLine) %>%
  summarise(
    mean_RelExp = mean(RelExp, na.rm = TRUE),
    sem_RelExp  = sd(RelExp, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )


colours <- c(
  "Control" = "grey50",  
  "COCSM"   = "mediumpurple")  

  
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
  facet_wrap(~ CellLine, scales = "free_y") +          # one panel per gene
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
  facet_wrap(~ CellLine, scales = "free_y") + 
  scale_fill_manual(values = colours) +
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  labs(y = "Fold Change", x = "Sample Group")
print(box)
ggsave(filename = "cell-proliferation2-box-xleg.png", plot = box, width = 6, height = 4, dpi = 800)

#loop over each cell line so that I have individual graphs:
# Get all unique genes
cell.lines <- unique(data$CellLine)

# Loop through each gene and save a plot
for (cell.line in cell.lines) {
  
  # Filter for this gene
  df_gene <- data %>% filter(CellLine == cell.line)
  
  # Create plot
  p <- ggplot(df_gene, aes(x = Condition, y = RelExp, fill = Condition)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, color = "black") +
    #geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_manual(values = colours) +
    theme_bw(base_size = 25) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    labs(y = "Fold-Change", x = "Sample Group",
         title = cell.line, )
  
  # Print (optional)
  print(p)
  
  # Save plot as PNG
  ggsave(filename = paste0("prolif-", cell.line, ".png"),
         plot = p, width = 6, height = 8, dpi = 800)
}


#statistical analysis
#Check normality assumptions
#histogram - I think not as good for relatively small n
library(car)
#q-q
ggplot(data, aes(sample = RelExp)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~CellLine, scales = "free") +
  theme_minimal()


# Output file path
output_file <- "cellproliferation-cell-line-stats.txt"

# Start redirecting output
sink(output_file)

#Shapiro Wilk Test for normality
by(data$RelExp, data$CellLine, shapiro.test)
#loop for levene's and anova
#Check homogeneity of variances(homoscedacity)
#groups should have roughly equal variances

#loop for levene's and anova
#Check homogeneity of variances(homoscedacity)
#groups should have roughly equal variances

cells <- unique(data$CellLine)

for (cell in cells) {
  
  cat("\n===== Cell Line", cell, "=====\n")
  
  # Subset the gene
 cell_data <- subset(data, CellLine == cell)
  
  # Levene's test
  print("Levene's Test:")
  print(leveneTest(RelExp ~ Condition, data = cell_data))
  
  # ANOVA
  anova_result <- aov(RelExp ~ Condition, data = cell_data)
  print("ANOVA Summary:")
  print(summary(anova_result))
  
  # Tukey post-hoc
  print("Tukey HSD:")
  print(TukeyHSD(anova_result))
  #nonparametric 
  print("Kruskal Wallis:")
  print(kruskal.test(RelExp ~ Condition, data = cell_data))
  
  # Post-hoc pairwise Wilcoxon (Holm correction)
  cat("\nPost-hoc Pairwise Wilcoxon (Holm):\n")
  print(pairwise.wilcox.test(
    cell_data$RelExp, 
    cell_data$Condition, 
    p.adjust.method = "holm"
  ))
  
}

# Stop redirecting output
sink()

library(ggplot2)
library(dplyr)
library(car)

# Check normality visually (QQ plots per cell line)

ggplot(data, aes(sample = RelExp)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~CellLine, scales = "free") +
  theme_minimal()

# Output file for all stats
output_file <- "cellproliferation-cell-line-stats.txt"
sink(output_file)

# Shapiro–Wilk test for each CellLine × Condition

cat("===== Shapiro–Wilk Normality Tests (by Cell Line × Condition) =====\n")
#Shapiro Wilk Test for normality
by(data$RelExp, data$CellLine, shapiro.test)


# Per-cell line loop: Levene + t-test or Wilcoxon
cells <- unique(data$CellLine)

for (cell in cells) {
  cat("\n\n============================\n")
  cat("Cell Line:", cell, "\n")
  cat("============================\n")
  
  # Subset data for this cell line
  cell_data <- subset(data, CellLine == cell)
  
  
  # Levene's test
  print("Levene's Test:")
  print(leveneTest(RelExp ~ Condition, data = cell_data))
  
  #parametric t-test
  cat("\nParametric t-test (two-sided, equal variance assumed):\n")
  print(t.test(RelExp ~ Condition, data = cell_data, var.equal = TRUE))
  
   # Nonparametric Wilcoxon rank-sum test
  cat("\nNonparametric Wilcoxon rank-sum test:\n")
  print(wilcox.test(RelExp ~ Condition, data = cell_data, exact = FALSE))
}

# Stop redirecting output
sink()
