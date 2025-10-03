
getwd()
setwd("/mnt/data/home/sarahsczelecki/qPCR-OSM")
library(tidyverse)
library(ggplot2)
library(tidyr)

data <- read.csv("/mnt/data/home/sarahsczelecki/input-files/cell-cycle-genes-combined-for-R.csv")

head(data)
str(data)

#relabel calibrators
qpcr <- data %>%
  mutate(Sample = ifelse(Sample %in% c("Cal1", "Cal2"), "Cal", Sample))

# Define your endogenous control and calibrator
endogenous <- "Rpl19-1"   # replace with your reference gene
calibrator <- "Cal"       # replace with your baseline sample


#average technical replicates
qpcr_avg <- qpcr %>%
  group_by(Sample, Gene) %>%
  summarise(mean_Ct = mean(AvgCt, na.rm = TRUE))

#Spread so each sample has columns for each gene
qpcr_wide <- qpcr_avg %>%
  pivot_wider(names_from = Gene, values_from = mean_Ct)

# 4. Calculate ΔCt (Target - Endogenous control)
# assume multiple target genes, keep endogenous control separately
qpcr_dCt <- qpcr_wide %>%
  mutate(across(all_of(setdiff(names(qpcr_wide), c("Sample", endogenous))),
                ~ .x - .data[[endogenous]],
                .names = "dCt_{.col}"))

#Gather ΔCt values into long format
qpcr_dCt_long <- qpcr_dCt %>%
  select(Sample, starts_with("dCt_")) %>%
  pivot_longer(-Sample, names_to = "Target", values_to = "dCt") %>%
  mutate(Target = sub("dCt_", "", Target))

# ΔΔCt relative to calibrator sample
qpcr_ddCt <- qpcr_dCt_long %>%
  group_by(Target) %>%
  mutate(calibrator_dCt = dCt[Sample == calibrator],
         ddCt = dCt - calibrator_dCt,
         RelExp = 2^(-ddCt)) %>%
  ungroup()

# Final output
print(qpcr_ddCt)

#add condition column
qpcr_ddCt <- qpcr_ddCt %>%
  mutate(Condition = case_when(
    Sample %in% c("Ctrl1", "Ctrl2", "Ctrl3", "Ctrl4") ~ "Control",
    Sample %in% c("rCOCSM1", "rCOCSM2", "rCOCSM3", "rCOCSM4") ~ "COCSM",
    Sample %in% c("rDNO1", "rDNO2", "rDNO3", "rDNO4") ~ "DNOSM",
    Sample %in% c("rCoCul1", "rCoCul2", "rCoCul3", "rCoCul4") ~ "CO-CUL",
    Sample %in% c("GDF9:BMP15 1:1 - 1", "GDF9:BMP15 1:1 - 2", "GDF9:BMP15 1:1 - 3", "GDF9:BMP15 1:1 - 4") ~ "OSF1",
    Sample %in% c("GDF9:BMP15 1:.2 - 1", "GDF9:BMP15 1:.2 - 2", "GDF9:BMP15 1:.2 - 3", "GDF9:BMP15 1:.2 - 4") ~ "OSF2",
    TRUE ~ Sample   # keep other names as-is
  ))

#now make relative to the control group
# compute mean dCt for Control group per Target
control_means <- qpcr_ddCt %>%
  filter(Condition == "Control") %>%
  group_by(Target) %>%
  summarise(mean_control_rel = mean(RelExp, na.rm = TRUE), .groups = "drop")

qpcr_ddCt_norm <- qpcr_ddCt %>%
  left_join(control_means, by = "Target") %>%
  mutate(RelExp_ctrl = RelExp / mean_control_rel)

qpcr_ddCt_norm <- qpcr_ddCt_norm %>%
  mutate(Condition = factor(Condition, 
                            levels = c("Control", "DNOSM", "COCSM", "CO-CUL", "OSF1", "OSF2")))
levels(qpcr_ddCt_norm$Condition)

qpcr_summary <- qpcr_ddCt_norm %>%
  filter(Condition != "Calibrator") %>%
  filter(Condition != "OSF1") %>%
  filter(Condition != "OSF2") %>%
  group_by(Target, Condition) %>%
  summarise(
    mean_RelExp = mean(RelExp_ctrl, na.rm = TRUE),
    sem_RelExp  = sd(RelExp_ctrl, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  )
colours <- c(
  "Control" = "#A09E9F", 
  "DNOSM"   = "#1E88E5",  
  "COCSM"   = "#FFC107",  
  "CO-CUL"  = "#004D40") 

#box and whisker
qpcr_ddCt_norm <- qpcr_ddCt_norm %>%
  filter(Condition != "Calibrator") %>%
  filter(Condition != "OSF1") %>%
  filter(Condition != "OSF2")

# Plot
ggplot() +
  # Mean bars
  geom_col(data = qpcr_summary, aes(x = Condition, y = mean_RelExp, fill = Condition), alpha = 1) +
  # SEM error bars
  geom_errorbar(data = qpcr_summary, 
                aes(x = Condition, ymin = mean_RelExp - sem_RelExp, ymax = mean_RelExp + sem_RelExp),
                width = 0.2, color = "black") +
  # Individual points
  geom_jitter(data = qpcr_ddCt_norm %>% filter(Condition != "Calibrator"), 
              aes(x = Condition, y = RelExp_ctrl), 
              width = 0.15, size = 1, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +  # baseline = 1
  facet_wrap(~ Target, scales = "free_y") +          # one panel per gene
  scale_fill_manual(values = colours) +
  theme_light(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(y = "Relative Expression to Rpl19", x = "Sample Group")

# Get all unique genes
genes <- unique(qpcr_ddCt_norm$Target)

# Loop through each gene and save a plot
for (gene in genes) {
  
  # Filter for this gene
  df_gene <- qpcr_ddCt_norm %>% filter(Target == gene)
  
  # Create plot
  p <- ggplot(df_gene, aes(x = Condition, y = RelExp_ctrl, fill = Condition)) +
    geom_boxplot(alpha = 1, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2, color = "black") +
    #geom_hline(yintercept = 1, linetype = "dashed") +
    scale_fill_manual(values = colours) +
    theme_bw(base_size = 25) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none") +
    labs(y = "Relative Expression to Rpl19", x = "Sample Group",
         title = gene, )
  
  # Print (optional)
  print(p)
  
  # Save plot as PNG
  ggsave(filename = paste0("genes-", gene, ".png"),
         plot = p, width = 6, height = 8, dpi = 800)
}

#statistical analysis
library(car)
#Check normality assumptions
#histogram - I think not as good for relatively small n
ggplot(qpcr_ddCt_norm, aes(x = RelExp_ctrl, fill = Target)) +
  geom_histogram(color = "black", bins = 10, alpha = 0.6) +
  facet_wrap(~Target, scales = "free") +
  theme_minimal()

#q-q
ggplot(qpcr_ddCt_norm, aes(sample = RelExp_ctrl)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~Target, scales = "free") +
  theme_minimal()

# Output file path
output_file <- "cellcycle-stats.txt"

# Start redirecting output
sink(output_file)

#Shapiro Wilk Test for normality
by(qpcr_ddCt_norm$RelExp_ctrl, qpcr_ddCt_norm$Target, shapiro.test)
#shapiro.test(ratio$Ratio_Bax_Bcl2)

#loop for levene's and anova
#Check homogeneity of variances(homoscedacity)
#groups should have roughly equal variances

genes <- c("Cdk2", "CyclinA2", "CyclinD1", "CyclinE1", "p21", "p27")

for (gene in genes) {
  
  cat("\n===== Gene:", gene, "=====\n")
  
  # Subset the gene
  gene_data <- subset(qpcr_ddCt_norm, Target == gene)
  
  # Levene's test
  print("Levene's Test:")
  print(leveneTest(RelExp_ctrl ~ Condition, data = gene_data))
  
  # ANOVA
  anova_result <- aov(RelExp_ctrl ~ Condition, data = gene_data)
  print("ANOVA Summary:")
  print(summary(anova_result))
  
  # Tukey post-hoc
  print("Tukey HSD:")
  print(TukeyHSD(anova_result))
  #nonparametric 
  print("Kruskal Wallis:")
  print(kruskal.test(RelExp_ctrl ~ Condition, data = gene_data))
  
  # Post-hoc pairwise Wilcoxon (Holm correction)
  cat("\nPost-hoc Pairwise Wilcoxon (Holm):\n")
  print(pairwise.wilcox.test(
    gene_data$RelExp_ctrl, 
    gene_data$Condition, 
    p.adjust.method = "holm"
  ))
  
}

# Stop redirecting output
sink()

#nonparametric 
#kruskal.test(RelExp_ctrl ~ Condition, data = bax)



