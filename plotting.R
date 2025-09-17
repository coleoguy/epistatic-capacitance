# Andres Barboza
# 
# Figures

library(dplyr)
library(tidyr)
library(ggplot2)
library(wesanderson)
library(purrr)

# Helper: tidy long
tidy_variance <- function(df, xvar) {
  df %>%
    pivot_longer(
      cols = matches("^(mean|sd)_"),
      names_to = c(".value","component"),
      names_pattern = "(mean|sd)_(.*)"
    ) %>%
    mutate(
      mean      = pmax(mean, 0),
      ymin      = pmax(mean - sd, 0),
      ymax      = pmin(mean + sd, 1),
      component = factor(component,
                         levels=c("add","dom","epi"),
                         labels=c("Additive","Dominance","Epistasis"))
    ) %>%
    rename(!!xvar := all_of(xvar))
}

# Palettes
pal_arch <- wes_palette("FantasticFox1", length(arches), type="continuous")
names(pal_arch) <- arches
# define color palette for the three components
pal_comp <- wes_palette("FantasticFox1", 3, type = "continuous")
names(pal_comp) <- c("Additive", "Dominance", "Epistasis")




loci_long <- tidy_variance(df_loci, "loci")
#loci_long <- readRDS("~/Documents/GitHub/epistatic-capacitance/results/loci")

p1 <- ggplot(loci_long, aes(x = loci, y = mean,
                            color = component, fill = component)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~architecture, nrow = 1) +
  labs(
    x     = "Number of Loci",
    y     = "Genetic Variance",
    title = "Effect of Number of Loci on Genetic Variance",
    color = "Variance\nComponent",
    fill  = "Variance\nComponent"
  ) +
  scale_color_manual(values = pal_comp) +
  scale_fill_manual(values = pal_comp) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(min(loci_range), max(loci_range), length.out = 5)) +
  theme_minimal() +
  theme(
    panel.border    = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "bottom",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    plot.title      = element_text(size = 20),
    axis.title      = element_text(size = 16),
    axis.text       = element_text(size = 14)
  )



sigma_long <- tidy_variance(df_sigma, "sigma")
#sigma_long <- readRDS("~/Documents/GitHub/epistatic-capacitance/results/sigma")

p2 <- ggplot(sigma_long, aes(x = sigma, y = mean,
                             color = component, fill = component)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~architecture, nrow = 1) +
  labs(
    x     = "Selection Strength (σ)",
    y     = "Genetic Variance",
    title = "Effect of Selection Strength on Genetic Variance",
    color = "Variance\nComponent",
    fill  = "Variance\nComponent"
  ) +
  scale_color_manual(values = pal_comp) +
  scale_fill_manual(values = pal_comp) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(1, 10, by = 2)) +
  theme_minimal() +
  theme(
    panel.border    = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "bottom",
    legend.title    = element_text(size = 14),
    legend.text     = element_text(size = 12),
    plot.title      = element_text(size = 20),
    axis.title      = element_text(size = 16),
    axis.text       = element_text(size = 14)
  )


# and tidy up for plotting
N_long <- tidy_variance(df_N, "N")


p3 <- ggplot(N_long, aes(x = N, y = mean,
                         color = component, fill = component)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~architecture, nrow = 1) +
  labs(
    x     = "Population Size (N)",
    y     = "Genetic Variance",
    title = "Effect of Population Size on Genetic Variance",
    color = "Variance\nComponent",
    fill  = "Variance\nComponent"
  ) +
  scale_color_manual(values = pal_comp) +
  scale_fill_manual(values = pal_comp) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = N_values) +
  theme_minimal() +
  theme(
    panel.border   = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "bottom",
    legend.title   = element_text(size = 14),
    legend.text    = element_text(size = 12),
    plot.title     = element_text(size = 20),
    axis.title     = element_text(size = 16),
    axis.text      = element_text(size = 14)
  )


# summarise for plot 4
df_time_sum <- df_time %>%
  group_by(architecture, bottleneck, time) %>%
  summarise(
    mean_diff = mean(diff),
    min_diff  = min(diff),
    max_diff  = max(diff),
    .groups   = "drop"
  )
pal_b <- c("FALSE" = "darkblue", "TRUE" = "darkred")

p4 <- ggplot(df_time_sum,
             aes(x = time, y = mean_diff,
                 color = factor(bottleneck),
                 fill  = factor(bottleneck))) +
  geom_ribbon(aes(ymin = min_diff, ymax = max_diff),
              alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~architecture, nrow = 1) +
  scale_color_manual(
    name   = "Bottleneck",
    values = pal_b,
    labels = c("No", "Yes")
  ) +
  scale_fill_manual(
    name   = "Bottleneck",
    values = pal_b,
    labels = c("No", "Yes")
  ) +
  labs(
    x     = "Generations after optimum shift",
    y     = expression("|mean phenotype – optimum|"),
    title = "Response to selection with vs. without bottleneck"
  ) +
  theme_minimal() +
  theme(
    panel.border    = element_rect(color = "darkgray", fill = NA, size = 1),
    legend.position = "bottom",
    plot.title      = element_text(size = 18),
    axis.title      = element_text(size = 14),
    axis.text       = element_text(size = 12)
  )

# Print plots
print(p1)
print(p2)
print(p3)
print(p4)



