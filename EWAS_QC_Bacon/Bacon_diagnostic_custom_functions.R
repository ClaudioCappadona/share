#### Function to plot Bacon mixture distribution ####

mixture_dist_plot <- function(bcobj){

bc <- bcobj
# mixture weights
p <- estimates(bc) %>% as.data.frame() %>% select(starts_with("p.")) %>% unlist() %>% as.numeric()

# means (bias)
mu <- estimates(bc) %>% as.data.frame() %>% select(starts_with("mu.")) %>% unlist() %>% as.numeric()

# standard deviations (inflation)
sigma <- estimates(bc) %>% as.data.frame() %>% select(starts_with("sigma.")) %>% unlist() %>% as.numeric()

# X grid spanning all distributions
xgrid <- seq(-5, 5, length.out = 500)

# Individual component densities
dens0 <- dnorm(xgrid, mu[1], sigma[1])
dens1 <- dnorm(xgrid, mu[2], sigma[2])
dens2 <- dnorm(xgrid, mu[3], sigma[3])

# Mixture density (weighted sum using your p values)
mixture <- p[1]*dens0 + p[2]*dens1 + p[3]*dens2

# build data frame for ggplot
plotdf <- data.frame(
  x = xgrid,
  mixture = mixture,
  comp0 = p[1] * dens0,
  comp1 = p[2] * dens1,
  comp2 = p[3] * dens2
) %>%
  pivot_longer(
    cols = c(mixture, comp0, comp1, comp2),
    names_to = "component",
    values_to = "density"
  )

# labels for legend (same as base‑R legend order)
lab_mixture <- "Mixture"
lab_comp0   <- paste0("Comp 0 (", round(100*p[1], digits = 3), "%)")
lab_comp1   <- paste0("Comp 1 (", round(100*p[2], digits = 3), "%)")
lab_comp2   <- paste0("Comp 2 (", round(100*p[3], digits = 3), "%)")

plotdf$label <- case_when(
  plotdf$component == "mixture" ~ lab_mixture,
  plotdf$component == "comp0"   ~ lab_comp0,
  plotdf$component == "comp1"   ~ lab_comp1,
  plotdf$component == "comp2"   ~ lab_comp2
)

# ensure legend order matches original legend order
plotdf$label <- factor(plotdf$label,
                   levels = c(lab_mixture, lab_comp0, lab_comp1, lab_comp2))

# colors and line widths
cols <- c(
  lab_mixture = "black",
  lab_comp0   = "#439161",
  lab_comp1   = "#E69F00",
  lab_comp2   = "#56B4E9"
)

cols <- c("black", "#439161", "#E69F00", "#56B4E9")

lwds <- c(
  lab_mixture = 4,
  lab_comp0   = 2,
  lab_comp1   = 2,
  lab_comp2   = 2
)

ggplot(plotdf, aes(x = x, y = density, colour = label, size = label)) +
  geom_line(size=1.2) +
  scale_color_manual(values = cols) +
  #scale_size_manual(values = lwds) +
  labs(
    x = "x",
    y = "Density",
    title = "Normal Mixture Distribution",
    colour = NULL,
    size = NULL
  ) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.position = c(0.95, 0.95),
        legend.justification = c(1, 1))

}



#### Function to make tstat correction comparison plot ####

tstat_corr_comparison <- function(bcobj) {

bc <- bcobj

# compute density objects
d_uc <- density(tstat(bc, corrected = FALSE))
d_c  <- density(tstat(bc, corrected = TRUE))

# make a data frame from the first density for plotting
df_uc <- data.frame(x = d_uc$x, y = d_uc$y, group = "Not corrected")
df_c  <- data.frame(x = d_c$x, y = d_c$y, group = "Corrected")

df <- rbind(df_uc, df_c)

# draw the ggplot version
p <- ggplot(df, aes(x = x, y = y, color = group)) +
  geom_line(lwd = 1.5) +
  labs(
    x = "tstat",
    y = "Density",
    title = paste0("T-stat correction"),
    color = NULL
  ) +
  scale_color_manual(
    values = c("Not corrected" = "blue", "Corrected" = "red")
  ) +
  
  annotate("text", x=0.90*(min(df$x)),y=0.88*(max(df$y)), hjust = 0, color="black",size=3.5,parse=TRUE,label= c(bquote(sigma[Bacon] ~ "=" ~.(as.numeric(round(inflation(bc), 2)))))) +
  annotate("text", x=0.90*(min(df$x)),y=0.82*(max(df$y)), hjust = 0, color="black",size=3.5,parse=TRUE,label= c(bquote(mu[Bacon] ~ "=" ~.(as.numeric(round(bias(bc), 2)))))) +        
         
  theme_minimal() +
  theme(legend.title = element_blank(), 
        legend.position = c(0.97, 0.97),
        legend.justification = c(1, 1))


p + labs(
  title = paste0("T-stat correction"))
}

