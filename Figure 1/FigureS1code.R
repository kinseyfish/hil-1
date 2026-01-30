df <- data.frame(
  time = factor(c("0 hr", "1 hr", "3 hr", "6 hr"),
                levels = c("0 hr", "1 hr", "3 hr", "6 hr")),
  mean = c(444.563, 7.3994, 7.80738, 7.56978),
  ci_low = c(211.811, 5.56984, 5.83379, 4.07466),
  ci_high = c(677.315, 9.22896, 9.78098, 11.0649)
)
ggplot(df, aes(x = time, y = mean)) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.2
  ) +
  labs(
    x = "Hours of recovery",
    y = "hil-1 expression"
  ) +
  theme_classic()+ 
  geom_line(group = 1)


