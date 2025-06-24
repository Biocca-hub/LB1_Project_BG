library(ggplot2)
data <- read.csv("MCC.csv", sep=";")
ggplot(data, aes(x = E.Value, y = MCC, color = factor(Set), group = Set)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = data$E.Value) + 
  scale_color_manual(values = c("paleturquoise", "darkturquoise")) +  
  labs(title = "MCC vs E-Value",
       x = "E-Value (log-scale)",
       y = "MCC",
       color = "Set") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Center and format title
    legend.position = c(0.85, 0.2),  # Bottom right (x=85%, y=20%)
    legend.background = element_rect(fill = "white", color = "gray", size = 0.3)
  )
