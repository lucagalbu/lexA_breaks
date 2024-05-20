library(ggplot2)
library(dplyr)

# Constitutive expression -------------------------------------------------
data <- read.table("results/constitutiveExpressers", sep = ";", header = TRUE) |> tibble()
rm <- unique(data$rm)
dm <- unique(data$dm)
rp <- unique(data$rp)
rf <- unique(data$rf)

expected <- data.frame(time = seq(0, max(data$time), length.out = 3000))
expected["mrna"] = rm/dm*(1-exp(-dm*expected$time))

ggplot(mapping = aes(time, mrna)) + 
  geom_line(data = data, aes(col = "Simulation")) + 
  geom_line(data = expected, aes(col = "Expected")) +
  theme_minimal() +
  xlab("Time [min]") +
  ggtitle("mRNA - Constitutive expresser")

data |> ggplot(aes(time, protein)) + 
  geom_line(aes(col = "Simulation")) +
  geom_hline(yintercept = rp*rm/(rf*dm), aes(col = "Expected mean")) +
  theme_minimal() +
  xlab("Time [min]") +
  ggtitle("Protein - Constitutive expresser")


rm(data, rm, dm, rp, rf, expected)

# Break repair distribution -----------------------------------------------
data <- read.table("results/breakRepairRate", sep = ";", header = TRUE) |> 
  as_tibble() |> rename(repairRate = rr) |> 
  group_by(iter, repairRate) |> filter(time == max(time)) |> 
  rename(repairTime = time)

assertthat::assert_that(unique(data$reaction) == "breakRepaired")

r_exp <- data |> 
  group_by(repairRate) |> 
  reframe(repairTime = rexp(n = n(), rate = unique(repairRate)))
ggplot(mapping = aes(repairTime, log(1-..y..))) + 
  stat_ecdf(data = data, aes(col = "Simulation")) + 
  stat_ecdf(data = r_exp, aes(col = "R exponential numbers")) + 
  geom_abline(data = data, aes(slope = -repairRate, intercept = 0)) + 
  facet_wrap(~repairRate, scales = "free_x", labeller = label_both, ncol = 2) + 
  theme_minimal() +
  xlab("Repair time [min]") + ylab("Log CDF")

rm(data, r_exp)

# LexA dynamics -----------------------------------------------------------
data <- read.table("results/lexaDynamics", sep = ";", header = TRUE) |>
  as_tibble() |> 
  filter(time < 10)

lss <- unique(data$lss)
ld <- unique(data$ld)

repairTime = data |> filter(reaction == "breakRepaired") |> pull(time)

expected <- bind_rows(
  tibble(time = seq(0, repairTime, length.out = 1500)) |> 
    mutate(lexa = lss*exp(-ld * time)),
  tibble(time = seq(repairTime, max(data$time), length.out = 1500)) |> 
    mutate(lexa = lss*(exp(-ld * repairTime)-1)*exp(-ld*(time-repairTime)) + lss)
)

ggplot(mapping = aes(time, lexa)) + 
  geom_line(data = expected, aes(col = "Expected"), linewidth = 1) + 
  geom_point(data = data, aes(col = "Simulation")) +
  geom_vline(aes(xintercept = repairTime, col = "Break repaired")) +
  theme_minimal() + xlab("Time [min]") + ylab("LexA binding [1/min]")

rm(data, repairTime, lss, ld, expected)


# LexA binding / unbinding ------------------------------------------------
data <- read.table("results/lexaRates", sep = ";", header = TRUE) |>
  as_tibble() |> 
  mutate(dt = time - lag(time)) |> 
  filter(!is.na(dt))

expected_slopes <- data.frame(
  reaction = c("lexaBinding", "lexaUnbinding"),
  slope = c(-unique(data$lexa), -unique(data$loff))
)

data |> ggplot(mapping = aes(dt, y = log(1-..y..))) + 
  geom_abline(data = expected_slopes, aes(slope = slope, intercept = 0)) +
  stat_ecdf(aes(col = "Simulation")) +
  facet_wrap(~reaction, scales = "free_x") +
  xlab("Time between two reactions [min]") +
  ylab("log reverse CDF") +
  theme_minimal() 

rm(data, expected_slopes)
