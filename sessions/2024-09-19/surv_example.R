library(haven)
library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(ggsurvfit)

# some of this is inspired by:
# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html

# Change this to match your disk layout
caerphilly_dat <- read_dta("D:/stovring/SDCA/EpiSpace_EpiStats/Survival_I/Caerphilly1.dta")

# Select relevant variables
caerphilly_dat <- caerphilly_dat %>%
  select(id:diabetes, smoking)

# Adding some basic variables on survival in study and smoking status
caerphilly_dat <- caerphilly_dat %>%
  mutate(
    os_dur = as.numeric(difftime(pmin(dthdate, emdate, eosdate, na.rm = TRUE), examdate) / 365.25),
    status = ifelse(is.na(dthdate), 0, 1),
    cursmoker = ifelse(smoking >= 3, 1, 0)
  )

# A survival function object
s1 <- survfit(Surv(os_dur, status) ~ 1, data = caerphilly_dat)

str(s1)

# Simple Kaplan-Meier plot, time in study as time scale

survfit2(Surv(os_dur, status) ~ 1, data = caerphilly_dat) %>%
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  ) +
  add_confidence_interval() +
  add_risktable()

# As above, but divided into current smokers (1) vs non/ex-smokers (0)
survfit2(Surv(os_dur, status) ~ cursmoker, data = caerphilly_dat) %>%
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  ) +
  add_confidence_interval() +
  add_risktable()

# Log-rank test for difference in survival
survdiff(Surv(os_dur, status) ~ cursmoker, data = caerphilly_dat)


# Age as time scale - creating variables
caerphilly_dat <- caerphilly_dat %>%
  mutate(
    agein = as.numeric(difftime(examdate, birthdate)) / 365.25,
    ageout = as.numeric(difftime(
      pmin(dthdate, emdate, eosdate, na.rm = TRUE), birthdate)) / 365.25
    )

# Kaplan-Meier plot, age as time scale, note how the x-axis is truncated
survfit2(Surv(agein, ageout, status) ~ 1, data = caerphilly_dat) %>%
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  ) +
  scale_x_continuous(expand = c(0,0), limits = c(45, 85),
                     breaks = seq(45, 85, 10)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1),
                     breaks = seq(0, 1, .2)) +
  add_confidence_interval() +
  add_risktable()

# Kaplan-Meier plot, age as time scale, by current smoking status
survfit2(Surv(agein, ageout, status) ~ cursmoker, data = caerphilly_dat) %>%
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  ) +
  scale_x_continuous(expand = c(0,0), limits = c(45, 80),
                     breaks = seq(45, 80, 10)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1),
                     breaks = seq(0, 1, .2)) +
  add_confidence_interval() +
  add_risktable()

# Example with competing risk
library(tidycmprsk)
# Time in study as time scale - time to MI (0 censored, 1 MI, 2 death)
caerphilly_dat <- caerphilly_dat %>%
  mutate(
   midur = as.numeric(difftime(
      pmin(dthdate, midate, emdate, eosdate, na.rm = TRUE), examdate)) / 365.25,
    mistatus = as.factor(ifelse(!is.na(midate), 1, 0)
      + 2 * ifelse(!is.na(dthdate) & is.na(midate), 1, 0)),
    mistatus_num = as.numeric(mistatus=="1")
  )

cuminc(Surv(midur, mistatus) ~ 1, data = caerphilly_dat) %>%
  ggcuminc(outcome = c("1")) +
  ylim(c(0, .3)) +
  labs(
    x = "Years"
  )

# Incorrect - death used as censoring
survfit2(Surv(midur, mistatus_num) ~ 1, data = caerphilly_dat) %>%
  ggsurvfit() +
  labs(
    x = "Years",
    y = "MI cumulative incidence (censored for death)"
  ) +
  add_confidence_interval() +
  add_risktable()

