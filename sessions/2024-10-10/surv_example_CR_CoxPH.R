library(haven)
library(knitr)
library(dplyr)
library(survival)
library(ggplot2)
library(tibble)
library(ggsurvfit)
library(splines)
library(Epi)
library(labelled)

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
    exitdate = as.Date(pmin(emdate, dthdate, eosdate, na.rm = TRUE)),
    midate = as.Date(ifelse(midate == dthdate & !is.na(dthdate), NA, midate)),
    os_dur = as.numeric(difftime(pmin(dthdate, emdate,
                                      eosdate, na.rm = TRUE), examdate) / 365.25),
    status = ifelse(is.na(dthdate), 0, 1),
    cursmoker = factor(ifelse(smoking >= 3, 1, 0), labels = c("No", "Yes")),
    binsmoker = ifelse(smoking >= 3, 1, 0),
    fact_smok = factor(smoking, labels = c("Never", "Ex >5yrs",
                                              "Ex 1-4yrs", "Cur <15/day", "Cur >15/day" )),
    agein = as.numeric(difftime(examdate, birthdate)) / 365.25,
    agecat = as.factor(cut(agein,
                           breaks = c(0, 50, 55, 200),
                           labels = c("40-49", "50-54", "55+")))
  )
caerphilly_dat <- caerphilly_dat %>%
  set_variable_labels(exitdate = "Date of exit")
summary(caerphilly_dat$agecat)

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
  ggsurvfit(type = "survival") +
  labs(
    x = "Years",
    y = "MI cumulative incidence (censored for death)"
  ) +
  add_confidence_interval() +
  add_risktable()


### Cox PH regression
# smoking status at entry as covariate

model_smoking <- coxph(Surv(os_dur, status) ~ fact_smok, data = caerphilly_dat)
summary(model_smoking)

# adjusted for age at entry (categorical)

model_smoking_agecat <- coxph(Surv(os_dur, status) ~ fact_smok + agecat,
                       data = caerphilly_dat)
summary(model_smoking_agecat)


# adjusted for age at entry (natural cubic spline)
model_smoking_splage <- coxph(Surv(os_dur, status) ~ fact_smok + ns(agein, df = 2),
                              data = caerphilly_dat)
summary(model_smoking_splage)


## Checking PH assumption wrt smoking
survfit2(Surv(os_dur, status) ~ strata(fact_smok), data = caerphilly_dat) %>%
  ggsurvfit(type = "cloglog") +
  labs(
    x = "Years",
    y = "log(-log(S(t)))"
  ) +
  scale_x_log10()

## Checking PH assumption wrt current smoking
coxph(Surv(os_dur, status) ~ cursmoker, data = caerphilly_dat)
survfit2(Surv(os_dur, status) ~ strata(cursmoker), data = caerphilly_dat) %>%
  ggsurvfit(type = "cloglog") +
  labs(
    x = "Years",
    y = "log(-log(S(t)))"
  ) +
  scale_x_log10()

# Note: We use a binary numeric variable for smoking, not a factor
# Different HR before or after 5 yrs?
coxph(Surv(os_dur, status) ~ binsmoker + tt(binsmoker) + ns(agein, df = 2),
      data = caerphilly_dat,
      tt =function(x, t, ...) x * (t + 5)
)

# Effect changes log-linearly over time
coxph(Surv(os_dur, status) ~ binsmoker + tt(binsmoker) + ns(agein, df = 2),
      data = caerphilly_dat,
      tt =function(x, t, ...) x * log(t)
)




# Splitting on time of MI using Epi package (BXC)
# Note how we change dates to fractions of years

lex_obj_data <- caerphilly_dat %>%
  mutate(
    examyr = cal.yr(examdate),
    birthyr = cal.yr(birthdate),
    exityr = cal.yr(exitdate),
    dthyr = cal.yr(dthdate),
    miyr = cal.yr(midate)
  ) %>%
  select(id, examyr, birthyr, exityr, dthyr, miyr, status, smoking,
         fact_smok, agein)

# Setting up Lexis object

dmL <- Lexis(entry = list(per = examyr,
                          age = examyr - birthyr,
                          tfIncl = 0),
             exit = list(per = exityr),
             exit.status = factor(!is.na(dthyr),
                                  labels = c("Cens", "Dead")),
             data = lex_obj_data)
str(dmL)


dmC <- cutLexis(data = dmL,
                cut = dmL$miyr,
                timescale = "per",
                new.state = "MI",
                new.scale = "tfMI")

# Same analysis using original data-set...
coxph(Surv(age, age+lex.dur, lex.Xst == "Dead") ~ fact_smok,
      data = dmL)
# ... or cut data-set
coxph(Surv(age, age+lex.dur, lex.Xst == "Dead") ~ fact_smok,
      data = dmC)

# In cut data-set, we can examine effect of having had an MI
coxph(Surv(age, age+lex.dur, lex.Xst == "Dead") ~ fact_smok + lex.Cst,
      data = dmC)



# Using standard splitting of survival follow-up time (harder and less
# flexible)

split_data <- survSplit(Surv(os_dur, status) ~ .,
                        data = caerphilly_dat,
                        cut = c(0, 5, 10, 20),
                        episode = "FUper")

split_data <- split_data %>%
  mutate(
    FUper = factor(FUper)
  )

period_smoker <- coxph(Surv(tstart, os_dur, status) ~ 
                         binsmoker : FUper + ns(agein, df = 2),
      data = split_data)

summary(period_smoker)

