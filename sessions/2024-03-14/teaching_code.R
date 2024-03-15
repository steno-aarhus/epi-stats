#Acknowledgement: This material is based on a workshop organized by CausalLab Karolinska Institutet with Anthony Matthews and Conor MacDonald  

library(data.table)
library(MASS)
#Loading the data
library(haven)
library(here)


### load in example data
targetpop <- read.csv2(here("data", "targetpop.csv")) 
View(targetpop)

studypop <- read.csv2(here("data", "studypop.csv"))
View(studypop)
summary(studypop)


table(studypop$treatment, studypop$Y)
prop.table(table(studypop$treatment, studypop$Y))
table(studypop$treatment, studypop$Y)
YA = 4698/(4698+85343)
YB = 6321/(6321+84129)
RD = YB - YA
RR = YB / YA

#risk difference according to hypertension status
# no hypertension
studypop_hta0 = subset(studypop, hypertension == 0)
prop.table(table(studypop_hta0$treatment, studypop_hta0$Y))
table(studypop_hta0$treatment, studypop_hta0$Y)
YA0 = 915/(915+31519)
YB0 = 1122/(1122+31453)
RD0 = YB0 - YA0
RR0 = YB0 / YA0


# with hypertension
studypop_hta1 = subset(studypop, hypertension == 1)
prop.table(table(studypop_hta1$treatment, studypop_hta1$Y))
table(studypop_hta1$treatment, studypop_hta1$Y)
YA1 = 3783/(3783+53824)
YB1 = 5199/(5199+52676)
RD1 = YB1 - YA1
RR1 = YB1 / YA1

# distribution of hypertension in the trial and target population
summary(studypop$hypertension)
summary(targetpop$hypertension)

#Fitting our model in the trial data, including interactions between treatment and effect modifiers
logistic_model <- glm('Y ~ treatment*(Sex + age + hypertension + heartdisease)', data = studypop, family = binomial())

#Making events occur in the target population under each treatment and assign them the correct treatment
### explain a little here
targetpopA = targetpop
targetpopA$treatment = "A"
targetpopA$predA <- predict(logistic_model, type = "response", newdata = targetpopA)

targetpopB = targetpop
targetpopB$treatment = "B"
targetpopB$predB <- predict(logistic_model, type = "response", newdata = targetpopB)

#estimate the average risk under each treatment 
YA = mean(targetpopA$predA)
YB = mean(targetpopB$predB)

# estimate the average treatment effects
RD = YB - YA
RR = YB/YA

RD
RR
