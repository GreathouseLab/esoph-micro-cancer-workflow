
```{r}

library(lme4) # for multilevel models
library(tidyverse) # for data manipulation and plots
library(effects) #for plotting parameter effects
library(jtools) #for transformaing model summaries
library(eha) #for the data sets used in this tutorial
library(discSurv) #discrete-time survival analysis tool kit
library(brms) #Bayesian discrete-time survival analysis


# generate long format dataset
full_dat_long <- data.frame(matrix(nrow=sum(full_dat$etime),ncol=ncol(full_dat)+3))
colnames(full_dat_long) <- c("enter", "exit", "event", colnames(full_dat))

i <- j <- 1; k <- 0
ids <- full_dat$ID
for(i in 1:length(ids)){
  
  svec <- full_dat[full_dat$ID == ids[i], , drop=T]
  k <- k + as.numeric(svec["etime"])
  full_dat_long$enter[j:k] <- seq(0, as.numeric(svec["etime"]-1), 1)
  full_dat_long$exit[j:k] <- full_dat_long$enter[j:k] + 1
  full_dat_long$event[j:k] <- c(rep(0, as.numeric(svec["etime"]-1)),1)
  full_dat_long[j:k,-c(1:3)] <- svec
  
  j <- k + 1
}

```


### BAseline Model

```{r}

# fit Gompertz Curve
full_dat_long %>%
  group_by(exit) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  filter(event > 0) %>%
  ggplot(aes(x = exit, y = log(-log(1-hazard)))) +
  geom_point() +
  geom_smooth()

Gompertz_Model_Baseline <- glm(formula = event ~ exit,
                               family = binomial(link = "cloglog"),
                               data = full_dat_long)

summary(Gompertz_Model_Baseline)
summ(Gompertz_Model_Baseline, exp = T) #exp = T means that we want exponentiated estimates


full_dat_long %>%
  group_by(exit, Abundance) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  filter(event > 0) %>%
  ggplot(aes(x = exit, 
             y = log(-log(1-hazard)),
             col = Abundance)) +
  geom_point() +
  geom_smooth()



```






library(lme4) # for multilevel models
library(tidyverse) # for data manipulation and plots
library(effects) #for plotting parameter effects
library(jtools) #for transformaing model summaries
library(eha) #for the data sets used in this tutorial
library(discSurv) #discrete-time survival analysis tool kit
library(brms) #Bayesian discrete-time survival analysis

data("scania")
head(scania)

Scania_Person <- scania %>%
  mutate(exit = ceiling(exit),
         birthdate = floor(birthdate),
         spell = exit - enter) %>% #spell refers to the observed duration of a person
  mutate(enter = enter - 50,
         exit = exit - 50)

head(Scania_Person)

set.seed(123)
Scania_Person_Train <- sample_frac(Scania_Person, 0.8)
Scania_Person_Test <- Scania_Person[!Scania_Person$id %in% Scania_Person_Train$id,]

#convert the training set
Scania_PersonPeriod_Train <- dataLong(
  dataSet = Scania_Person_Train, 
  timeColumn = "spell", 
  censColumn = "event",
  timeAsFactor = F) %>%
  as_tibble() %>%
  mutate(enter = timeInt - 1,
         age = enter + 50) %>%
  dplyr::select(-obj, -event, -exit) %>%
  rename(event = y,
         exit = timeInt) %>%
  mutate(year = age + birthdate) %>%
  dplyr::select(id, enter, exit, event, everything()) %>%
  left_join(logrye, by = "year") #joined with the `logrye` data for a variable on yearly food prices

head(Scania_PersonPeriod_Train)
#convert the test set
Scania_PersonPeriod_Test <- dataLong(
  dataSet = Scania_Person_Test, 
  timeColumn = "spell", 
  censColumn = "event",
  timeAsFactor = F) %>%
  as_tibble() %>%
  mutate(enter = timeInt - 1,
         age = enter + 50) %>%
  dplyr::select(-obj, -event, -exit) %>%
  rename(event = y,
         exit = timeInt) %>%
  mutate(year = age + birthdate) %>%
  dplyr::select(id, enter, exit, event, everything()) %>%
  left_join(logrye, by = "year")

head(Scania_PersonPeriod_Test)


Scania_PersonPeriod_Train %>%
  group_by(exit) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  ggplot(aes(x = exit, y = log(-log(1-hazard)))) +
    geom_point() +
    geom_smooth()

Gompertz_Model_Baseline <- glm(
  formula = event ~ exit,
  family = binomial(link = "cloglog"),
  data = Scania_PersonPeriod_Train)

summary(Gompertz_Model_Baseline)

summ(Gompertz_Model_Baseline, exp = T) #exp = T means that we want exponentiated estimates


Scania_PersonPeriod_Train %>%
  group_by(exit, sex) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  ggplot(aes(x = exit, 
             y = log(-log(1-hazard)),
             col = sex)) +
  geom_point() +
  geom_smooth()


Scania_PersonPeriod_Train %>%
  group_by(exit, ses) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  ggplot(aes(x = exit, 
             y = log(-log(1-hazard)),
             col = ses)) +
  geom_point() +
  geom_smooth()


Scania_PersonPeriod_Train %>%
  group_by(exit, immigrant) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  ggplot(aes(x = exit, 
             y = log(-log(1-hazard)),
             col = immigrant)) +
  geom_point() +
  geom_smooth()


Scania_PersonPeriod_Train %>%
  mutate(foodprices = cut_interval(foodprices, 10, labels = F)) %>%
  group_by(foodprices) %>%
  summarise(event = sum(event),
            total = n()) %>%
  mutate(hazard = event/total) %>%
  ggplot(aes(x = foodprices, y = log(-log(1-hazard)))) +
  geom_point() +
  geom_smooth(method = "lm")


Gompertz_Model_Full <- glm(
  formula = event ~ exit + sex + ses + immigrant + foodprices,
  family = binomial(link = "cloglog"),
  data = Scania_PersonPeriod_Train)

summary(Gompertz_Model_Full)

summ(Gompertz_Model_Full, exp = T)

summ(Gompertz_Model_Full, exp = T, scale = T)


plot_summs(Gompertz_Model_Full, exp = T, scale = T)







