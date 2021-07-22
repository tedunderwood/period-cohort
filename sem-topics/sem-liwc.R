
#Load in data
liwc <- read_tsv("./liwc/liwc_w_meta.tsv")

#Grouping by LIWC categories
liwc_categories <- tibble(liwc_var = names(liwc)[9:89]) %>%
  mutate(category = c("word count", rep("summary language", 7),
                      rep("linguistic dimensions", 15),
                      rep("other grammar", 6), 
                      rep("affective processes", 6),
                      rep("social processes", 5), 
                      rep("cognitive processes", 7),
                      rep("perceptual processes", 4),
                      rep("biological processes", 5),
                      rep("drives", 6),
                      rep("time orientations", 4),
                      rep("relativity", 4),
                      rep("personal concerns", 6),
                      rep("informal language", 5)))

for (varnum in seq(9, 89)){
  liwc[ , varnum] <- scale(liwc[ , varnum])[, 1]
}

#Selecting three works by each author who has
#at least three works
three_works <- liwc %>% 
  filter(firstpub < 1990 & firstpub > 1889) %>%
  group_by(hathi_author) %>% 
  filter(n() >= 3) %>%
  arrange(hathi_author, firstpub) 

#Arranging data
d <- three_works %>% select(hathi_author, firstpub, birthyear, WC:filler) %>%
  pivot_longer(WC:filler, names_to = "variable", values_to = "y") %>% 
  group_by(hathi_author, variable) %>%
  mutate(y1 = lag(y, 2),
         y2 = lag(y, 1),
         y3 = y,
         d1 = lag(firstpub, 2),
         d2 = lag(firstpub, 1),
         d3 = firstpub,
         time12 = d2 - d1, 
         time23 = d3 - d2) %>%
  select(-c(y, firstpub)) %>%
  filter(!is.na(y1), !is.na(y2)) %>%
  drop_na() %>% 
  group_by(variable) %>%
  mutate(age = d1 - birthyear) %>% 
  group_by(hathi_author, variable) %>% mutate(triples = n()) 


#SEM syntax
aum1_mod <- "

  # structural part 
    y3 ~ alpha*1 + rho*y2 
    y2 ~ alpha*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (constrained)
    U ~~ 0*y1
    
  # constraint var(U) = 0
    U ~~ 0*U
    
"

## AUM2
aum2_mod <- "

  # structural part 
    y3 ~ alpha3*1 + rho*y2 
    y2 ~ alpha2*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (constrained)
    U ~~ 0*y1
    
  # constraint var(U) = 0
    U ~~ 0*U
    
"

## AUM3
aum3_mod <- "

  # structural part 
    y3 ~ alpha*1 + rho*y2 
    y2 ~ alpha*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance
    U ~~ tau*y1
    
  # var(U)
    U ~~ U
  
"

## AUM4
aum4_mod <- "

  # structural part 
    y3 ~ alpha3*1 + rho*y2 
    y2 ~ alpha2*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (tau)
    U ~~ tau*y1
    
  # var(U)
    U ~~ U
  
"

## SDM1
sdm1_mod <- "

  # structural part 
    y3 ~ alpha*1 + rho*y2 
    y2 ~ alpha*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (tau)
    U ~~ tau*y1
    
  # var(U)
    U ~~ U
    
  # constraint
    rho == 0
  
"

## SDM2
sdm2_mod <- "

  # structural part 
    y3 ~ alpha3*1 + rho*y2 
    y2 ~ alpha2*1 + rho*y1
    
  # measurement part
    U =~ 1*y3 + 1*y2

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    
  # covariance (tau)
    U ~~ tau*y1
    
  # var(U)
    U ~~ U
    
  # constraint
    rho == 0
  
"

## CFA
cfa_mod <- "

  # measurement part
    U =~ 1*y3 + 1*y2 + 1*y1

  # variances
    y3 ~~ v*y3
    y2 ~~ v*y2
    y1 ~~ v*y1
    
  # var(U)
    U ~~ U

  # intercepts
    y1 ~ int*1
    y2 ~ int*1
    y3 ~ int*1

"

# fitting functions
aum1_fit <- function (x) {
  sem(aum1_mod,
      check.post = TRUE,
      se = "none",
      data = x)
}

aum2_fit <- function (x) {
  sem(aum2_mod,
      check.post = TRUE,
      se = "none",
      data = x)
}

aum3_fit <- function (x) {
  sem(aum3_mod,
      check.post = TRUE,
      se = "none",
      data = x)
}

aum4_fit <- function (x) {
  sem(aum4_mod,
      check.post = TRUE,
      se = "none",
      data = x)
}

sdm1_fit <- function (x) {
  sem(sdm1_mod,
      check.post = TRUE,
      se = "none",
      data = x)
}

sdm2_fit <- function (x) {
  sem(sdm2_mod,
      check.post = TRUE,
      se = "none",
      data = x)
}

cfa_fit <- function (x) {
  sem(cfa_mod,
      check.post = TRUE,
      se = "none",
      data = x)
}

p_aum1_fit <- possibly(aum1_fit, otherwise = "NOPE")
p_aum2_fit <- possibly(aum2_fit, otherwise = "NOPE")
p_aum3_fit <- possibly(aum3_fit, otherwise = "NOPE")
p_aum4_fit <- possibly(aum4_fit, otherwise = "NOPE")
p_sdm1_fit <- possibly(sdm1_fit, otherwise = "NOPE")
p_sdm2_fit <- possibly(sdm2_fit, otherwise = "NOPE")
p_cfa_fit <- possibly(cfa_fit, otherwise = "NOPE")

#Running the models
d <- d %>% 
  group_by(variable) %>% nest() %>%
  mutate(AUM1 = map(data, p_aum1_fit),
         AUM2 = map(data, p_aum2_fit),
         AUM3 = map(data, p_aum3_fit),
         AUM4 = map(data, p_aum4_fit),
         SDM1 = map(data, p_sdm1_fit),
         SDM2 = map(data, p_sdm2_fit),
         CFA = map(data, p_cfa_fit))


p_glance = possibly(glance, otherwise = "NOPE")

#Extract model symmary information
results <- d %>%
  pivot_longer(AUM1:CFA, names_to = "mod_spec", values_to = "mod_object") %>% 
  mutate(glanced = map(mod_object, p_glance)) %>% 
  unnest(glanced) %>% 
  select(variable, mod_spec, BIC, chisq, npar, rmsea, converged, nobs) %>% 
  mutate(type = case_when(
    mod_spec == "SDM1" | mod_spec == "SDM2" ~ "SDM",
    mod_spec == "CFA" ~ "CFA",
    TRUE ~ "AUM")) %>% 
  ungroup()

#Determine best fit of each, then overall best fit
winners <- results %>%
  filter(type != "CFA") %>% 
  group_by(variable, type) %>% 
  filter(BIC == min(BIC)) %>% 
  group_by(variable) %>% 
  mutate(BIC_diff = BIC - max(BIC)) %>% 
  filter(BIC == min(BIC)) %>% 
  mutate(verdict = if_else(BIC_diff >= -2, "Inconclusive", type)) %>% 
  ungroup()

table(winners$verdict)


winners %>%
  left_join(liwc_categories, by = c("variable"="liwc_var")) %>%
  group_by(category, verdict) %>% summarise(count = n()) %>%
  mutate(pct = count/sum(count), total_count = sum(count),
         category = paste(category, " (", total_count, ")", sep = "")) %>%
  ggplot(aes(x = category, y = pct, fill = verdict)) + 
  geom_bar(stat = "identity", color = "black") + 
  coord_flip() + 
  theme_classic() + 
  theme(legend.position = "top") + 
  labs(x = "", y = "Count", fill = "Models:") + 
  scale_fill_brewer(type = "qual", palette = "Set2")

  
  
  