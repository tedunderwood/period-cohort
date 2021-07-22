
#Load in data
liwc <- read_tsv("~/Dropbox/period-cohort/liwc/liwc_w_meta.tsv")

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

#Selecting three works by each author who has
#at least three works
three_works <- liwc %>% 
  mutate(WC = WC/10000, WPS = WPS/10000) %>%
  group_by(hathi_author) %>% 
  filter(n() >= 3) %>%
  slice_sample(n = 3) %>%
  arrange(hathi_author, firstpub) %>%
  mutate(wave = 1:3)

#Arranging data
d <- three_works %>% select(hathi_author, wave, WC:filler) %>%
  pivot_longer(WC:filler, names_to = "variable", values_to = "y") %>% 
  pivot_wider(names_from = wave, values_from = y, names_prefix = "y") %>%
  arrange(variable) %>%
  drop_na() %>% 
  group_by(variable) %>% 
  nest()

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
  select(variable, mod_spec, bic, chisq, npar, rmsea, converged, nobs) %>% 
  mutate(type = case_when(
    mod_spec == "SDM1" | mod_spec == "SDM2" ~ "SDM",
    mod_spec == "CFA" ~ "CFA",
    TRUE ~ "AUM")) %>% 
  ungroup()

#Determine best fit of each, then overall best fit
winners <- results %>%
  filter(type != "CFA") %>% 
  group_by(variable, type) %>% 
  filter(bic == min(bic)) %>% 
  group_by(variable) %>% 
  mutate(BIC_diff = bic - max(bic)) %>% 
  filter(bic == min(bic)) %>% 
  mutate(verdict = if_else(BIC_diff >= -2, "Inconclusive", type)) %>% 
  ungroup()



preferred_models <- vector(mode = "list", length = 1000)
for (i in 1:100) {
  #Extract three works from each author
  three_works <- liwc %>% 
    mutate(WC = WC/10000, WPS = WPS/10000) %>%
    group_by(hathi_author) %>% 
    filter(n() >= 3) %>%
    slice_sample(n = 3) %>%
    arrange(hathi_author, firstpub) %>%
    mutate(wave = 1:3)
  
  #Arranging data
  d <- three_works %>% select(hathi_author, wave, WC:filler) %>%
    pivot_longer(WC:filler, names_to = "variable", values_to = "y") %>% 
    pivot_wider(names_from = wave, values_from = y, names_prefix = "y") %>%
    arrange(variable) %>%
    drop_na() %>% 
    group_by(variable) %>% 
    nest()
  
  d <- d %>% 
    mutate(AUM1 = map(data, p_aum1_fit),
           AUM2 = map(data, p_aum2_fit),
           AUM3 = map(data, p_aum3_fit),
           AUM4 = map(data, p_aum4_fit),
           SDM1 = map(data, p_sdm1_fit),
           SDM2 = map(data, p_sdm2_fit))
  
  results <- d %>%
    pivot_longer(AUM1:SDM2, names_to = "mod_spec", values_to = "mod_object") %>% 
    mutate(glanced = map(mod_object, p_glance)) %>% 
    unnest(glanced) %>% 
    select(variable, mod_spec, bic, chisq, npar, rmsea, converged, nobs) %>% 
    mutate(type = case_when(
      mod_spec == "SDM1" | mod_spec == "SDM2" ~ "SDM",
      TRUE ~ "AUM")) %>% 
    ungroup()
  
  preferred_models[[i]] <- results %>%
    group_by(variable, type) %>% 
    filter(bic == min(bic)) %>% 
    group_by(variable) %>% 
    mutate(BIC_diff = bic - max(bic)) %>% 
    filter(bic == min(bic)) %>% 
    mutate(verdict = if_else(BIC_diff >= -2, "Inconclusive", type)) %>% 
    ungroup() %>%
    mutate(iter = i)
  
  print(i)
  
}



bind_rows(preferred_models) %>%
  arrange(variable, mod_spec) %>%
  group_by(variable, verdict) %>%
  summarise(n = n()) %>%
  left_join(liwc_categories, by = c("variable" = "liwc_var")) %>%
  group_by(variable) %>% 
  ggplot(aes(x = variable, y = n, fill = verdict)) + 
  geom_bar(stat = "identity", color = "black") + 
  coord_flip() +
  scale_fill_brewer(palette = "Set2") + 
  facet_wrap(~category, scales = "free") +
  labs(fill = "Preferred\nModel",
       y = "Percent",
       x = "")

preferred_model_df <- bind_rows(preferred_models) %>%
  arrange(variable, mod_spec) 

save(preferred_model_df, file = "~/Dropbox/period-cohort/analysis/sem-results.Rdata")

