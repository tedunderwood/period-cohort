Finding the best models through a grid search
================

Needed packages.

``` r
library(tidyverse)
library(rstatix)
```

Bring in the data. Limiting to 1890-1989. Not requiring multiple works,
but requiring us\_national. Once this chunk is done executing, there
should be 5,573 books in d.

``` r
d <- read_tsv("../liwc/liwc_w_meta.tsv")

d <- d %>%
  filter(firstpub < 1990 & firstpub > 1889) %>% 
  filter(us_national == TRUE) %>%
  rename(author = hathi_author)
```

Now we convert the numeric columns for birthyear and firstpub into
factors. Six different factors are created for each numeric variable,
using different binwidths.

To avoid discarding data from the tails, we use a function that creates
ad-hoc breaks and requires each new factor level to have at least 20
books. So the breaks toward the end, and especially at the beginning of
birthyear, may be more widely spaced than in the middle of the timeline.

``` r
discretize <- function(numericvar, width) {
  minval = min(numericvar) - 1
  maxval = max(numericvar) + 1
  breaks = c()
  for (i in seq(minval, maxval, width)) {
    breaks <- c(breaks, i)
  }
  breaks <- c(breaks, maxval)
  
  neededbreaks = c(breaks[1])
  previous = breaks[1]
  for (idx in seq(2, length(breaks) - 1)) {
    prevct = sum(numericvar >= previous & numericvar < breaks[idx])
    nextct = sum(numericvar >= breaks[idx] & numericvar < breaks[idx + 1])
    if (prevct > 20 & nextct > 20 & breaks[idx] != 1989 ){
      previous = breaks[idx]
      neededbreaks <- c(neededbreaks, breaks[idx])
    }
  }
breaks <- c(neededbreaks, maxval)

# The algorithm above ensures that breaks are less tightly packed in tails of 
# distributions. But when the list of breaks gets quite short it can be
# tricky to arrange the breaks to compromise between even chronological
# spacing and (roughly!) even numbers of volumes in factor levels.
# The following cases are therefore addressed manually:

if (minval == 1889 & width == 24) breaks = c(1889, 1914, 1938, 1961, 1990)
if (minval == 1808 & width == 24) breaks = c(1808, 1854, 1878, 1902, 1926, 1990)
if (minval == 1808 & width == 20) breaks = c(1808, 1850, 1870, 1890, 1910, 1930, 1990)
if (minval == 1889 & width == 16) breaks = c(1889, 1905, 1921, 1937, 1953, 1968, 1982, 1990)

print(breaks)
result <- cut(numericvar, breaks = as.integer(breaks), labels = as.character(breaks[1: length(breaks) -1]))
result
}

for (width in seq(4, 24, 4)){
  bylabel = paste('by_', as.character(width), sep = '')
  fplabel = paste('fp_', as.character(width), sep = '')
  d[bylabel] <- discretize(as.integer(d$birthyear), width)
  d[fplabel] <- discretize(as.integer(d$firstpub), width)
}
##  [1] 1808 1836 1840 1844 1848 1852 1856 1860 1864 1868 1872 1876 1880 1884
## [15] 1888 1892 1896 1900 1904 1908 1912 1916 1920 1924 1928 1932 1936 1940
## [29] 1944 1948 1952 1963
##  [1] 1889 1893 1897 1901 1905 1909 1913 1917 1921 1925 1929 1933 1937 1941
## [15] 1945 1949 1953 1957 1961 1965 1969 1973 1977 1981 1985 1990
##  [1] 1808 1832 1840 1848 1856 1864 1872 1880 1888 1896 1904 1912 1920 1928
## [15] 1936 1944 1952 1963
##  [1] 1889 1897 1905 1913 1921 1929 1937 1945 1953 1961 1969 1977 1985 1990
##  [1] 1808 1832 1844 1856 1868 1880 1892 1904 1916 1928 1940 1952 1963
##  [1] 1889 1901 1913 1925 1937 1949 1961 1973 1985 1990
##  [1] 1808 1840 1856 1872 1888 1904 1920 1936 1952 1963
## [1] 1889 1905 1921 1937 1953 1968 1982 1990
## [1] 1808 1850 1870 1890 1910 1930 1990
## [1] 1889 1909 1929 1949 1969 1990
## [1] 1808 1854 1878 1902 1926 1990
## [1] 1889 1914 1938 1961 1990
```

SCALE LIWC COLUMNS WITH BIZARRELY UNEVEN SCALES

All the dependent variables are converted to zscores, so they will have
comparable scales.

``` r
for (varnum in seq(10, 89)){
  d[ , varnum] <- scale(d[ , varnum])[, 1]
}
```

Now we actually do a grid search to find the best binwdiths for firstpub
and birthyear. We select the model with highest overall r2, and in doing
that use fivefold cross-validation *on unseen authors*.

After selecting the best model we find out how variance is apportioned
across variables.

``` r
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

r_squared <- function(vals, preds) {
  1 - (sum((vals - preds)^2) / sum((vals - mean(preds))^2))
}

cross_validate <- function(modelstring, data, depvar) {
  
  # Testing a model specified by modelstring.
  # We test it on out-of-sample authors.
  
  authors <- unique(d$author)
  authors <- sample(authors)
  tenauthsets <- chunk2(authors, 10)
  rsquaredvals <- c()
  
  for (authset in tenauthsets){
    dtest <- d[d$author %in% authset,  ]
    dtrain <- d[!d$author %in% authset,  ]
    model <- lm(as.formula(modelstring), data = dtrain)
    oos_predictions <- predict(model, newdata = dtest)
    r2 <- r_squared(dtest[depvar], oos_predictions)
    rsquaredvals <- c(rsquaredvals, r2)
  }
  mean(rsquaredvals)
}

partial_r2s_oos <- function(modelstring, data, depvar, cohortvar, periodvar) {
  
  # See how much difference variables make out of sample.
  
  authors <- unique(d$author)
  authors <- sample(authors)
  tenauthsets <- chunk2(authors, 10)
  
  periodr2s <- c()
  cohortr2s <- c()
  
  for (authset in tenauthsets){
    dtest <- d[d$author %in% authset,  ]
    dtrain <- d[!d$author %in% authset,  ]
    model <- lm(as.formula(modelstring), data = dtrain)
    oos_predictions <- predict(model, newdata = dtest)
    baser2 <- r_squared(dtest[depvar], oos_predictions)
    
    cohortdata <- data.frame(dtest)
    cohortdata[cohortvar] <- sample(cohortdata[cohortvar][[1]])
    oos_predictions <- predict(model, newdata = cohortdata)
    losecohortr2 <- baser2 - r_squared(cohortdata[depvar], oos_predictions)
    # cat(baser2, ' ', r_squared(cohortdata[depvar], oos_predictions), '\n')
    
    perioddata <- data.frame(dtest)
    perioddata[periodvar] <- sample(perioddata[periodvar][[1]])
    oos_predictions <- predict(model, newdata = perioddata)
    loseperiodr2 <- baser2 - r_squared(cohortdata[depvar], oos_predictions)
    
    periodr2s <- c(periodr2s, loseperiodr2)
    cohortr2s <- c(cohortr2s, losecohortr2)
    
    # oldcor = sum(dtest[cohortvar] == cohortdata[cohortvar])
    # cat(oldcor, '\n')
  }
  returnvalue <- list(periodr2 = mean(periodr2s), cohortr2  = mean(cohortr2s))
  returnvalue
}

periodnames = c('fp_4', 'fp_8', 'fp_12', 'fp_16', 'fp_20', 'fp_24')
cohortnames = c('by_4', 'by_8', 'by_12', 'by_16', 'by_20', 'by_24')

varname.col <- c()
cmse.col <- c()
pmse.col <- c()
r2.col <- c()
delta.col <- c()
adjdelta.col <- c()
bywidth.col <- c()
fpwidth.col <- c()
bydf.col <- c()
fpdf.col <- c()
pmse_oos.col <- c()
cmse_oos.col <- c()
delta_oos.col <- c()
r2_oos.col <- c()

# Iterate across dependent variables:

for (varnum in seq(10, 89)){
  depvar <- colnames(d)[varnum]
  if (depvar == 'function') next  # that word breaks my as.formula!     
  bestr2 <- -0.5
  bestmodel <- 'the unknown ideal'
  bestbywidth <- 50
  bestfpwidth <- 50
  
  # For each dependent variables, do a grid search to find the
  # granularity of factors that produce the highest r2.
  # This is assessed using the cross_validate function above,
  # which tests on out-of-sample authors.
  
  for (bywidth in seq(4, 24, 4)) {
    for (fpwidth in seq(4, 24, 4)) {
      
      bylabel = paste('by_', as.character(bywidth), sep = '')
      fplabel = paste('fp_', as.character(fpwidth), sep = '')
      modelstring <- paste(depvar, '~', 'authorage + I(authorage^2) + ', bylabel, '+', fplabel)
      thisr2 <- cross_validate(modelstring, d, depvar)
      
      if (thisr2 > bestr2) {
        bestr2 <- thisr2
        bestmodel <- modelstring
        bestbywidth <- bylabel
        bestfpwidth <- fplabel
      }
    }
  }
  
  # We're done finding the best model specification.
  # Now we train it on all data and examine sums of squares.
  
  r2_oos = bestr2
  variances_oos <- partial_r2s_oos(bestmodel, d, depvar, bestbywidth, bestfpwidth)
  
  cmse_oos = variances_oos$cohortr2
  pmse_oos = variances_oos$periodr2
  
  if (cmse_oos < 0) cmse_oos = .00000001    # debatable choice here
  if (pmse_oos < 0) pmse_oos = .00000001
  
  model <- lm(as.formula(bestmodel), data = d)
  thisr2 <- summary(model)$r.squared
  at <- anova_test(model, detailed = TRUE)
  
  cmse = 0
  pmse = 0
  cdf = 0
  pdf = 0
  
  for (rownum in seq(3, dim(at)[1])) {
    variable_name = at[rownum, 1]
    if (variable_name %in% cohortnames){
      cmse = cmse + at[rownum, 2]
      cdf = cdf + at[rownum, 4]
    }
    else if (variable_name %in% periodnames) {
      pmse = pmse + at[rownum, 2]
      pdf = pdf + at[rownum, 4]
    }
    else {
      print('Model error.')
    }
  }
  
  delta <- cmse / (cmse + pmse)
  delta_oos <- cmse_oos / (cmse_oos + pmse_oos)
  adjdelta = (cmse/cdf) / ((cmse/cdf) + (pmse / pdf))
  cat(depvar, bestbywidth, bestfpwidth, delta, adjdelta, thisr2, '\n')
  varname.col <- c(varname.col, depvar)
  cmse.col <- c(cmse.col, cmse)
  pmse.col <- c(pmse.col, pmse)
  pmse_oos.col <- c(pmse_oos.col, pmse_oos)
  cmse_oos.col <- c(cmse_oos.col, cmse_oos)
  delta_oos.col <- c(delta_oos.col, delta_oos)
  r2.col <- c(r2.col, thisr2)
  r2_oos.col <- c(r2_oos.col, r2_oos)
  delta.col <- c(delta.col, delta)
  adjdelta.col <- c(adjdelta.col, adjdelta)
  bywidth.col <- c(bywidth.col, bestbywidth)
  fpwidth.col <- c(fpwidth.col, bestfpwidth)
  bydf.col <- c(bydf.col, cdf)
  fpdf.col <- c(fpdf.col, pdf)
}
## Analytic by_20 fp_16 0.764659 0.795876 0.03284901 
## Clout by_24 fp_20 0.2792976 0.2792976 0.0195134 
## Authentic by_24 fp_24 0.06147027 0.04682225 0.0270067 
## Tone by_20 fp_8 0.5371545 0.7358212 0.1413012 
## WPS by_20 fp_16 0.5564329 0.6008527 0.01042196 
## Sixltr by_20 fp_4 0.3690147 0.7373363 0.03755357 
## Dic by_20 fp_8 0.4613898 0.6727655 0.03607194 
## was_function by_24 fp_8 0.3663573 0.6343064 0.02543075 
## pronoun by_24 fp_8 0.4528954 0.7129252 0.01866339 
## ppron by_24 fp_16 0.4645593 0.5654876 0.01780649 
## i by_12 fp_24 0.7898758 0.5062232 0.04836715 
## we by_20 fp_16 0.3970308 0.441388 0.02076364 
## you by_20 fp_20 0.811188 0.7746234 0.02231925 
## shehe by_24 fp_8 0.14862 0.3436989 0.02445865 
## they by_16 fp_24 0.8394206 0.6621957 0.02275039 
## ipron by_20 fp_24 0.9310704 0.8901646 0.03427086 
## article by_16 fp_16 0.9433552 0.9258732 0.02894608 
## prep by_20 fp_16 0.5958953 0.6389276 0.06782346 
## auxverb by_12 fp_4 0.4938075 0.6803514 0.08018414 
## adverb by_20 fp_8 0.2662042 0.4654307 0.05704687 
## conj by_20 fp_20 0.8212386 0.7861072 0.03537081 
## negate by_12 fp_20 0.7480399 0.5191368 0.04346786 
## verb by_20 fp_24 0.6917732 0.5738549 0.06912026 
## adj by_24 fp_16 0.7587389 0.8250932 0.0361977 
## compare by_24 fp_16 0.7410509 0.8110586 0.01501431 
## interrog by_20 fp_16 0.7300701 0.7644619 0.04779582 
## number by_16 fp_24 0.821529 0.6331866 0.02676859 
## quant by_20 fp_12 0.4435272 0.5604884 0.04761438 
## affect by_20 fp_4 0.2368163 0.5983036 0.1257349 
## posemo by_12 fp_16 0.8780591 0.7970635 0.1845748 
## negemo by_20 fp_20 0.4741982 0.4191065 0.02744558 
## anx by_20 fp_20 0.4064345 0.3539163 0.01963762 
## anger by_20 fp_8 0.5778579 0.7666435 0.0563979 
## sad by_20 fp_4 0.1158634 0.3861359 0.1121442 
## social by_24 fp_20 0.200333 0.200333 0.012627 
## family by_20 fp_20 0.6000416 0.5454975 0.01686534 
## friend by_20 fp_16 0.3261272 0.3673894 0.1180271 
## female by_20 fp_16 0.4109245 0.4556609 0.02205012 
## male by_20 fp_24 0.7646045 0.6608907 0.02009567 
## cogproc by_24 fp_24 0.7019095 0.638469 0.03837287 
## insight by_20 fp_24 0.8525951 0.7763074 0.01888291 
## cause by_20 fp_20 0.684138 0.6340685 0.05834792 
## discrep by_24 fp_24 0.5008953 0.4294487 0.0393135 
## tentat by_24 fp_24 0.3517581 0.2892556 0.04469576 
## certain by_12 fp_16 0.9364839 0.8894074 0.04960056 
## differ by_24 fp_8 0.3111936 0.5754365 0.0420562 
## percept by_20 fp_8 0.6087225 0.7887512 0.1488818 
## see by_20 fp_24 0.8343484 0.7513712 0.05674099 
## hear by_20 fp_12 0.672359 0.7665402 0.1174244 
## feel by_20 fp_12 0.5892665 0.6965532 0.1192195 
## bio by_20 fp_12 0.5733416 0.682547 0.1519411 
## body by_20 fp_12 0.4973265 0.6128504 0.129908 
## health by_24 fp_16 0.2503989 0.333806 0.02990671 
## sexual by_16 fp_4 0.1080393 0.266527 0.1193694 
## ingest by_24 fp_20 0.5516243 0.5516243 0.0864737 
## drives by_20 fp_16 0.3874002 0.4314514 0.04170832 
## affiliation by_20 fp_20 0.6756517 0.6249744 0.0372291 
## achieve by_20 fp_4 0.2062604 0.5550262 0.03904032 
## power by_20 fp_16 0.3135485 0.3540555 0.01847998 
## reward by_24 fp_20 0.2895897 0.2895897 0.02494999 
## risk by_24 fp_20 0.4061171 0.4061171 0.04721931 
## focuspast by_20 fp_16 0.6394604 0.6803422 0.06976837 
## focuspresent by_20 fp_20 0.8401065 0.8078154 0.03488815 
## focusfuture by_12 fp_12 0.6198644 0.5425267 0.07331195 
## relativ by_20 fp_20 0.4382149 0.3842486 0.04709023 
## motion by_20 fp_20 0.6135716 0.5595182 0.04424091 
## space by_20 fp_24 0.54825 0.4213526 0.07654873 
## time by_16 fp_16 0.6556603 0.5881523 0.01777313 
## work by_20 fp_8 0.5449647 0.7418899 0.02971034 
## leisure by_20 fp_8 0.150461 0.2982759 0.07979021 
## home by_24 fp_12 0.406328 0.5778566 0.04309691 
## money by_12 fp_20 0.6838588 0.4402769 0.0273939 
## relig by_24 fp_16 0.3313876 0.4264253 0.01404888 
## death by_20 fp_20 0.6864327 0.6365336 0.01962599 
## informal by_20 fp_12 0.4501533 0.5670812 0.02705442 
## swear by_20 fp_16 0.4782721 0.5238208 0.2145503 
## netspeak by_8 fp_12 0.7315711 0.5767537 0.02747002 
## assent by_16 fp_8 0.5847278 0.6786724 0.02618648 
## nonflu by_20 fp_24 0.6204515 0.49516 0.02993087 
## filler by_20 fp_24 0.5778458 0.4509361 0.07932108

res.df <- data.frame(depvar = varname.col, cmse = cmse.col, pmse = pmse.col, totalr2 = r2.col, delta = delta.col, adjdelta = adjdelta.col, bywidth = bywidth.col, fpwidth = fpwidth.col, bydf = bydf.col, fpdf = fpdf.col, pmse_oos = pmse_oos.col, cmse_oos = cmse_oos.col, delta_oos = delta_oos.col, r2_oos = r2_oos.col)
write.csv(res.df, file = 'crossvalidated_delta_knit.tsv', quote = FALSE, row.names = FALSE)
```

RESULTS:

    ## Mean delta is  0.53977
    ## If we adjust for df it is  0.57346
    ## The average r2 is  0.05315
    ## The weighted average, sum(cmse.col) / (sum(cmse.col) + sum(pmse.col)):
    ## 0.53426
    ## 
    ## Mean delta measured oos is  0.46374
    ## but as a weighted average it is: 0.48683

<img src="gridsearch_oos_files/figure-gfm/unnamed-chunk-7-1.svg" width="100%" style="display: block; margin: auto;" />
