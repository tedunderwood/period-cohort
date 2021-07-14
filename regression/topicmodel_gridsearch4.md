Estimating delta through a grid search on the topic model
================

Needed packages.

``` r
library(tidyverse)
library(rstatix)
```

Bring in the data. Limiting to 1890-1989. Not requiring multiple works,
but requiring us\_national. Once this chunk is done executing, there
should be 5,572 books in
d.

``` r
book <- read.csv('/Users/tunder/Dropbox/python/cohort/regression/bookleveltopicdata.tsv', sep = '\t')
book[['us_national']] <- as.logical(book$us_national)
d <- book %>%
  filter(firstpub < 1990 & firstpub > 1889) %>% 
  filter(us_national == TRUE) %>%
  rename(author = hathi_author, authorage = age)
```

All the dependent variables are converted to zscores, so they will have
comparable scales.

``` r
for (varnum in seq(1, 200)){
  d[ , varnum] <- scale(d[ , varnum])[, 1]
}
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

Now we actually do a grid search to find the best binwdiths for firstpub
and birthyear. We select the model with highest overall r2, and in doing
that use tenfold cross-validation *on unseen authors*.

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
  correcty <- c()
  ypredictions <- c()
  
  for (authset in tenauthsets){
    dtest <- d[d$author %in% authset,  ]
    dtrain <- d[!d$author %in% authset,  ]
    model <- lm(as.formula(modelstring), data = dtrain)
    oos_predictions <- predict(model, newdata = dtest)
    ypredictions <- c(ypredictions, oos_predictions)
    correcty <- c(correcty, dtest[[depvar]])
    # r2 <- r_squared(dtest[depvar], oos_predictions)
    # rsquaredvals <- c(rsquaredvals, r2)
  }
  accurater2 <- r_squared(correcty, ypredictions)
  # badr2 <-mean(rsquaredvals)
  # print(abs(badr2-accurater2) / accurater2)
  accurater2
}

partial_r2s_oos <- function(modelstring, data, depvar, cohortvar, periodvar) {
  
  # See how much difference variables make out of sample.
  # Our strategy is to make predictions out of sample:
  # 1) using actual predictive variables
  # 2) with the cohort variable permuted
  # 3) with the period variable permuted
  #
  # The differences between 1 and 2 or 3 reflect lost predictive power.
  
  authors <- unique(d$author)
  authors <- sample(authors)
  tenauthsets <- chunk2(authors, 10)
  
  baseprediction <- c()
  predictionsanscohort <- c()
  predictionsansperiod <- c()
  correctyvals <- c()
  
  for (authset in tenauthsets){
    dtest <- d[d$author %in% authset,  ]
    dtrain <- d[!d$author %in% authset,  ]
    model <- lm(as.formula(modelstring), data = dtrain)
    oos_predictions <- predict(model, newdata = dtest)
    baseprediction <- c(baseprediction, oos_predictions)
    correctyvals <- c(correctyvals, dtest[depvar][[1]])
    
    cohortdata <- data.frame(dtest)
    cohortdata[cohortvar] <- sample(cohortdata[cohortvar][[1]])
    oos_predictions <- predict(model, newdata = cohortdata)
    predictionsanscohort <- c(predictionsanscohort, oos_predictions)
    
    perioddata <- data.frame(dtest)
    perioddata[periodvar] <- sample(perioddata[periodvar][[1]])
    oos_predictions <- predict(model, newdata = perioddata)
    predictionsansperiod <- c(predictionsansperiod, oos_predictions)
  }
  #cat(correctyvals, '\n')
  # cat(baseprediction, '\n')
  baser2 <- r_squared(correctyvals, baseprediction)
  sanscohortr2 <- r_squared(correctyvals, predictionsanscohort)
  sansperiodr2 <- r_squared(correctyvals, predictionsansperiod)
  
  returnvalue <- list(periodr2 = baser2 - sansperiodr2, 
                      cohortr2  = baser2 - sanscohortr2,
                      baser2 = baser2)
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
age.col <- c()

# Iterate across dependent variables:

for (varnum in seq(1, 200)){
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
  
  variances_oos <- partial_r2s_oos(bestmodel, d, depvar, bestbywidth, bestfpwidth)
  
  cmse_oos = variances_oos$cohortr2
  pmse_oos = variances_oos$periodr2
  r2_oos = variances_oos$baser2
  
  if (cmse_oos < 0) cmse_oos = .00000001    # debatable choice here
  if (pmse_oos < 0) pmse_oos = .00000001
  
  model <- lm(as.formula(bestmodel), data = d)
  thisr2 <- summary(model)$r.squared
  at <- anova_test(model, detailed = TRUE)
  
  cmse = 0
  pmse = 0
  cdf = 0
  pdf = 0
  
  agemse = at[1, 2] + at[2, 2]
  
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
  age.col <- c(age.col, agemse)
}
## t0 by_20 fp_4 0.04851111 0.1966098 0.04767901 
## t1 by_20 fp_8 0.5610876 0.7541827 0.01885585 
## t2 by_16 fp_16 0.6580685 0.5907378 0.01027996 
## t3 by_20 fp_8 0.2419042 0.4336931 0.1351399 
## t4 by_12 fp_16 0.7965398 0.6810651 0.02772009 
## t5 by_12 fp_4 0.6930279 0.8312444 0.02242733 
## t6 by_20 fp_4 0.490619 0.8221654 0.05882429 
## t7 by_24 fp_8 0.1836041 0.4028737 0.1459044 
## t8 by_12 fp_16 0.8858069 0.8088372 0.01328326 
## t9 by_20 fp_20 0.6328506 0.5796464 0.02795295 
## t10 by_20 fp_12 0.4487588 0.5656971 0.01605452 
## t11 by_20 fp_16 0.3684162 0.4117596 0.01087445 
## t12 by_8 fp_4 0.6882432 0.7680591 0.1218757 
## t13 by_20 fp_8 0.4482599 0.6610025 0.1963967 
## t14 by_24 fp_24 0.9553213 0.9413026 0.02094515 
## t15 by_24 fp_24 0.8717756 0.8360419 0.01183742 
## t16 by_24 fp_20 0.123179 0.123179 0.02345371 
## t17 by_12 fp_16 0.838479 0.739008 0.03034723 
## t18 by_24 fp_8 0.04411819 0.121623 0.06828292 
## t19 by_24 fp_12 0.4455424 0.6164363 0.1051871 
## t20 by_20 fp_16 0.5533999 0.5979039 0.04450898 
## t21 by_12 fp_12 0.6112996 0.5335307 0.04343828 
## t22 by_16 fp_24 0.4658049 0.2464152 0.04970611 
## t23 by_16 fp_24 0.9441672 0.8637877 0.01893524 
## t24 by_24 fp_8 0.2387142 0.4847224 0.02361904 
## t25 by_20 fp_24 0.5026422 0.3774803 0.01076076 
## t26 by_20 fp_16 0.4542776 0.4997299 0.02705128 
## t27 by_20 fp_12 0.4029823 0.5192278 0.05482702 
## t28 by_12 fp_20 0.693734 0.4516599 0.04769548 
## t29 by_24 fp_24 0.6979547 0.6341113 0.01486523 
## t30 by_24 fp_24 0.5820607 0.5108874 0.004804202 
## t31 by_24 fp_16 0.8132341 0.8672233 0.01072806 
## t32 by_20 fp_24 0.8250926 0.7389297 0.01688626 
## t33 by_24 fp_12 0.3856793 0.5566646 0.0415734 
## t34 by_12 fp_24 0.8054395 0.5303034 0.1292351 
## t35 by_16 fp_8 0.1385654 0.1943808 0.2448881 
## t36 by_24 fp_20 0.2388103 0.2388103 0.01627368 
## t37 by_12 fp_12 0.6819472 0.6092787 0.007752665 
## t38 by_12 fp_12 0.8926443 0.8580986 0.118325 
## t39 by_20 fp_24 0.8070207 0.71503 0.035388 
## t40 by_12 fp_8 0.4244928 0.4458766 0.08048956 
## t41 by_20 fp_20 0.5323728 0.4766492 0.01617686 
## t42 by_16 fp_12 0.5749369 0.5749369 0.01315807 
## t43 by_4 fp_4 0.7892892 0.7497919 0.1158185 
## t44 by_12 fp_16 0.4093032 0.2742864 0.01887089 
## t45 by_12 fp_8 0.6189165 0.6392161 0.009959986 
## t46 by_16 fp_4 0.1178573 0.2861274 0.165044 
## t47 by_20 fp_24 0.5076946 0.3822416 0.0236695 
## t48 by_8 fp_16 0.8134392 0.6205032 0.07587401 
## t49 by_24 fp_8 0.2882683 0.5485473 0.01311753 
## t50 by_8 fp_4 0.4884658 0.5888759 0.1734154 
## t51 by_4 fp_16 0.8451501 0.5218911 0.09405237 
## t52 by_20 fp_24 0.8830605 0.8191963 0.0184368 
## t53 by_20 fp_20 0.5378849 0.4821793 0.07784039 
## t54 by_12 fp_12 0.3842954 0.3122096 0.03863586 
## t55 by_16 fp_4 0.7471491 0.8986284 0.1354064 
## t56 by_20 fp_16 0.4380413 0.4833078 0.009061902 
## t57 by_24 fp_24 0.7210375 0.6596942 0.04243867 
## t58 by_20 fp_20 0.5285357 0.4728077 0.01918478 
## t59 by_20 fp_16 0.1747663 0.2026367 0.04880672 
## t60 by_24 fp_24 0.5769398 0.5056351 0.002718322 
## t61 by_20 fp_12 0.360759 0.4745052 0.02688834 
## t62 by_24 fp_24 0.6601287 0.5929528 0.03380776 
## t63 by_20 fp_8 0.2518794 0.4469146 0.01384207 
## t64 by_20 fp_4 0.2454458 0.6095843 0.04523657 
## t65 by_20 fp_8 0.5675768 0.759043 0.1333161 
## t66 by_12 fp_20 0.8434251 0.6620259 0.02290221 
## t67 by_20 fp_12 0.3199351 0.4294571 0.06031111 
## t68 by_24 fp_24 0.7700065 0.7151778 0.001148831 
## t69 by_20 fp_20 0.7614446 0.7185888 0.009023558 
## t70 by_20 fp_20 0.7368217 0.6913356 0.1444666 
## t71 by_24 fp_16 0.7382551 0.8088236 0.02887926 
## t72 by_24 fp_24 0.5430498 0.4712679 0.0214051 
## t73 by_24 fp_8 0.07154877 0.187776 0.01678013 
## t74 by_20 fp_12 0.2970774 0.4034163 0.015021 
## t75 by_24 fp_16 0.6633433 0.7471924 0.02334792 
## t76 by_20 fp_20 0.6927155 0.6432967 0.04626555 
## t77 by_24 fp_20 0.3701185 0.3701185 0.005010808 
## t78 by_12 fp_12 0.8860603 0.8497526 0.1861227 
## t79 by_24 fp_20 0.7671758 0.7671758 0.01433085 
## t80 by_4 fp_16 0.8109439 0.4617542 0.04670697 
## t81 by_24 fp_24 0.16636 0.1301843 0.00238839 
## t82 by_20 fp_12 0.5859818 0.6936805 0.004338471 
## t83 by_20 fp_24 0.8754485 0.8083294 0.1143114 
## t84 by_12 fp_12 0.5387848 0.4593395 0.0404373 
## t85 by_4 fp_20 0.9164653 0.5939598 0.03360796 
## t86 by_20 fp_16 0.4655518 0.5110757 0.01803582 
## t87 by_8 fp_20 0.9884726 0.9554315 0.2335646 
## t88 by_20 fp_24 0.3654115 0.2567789 0.05739073 
## t89 by_4 fp_16 0.8805704 0.5958985 0.1540655 
## t90 by_24 fp_12 0.390419 0.5615847 0.02303312 
## t91 by_24 fp_20 0.7122989 0.7122989 0.01079276 
## t92 by_24 fp_24 0.9151492 0.8899774 0.01015285 
## t93 by_24 fp_8 0.0803428 0.2076604 0.01863703 
## t94 by_20 fp_24 0.6510857 0.5282173 0.007134845 
## t95 by_16 fp_20 0.4650078 0.3029382 0.007958039 
## t96 by_24 fp_20 0.125191 0.125191 0.004141965 
## t97 by_24 fp_16 0.2471552 0.3299575 0.1729207 
## t98 by_24 fp_24 0.1449308 0.1127846 0.02083603 
## t99 by_24 fp_12 0.2483133 0.3978382 0.01512726 
## t100 by_20 fp_20 0.4700827 0.4150916 0.04749334 
## t101 by_20 fp_20 0.562851 0.5073988 0.05479139 
## t102 by_24 fp_24 0.542347 0.4705623 0.01038248 
## t103 by_16 fp_24 0.7680236 0.5538789 0.02334998 
## t104 by_20 fp_16 0.7528534 0.7851965 0.01105828 
## t105 by_20 fp_24 0.8118401 0.7213535 0.009594842 
## t106 by_20 fp_8 0.3821096 0.5974534 0.02161384 
## t107 by_24 fp_12 0.3014682 0.4632741 0.03272254 
## t108 by_12 fp_8 0.6244745 0.644648 0.0381548 
## t109 by_16 fp_24 0.9373344 0.8486946 0.07139357 
## t110 by_24 fp_24 0.4024523 0.3356056 0.01817566 
## t111 by_24 fp_24 0.3540557 0.2913284 0.0154562 
## t112 by_16 fp_24 0.9373142 0.8486503 0.003189668 
## t113 by_16 fp_8 0.4850307 0.5855429 0.137029 
## t114 by_24 fp_20 0.2985214 0.2985214 0.04606412 
## t115 by_20 fp_8 0.2571696 0.4538163 0.1057102 
## t116 by_16 fp_20 0.6424755 0.4732699 0.07233299 
## t117 by_16 fp_24 0.3408513 0.1624199 0.02604788 
## t118 by_24 fp_20 0.3113952 0.3113952 0.06326332 
## t119 by_20 fp_16 0.7326868 0.7668517 0.04084677 
## t120 by_20 fp_8 0.4928875 0.699941 0.01480073 
## t121 by_24 fp_20 0.4581035 0.4581035 0.005902859 
## t122 by_16 fp_24 0.5522561 0.3162544 0.0234387 
## t123 by_20 fp_12 0.56051 0.6711157 0.02328764 
## t124 by_20 fp_8 0.2213473 0.405557 0.1133474 
## t125 by_24 fp_20 0.7319665 0.7319665 0.0247748 
## t126 by_12 fp_20 0.7767218 0.5584966 0.1340342 
## t127 by_20 fp_8 0.2229879 0.4078478 0.2385047 
## t128 by_24 fp_20 0.7092722 0.7092722 0.01424119 
## t129 by_24 fp_24 0.7107874 0.6482896 0.1539451 
## t130 by_12 fp_8 0.8783378 0.8873341 0.1464284 
## t131 by_16 fp_16 0.2808852 0.2265743 0.3632383 
## t132 by_16 fp_8 0.5825304 0.6766973 0.02521979 
## t133 by_24 fp_24 0.7409962 0.6821067 0.004198461 
## t134 by_24 fp_20 0.7989193 0.7989193 0.0069483 
## t135 by_24 fp_12 0.674706 0.8057605 0.05329944 
## t136 by_16 fp_20 0.905165 0.8267593 0.09586584 
## t137 by_24 fp_24 0.5835395 0.5124071 0.0113247 
## t138 by_24 fp_20 0.4465982 0.4465982 0.003996356 
## t139 by_12 fp_4 0.825238 0.9115255 0.1221431 
## t140 by_12 fp_16 0.8791412 0.7986997 0.04836262 
## t141 by_16 fp_16 0.7484522 0.6905501 0.1202303 
## t142 by_20 fp_24 0.8302395 0.7458311 0.219393 
## t143 by_20 fp_8 0.5403802 0.7383368 0.04951988 
## t144 by_20 fp_16 0.5619114 0.6061708 0.008521986 
## t145 by_8 fp_12 0.7654232 0.6199884 0.0285369 
## t146 by_20 fp_16 0.5297468 0.5747969 0.009204492 
## t147 by_12 fp_12 0.7733074 0.7127193 0.03365675 
## t148 by_12 fp_12 0.6106679 0.5328692 0.03318607 
## t149 by_24 fp_24 0.7371338 0.6777482 0.005421159 
## t150 by_24 fp_12 0.2374832 0.3838165 0.008675597 
## t151 by_12 fp_16 0.5729726 0.422591 0.03233464 
## t152 by_16 fp_12 0.4921725 0.4921725 0.02106447 
## t153 by_20 fp_20 0.9292467 0.9130955 0.007559265 
## t154 by_20 fp_24 0.8555655 0.7804192 0.0399142 
## t155 by_24 fp_8 0.01205856 0.03532377 0.147405 
## t156 by_20 fp_24 0.8877497 0.8259416 0.005174927 
## t157 by_16 fp_4 0.2263339 0.4674171 0.2329861 
## t158 by_24 fp_8 0.1837757 0.403149 0.02236869 
## t159 by_20 fp_20 0.5449749 0.4893126 0.01184148 
## t160 by_16 fp_24 0.9664806 0.9153443 0.01631543 
## t161 by_24 fp_20 0.8685325 0.8685325 0.007489154 
## t162 by_16 fp_20 0.3771119 0.2323709 0.008634481 
## t163 by_24 fp_16 0.4808623 0.5814861 0.004614095 
## t164 by_24 fp_8 0.09662028 0.242919 0.04491128 
## t165 by_24 fp_24 0.3275842 0.2676039 0.02348691 
## t166 by_20 fp_12 0.04623116 0.0719734 0.006093831 
## t167 by_20 fp_4 0.1652672 0.4872697 0.05009026 
## t168 by_20 fp_8 0.374187 0.5893243 0.02012077 
## t169 by_12 fp_20 0.8786015 0.7246512 0.002203251 
## t170 by_16 fp_24 0.9088798 0.7890493 0.006812454 
## t171 by_16 fp_24 0.9304277 0.8337511 0.02143603 
## t172 by_20 fp_8 0.21371 0.3947863 0.02647118 
## t173 by_24 fp_20 0.5972994 0.5972994 0.01983763 
## t174 by_24 fp_20 0.4936407 0.4936407 0.007188798 
## t175 by_16 fp_4 0.3243085 0.5901465 0.03683703 
## t176 by_12 fp_24 0.9580473 0.8616511 0.01229351 
## t177 by_16 fp_16 0.4740201 0.4033093 0.01300839 
## t178 by_20 fp_12 0.7560619 0.8321878 0.02884707 
## t179 by_16 fp_20 0.8600352 0.7544402 0.0525245 
## t180 by_20 fp_24 0.9107013 0.8595313 0.00958726 
## t181 by_20 fp_20 0.7470553 0.702624 0.0452843 
## t182 by_24 fp_24 0.2545211 0.2038627 0.0125026 
## t183 by_20 fp_20 0.4784466 0.4232587 0.004873196 
## t184 by_20 fp_8 0.3755014 0.5906811 0.01557065 
## t185 by_16 fp_20 0.7708098 0.6270874 0.00658974 
## t186 by_20 fp_8 0.4521481 0.6645135 0.007768495 
## t187 by_24 fp_24 0.6375919 0.5688708 0.02085719 
## t188 by_24 fp_24 0.4104308 0.3430196 0.01445605 
## t189 by_20 fp_24 0.8818834 0.8175091 0.01486525 
## t190 by_24 fp_20 0.7436924 0.7436924 0.1358586 
## t191 by_16 fp_12 0.4266438 0.4266438 0.04803308 
## t192 by_20 fp_8 0.3658438 0.5806348 0.02013446 
## t193 by_20 fp_4 0.4117604 0.770639 0.2095939 
## t194 by_24 fp_20 0.4744003 0.4744003 0.05875816 
## t195 by_16 fp_24 0.523247 0.2915694 0.01058919 
## t196 by_20 fp_16 0.5773146 0.6210673 0.01708944 
## t197 by_24 fp_20 0.259134 0.259134 0.01460659 
## t198 by_20 fp_8 0.2360995 0.4258716 0.1333206 
## t199 by_24 fp_20 0.3620791 0.3620791 0.006263922

res.df <- data.frame(depvar = varname.col, cmse = cmse.col, pmse = pmse.col, totalr2 = r2.col, delta = delta.col, adjdelta = adjdelta.col, bywidth = bywidth.col, fpwidth = fpwidth.col, bydf = bydf.col, fpdf = fpdf.col, pmse_oos = pmse_oos.col, cmse_oos = cmse_oos.col, delta_oos = delta_oos.col, r2_oos = r2_oos.col, agemse = age.col)
write.csv(res.df, file = 'topicmodel_deltas4.csv', quote = FALSE, row.names = FALSE)
```

RESULTS:

    ## Mean delta is  0.55534
    ## If we adjust for df it is  0.55979
    ## The average r2 is  0.0492
    ## The weighted average, sum(cmse.col) / (sum(cmse.col) + sum(pmse.col)):
    ## 0.55206
    ## 
    ## Mean delta measured oos is  0.51125
    ## but as a weighted average it is: 0.5003

<img src="topicmodel_gridsearch_files/figure-gfm/unnamed-chunk-7-1.svg" width="100%" style="display: block; margin: auto;" />
