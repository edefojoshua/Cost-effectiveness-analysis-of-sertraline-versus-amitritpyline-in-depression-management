Cost effectiveness analysis of sertraline versus Amitriptyline
================
Joshua Edefo
2024-01-11

libraries

``` r
# Markov cohort simulation with probability Sensitivity analysis of Sertraline
library(markovchain)
```

    ## Warning: package 'markovchain' was built under R version 4.3.2

``` r
library(cartography)
```

    ## Warning: package 'cartography' was built under R version 4.3.2

``` r
library(rmarkdown)
```

    ## Warning: package 'rmarkdown' was built under R version 4.3.2

``` r
library(sf)
```

    ## Warning: package 'sf' was built under R version 4.3.2

``` r
library(sp)
```

    ## Warning: package 'sp' was built under R version 4.3.2

Sertraline

    Diagrammatic presentation of transition matrix

    ```r
    Set.matrix<-new("markovchain", states=c("Depressed", "Borderline", "Depression_free"), 
             transitionMatrix=matrix(data=c(0.1, 0.3, 0.6, 0.3, 0.1, 0.6, 0.1, 0.1, 0.8), 
                                       byrow = TRUE, nrow =  3), name = "Sertraline") 
    plot(Set.matrix, package = "daigram")

![](Cost-effectiveness-analysis-of-sertraline-versus-amitritpyline-in-depression-management._files/figure-gfm/ab-1.png)<!-- -->

``` r
# R code inside a function
N_psa <- 10000

# distribution of transition probability matrix using beta distribution with SD= 25% of the mean value
params<-data.frame(
  p_Depressed_Depressed   = rbeta(N_psa, 10, 90),
  p_Depressed_Borderline  = rbeta(N_psa,  30, 70),
  p_Depressed_Depressionfree  = rbeta(N_psa,  60, 40),
  p_Bordeline_Depressed  = rbeta(N_psa,  30, 70),
  p_Borderline_Borderline = rbeta(N_psa, 10, 90),
  p_Borderline_Depressionfree = rbeta(N_psa,  60, 40),
  p_Depressionfree_Depressed  = rbeta(N_psa,  10, 90),
  p_Depressionfree_Borderline = rbeta(N_psa,  10,  90),
  p_Depressionfree_Depressionfree = rbeta(N_psa, 80, 20),
  
  # distribution of cost using gamma distribution with SD= 25% of the mean value
  c_Depressed = rgamma(N_psa, shape =16, scale = 1432/16),
  c_Borderline = rgamma(N_psa, shape =16, scale =4296/16),
  c_Depressionfree = rgamma(N_psa, shape =16, scale =8593/16),
  
  # distribution of effect using gamma distribution with SD= 25% of the mean value
  q_Depressed = rgamma (N_psa, shape =16, scale = 2/16),
  q_Borderline = rgamma(N_psa, shape =16, scale = 6/16),
  q_Depressionfree = rgamma(N_psa, shape =16, scale = 9/16)
)

model.se<- function(.params) { 
  with(.params, {
    
    n_t<-40     # no of cycles
    n_s<-3      # no of states
    n_c<-1000   #  no of persons in cohort
    
    v_state_names<-c("Depressed", "Borderline", "Depressionfree") # names of states
    
    # probability tranistion matrix
    m_P<-matrix(0, nrow = 3, ncol =3, 
                dimnames = list(from = v_state_names, to = v_state_names ))
    m_P["Depressed", "Depressed"]<-p_Depressed_Depressed 
    m_P["Depressed", "Borderline"]<- p_Depressed_Borderline
    m_P["Depressed", "Depressionfree"]<-  p_Depressed_Depressionfree
    m_P["Borderline", "Depressed"]<-p_Bordeline_Depressed
    m_P["Borderline", "Borderline"]<-p_Borderline_Borderline 
    m_P["Borderline", "Depressionfree"]<-  p_Borderline_Depressionfree
    m_P["Depressionfree", "Depressed"]<-p_Depressionfree_Depressed
    m_P["Depressionfree", "Borderline"]<-p_Depressionfree_Borderline
    m_P["Depressionfree", "Depressionfree"]<-p_Depressionfree_Depressionfree
    
    # State membership
    state_membership<-array(NA_real_,  dim= c(n_t, n_s), dimnames = list (cycle =1:n_t, state = v_state_names ))
    state_membership[1, ] <- c(n_c, 0, 0)
    for (i in 2:n_t) {state_membership[i, ] <-state_membership[i-1, ] %*% m_P}
    
    # payoffs
    m_payoffs <- matrix(0, nrow = 3, ncol = 2, 
                        dimnames = list(v_state_names, payoff= c("Cost", "HDS")))
    m_payoffs["Depressed", "Cost"]<- c_Depressed
    m_payoffs["Depressed", "HDS"]<- q_Depressed
    m_payoffs["Borderline", "Cost"]<- c_Borderline
    m_payoffs["Borderline", "HDS"]<- q_Borderline
    m_payoffs["Depressionfree", "Cost"]<- c_Depressionfree
    m_payoffs["Depressionfree", "HDS"]<- q_Depressionfree
    
    payoff_trace.se<-state_membership %*% m_payoffs
    payoff_trace.se
    
    summary_results.se = colSums(payoff_trace.se)/n_c 
    
  })
}
psa_results.se <-
  t(sapply(
    X = split(params, 1:N_psa), 
    FUN = model.se,
    simplify = TRUE ))

write.csv(psa_results.se, file="psa_coef3_se.csv")

plot(
  psa_results.se[, 2],
  psa_results.se[, 1],
  type = 'p',
  xlab = "HDS",
  ylab = "Cost"
)
```

![](Cost-effectiveness-analysis-of-sertraline-versus-amitritpyline-in-depression-management._files/figure-gfm/b-1.png)<!-- -->
Amitriptyline

Diagrammatic presentation of transition matrix

``` r
Ami.matrix<-new("markovchain", states=c("Depressed", "Borderline", "Depression_free"), 
         transitionMatrix=matrix(data=c(0.2, 0.3, 0.5, 0.4, 0.2, 0.4, 0.2, 0.2, 0.6), 
                                   byrow = TRUE, nrow =  3), name = "Amitriptyline") 
plot(Ami.matrix, package = "daigram")
```

![](Cost-effectiveness-analysis-of-sertraline-versus-amitritpyline-in-depression-management._files/figure-gfm/ca-1.png)<!-- -->
Markov cohort simulation with probability Sensitivity analysis of
amitriptlyine

``` r
N_psa <- 10000

# R code inside a function

# distribution of transition probability matrix using beta distribution with SD= 25% of the mean value
params<-data.frame(
  p_Depressed_Depressed   = rbeta(N_psa, 20, 80),
  p_Depressed_Borderline  = rbeta(N_psa,  30, 70),
  p_Depressed_Depressionfree  = rbeta(N_psa,  50, 50),
  p_Bordeline_Depressed  = rbeta(N_psa,  40, 60),
  p_Borderline_Borderline = rbeta(N_psa, 20, 80),
  p_Borderline_Depressionfree = rbeta(N_psa,  40, 80),
  p_Depressionfree_Depressed  = rbeta(N_psa,  20, 80),
  p_Depressionfree_Borderline = rbeta(N_psa,  20,  80),
  p_Depressionfree_Depressionfree = rbeta(N_psa, 60, 40),
 
  # distribution of cost using gamma distribution with SD= 25% of the mean value 
  c_Depressed = rgamma(N_psa, shape = 16.05, scale = 41.19),
  c_Borderline = rgamma(N_psa, shape = 16, scale = 62), 
  c_Depressionfree = rgamma(N_psa, shape = 16.01, scale = 103.19),
  
  # distribution of cost using gamma distribution with SD= 25% of the mean value
  q_Depressed = rgamma(N_psa, shape = 16., scale = 0.125),
  q_Borderline = rgamma(N_psa, shape = 16.13, scale = 0.31),
  q_Depressionfree = rgamma(N_psa, shape = 16.05, scale = 0.5)
)

model.ami<- function(.params) { 
  with(.params, {
    
    n_t<-40
    n_s<-3
    n_c<-1000
    
    v_state_names<-c("Depressed", "Borderline", "Depressionfree")
    m_P<-matrix(0, nrow = 3, ncol =3, 
                dimnames = list(from = v_state_names, to = v_state_names ))
    m_P["Depressed", "Depressed"]<-p_Depressed_Depressed 
    m_P["Depressed", "Borderline"]<- p_Depressed_Borderline
    m_P["Depressed", "Depressionfree"]<-  p_Depressed_Depressionfree
    m_P["Borderline", "Depressed"]<-p_Bordeline_Depressed
    m_P["Borderline", "Borderline"]<-p_Borderline_Borderline 
    m_P["Borderline", "Depressionfree"]<-  p_Borderline_Depressionfree
    m_P["Depressionfree", "Depressed"]<-p_Depressionfree_Depressed
    m_P["Depressionfree", "Borderline"]<-p_Depressionfree_Borderline
    m_P["Depressionfree", "Depressionfree"]<-p_Depressionfree_Depressionfree
    
    # State membership
    state_membership<-array(NA_real_,  dim= c(n_t, n_s), dimnames = list (cycle =1:n_t, state = v_state_names ))
    state_membership[1, ] <- c(n_c, 0, 0)
    for (i in 2:n_t) {state_membership[i, ] <-state_membership[i-1, ] %*% m_P}
    
    # payoffs
    m_payoffs <- matrix(0, nrow = 3, ncol = 2, 
                        dimnames = list(v_state_names, payoff= c("Cost", "HDS")))
    m_payoffs["Depressed", "Cost"]<- c_Depressed
    m_payoffs["Depressed", "HDS"]<- q_Depressed
    m_payoffs["Borderline", "Cost"]<- c_Borderline
    m_payoffs["Borderline", "HDS"]<- q_Borderline
    m_payoffs["Depressionfree", "Cost"]<- c_Depressionfree
    m_payoffs["Depressionfree", "HDS"]<- q_Depressionfree
    
    payoff_trace.ami<-state_membership %*% m_payoffs
    payoff_trace.ami
    
    summary_results.ami = colSums(payoff_trace.ami)/n_c 
    
  })
}
psa_results.ami <-
  t(sapply(
    X = split(params, 1:N_psa), 
    FUN = model.ami,
    simplify = TRUE ))

write.csv(psa_results.ami, file="psa_coef3.csv")

plot(
  psa_results.ami[, 2],
  psa_results.ami[, 1],
  type = 'p',
  xlab = "HDS",
  ylab = "Cost"

)
```

![](Cost-effectiveness-analysis-of-sertraline-versus-amitritpyline-in-depression-management._files/figure-gfm/c-1.png)<!-- -->
Ecconomic evaluation

``` r
# incremental cost
data.result<- data.frame(psa_results.se, psa_results.ami)


inc_cost<- (data.result$Cost - data.result$Cost.1) 


# incremental effect
inc_effect<- (data.result$HDS - data.result$HDS.1)

# incremental cost effectiveness ratio
inc.CER <- inc_cost/inc_effect


# using a willingness to pay threshold of NGN18,000

wtp =  18000

psa_wtp <- sum(inc.CER <= wtp)

# percentage of cost effectiveness of sertraline over amitriptyline at wtp

result <- (psa_wtp /N_psa) * 100
 
result
```

    ## [1] 98.53

session information

``` r
sessionInfo()
```

    ## R version 4.3.1 (2023-06-16 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 11 x64 (build 22631)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.utf8 
    ## [2] LC_CTYPE=English_United Kingdom.utf8   
    ## [3] LC_MONETARY=English_United Kingdom.utf8
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.utf8    
    ## 
    ## time zone: Europe/London
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] sp_2.1-2          sf_1.0-15         rmarkdown_2.25    cartography_3.1.4
    ## [5] markovchain_0.9.5
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Matrix_1.6-1.1     expm_0.999-7       dplyr_1.1.3        compiler_4.3.1    
    ##  [5] tidyselect_1.2.0   Rcpp_1.0.11        parallel_4.3.1     yaml_2.3.7        
    ##  [9] fastmap_1.1.1      lattice_0.21-8     R6_2.5.1           generics_0.1.3    
    ## [13] igraph_1.5.1       classInt_0.4-10    knitr_1.44         tibble_3.2.1      
    ## [17] units_0.8-5        DBI_1.2.0          pillar_1.9.0       rlang_1.1.1       
    ## [21] utf8_1.2.3         xfun_0.40          RcppParallel_5.1.7 cli_3.6.1         
    ## [25] magrittr_2.0.3     class_7.3-22       digest_0.6.33      grid_4.3.1        
    ## [29] rstudioapi_0.15.0  lifecycle_1.0.3    vctrs_0.6.3        KernSmooth_2.23-21
    ## [33] proxy_0.4-27       evaluate_0.21      glue_1.6.2         stats4_4.3.1      
    ## [37] fansi_1.0.4        e1071_1.7-14       tools_4.3.1        pkgconfig_2.0.3   
    ## [41] htmltools_0.5.6
