library(tictoc)
library(data.table)
library(brms)
library(progressr)
library(tidyverse)
library(plotly)
library(furrr)
library(stringr)
library(dials)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

memory.limit(80000)


# future::plan(multicore,
#              workers = 4)

plan(multisession)

sample_specification = list_to_sim[[1]]
# creates the sample
build_sample <- function(
  sample_specification
){
  
  
  u1 = sample_specification["u"]
  effect_size_u = sample_specification["effect_size_u"]
  allocate_n1 = sample_specification["allocate_n1"]
  allocate_n2 = sample_specification["allocate_n2"]
  sample_size = sample_specification["sample_size"]
  sigma_u = sample_specification['sigma_u']
  effect_size_sigma = sample_specification['effect_size_sigma']
  
  # get the sample size right for the two groups while respecting the allocations 
  number_batches = floor(sample_size/(allocate_n1 + allocate_n2))
  n1 = number_batches * allocate_n1
  n2 = number_batches * allocate_n2
  total_sample_size = n1 + n2
  
  
  # create the coefficient matrix
  treatment_vector = c(rep(0,n1),
                       rep(1,n2))
  
  intercept_vector = rep(1,total_sample_size)
  
  x_matrix <-   matrix(c(intercept_vector,
                         treatment_vector),
                       nrow = total_sample_size)
  
  
  # create the mean regression vector
  coeff_mean_vector <- matrix(c(u1,effect_size_u))
  u_regression <- x_matrix %*% coeff_mean_vector
  
  
  # create the variance regression vector
  coeff_sigma_vector <- matrix(c(sigma_u,effect_size_sigma))
  sigma_regression <- x_matrix %*% coeff_sigma_vector
  
  
  y <- rnorm(total_sample_size,
             mean = u_regression,
             sd = sigma_regression)
  
  x <- c(rep("group one",n1),rep("group two",n2))
  
  
  df <- data.frame(x = factor(x),
                   y = y)
  
  
  model <- update(regress,newdata = df,formula. = bf(y ~ x,sigma ~ x))
  
  
  values <-   as.vector(t(fixef(model)[c(1,3),c(1,3,4)])) 
  

  hyp <- c("exp(sigma_Intercept) = 0",
           "exp(sigma_Intercept + sigma_xgrouptwo) = 0")

  values <-  c(values,as.vector(t(hypothesis(model, hyp)$hypothesis[1,c(2,4,5)])))
  
  hyp <- "exp(sigma_Intercept + sigma_xgrouptwo) > exp(sigma_Intercept)"
  
  values <-  c(values,as.vector(t(hypothesis(model, hyp)$hypothesis[1,c(2,4,5)])))
  
  
  return(values)
}




# set up the specs for the data generation 

parameters_grid <- expand.grid(
  allocation = c('1-1','1-3','1-5'),
  sample_size = seq(30,400,10),
  u = c(2),
  effect_size_u = c(0.5,0.8), 
  sigma_u = c(1),
  effect_size_sigma = c(0,1,2),
  nsim = 1:20 
)

time_est <- nrow(parameters_grid) / 360 * 3.6 / 60

# add on a unique counter for split into a list
parameters_grid$unique = 1:nrow(parameters_grid)


# split the allocation to get the allocate parameters
parameters_grid <- parameters_grid %>% 
  mutate(
    allocate_n1 = as.numeric(str_split(allocation,"-",simplify = T)[,1]),
    allocate_n2 = as.numeric(str_split(allocation,"-",simplify = T)[,2]),
  ) 



#  create the list to feed into the future
list_to_sim <- parameters_grid %>%
  dplyr::select(!allocation) %>% 
  group_split(unique) %>%
  map(unlist)



# initiate the regression 
# dataset <- data.frame(
#   x = factor(c(rep("group one",20),rep("group two",20))),
#   y = rnorm(40)
# )
# 
# regress <- brm(bf(y ~ x,sigma ~ x), dataset)



# run the sims
tic()
info_list <- future_map( 
  .x = list_to_sim,
  .f =  build_sample
)
toc()

saveRDS(info_list,
        paste0(getwd(),"/documents/two_sample_power_files/info_list_table_problem.rds"))


#16294/13680 for one simulation 

info <- data.table(matrix(unlist(info_list), nrow=length(info_list), byrow=TRUE))




names(info) <- c(
  'Intercept',
  'Intercept_q2.5',
  'Intercept_q97.5',
  
  'xgrouptwo',
  'xgrouptwo_q2.5',
  'xgrouptwo_q97.5',
  
  
  'sigma_Intercept',
  'sigma_Intercept_q2.5',
  'sigma_Intercept_q97.5',
  
  'sigma_xgrouptwo',
  'sigma_xgrouptwo_q2.5',
  'sigma_xgrouptwo_q97.5'
  
)


# deal with after simulation 

info <- cbind(parameters_grid,info)

Intercept_wantcover <- 0.7
xgrouptwo_wantcover <- 0.7
sigma_Intercept_wantcover <- 0.7
sigma_xgrouptwo_wantcover <- 0.7

info <- lazy_dt(info)


info <- info %>%
  
  mutate(
    
    Intercept_cover = Intercept_q97.5 - Intercept_q2.5 < Intercept_wantcover,
    Intercept_detect =  u > Intercept_q2.5 &  u < Intercept_q97.5 ,
    Intercept_length =  Intercept_q97.5 -  Intercept_q2.5 ,
    
    
    
    xgrouptwo_cover = xgrouptwo_q97.5 - xgrouptwo_q2.5  < xgrouptwo_wantcover,
    xgrouptwo_detect =  effect_size_u > xgrouptwo_q2.5 &  effect_size_u < xgrouptwo_q97.5, 
    xgrouptwo_length = xgrouptwo_q97.5 - xgrouptwo_q2.5,
    
    
    sigma_Intercept_cover = sigma_Intercept_q97.5 - sigma_Intercept_q2.5 < sigma_Intercept_wantcover,
    sigma_Intercept_detect =   sigma_u > sigma_Intercept_q2.5 &  sigma_u < sigma_Intercept_q97.5, 
    sigma_Intercept_length = sigma_Intercept_q97.5 - sigma_Intercept_q2.5 ,
    
    
    sigma_xgrouptwo_cover = sigma_xgrouptwo_q97.5 - sigma_xgrouptwo_q2.5 < sigma_xgrouptwo_wantcover,
    sigma_xgrouptwo_detect = effect_size_sigma > sigma_xgrouptwo_q2.5 &  effect_size_sigma < sigma_xgrouptwo_q97.5, 
    sigma_xgrouptwo_length = sigma_xgrouptwo_q97.5 - sigma_xgrouptwo_q2.5 
    
    
  ) %>% 
  
  mutate(para_com = paste( 
    as.character(u),
    "-",
    as.character(effect_size_u),
    "-",
    as.character(sigma_u),
    "-",
    as.character(effect_size_sigma)
                           ),
    .keep = "unused") %>% 
  
  dplyr::select(!c(allocate_n1,allocate_n2,nsim,unique,.keep,u,effect_size_u,effect_size_sigma,sigma_u)) %>% 
  
  mutate(
    across(ends_with("cover"),as.numeric),
    across(ends_with("detect"),as.numeric),
    
  ) %>% 
  
  group_by(
    allocation,
    para_com,
    sample_size
  ) %>%
  
  mutate(
     across(.cols = everything(),
             mean) 
  ) %>% 
  distinct() %>% 
  arrange(allocation,
          para_com,
          sample_size) %>% 
  
  as_tibble()

info$n1 <- 0
info$n2 <- 0

for(i in 1:nrow(info)){
  allocate_adjust <- as.numeric(unlist(str_split(info$allocation[i],"-")))
  number_batches = floor(info$sample_size[i]/(allocate_adjust[1] + allocate_adjust[2]))
  n1 = number_batches * allocate_adjust[1]
  n2 = number_batches * allocate_adjust[2]
  info$n1[i] = n1
  info$n2[i] = n2
  info$sample_size[i] = n1 + n2
}




# graphs

groupone_info <- info %>%
  filter(para_com == "2 - 0.5 - 1 - 0") %>% 
  ggplot(aes(x = sample_size,
             y = Intercept,
             color = allocation,
             ymin = Intercept_q2.5,
             ymax = Intercept_q97.5,
             text = paste('interval length: ', Intercept_length,
                          '</br>interval length criteria satisfied: ', Intercept_cover,
                          '</br>detect probability: ', Intercept_detect,
                          '</br>group one sample size: ',n1,
                          '</br>group two sample size: ',n2
             ))) +
  geom_hline(yintercept = c(0, .5), color = "white") +
  geom_pointrange(fatten = 1/2) +
  # scale_x_discrete("reordered by the lower level of the 95% intervals", breaks = NULL) + 
  geom_point() + 
  ylim(min(info$Intercept_q2.5) - 0.5,max(info$Intercept_q97.5) + 0.5) + 
  geom_hline(yintercept = 2) + 
  labs(x = "Sample Size (Total)",
       y = "Location Effect Size",
       title = "Properties for Main Effect across Sample Sizes") + 
  theme_bw() 


ggplotly(groupone_info) %>% 
  layout(hovermode = "x")


grouptwo_info <- info %>%
  filter(para_com == "2 - 0.5 - 1 - 0") %>% 
  ggplot(aes(x = sample_size,
             y = xgrouptwo,
             color = allocation,
             ymin = xgrouptwo_q2.5,
             ymax = xgrouptwo_q97.5,
             text = paste('interval length: ', xgrouptwo_length,
                          '</br>interval length criteria satisfied: ', xgrouptwo_cover,
                          '</br>detect probability: ', xgrouptwo_detect,
                          '</br>group one sample size: ',n1,
                          '</br>group two sample size: ',n2
             ))) +
  geom_hline(yintercept = c(0, .5), color = "white") +
  geom_pointrange(fatten = 1/2) +
  # scale_x_discrete("reordered by the lower level of the 95% intervals", breaks = NULL) + 
  geom_point() + 
  ylim(min(info$xgrouptwo_q2.5) - 0.5,max(info$xgrouptwo_q97.5) + 0.5) + 
  geom_hline(yintercept = 0.5) + 
  labs(x = "Sample Size (Total)",
       y = "Location Effect Size",
       title = "Properties of the Location Effect across Sample Sizes") + 
  theme_bw() 


ggplotly(grouptwo_info) %>% 
  layout(hovermode = "x")



groupone_sd_info <- info %>%
  filter(para_com == "2 - 0.5 - 1 - 0") %>% 
  ggplot(aes(x = sample_size,
             y = sigma_Intercept,
             color = allocation,
             ymin = sigma_Intercept_q2.5,
             ymax = sigma_Intercept_q97.5,
             text = paste('interval length: ', sigma_Intercept_length,
                          '</br>interval length criteria satisfied: ', sigma_Intercept_cover,
                          '</br>detect probability: ', sigma_Intercept_detect,
                          '</br>group one sample size: ',n1,
                          '</br>group two sample size: ',n2
             ))) +
  # geom_hline(yintercept = c(0, .5), color = "white") +
  geom_pointrange(fatten = 1/2) +
  # scale_x_discrete("reordered by the lower level of the 95% intervals", breaks = NULL) + 
  geom_point() + 
  ylim(min(info$sigma_Intercept_q2.5) - 0.5,max(info$sigma_Intercept_q97.5) + 0.5) + 
  geom_hline(yintercept = 1) + 
  labs(x = "Sample Size (Total)",
       y = "Location Effect Size",
       title = "Properties for Main Effect across Sample Sizes") + 
  theme_bw() 


ggplotly(groupone_sd_info) %>% 
  layout(hovermode = "x")



grouptwo_sd_info <- info %>%
  filter(para_com == "2 - 0.5 - 1 - 0") %>% 
  ggplot(aes(x = sample_size,
             y = sigma_xgrouptwo,
             color = allocation,
             ymin = sigma_xgrouptwo_q2.5,
             ymax = sigma_xgrouptwo_q97.5,
             text = paste('interval length: ', sigma_xgrouptwo_length,
                          '</br>interval length criteria satisfied: ', sigma_xgrouptwo_cover,
                          '</br>detect probability: ', sigma_xgrouptwo_detect,
                          '</br>group one sample size: ',n1,
                          '</br>group two sample size: ',n2
             ))) +
  geom_hline(yintercept = c(0, .5), color = "white") +
  geom_pointrange(fatten = 1/2) +
  # scale_x_discrete("reordered by the lower level of the 95% intervals", breaks = NULL) + 
  geom_point() + 
  ylim(min(info$sigma_xgrouptwo_q2.5) - 0.5,max(info$sigma_xgrouptwo_q97.5) + 0.5) + 
  geom_hline(yintercept = 0) + 
  labs(x = "Sample Size (Total)",
       y = "Location Effect Size",
       title = "Properties of the Scale Effect across Sample Sizes") + 
  theme_bw() 


ggplotly(grouptwo_sd_info) %>% 
  layout(hovermode = "x")




saveRDS(info,
        paste0(getwd(),"/documents/two_sample_power_files/all_problems.rds"))




# intercep intercep2.5 intercept 97.5 effect effect2.5 effect 97.5 sigmaintercept sigmaintere2.5 siginter97.5 sigmaeffec sigmeffec2.5 sigmaeffect97.5



# runs the regression and extracts the info we want
model_extract <- function(dataframeobject){
  model <- update(regress,newdata = dataframeobject,formula. = bf(y ~ x,sigma ~ x))
  
  hyp <- c("exp(sigma_Intercept) = 0",
           "exp(sigma_Intercept + sigma_xgrouptwo) = 0")
  sigma_info <- hypothesis(model, hyp)$hypothesis[,2:5]
  
  hyp <- "exp(sigma_Intercept + sigma_xgrouptwo) > exp(sigma_Intercept)"
  sigma_info2 <- hypothesis(model, hyp)$hypothesis[,2:5]
  
  
  
  estimates <- fixef(model)[,c(1,3,4)]
  
  return(c(
    estimates
  ))
}





x <- c(rep("group one",n1),rep("group two",n2))


df <- data.frame(x = factor(x),
                 y = y)



hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_xgrouptwo) = 0")
(hyp <- hypothesis(regress, hyp))
plot(hyp)
# 
hyp <- "exp(sigma_Intercept + sigma_xgrouptwo) > exp(sigma_Intercept)"
(hyp <- hypothesis(regress, hyp))
# 
plot(hyp, chars = NULL)


a <-  model_extract(dataset)
# regress <- brm(y ~ x, dataset)
# regress <- update(regress,dataset)

# saveRDS(regress,
#         paste0(getwd(),"/documents/two_sample_power_files/regression1.rds"))

# summary(regress)
# plot(regress, ask = FALSE)
# plot(conditional_effects(regress), points = TRUE)

# b0pop <- 2
# b1pop <- 5
# y <- list()
# xdat <- list()
# nsims <- 10
Intercept_wantcover <- 0.7
xgrouptwo_wantcover <- 0.7
sigma_Intercept_wantcover <- 0.7
sigma_xgrouptwo_wantcover <- 0.7


samplerange <- seq(30,400,10)
nsims = 20
allocation <- list(
  c(1,1),
  c(1,3),
  c(1,5)
)


parameters_grid <- expand.grid(
  pop_u = c(2),
  pop_effect_size_u = c(0.5,0.3,0.8), 
  pop_sigma_u1 = c(1),
  pop_effect_size_sigma = c(1,2)
)

parameters <- parameters_grid %>% purrr::transpose()




allocation_data <- tibble()

tic()

for(allocate in allocation){
  
  
  for(params in parameters){
  
  toplot <- data.frame(
    samplesize  = samplerange,
    Intercept = 0,
    Intercept_q2.5 = 0,
    Intercept_q97.5= 0,
    sigma_Intercept= 0,
    sigma_Intercept_q2.5= 0,
    sigma_Intercept_q97.5= 0,
    xgrouptwo= 0,
    xgrouptwo_q2.5= 0,
    xgrouptwo_q97.5= 0,
    sigma_xgrouptwo= 0,
    sigma_xgrouptwo_q2.5= 0,
    sigma_xgrouptwo_q97.5= 0,
    Intercept_cover = 0,
    Intercept_detect = 0,
    xgrouptwo_cover = 0,
    xgrouptwo_detect = 0,
    sigma_Intercept_cover = 0,
    sigma_Intercept_detect = 0,
    sigma_xgrouptwo_cover = 0,
    sigma_xgrouptwo_detect =0
  )
  
  
  
  
  
  
  
  for(n in samplerange){
    # generate the datasets
    dataframes <- future_map( 
      .x = 1:nsims,
      .f =  build_sample,
      u1 = params[["pop_u"]],
      effect_size_u = params[["pop_effect_size_u"]],
      allocate_n1 = allocate[1],
      allocate_n2 = allocate[2],
      sample_size = n,
      sigma_u1 = params[["pop_sigma_u1"]],
      effect_size_sigma = params[["pop_effect_size_sigma"]],
      .options = furrr_options(seed = T)
    )
  
    
    
    
    # extract the values from each dataset after running regression  
    values <- future_map(
      .x = dataframes,
      .f = model_extract
    )
    
  
    # curate the results
    results <- data.frame(matrix(unlist(values), nrow =nsims, byrow=TRUE))
    names(results) <- c(
      'Intercept',
      'sigma_Intercept',
      'xgrouptwo',
      'sigma_xgrouptwo',
      
      'Intercept_q2.5',
      'sigma_Intercept_q2.5',
      'xgrouptwo_q2.5',
      'sigma_xgrouptwo_q2.5',
      
      'Intercept_q97.5',
      'sigma_Intercept_q97.5',
      'xgrouptwo_q97.5',
      'sigma_xgrouptwo_q97.5'
    )
  

    
    
    
    # check for interval coverage and whether the parameter is inside the credible interval 
   
  
    # extract the mean of the parameters descriptors for the simulations for each n 
    
    toplot[[toplot$samplesize == n,'Intercept']] <- mean(results$Intercept)
    toplot[[toplot$samplesize == n,'Intercept_q2.5']] <- mean(results$Intercept_q2.5)
    toplot[[toplot$samplesize == n,'Intercept_q97.5']] <- mean(results$Intercept_q97.5)
    toplot[[toplot$samplesize == n,'sigma_Intercept']] <- mean(results$sigma_Intercept)
    toplot[[toplot$samplesize == n,'sigma_Intercept_q2.5']] <- mean(results$sigma_Intercept_q2.5)
    toplot[[toplot$samplesize == n,'sigma_Intercept_q97.5']] <- mean(results$sigma_Intercept_q97.5)
    toplot[[toplot$samplesize == n,'xgrouptwo']] <- mean(results$xgrouptwo)
    toplot[[toplot$samplesize == n,'xgrouptwo_q2.5']] <- mean(results$xgrouptwo_q2.5)
    toplot[[toplot$samplesize == n,'xgrouptwo_q97.5']] <- mean(results$xgrouptwo_q97.5)
    toplot[[toplot$samplesize == n,'sigma_xgrouptwo']] <- mean(results$sigma_xgrouptwo)
    toplot[[toplot$samplesize == n,'sigma_xgrouptwo_q2.5']] <- mean(results$sigma_xgrouptwo_q2.5)
    toplot[[toplot$samplesize == n,'sigma_xgrouptwo_q97.5']] <- mean(results$sigma_xgrouptwo_q97.5)
    
    toplot[[toplot$samplesize == n,'Intercept_cover']] <- mean(results$Intercept_cover)
    toplot[[toplot$samplesize == n,'Intercept_detect']] <- mean(results$Intercept_detect)
    
    
    toplot[[toplot$samplesize == n,'xgrouptwo_cover']] <- mean(results$xgrouptwo_cover)
    toplot[[toplot$samplesize == n,'xgrouptwo_detect']] <- mean(results$xgrouptwo_detect)
    
    toplot[[toplot$samplesize == n,'sigma_Intercept_cover']] <- mean(results$sigma_Intercept_cover)
    toplot[[toplot$samplesize == n,'sigma_Intercept_detect']] <- mean(results$sigma_Intercept_detect)
    
    toplot[[toplot$samplesize == n,'sigma_xgrouptwo_cover']] <- mean(results$sigma_xgrouptwo_cover)
    toplot[[toplot$samplesize == n,'sigma_xgrouptwo_detect']] <- mean(results$sigma_xgrouptwo_detect)
    
    
    toplot <- as.data.frame(apply(toplot,2,round,digits = 2))
    
  }
  
  
  allocation_data <- rbind(allocation_data,
                           toplot)
}
  }




toc()




# make adjustments to the sample size to respect the allocation and max sample size limit
allo <- paste(allocation)
allo <- str_replace_all(allo,c("c" = "",
                               "," = "",
                               "[//(]"="",
                               "[//)]" = "",
                               " " = "-"))

para <- paste(unlist(parameters))

para <- split(para, ceiling(seq_along(para)/4))

allocation_data$params <- rep(para,each = length(samplerange),times = length(allocation) )



allocation_data$allocation <- factor(
  rep(allo,each = nrow(parameters_grid) * length(samplerange))
)






allocation_data <- allocation_data %>% 
  mutate(
    Intercept_length = Intercept_q97.5 - Intercept_q2.5,
    sigma_Intercept_length = sigma_Intercept_q97.5 - sigma_Intercept_q2.5,
    xgrouptwo_length = xgrouptwo_q97.5 - xgrouptwo_q2.5,
    sigma_xgrouptwo_length = sigma_xgrouptwo_q97.5 - sigma_xgrouptwo_q2.5
  ) 


show_table <- allocation_data %>% 
  dplyr::select(
    allocation,
    params,
    samplesize,
    Intercept_detect,
    Intercept_length,
    xgrouptwo_detect,
    xgrouptwo_length,
    sigma_Intercept_detect,
    sigma_Intercept_length,
    sigma_xgrouptwo_detect,
    sigma_xgrouptwo_length
  )


saveRDS(show_table,
        paste0(getwd(),"/documents/two_sample_power_files/show_table_problem.rds"))

saveRDS(allocation_data,
        paste0(getwd(),"/documents/two_sample_power_files/allocation_data_problem_all.rds"))


