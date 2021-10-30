library(tictoc)
library(data.table)
library(brms)
library(progressr)
library(tidyverse)
library(plotly)
library(furrr)
library(stringr)


future::plan(multicore,
             workers = 4)



# creates the sample
build_sample <- function(
  nsims,
  u1,
  effect_size_u,
  allocate_n1,
  allocate_n2,
  sample_size,
  sigma_u1,
  effect_size_sigma
){
  
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
  coeff_sigma_vector <- matrix(c(sigma_u1,effect_size_sigma))
  sigma_regression <- x_matrix %*% coeff_sigma_vector
  
  
  y <- rnorm(total_sample_size,
             mean = u_regression,
             sd = sigma_regression)
  
  x <- c(rep("group one",n1),rep("group two",n2))
  
  
  df <- data.frame(x = factor(x),
                   y = y)
  # hist(c(y1,y2))
  
  
  return(df)
}







dataset <- build_sample(
  nsims = 1,
  u1 = 2,
  effect_size_u = 5,
  allocate_n1 = 1,
  allocate_n2 = 1,
  sample_size = 50,
  sigma_u1 = 1,
  effect_size_sigma = 0
)


ggplot(dataset,
       aes(y,
           fill = x)) +
  geom_histogram()

dataset %>%
  group_by(x) %>%
  summarise(
    mean = mean(y)
  )


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

a <-  model_extract(dataset)

# initial regress
# u1 = 2
# pop_u2 = 2.8
# allocate_n1 = 1
# allocate_n2 = 1
# n = 50
# sigma_u1 = 1
# sigma_u2 = 1
# dataset <- build_sample(1,par_sample_size = 20)
regress <- brm(bf(y ~ x,sigma ~ x), dataset)
# regress <- update(regress,dataset,formula. = bf(y ~ x,sigma ~ x))
# summary(regress)
# plot(regress, ask = FALSE)
# plot(conditional_effects(regress), points = TRUE)

hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_xgrouptwo) = 0")
(hyp <- hypothesis(regress, hyp))
plot(hyp)
# 
hyp <- "exp(sigma_Intercept + sigma_xgrouptwo) > exp(sigma_Intercept)"
(hyp <- hypothesis(regress, hyp))
# 
plot(hyp, chars = NULL)

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


samplerange <- seq(30,140,5)
nsims = 5
allocation <- list(
  c(1,1),
  c(1,3),
  c(1,5)
)
pop_u = 2
pop_effect_size_u = 0.8 
pop_sigma_u1 = 1
pop_effect_size_sigma = 2



allocation_data <- tibble()

tic()

for(allocate in allocation){
  
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
    dataframes <- future_map( .x = 1:nsims,
                              .f =  build_sample,
                              u1 = pop_u,
                              effect_size_u = pop_effect_size_u,
                              allocate_n1 = allocate[1],
                              allocate_n2 = allocate[2],
                              sample_size = n,
                              sigma_u1 = pop_sigma_u1,
                              effect_size_sigma = pop_effect_size_sigma,
                              .options = furrr_options(seed = T)
    )
    
    
    
    
    # extract the values from each dataset after running regression  
    values <- future_map(.x = dataframes,
                         .f = model_extract)
    
    # curate the results
    results <- data.frame(matrix(unlist(values), nrow =nsims, byrow=TRUE))
    names(results) <- c(
      'Intercept',
      'Intercept_q2.5',
      'Intercept_q97.5',
      'sigma_Intercept',
      'sigma_Intercept_q2.5',
      'sigma_Intercept_q97.5',
      'xgrouptwo',
      'xgrouptwo_q2.5',
      'xgrouptwo_q97.5',
      'sigma_xgrouptwo',
      'sigma_xgrouptwo_q2.5',
      'sigma_xgrouptwo_q97.5'
    )
    
    
    # exp the coefficitents for the sigmas intercept because its on the log scale
    # results$sigma_xgrouptwo <- exp(results$sigma_xgrouptwo  + results$sigma_Intercept)
    # results$sigma_xgrouptwo_q2.5 <- exp(results$sigma_xgrouptwo_q2.5)
    # results$sigma_xgrouptwo_q97.5 <- exp(results$sigma_xgrouptwo_q97.5)
    # 
    # results$sigma_Intercept <- exp(results$sigma_Intercept)
    # results$sigma_Intercept_q2.5 <- exp(results$sigma_Intercept_q2.5)
    # results$sigma_Intercept_q97.5 <- exp(results$sigma_Intercept_q97.5)
    
    
    
    # check for interval coverage and whether the parameter is inside the credible interval 
    results$Intercept_cover <- results$Intercept_q97.5 - results$Intercept_q2.5 < Intercept_wantcover
    results$Intercept_detect <-  pop_u > results$Intercept_q2.5 &  pop_u < results$Intercept_q97.5 
    
    
    results$xgrouptwo_cover <- results$xgrouptwo_q97.5 - results$xgrouptwo_q2.5  < xgrouptwo_wantcover
    results$xgrouptwo_detect <-  pop_effect_size_u > results$xgrouptwo_q2.5 &  pop_effect_size_u < results$xgrouptwo_q97.5 
    
    
    results$sigma_Intercept_cover <- results$sigma_Intercept_q97.5 - results$sigma_Intercept_q2.5 < sigma_Intercept_wantcover
    results$sigma_Intercept_detect <-   pop_sigma_u1 > results$sigma_Intercept_q2.5 &  pop_sigma_u1 < results$sigma_Intercept_q97.5 
    
    
    results$sigma_xgrouptwo_cover <- results$sigma_xgrouptwo_q97.5 - results$sigma_xgrouptwo_q2.5 < sigma_xgrouptwo_wantcover
    results$sigma_xgrouptwo_detect <- pop_effect_size_sigma > results$sigma_xgrouptwo_q2.5 &  pop_effect_size_sigma < results$sigma_xgrouptwo_q97.5 
    
    
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
    
    
    toplot <- as.data.frame(apply(toplot,2,round,digits = 3))
    
  }
  
  allocation_data <- rbind(allocation_data,
                           toplot)
  
}

toc()




# make adjustments to the sample size to respect the allocation and max sample size limit
allo <- paste(allocation)
allo <- str_replace_all(allo,c("c" = "",
                               "," = "",
                               "[//(]"="",
                               "[//)]" = "",
                               " " = "-"))
allocation_data$allocation <- factor(
  rep(allo,each = length(samplerange))
)

for(i in 1:nrow(allocation_data)){
  allocate_adjust <- as.numeric(unlist(str_split(allocation_data$allocation[i],"-")))
  number_batches = floor(allocation_data$samplesize[i]/(allocate_adjust[1] + allocate_adjust[2]))
  n1 = number_batches * allocate_adjust[1]
  n2 = number_batches * allocate_adjust[2]
  allocation_data$samplesize[i] = n1 + n2
}





allocation_data <- allocation_data %>% 
  mutate(
    Intercept_length = Intercept_q97.5 - Intercept_q2.5,
    sigma_Intercept_length = sigma_Intercept_q97.5 - sigma_Intercept_q2.5,
    xgrouptwo_length = xgrouptwo_q97.5 - xgrouptwo_q2.5,
    sigma_xgrouptwo_length = sigma_xgrouptwo_q97.5 - sigma_xgrouptwo_q2.5
  )

saveRDS(allocation_data,
        paste0(getwd(),"/documents/two_sample_power_files/allocation_data_problem_two.rds"))


groupone_info <- allocation_data %>% 
  ggplot(aes(x = samplesize,
             y = Intercept,
             color = allocation,
             ymin = Intercept_q2.5,
             ymax = Intercept_q97.5,
             text = paste('interval length: ', Intercept_length,
                          '</br>interval length criteria satisfied: ', Intercept_cover,
                          '</br>detect probability: ', Intercept_detect
             ))) +
  geom_hline(yintercept = c(0, .5), color = "white") +
  geom_pointrange(fatten = 1/2) +
  scale_x_discrete("reordered by the lower level of the 95% intervals", breaks = NULL) + 
  geom_point() + 
  ylim(min(allocation_data$Intercept_q2.5) - 0.5,max(allocation_data$Intercept_q97.5) + 0.5) + 
  geom_hline(yintercept = pop_u) + 
  theme_bw()

ggplotly(groupone_info)


grouptwo_info <- allocation_data %>% 
  ggplot(aes(x = samplesize,
             y = xgrouptwo,
             color = allocation,
             ymin = xgrouptwo_q2.5,
             ymax = xgrouptwo_q97.5,
             text = paste('interval length: ', xgrouptwo_length,
                          '</br>interval length criteria satisfied: ', xgrouptwo_cover,
                          '</br>detect probability: ', xgrouptwo_detect
             ))) +
  geom_hline(yintercept = c(0, .5), color = "white") +
  geom_pointrange(fatten = 1/2) +
  # scale_x_discrete("reordered by the lower level of the 95% intervals", breaks = NULL) + 
  geom_point() + 
  ylim(min(allocation_data$xgrouptwo_q2.5) - 0.5,max(allocation_data$xgrouptwo_q97.5) + 0.5) + 
  geom_hline(yintercept = pop_effect_size_u) + 
  labs(x = "Sample Size (Total)",
       y = "Location Effect Size",
       title = "Properties for Main Effect across Sample Sizes") + 
  theme_bw() 


ggplotly(grouptwo_info) %>% 
  layout(hovermode = "x")



groupone_sd_info <- allocation_data %>% 
  ggplot(aes(x = samplesize,
             y = sigma_Intercept,
             color = allocation,
             ymin = sigma_Intercept_q2.5,
             ymax = sigma_Intercept_q97.5,
             text = paste('interval length: ', sigma_Intercept_length,
                          '</br>interval length criteria satisfied: ', sigma_Intercept_cover,
                          '</br>detect probability: ', sigma_Intercept_detect
             )
  )) +
  geom_hline(yintercept = c(0, .5), color = "white") +
  geom_pointrange(fatten = 1/2) +
  scale_x_discrete("reordered by the lower level of the 95% intervals", breaks = NULL) + 
  geom_point() + 
  ylim(min(allocation_data$sigma_Intercept_q2.5) - 0.5,max(allocation_data$sigma_Intercept_q97.5) + 0.5) + 
  geom_hline(yintercept = pop_sigma_u1) + 
  theme_bw()

ggplotly(groupone_sd_info) %>% 
  layout(hovermode = "x")



grouptwo_sd_info <- allocation_data %>% 
  ggplot(aes(x = samplesize,
             y = sigma_xgrouptwo,
             color = allocation,
             ymin = sigma_xgrouptwo_q2.5,
             ymax = sigma_xgrouptwo_q97.5,
             text = paste('interval length: ', sigma_xgrouptwo_length,
                          '</br>interval length criteria satisfied: ', sigma_xgrouptwo_cover,
                          '</br>detect probability: ', sigma_xgrouptwo_detect
             )
  )) +
  geom_hline(yintercept = c(0, .5), color = "white") +
  geom_pointrange(fatten = 1/2) +
  # scale_x_discrete("reordered by the lower level of the 95% intervals", breaks = NULL) + 
  geom_point() + 
  ylim(min(allocation_data$sigma_xgrouptwo_q2.5) - 0.5,max(allocation_data$sigma_xgrouptwo_q97.5) + 0.5) + 
  geom_hline(yintercept = pop_effect_size_sigma) + 
  labs(x = "Sample Size (Total)",
       y = "Scale Effect Size",
       title = "Properties for Scale Effect for Group Two across Sample Sizes") + 
  theme_bw() 

ggplotly(grouptwo_sd_info) %>% 
  layout(hovermode = "x")

