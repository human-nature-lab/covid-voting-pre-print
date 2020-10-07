`%notin%` <- Negate(`%in%`)

library("tibble")
library("magrittr")

library(rstan)
library(epidemia)

options(mc.cores = parallel::detectCores()) 
rstan_options(auto_write = TRUE)

args_us <- EuropeCovid

load("epidemia_data.Rdat")

# general state filter
data_us_f <- data_us %>% 
  dplyr::filter(
    country != "Northern Mariana Islands", 
    country != "Virgin Islands",
    country != "Puerto Rico",
    country != "Nebraska")

pop_df_f <- pop_df %>% 
  dplyr::filter(
    country != "Northern Mariana Islands", 
    country != "Virgin Islands",
    country != "Puerto Rico",
    country != "Nebraska")

# time filter
data_us_f <- data_us_f %>% 
  dplyr::filter(
    date >= "2020-03-15",
    date <= "2020-05-21")


exclude <- c(
	"Alaska", "Hawaii",
	"Montana", "South Dakota", "Vermont", "Wyoming",
	"Rhode Island", "Utah", "Alabama", "Colorado")

states <- setdiff(as.character(data_us_f$country), exclude)

# subset to fit
args_us$group_subset <- states

# specify for epidemia
args_us$data <- as.data.frame(data_us_f)
args_us$pops <- as.data.frame(pop_df_f)

args_us$rt <- epirt(
  formula = R(country,date) ~ 0 + (1 | country) + stay + primary,
  prior = rstanarm::normal(scale = 0.5))

args_us$algorithm <- "sampling"
args_us$sampling_args <- list(
  iter = 12000,
  control = list(adapt_delta = 0.95, max_treedepth = 15))

args_us$init_run <- TRUE
fit <- do.call(epim, args_us)

save(fit, file = "stay_primary.Rdat")