# epidemia new analysis
# separate fits for early and late time periods

`%notin%` <- Negate(`%in%`)

library("tibble")
library("magrittr")

library(rstan)
library(epidemia)

options(mc.cores = parallel::detectCores() - 5) 
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
    date >= "2020-03-01",
    date <= "2020-06-01")


exclude <- c(
	"Alaska", "Hawaii",
	"Montana", "South Dakota", "Vermont", "Wyoming",
	"Rhode Island", "Utah", "Alabama", "Colorado", "Massachusetts", "Tennessee", "Arkansas", "Maine", "Idaho", "Michigan", "California", "Texas", "Missouri", "Mississippi", "North Carolina", "North Dakota", "Oklahoma", "Virginia")


# AL,
# AR, ID, ME, MI, MN, MO, MS, NC, ND, OK, TN, TX, VT, and VA
# based on little to no activity and 
# presence of negative daily counts
# this exclusion list will need to be adjusted for late-phase counts

states <- setdiff(as.character(data_us_f$country), exclude)

# subset to fit
args_us$group_subset <- states
# , "Rhode Island",
# "Maryland", "Delaware", "Pennsylvania",
# "Colorado", "Washington", "Oregon", "Wyoming")

# specify for epidemia
args_us$data <- as.data.frame(data_us_f)
args_us$pops <- as.data.frame(pop_df_f)

args_us$rt <- epirt(
  formula = R(country,date) ~ 0 + (1 | country) + stay + primary,
  prior = rstanarm::normal(scale = 0.5))
  # prior = shifted_gamma(shape=1/6, scale=1, shift = log(1.05)/6),
  #   - alternate prior
  # prior_tau = rstanarm::exponential(0.03)
  #   - default assumption that seeds begin 30 days before 10 cum deaths 
  #     are observed
  #   - how is 10 set?

args_us$algorithm <- "sampling"
args_us$sampling_args <- list(
  iter = 8000,
  control = list(adapt_delta = 0.8, max_treedepth = 12))

args_us$init_run <- TRUE
fit <- do.call(epim, args_us)

save(fit, file = "stay_primary_new_analsis.Rdat")