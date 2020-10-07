# make plots

library("epidemia")
library("tibble")
library("lubridate")
library("ggplot2")

load("stay_primary.RData")
load("epidemia_data.RData")

print(fit)

theme_fivethirtyeight_ALT <- function(base_size = 12, base_family = "sans") {
  colors <- deframe(ggthemes::ggthemes_data[["fivethirtyeight"]])
  (ggthemes::theme_foundation(base_size = base_size, base_family = base_family)
   + theme(
     line = element_line(colour = "black"),
     rect = element_rect(fill = colors["Light Gray"],
                         linetype = 0, colour = NA),
     text = element_text(colour = colors["Dark Gray"]),
     axis.ticks = element_blank(),
     axis.line = element_blank(),
     legend.background = element_rect(),
     legend.position = "bottom",
     legend.direction = "horizontal",
     legend.box = "vertical",
     panel.grid = element_line(colour = NULL),
     panel.grid.major =
       element_line(colour = colors["Medium Gray"]),
       panel.grid.major.x = element_blank(),
     panel.grid.minor = element_blank(),
     plot.title = element_text(hjust = 0, size = rel(1.5), face = "bold"),
     plot.margin = unit(c(1, 1, 1, 1), "lines"),
     strip.background = element_rect()))
}

# construct data frame of model estimates and observations

pp <- posterior_predict(fit)

p <- as.matrix(fit)

qs <- c(0.025,0.1,0.25,0.5,0.75,0.9,0.975)
qs_name <- c(
  "lwr_95", "lwr_90", "lwr_50", "est", "upr_50", "upr_90", "upr_95")


pp_df <- apply(pp$draws, 2, quantile, probs = qs) %>%
  t() %>% as.data.frame() %>% tibble()
colnames(pp_df) <- qs_name

pp_df$t <- pp$time
pp_df$group <- pp$group
pp_df$obs <- fit$data$deaths

level_key_l <- c("lwr_50" = "50", "lwr_95" = "95", "lwr_90" = "90")
level_key_u <- c("upr_50" = "50", "upr_95" = "95", "upr_90" = "90")

pp_df_l <- pp_df %>% tidyr::pivot_longer(
  cols = lwr_95:lwr_50, names_to = ".width", values_to = ".lower") %>%
  dplyr::select(t, group, obs, est, .width, .lower)
pp_df_l$.width <- dplyr::recode(pp_df_l$.width, !!!level_key_l)

pp_df_u <- pp_df %>% tidyr::pivot_longer(
  cols = upr_50:upr_95, names_to = ".width", values_to = ".upper") %>%
  dplyr::select(t, group, .width, .upper)
pp_df_u$.width <- dplyr::recode(pp_df_u$.width, !!!level_key_u)

pp_df_g <- dplyr::left_join(pp_df_l, pp_df_u, by = c("t", "group", ".width"))


# dataframe for R-value
Rt <- posterior_rt(fit)

rt_df <- apply(Rt$draws, 2, quantile, probs = qs) %>%
  t() %>% as.data.frame() %>% tibble()

colnames(rt_df) <- qs_name
rt_df$t <- Rt$time
rt_df$group <- Rt$group

rt_df_l <- rt_df %>% tidyr::pivot_longer(
  cols = lwr_95:lwr_50, names_to = ".width", values_to = ".lower") %>%
  dplyr::select(t, group, est, .width, .lower)
rt_df_l$.width <- dplyr::recode(rt_df_l$.width, !!!level_key_l)

rt_df_u <- rt_df %>% tidyr::pivot_longer(
  cols = upr_50:upr_95, names_to = ".width", values_to = ".upper") %>%
  dplyr::select(t, group, .width, .upper)
rt_df_u$.width <- dplyr::recode(rt_df_u$.width, !!!level_key_u)

rt_df_g <- dplyr::left_join(rt_df_l, rt_df_u, by = c("t", "group", ".width"))

data_us <- data_us %>% dplyr::select(date, country, deaths, primary, stay)

# one state per plot
p_list <- as.character(unique(pp$group))
plte <- "Blues"

for (u in 1:length(p_list)) {
  p_dte <- data_us %>%
    dplyr::filter(country == p_list[u], primary == 1) %>%
    dplyr::pull(date)

  stay_dtes <- data_us %>%
    dplyr::filter(country == p_list[u], stay == 1) %>%
    dplyr::pull(date)

  if (length(p_dte) == 0) {
    p_dte <- as.Date("2020-01-01")
  }

  if(length(stay_dtes) == 0) {
    stay_dte <- as.Date("2021-01-01")
  } else {
    stay_dte <- min(stay_dtes[[1]])
  }

  # table of interventions
  interv <- tibble(
    intervention = c(
      "primary", "median time-to-death", "stay at home order"),
    value = c(p_dte[[1]], p_dte[[1]] + 21, stay_dte),
    clr = c("goldenrod4", "firebrick", "dodgerblue3"))

  pl_obs <- pp_df_g %>% dplyr::filter(group == p_list[u]) %>%
      ggplot(aes(x=t)) +
      ggdist::geom_lineribbon(
        aes(y = est, ymin = .lower, ymax = .upper), alpha = 0.7) +
      scale_fill_brewer(type = "seq", palette = plte) +
      geom_point(aes(y=obs)) +
      geom_vline(
        data = interv, aes(xintercept = value, color = intervention)) +
      scale_color_manual(values = interv$clr) +
      ylab("Daily deaths") + xlab("Date") +
      ggtitle(
      paste0(p_list[u], " Mortality")) +
      theme_fivethirtyeight_ALT()

  interv_rt <- tibble(
    intervention = c(
      "primary", "stay at home order"), 
    value = c(p_dte[[1]], stay_dte),
    clr = c("goldenrod4", "blue"))


  pl_rt <- rt_df_g %>% dplyr::filter(group == p_list[u]) %>%
    ggplot(aes(x=t)) +
    ggdist::geom_lineribbon(
      aes(
        y = est, ymin = .lower, ymax =.upper), 
        alpha = 0.7, show.legend = TRUE) +
    scale_fill_brewer(type = "seq", palette = plte) +
    geom_vline(
      data = interv_rt, aes(xintercept = value, color = intervention)) +
    scale_colour_brewer(type = "qual", palette = "Set1") +
    ylab(expression("R"[t])) + xlab("Date") +
    ggtitle(
      paste0(p_list[u], " Transmission")) +
    theme_fivethirtyeight_ALT()

  cmb <- ggpubr::ggarrange(pl_obs, pl_rt, common.legend = TRUE)


  ggsave(
    file = paste0("epi_plt/" ,p_list[u], "_comb.png"),
    cmb, width = 10, height = 7, bg = "gray94")

}