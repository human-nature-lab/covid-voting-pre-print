# functions to work with the Imai method, 
# data organization and results plotting

### functions to organize data
# add_pch <- function(df, cty_df) {
#   # if we want to assume the outcome is log-normal
#   # not reasonable!
#   counties <- unique(cty_df$fips)
#   df$pch_death <- NA
#   df$pch_cases <- NA
#   for (c in counties) {
#     df$pch_death[df$fips == c] <- c(
#       log(df$deathsCum[df$fips == c] + 0.001)[1],
#       diff(log(df$deathsCum[df$fips == c] + 0.001))
#     )
#     df$pch_cases[df$fips == c] <- c(
#       log(df$casesCum[df$fips == c] + 0.001)[1],
#       diff(log(df$casesCum[df$fips == c] + 0.001))
#     )
#   }
#   return(df)
# }

# mk_prim <- function(df, cty_df, lwr = 1, upr = 30, mail_only) {
#   # make primary indicator variable
#   # relevant range after the primary takes place, 
#   # guided by the incubation (5.2) and infection-to-death (17.8) periods
#   df$primary_on <- 0
#   s_list <- unique(df$state)
#   s_list <- setdiff(s_list, mail_only)

#   for (f in 1:length(s_list)) {
#     p <- cty_df$primary[cty_df$state == s_list[f]][1]
#     p_ste <- s_list[f]
#     if (p_ste %in% mail_only) {
#       p <- NA
#     }
#     if (is.na(p) == FALSE) {
#       rnge <- seq(as.Date(p + lwr)[1], as.Date(p + upr)[1], by = 1)
#     }
#     df$primary_on[df$date %in% rnge & df$state == p_ste] <- 1
#   }
#   return(df)
# }

# mk_prim_mk <- function(df, cty_df, delay = 0, mail_only) {
#   # make primary indicator variable for single day
#   # default primary day, optional to specify delay (e.g. 30 days after)
#   df$primary_mk <- 0
#   df$primary_s <- 0
#   s_list <- unique(cty_df$state)
#   s_list <- setdiff(s_list, mail_only)

#   for (f in 1:length(s_list)) {
#     p <- cty_df$primary[cty_df$state == s_list[f]][1] + delay
#     p_ste <- s_list[f]
#     if (p_ste %in% mail_only) {
#       p <- NA
#     }
#     if (is.na(p) == FALSE) {
#       df$primary_mk[df$date == p & df$state == p_ste] <- 1
#       df$primary_s[df$date > p & df$state == p_ste] <- 1
#     }
#   }
#   return(df)
# }

# mk_lck <- function(s_df, df, delay = 0) {
#   # stay at home order
#   # mk: date of + delay
#   # s: stable, on after + delay
#   df$lock_mk <- 0
#   df$lock_s <- 0
#   s_list <- unique(s_df$fips)

#   for (f in 1:length(s_list)) {
#     full <- usmap::fips_info(s_list[f])$full
#     ste <- s_df[s_df$fips == s_list[f],]
#     stay <- ste$StayAtHome + delay

#     if (is.na(stay) == FALSE) {
#       df$lock_mk[df$date == stay & df$state == full] <- 1
#       df$lock_s[df$date > stay & df$state == full] <- 1
#     }
#   }
#   return(df)
# }

indicator_mk <- function(ste_df, 
                         time_df, 
                         delay = 0, 
                         window_open = 13, 
                         window_close = 32, 
                         exclude_primary) {
  `%notin%` <- Negate(`%in%`)

  sl <- unique(ste_df$fips)*1000
  time_df$primary_day <- 0
  time_df$primary_period <- 0
  time_df$primary_window <- 0

  time_df$lock_day <- 0
  time_df$lock_period <- 0
  time_df$lock_window <- 0
  excl_fips <- as.integer(usmap::fips(exclude_primary))

  for (s in 1:length(sl)) {
    slw <- sl[s]
    sup <- slw + 999
    
    # primary
    if (as.integer(slw/1000) %in% excl_fips == TRUE) {
      p <- NA
    } else {
      p <- ste_df$primary_dte[ste_df$fips == slw/1000]
    }

    if(is.na(p) == FALSE) {
      p_d <- p + delay # 20% of time-to-death distribution or < that value by a few days (10th day)
      p_wo <- p + window_open # 20% of time-to-death distribution
      p_wc <- p + window_close # 80% of time-to-death distribution

      ind <- time_df$fips %in% (slw:sup)
      time_df$primary_day[time_df$fips %in% (slw:sup) & time_df$date == p] <- 1
      time_df$primary_period[time_df$fips %in% (slw:sup) & time_df$date >= p_d] <- 1
      time_df$primary_window[time_df$fips %in% (slw:sup) & time_df$date >= p_wo & time_df$date <= p_wc] <- 1
    }

    # stay at home
    l <- ste_df$StayAtHome[ste_df$fips == as.integer(slw/1000)]
    if (is.na(l) == FALSE) {
      l_d <- l + delay # 20% of time-to-death distribution or < that value by a few days (10th day)
      l_wo <- l + window_open # 20% of time-to-death distribution
      l_wc <- l + window_close # 80% of time-to-death distribution

      ind <- time_df$fips %in% (slw:sup)
      time_df$lock_day[time_df$fips %in% (slw:sup) & time_df$date == l] <- 1
      time_df$lock_period[time_df$fips %in% (slw:sup) & time_df$date >= l_d] <- 1
      time_df$lock_window[time_df$fips %in% (slw:sup) & time_df$date >= l_wo & time_df$date <= l_wc] <- 1
    }
  }
  return(time_df)
}

other_indicators <- function(df, orig) {
  # date of first case, death
  # not sure how to handle lockdown
  # assign very large positive value?
  # converts to running scale

  df$first_case <- as.Date(NA)
  df$first_death <- as.Date(NA)
  c_list <- unique(df$fips)
  for (s in 1:length(c_list)) {
    co <- c_list[s]
    ## cases
    cm <- min(df$date[df$fips == co & df$casesCum > 0], na.rm = TRUE)
    df$first_case[df$fips == co] <- as.Date(cm)
    ## deaths
    cd <- min(df$date[df$fips == co & df$deathsCum > 0], na.rm = TRUE)
    df$first_death[df$fips == co] <- as.Date(cd)
  }

  df$first_case <- as.integer(df$first_case - orig)
  df$first_death <- as.integer(df$first_death - orig)

  return(df)
}


# check if this was right
#   for (f in 1:length(f_list)) {
#     C <- cty_df[cty_df$state == s_list[f], ]

#     # range from 1 to 30 days
#     if (is.na(as.Date(C$primary)[1]) == FALSE) {
#       rnge <- seq(as.Date(C$primary + lwr)[1], as.Date(C$primary + upr)[1], by = 1)
#       s_nme <- usmap::fips_info(usmap::fips(s_list[f]))$full
#     } else {
#       rnge <- NA
#     }
#       df$primary_on[df$date %in% rnge & df$state == s_nme] <- 1
#   }
#   return(df)
# }

### functions to process and plot results
procres <- function(PE.r, death_window, model_name) {
  summ <- as.data.frame(summary(PE.r)$summary)
  summ$time <- row.names(summ)
  summ$t <- 1:nrow(summ)
  dw <- death_window - 10
  summ$window <- (summ$t >= dw[1] & summ$t <= dw[2])*1
  summ$model <- model_name
  return(summ)
}


plotres <- function(pe.summ, 
                    fname = "",
                    title = "Average treatment effect (on the treated)",
                    subtitle = "", 
                    ylab = "death rate per 1000 people",
                    caption = "The shaded region corresponds to middle 80% of the time-to-death distribution.",
                    death_window = 0) {
  if(length(death_window) < 3) {
    death_window <- c(13, 21, 32)
  }

  mn <- death_window[3] - 10
  mx <- death_window[1] - 10

  p_deaths_add <- pe.summ %>% 
    ggplot(aes(x = t)) +
      geom_rect(
          aes(
            xmin = mn, 
            xmax = mx, 
            ymin = -Inf, 
            ymax = Inf), 
          fill = "pink",
          alpha = 0.02) +
      geom_point(aes(y = estimate)) +
      geom_linerange(
        aes(ymin = `2.5%`, ymax = `97.5%`), alpha = 0.5) +
      scale_color_brewer(type = "div", palette = "Spectral") +
      geom_hline(yintercept = 0, color = 'darkblue') +
      #theme_classic() +
      theme_fivethirtyeight_ALT() +
      scale_x_continuous(
        breaks = seq(
          min(pe.summ$t), 
          max(pe.summ$t), 
          by = 1)) +
      ylab(ylab) + xlab("Days since 10 days after the primary") +
      labs(
        title = title, 
        subtitle = subtitle,
        caption = caption
    )
  if (fname != "") {
    ggsave(filename = fname,
      plot = p_deaths_add,
      height = 6,
      width = 6
      #bg = "transparent"
    )
  }
  return(p_deaths_add)
}


plotres_m <- function(pe.summ, 
                      fname = "results_pch_deaths_add.png",
                      title = "Average treatment effect (on the treated)",
                      subtitle = "", 
                      ylab = "death rate per 1000 peopls",
                      caption = "The pink region corresponds to middle 80% of the time-to-death distribution.",
                      death_window) {

  mn <- death_window[3] - 10
  mx <- death_window[1] - 10

  p_deaths_add <- pe.summ %>% 
    ggplot(aes(x = t)) +
      geom_rect(
          aes(
            xmin = mn, 
            xmax = mx, 
            ymin = -Inf, 
            ymax = Inf), 
          fill = "pink",
          alpha = 0.008) +
      geom_point(aes(y = estimate, color = model)) +
      geom_linerange(
        aes(ymin = `2.5%`, ymax = `97.5%`, color = model), alpha = 0.5) +
      scale_color_brewer(type = "div", palette = "Spectral") +
      geom_hline(yintercept = 0, color = 'darkblue') +
      #theme_classic() +
      theme_fivethirtyeight_ALT() +
      # theme(
      #   panel.background = element_rect(fill = "transparent", colour = NA),
      #   panel.grid.minor = element_blank(), 
      #   panel.grid.major = element_blank(),
      #   plot.background = element_rect(fill = "transparent", colour = NA),
      #   axis.line = element_line(colour = "black")
      # ) +
      scale_x_continuous(
        breaks = seq(
          min(pe.summ$t), 
          max(pe.summ$t), 
          by = 1)) +
      ylab(ylab) + xlab("Days since 10 days after the primary") +
      labs(
        title = title, 
        subtitle = subtitle,
        caption = caption
    )

  if (fname != "") {
    ggsave(filename = fname,
      plot = p_deaths_add,
      height = 6,
      width = 6
      #bg = "transparent"
    )
  }
  return(p_deaths_add)
}

CB_plot <- function(PM.r,
                    data,
                    covariates,
                    fname = "cb_s.png",
                    tt = "Covariate Balance") {
  CB.add <- get_covariate_balance(
  PM.r$att, 
  as.data.frame(data),
  covariates = covariates)

  CB.out <- tibble(as.data.frame(CB.add))

  t <- as.integer(
    stringr::str_split(
      rownames(CB.add), 
      pattern = "_", 
      simplify = TRUE)[,2]) * -1

  CB.out$t <- t

  CB.out %<>% tidyr::pivot_longer(
    -t, 
    names_to = "metric",
    values_to = "value")

  CB.out %<>% mutate(
    lt = case_when(
      metric == "cum_death_rte" ~ "solid",
      metric == "cum_case_rte" ~ "solid",
      metric != "cum_death_rte" ~ "dashed"),
    type = case_when(
      metric == "cum_case_rte" ~ "time-varying", 
      metric == "cum_death_rte" ~ "time-varying", 
      metric != "cum_death_rte" ~ "constant")
  )

  # cb.pl <- CB.out %>% 
  #   ggplot(aes(x = t)) +
  #     geom_line(aes(y = value, color = metric, linetype = type)) +
  #     scale_linetype_manual(values = CB.out$lt) +
  #     #ggthemes::scale_color_fivethirtyeight() +
  #     scale_color_brewer(palette = "Set3", type = "qual") +
  #     ggtitle("Refined covariate balance") +
  #     xlab("Days before treatment onset") +
  #     ylab("Standard difference") +
  #     scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
  #     ggthemes::theme_fivethirtyeight()

  CB.add.prerefine <- get_covariate_balance(
    PM.r$att, 
    as.data.frame(data),
    covariates = covariates,
    use.equal.weights = TRUE)

  CB.out.pr <- tibble(as.data.frame(CB.add.prerefine))

  CB.out.pr$t <- as.integer(
    stringr::str_split(
      rownames(CB.add.prerefine), 
      pattern = "_", 
      simplify = TRUE)[,2]) * -1

  CB.out.pr %<>% tidyr::pivot_longer(
    -t, 
    names_to = "metric",
    values_to = "value")

  CB.out.pr %<>% mutate(
    lt = case_when(
      metric == "cum_case_rte" ~ "solid",
      metric == "cum_death_rte" ~ "solid", 
      metric != "cum_death_rte" ~ "dashed"),
    type = case_when(
      metric == "cum_case_rte" ~ "time-varying", 
      metric == "cum_death_rte" ~ "time-varying", 
      metric != "cum_death_rte" ~ "constant")
  )

  # cb.pl.pr <- CB.out.pr %>% 
  #   ggplot(aes(x = t)) +
  #     geom_line(aes(y = value, color = metric, linetype = type)) +
  #     scale_linetype_manual(values = CB.out.pr$lt) +
  #     #ggthemes::scale_color_fivethirtyeight() +
  #     scale_color_brewer(palette = "Set3", type = "qual") +
  #     ggtitle("Pre-refinement covariate balance") +
  #     xlab("Days before treatment onset") +
  #     ylab("Standard difference") +
  #     scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
  #     ggthemes::theme_fivethirtyeight()

  # PL <- ggpubr::ggarrange(cb.pl.pr, cb.pl, 
  #         labels = c("A", "B"),
  #         ncol = 2, nrow = 1)

  CB.out$Refinement<- "Post"
  CB.out.pr$Refinement <- "Pre"

  CB.out.comb <- CB.out.pr %>% add_row(CB.out)
  CB.out.comb$Refinement <- factor(
    CB.out.comb$Refinement,
    levels = c("Pre", "Post"))

  cb.pl.comb <- CB.out.comb %>% 
    ggplot(aes(x = t)) +
      geom_line(aes(y = value, color = metric, linetype = type)) +
      scale_linetype_manual(values = CB.out.comb$lt) +
      #ggthemes::scale_color_fivethirtyeight() +
      scale_color_brewer(palette = "Set3", type = "qual") +
      scale_y_continuous(breaks = seq(-1, 1, 0.1)) +
      ggthemes::theme_fivethirtyeight() +
      facet_wrap(Refinement ~ .) +
      ggtitle(tt) +
      xlab("Days before treatment onset") +
      ylab("Standard difference")

  if (fname != "") {
    ggsave(filename = fname,
      plot = cb.pl.comb
      #bg = "transparent"
    )
  }

  return(cb.pl.comb)
}


match_frame <- function(pm = PM.s$att) {
  # create data-frame of matches
  # then filter it on mtch to get the candidates for county-of-interest
  for (k in 1:length(pm)) {
    nme <- stringr::str_split(names(pm[k]), "\\.", simplify = TRUE)
    loc <- tibble(usmap::fips_info(nme[1]))
    mt <- tibble(usmap::fips_info(pm[[k]]))
    if (length(mt) == 0 | nrow(mt) == 0) {
      next
    } else {
      mt$weight <- attr(pm[[k]], "weights")
      mt$trt <- 0
      loc$weight <- NA
      loc$trt <- 1
      mt %<>% add_row(loc)
      mt$match <- loc$fips
      mt$suffix <- nme[2]
      if (k == 1) {
        mf <- mt
      } else {
        mf %<>% add_row(mt)
      }
    }
  }
  return(mf)
}


mk_toss <- function(dFr.T, tr_q, quant = 3) {
  # toss out treatments based on turnout quantiles
  # toss indicates that the treated county data starting at 
  # the treatment window should be removed
  # tossL indicates the primary-held-counties that are below the threshold
  # so, filtering out toss=1 yields control and treated-high turnout
  # filtering out tossL=1 yields control and treated-low turnout
  s_list <- unique(dFr.T$state)
  dFr.T$toss <- 0
  dFr.T$tossL <- 0

  for (s in 1:length(s_list)) {
    dFr.T[dFr.T$state == s_list[s],]
    p_held <- any(dFr.T$primary_window[dFr.T$state == s_list[s]] == 1)
    if (p_held == TRUE) {
      c_list <- unique(dFr.T$fips[dFr.T$state == s_list[s]])
      for (co in 1:length(c_list)) {
        # toss if the primary-county is below p-tile threshold
        if (dFr.T$trn_rte[dFr.T$fips == c_list[co]][1] < tr_q[quant] | is.na(dFr.T$trn_rte[dFr.T$fips == c_list[co]][1])) {
          mw <- min(dFr.T$running[dFr.T$fips == c_list[co] & dFr.T$primary_window == 1])
          dFr.T$toss[dFr.T$fips == c_list[co] & dFr.T$running >= mw] <- 1
        }
        # toss if the primary-county is above or equal to p-tile threshold
        if (dFr.T$trn_rte[dFr.T$fips == c_list[co]][1] >= tr_q[quant] | is.na(dFr.T$trn_rte[dFr.T$fips == c_list[co]][1])) {
          mw <- min(dFr.T$running[dFr.T$fips == c_list[co] & dFr.T$primary_window == 1])
          dFr.T$tossL[dFr.T$fips == c_list[co] & dFr.T$running >= mw] <- 1
        }
      }
    }
  }
  return(dFr.T)
}



theme_fivethirtyeight_ALT <- function(base_size = 12, base_family = "sans") {
  colors <- deframe(ggthemes::ggthemes_data[["fivethirtyeight"]])
  (ggthemes::theme_foundation(base_size = base_size, base_family = base_family)
   + theme(
     line = element_line(colour = "black"),
     rect = element_rect(fill = colors["Light Gray"],
                         linetype = 0, colour = NA),
     text = element_text(colour = colors["Dark Gray"]),
     #axis.title = element_blank(),
     #axis.text = element_text(),
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
     # unfortunately, can't mimic subtitles TODO!
     plot.title = element_text(hjust = 0, size = rel(1.5), face = "bold"),
     plot.margin = unit(c(1, 1, 1, 1), "lines"),
     strip.background = element_rect()))
}