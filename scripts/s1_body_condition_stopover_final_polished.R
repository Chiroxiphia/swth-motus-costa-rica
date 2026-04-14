# S1 body condition vs initial minimum stopover duration at tagging sites
# Final polished script
#
# What this script does:
# - lets you select the banding Excel and the Motus Excel
# - standardizes site names to Finca Cántaros and Rancho Quemado
# - calculates body condition at capture as residual body mass relative to wing chord
#   (primary analysis) and tarsus length (sensitivity analysis)
# - calculates initial minimum stopover duration from the FIRST continuous post-capture
#   detection bout at the focal receiver associated with the tagging site
# - shows the figures in the RStudio Plots pane
# - saves the PNGs, summary table, and model summary
# - opens the main PNG automatically at the end
#
# Recommended supplementary figure:
#   s1_primary_wing_initial_stopover.png

required_pkgs <- c("readxl", "dplyr", "ggplot2", "lubridate", "stringr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)
library(stringr)

# -------------------------------------------------------------------
# OPTIONAL: set these manually if you do not want interactive selection
# -------------------------------------------------------------------
banding_file <- NULL
motus_file <- NULL
out_dir <- NULL

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
safe_character <- function(x) if (is.factor(x)) as.character(x) else as.character(x)

find_col <- function(df, choices) {
  nm <- names(df)
  hit <- nm[nm %in% choices]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

pick_file <- function(existing_path = NULL, label = "Select a file") {
  if (!is.null(existing_path) && nzchar(existing_path) && file.exists(existing_path)) {
    return(existing_path)
  }
  if (interactive()) {
    cat("\n", label, "\n", sep = "")
    return(file.choose())
  }
  stop(paste(label, "- set the path manually at the top of the script."))
}

pick_out_dir <- function(existing_dir = NULL, fallback_file = NULL) {
  if (!is.null(existing_dir) && nzchar(existing_dir) && dir.exists(existing_dir)) {
    return(existing_dir)
  }
  if (interactive()) {
    typed <- readline(
      prompt = paste0(
        "Output folder path (press Enter to use the same folder as ",
        basename(fallback_file), "): "
      )
    )
    if (nzchar(typed) && dir.exists(typed)) return(typed)
  }
  dirname(fallback_file)
}

parse_capture_datetime <- function(date_vec, time_vec) {
  x <- paste(safe_character(date_vec), safe_character(time_vec))
  out <- suppressWarnings(mdy_hm(x, tz = "America/Costa_Rica"))
  if (all(is.na(out))) {
    out <- suppressWarnings(parse_date_time(
      x,
      orders = c(
        "mdY HM", "mdY HMS", "Ymd HM", "Ymd HMS",
        "dmY HM", "dmY HMS", "mdy HM", "mdy HMS"
      ),
      tz = "America/Costa_Rica"
    ))
  }
  out
}

standardize_site <- function(site, location = NA_character_) {
  x <- str_to_lower(str_squish(paste(site, location)))
  case_when(
    str_detect(x, "finca") ~ "Finca Cántaros",
    str_detect(x, "cataros") ~ "Finca Cántaros",
    str_detect(x, "cantaros") ~ "Finca Cántaros",
    str_detect(x, "osa") ~ "Rancho Quemado",
    str_detect(x, "casa bontanica") ~ "Rancho Quemado",
    str_detect(x, "rancho") ~ "Rancho Quemado",
    str_detect(x, "quemado") ~ "Rancho Quemado",
    TRUE ~ str_squish(safe_character(site))
  )
}

receiver_pattern_for_site <- function(site) {
  dplyr::case_when(
    site == "Finca Cántaros" ~ "finca|c[aá]ntaros|cataros",
    site == "Rancho Quemado" ~ "rancho|quemado|osa|birds",
    TRUE ~ ""
  )
}

haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  6371 * 2 * atan2(sqrt(a), sqrt(1 - a))
}

build_condition_residuals <- function(df, size_var) {
  keep <- df %>% filter(!is.na(weight_g), !is.na(.data[[size_var]]))
  fit <- lm(reformulate(size_var, response = "weight_g"), data = keep)

  keep %>%
    mutate(
      condition_metric = paste0("residual mass ~ ", size_var),
      body_condition = resid(fit)
    )
}

choose_focal_receiver <- function(det_df, bird_row) {
  bird_row <- as.list(bird_row)
  site_regex <- receiver_pattern_for_site(bird_row$capture_site)
  bird_lat <- bird_row$tagDepLat
  bird_lon <- bird_row$tagDepLon

  cand <- det_df %>%
    filter(tagDeployID == bird_row$tagDeployID, motusTagID == bird_row$motusTagID)

  if (nrow(cand) == 0) return(NA_character_)

  if (!is.na(site_regex) && nzchar(site_regex)) {
    by_name <- cand %>%
      filter(str_detect(str_to_lower(recvDeployName), site_regex))
    if (nrow(by_name) > 0) {
      return(
        by_name %>%
          count(recvDeployName, sort = TRUE) %>%
          slice(1) %>%
          pull(recvDeployName)
      )
    }
  }

  if (!is.na(bird_lat) && !is.na(bird_lon) &&
      "recvDeployLat" %in% names(cand) && "recvDeployLon" %in% names(cand)) {

    cand2 <- cand %>%
      mutate(dist_km = haversine_km(bird_lat, bird_lon, recvDeployLat, recvDeployLon)) %>%
      filter(!is.na(dist_km))

    if (nrow(cand2) > 0) {
      return(
        cand2 %>%
          group_by(recvDeployName) %>%
          summarise(
            min_dist_km = min(dist_km, na.rm = TRUE),
            n_hits = n(),
            .groups = "drop"
          ) %>%
          arrange(min_dist_km, desc(n_hits)) %>%
          slice(1) %>%
          pull(recvDeployName)
      )
    }
  }

  cand %>%
    count(recvDeployName, sort = TRUE) %>%
    slice(1) %>%
    pull(recvDeployName)
}

calculate_tagging_site_stopover <- function(birds_df, det_df,
                                            max_gap_days = 5,
                                            max_bout_days = 60) {

  birds_with_receiver <- birds_df %>%
    rowwise() %>%
    mutate(focal_receiver = choose_focal_receiver(det_df, cur_data())) %>%
    ungroup()

  det_joined <- det_df %>%
    inner_join(
      birds_with_receiver %>%
        select(tagDeployID, motusTagID, capture_site, capture_datetime, focal_receiver),
      by = c("tagDeployID", "motusTagID")
    ) %>%
    filter(!is.na(focal_receiver)) %>%
    filter(recvDeployName == focal_receiver) %>%
    filter(ts >= capture_datetime) %>%
    mutate(detect_date = as.Date(with_tz(ts, tzone = "America/Costa_Rica"))) %>%
    distinct(tagDeployID, motusTagID, detect_date, .keep_all = TRUE) %>%
    arrange(tagDeployID, motusTagID, detect_date)

  out <- det_joined %>%
    group_by(tagDeployID, motusTagID) %>%
    group_modify(~{
      x <- .x %>% arrange(detect_date)

      if (nrow(x) == 0) {
        return(tibble(
          focal_receiver = NA_character_,
          first_detect = as.Date(NA),
          last_detect = as.Date(NA),
          n_detection_days = 0,
          stopover_days = NA_real_
        ))
      }

      first_date <- min(x$detect_date, na.rm = TRUE)

      x <- x %>%
        filter(detect_date <= first_date + max_bout_days)

      gaps <- c(0, as.numeric(diff(x$detect_date)))
      break_idx <- which(gaps > max_gap_days)[1]

      if (!is.na(break_idx)) {
        x <- x[1:(break_idx - 1), ]
      }

      tibble(
        focal_receiver = x$focal_receiver[1],
        first_detect = min(x$detect_date, na.rm = TRUE),
        last_detect = max(x$detect_date, na.rm = TRUE),
        n_detection_days = n_distinct(x$detect_date),
        stopover_days = as.numeric(max(x$detect_date, na.rm = TRUE) - min(x$detect_date, na.rm = TRUE))
      )
    }) %>%
    ungroup() %>%
    filter(n_detection_days >= 2)

  out
}

make_plot <- function(plot_df, x_label, out_file, show_plot = TRUE) {

  if (nrow(plot_df) < 3) {
    stop("Too few observations to draw the plot.")
  }

  ymax <- max(plot_df$stopover_days, na.rm = TRUE)
  ymax <- ceiling(ymax * 1.10)

  p <- ggplot(plot_df, aes(x = body_condition, y = stopover_days, color = capture_site)) +
    geom_point(size = 3, alpha = 0.9) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      color = "black",
      linewidth = 0.9
    ) +
    scale_color_manual(values = c(
      "Finca Cántaros" = "#1b9e77",
      "Rancho Quemado" = "#d95f02"
    )) +
    scale_y_continuous(
      limits = c(0, ymax),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(
      title = "Body condition at capture vs. initial minimum stopover duration",
      x = x_label,
      y = "Initial minimum stopover duration at tagging site (days)",
      color = "Capture site"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )

  if (show_plot) print(p)
  ggsave(out_file, p, width = 8.2, height = 5.6, dpi = 320)
  invisible(p)
}

# -------------------------------------------------------------------
# Select files
# -------------------------------------------------------------------
banding_file <- pick_file(banding_file, "Select the BANDING Excel file")
motus_file <- pick_file(motus_file, "Select the MOTUS Excel file")
out_dir <- pick_out_dir(out_dir, banding_file)

cat("\nBanding file:", banding_file, "\n")
cat("Motus file:", motus_file, "\n")
cat("Output folder:", out_dir, "\n\n")

# -------------------------------------------------------------------
# Read banding data
# -------------------------------------------------------------------
banding_raw <- read_excel(banding_file, skip = 1)

band_cols <- list(
  band = find_col(banding_raw, c("Band #")),
  date = find_col(banding_raw, c("Date")),
  time = find_col(banding_raw, c("Time")),
  site = find_col(banding_raw, c("Site")),
  location = find_col(banding_raw, c("Location")),
  sex = find_col(banding_raw, c("Sex")),
  age = find_col(banding_raw, c("Age")),
  fat = find_col(banding_raw, c("Fat score")),
  wing = find_col(banding_raw, c("Wing cord")),
  tarsus = find_col(banding_raw, c("Tarsus length")),
  weight = find_col(banding_raw, c("Weight")),
  motusTagID = find_col(banding_raw, c("motusTagID")),
  tagDeployID = find_col(banding_raw, c("tagDeployID"))
)

needed_band <- c("date", "time", "site", "weight", "wing", "tarsus", "motusTagID", "tagDeployID")
missing_band <- names(band_cols)[needed_band %in% names(band_cols) & is.na(unlist(band_cols[needed_band]))]
if (length(missing_band) > 0) {
  stop(paste("Missing required banding columns:", paste(missing_band, collapse = ", ")))
}

birds <- tibble(
  band_number  = if (!is.na(band_cols$band)) safe_character(banding_raw[[band_cols$band]]) else NA_character_,
  date_raw     = safe_character(banding_raw[[band_cols$date]]),
  time_raw     = safe_character(banding_raw[[band_cols$time]]),
  site_raw     = safe_character(banding_raw[[band_cols$site]]),
  location_raw = if (!is.na(band_cols$location)) safe_character(banding_raw[[band_cols$location]]) else NA_character_,
  sex          = if (!is.na(band_cols$sex)) safe_character(banding_raw[[band_cols$sex]]) else NA_character_,
  age          = if (!is.na(band_cols$age)) safe_character(banding_raw[[band_cols$age]]) else NA_character_,
  fat_score    = if (!is.na(band_cols$fat)) safe_numeric(banding_raw[[band_cols$fat]]) else NA_real_,
  wing_mm      = safe_numeric(banding_raw[[band_cols$wing]]),
  tarsus_mm    = safe_numeric(banding_raw[[band_cols$tarsus]]),
  weight_g     = safe_numeric(banding_raw[[band_cols$weight]]),
  motusTagID   = safe_numeric(banding_raw[[band_cols$motusTagID]]),
  tagDeployID  = safe_numeric(banding_raw[[band_cols$tagDeployID]])
) %>%
  mutate(
    capture_site = standardize_site(site_raw, location_raw),
    capture_datetime = parse_capture_datetime(date_raw, time_raw),
    julian_capture = yday(capture_datetime)
  ) %>%
  filter(!is.na(motusTagID), !is.na(tagDeployID)) %>%
  distinct(tagDeployID, motusTagID, .keep_all = TRUE)

cat("Banding rows with valid Motus IDs:", nrow(birds), "\n")
cat("Standardized sites:", paste(sort(unique(birds$capture_site)), collapse = ", "), "\n")
cat("Non-missing wing:", sum(!is.na(birds$wing_mm)),
    "| tarsus:", sum(!is.na(birds$tarsus_mm)),
    "| weight:", sum(!is.na(birds$weight_g)),
    "| fat score:", sum(!is.na(birds$fat_score)), "\n\n")

# -------------------------------------------------------------------
# Read Motus data
# -------------------------------------------------------------------
motus_raw <- read_excel(motus_file)

motus_cols <- list(
  ts = find_col(motus_raw, c("ts")),
  motusTagID = find_col(motus_raw, c("motusTagID")),
  tagDeployID = find_col(motus_raw, c("tagDeployID")),
  recvDeployName = find_col(motus_raw, c("recvDeployName")),
  recvDeployLat = find_col(motus_raw, c("recvDeployLat")),
  recvDeployLon = find_col(motus_raw, c("recvDeployLon")),
  tagDepLat = find_col(motus_raw, c("tagDepLat")),
  tagDepLon = find_col(motus_raw, c("tagDepLon")),
  motusFilter = find_col(motus_raw, c("motusFilter"))
)

needed_motus <- c("ts", "motusTagID", "tagDeployID", "recvDeployName")
missing_motus <- names(motus_cols)[needed_motus %in% names(motus_cols) & is.na(unlist(motus_cols[needed_motus]))]
if (length(missing_motus) > 0) {
  stop(paste("Missing required Motus columns:", paste(missing_motus, collapse = ", ")))
}

det <- tibble(
  ts = ymd_hms(safe_character(motus_raw[[motus_cols$ts]]), tz = "UTC", quiet = TRUE),
  motusTagID = safe_numeric(motus_raw[[motus_cols$motusTagID]]),
  tagDeployID = safe_numeric(motus_raw[[motus_cols$tagDeployID]]),
  recvDeployName = safe_character(motus_raw[[motus_cols$recvDeployName]]),
  recvDeployLat = if (!is.na(motus_cols$recvDeployLat)) safe_numeric(motus_raw[[motus_cols$recvDeployLat]]) else NA_real_,
  recvDeployLon = if (!is.na(motus_cols$recvDeployLon)) safe_numeric(motus_raw[[motus_cols$recvDeployLon]]) else NA_real_,
  tagDepLat = if (!is.na(motus_cols$tagDepLat)) safe_numeric(motus_raw[[motus_cols$tagDepLat]]) else NA_real_,
  tagDepLon = if (!is.na(motus_cols$tagDepLon)) safe_numeric(motus_raw[[motus_cols$tagDepLon]]) else NA_real_,
  motusFilter = if (!is.na(motus_cols$motusFilter)) safe_numeric(motus_raw[[motus_cols$motusFilter]]) else NA_real_
) %>%
  filter(!is.na(ts), !is.na(motusTagID), !is.na(tagDeployID), !is.na(recvDeployName))

if ("motusFilter" %in% names(det) && any(!is.na(det$motusFilter))) {
  det <- det %>% filter(is.na(motusFilter) | motusFilter == 1)
}

tag_coords <- det %>%
  group_by(tagDeployID, motusTagID) %>%
  summarise(
    tagDepLat = if (all(is.na(tagDepLat))) NA_real_ else tagDepLat[which(!is.na(tagDepLat))[1]],
    tagDepLon = if (all(is.na(tagDepLon))) NA_real_ else tagDepLon[which(!is.na(tagDepLon))[1]],
    .groups = "drop"
  )

birds <- birds %>% left_join(tag_coords, by = c("tagDeployID", "motusTagID"))

cat("Motus detection rows after filtering:", nrow(det), "\n")
cat("Detection span:", format(min(det$ts, na.rm = TRUE)), "to", format(max(det$ts, na.rm = TRUE)), "\n\n")

# -------------------------------------------------------------------
# Primary analysis
# -------------------------------------------------------------------
wing_df <- build_condition_residuals(birds, size_var = "wing_mm")
wing_stopover <- calculate_tagging_site_stopover(
  wing_df, det,
  max_gap_days = 5,
  max_bout_days = 60
)

wing_analysis <- wing_df %>%
  inner_join(wing_stopover, by = c("tagDeployID", "motusTagID")) %>%
  mutate(metric_version = "primary_wing")

cat("Primary wing-based analysis n =", nrow(wing_analysis), "\n")

# -------------------------------------------------------------------
# Sensitivity analysis
# -------------------------------------------------------------------
tarsus_df <- build_condition_residuals(birds, size_var = "tarsus_mm")
tarsus_stopover <- calculate_tagging_site_stopover(
  tarsus_df, det,
  max_gap_days = 5,
  max_bout_days = 60
)

tarsus_analysis <- tarsus_df %>%
  inner_join(tarsus_stopover, by = c("tagDeployID", "motusTagID")) %>%
  mutate(metric_version = "sensitivity_tarsus")

cat("Sensitivity tarsus-based analysis n =", nrow(tarsus_analysis), "\n\n")

# -------------------------------------------------------------------
# Models
# -------------------------------------------------------------------
fit_primary_simple <- lm(stopover_days ~ body_condition, data = wing_analysis)
fit_primary_full <- lm(stopover_days ~ body_condition + capture_site + age + julian_capture, data = wing_analysis)

fit_tarsus_simple <- lm(stopover_days ~ body_condition, data = tarsus_analysis)
fit_tarsus_full <- lm(stopover_days ~ body_condition + capture_site + age + julian_capture, data = tarsus_analysis)

# -------------------------------------------------------------------
# Export tables and model summary
# -------------------------------------------------------------------
summary_table <- bind_rows(
  wing_analysis %>%
    mutate(size_metric = "wing chord") %>%
    select(metric_version, size_metric, band_number, motusTagID, tagDeployID, capture_site,
           age, sex, date_raw, time_raw, julian_capture, weight_g, wing_mm, tarsus_mm,
           fat_score, body_condition, focal_receiver, first_detect, last_detect,
           n_detection_days, stopover_days),
  tarsus_analysis %>%
    mutate(size_metric = "tarsus") %>%
    select(metric_version, size_metric, band_number, motusTagID, tagDeployID, capture_site,
           age, sex, date_raw, time_raw, julian_capture, weight_g, wing_mm, tarsus_mm,
           fat_score, body_condition, focal_receiver, first_detect, last_detect,
           n_detection_days, stopover_days)
)

write.csv(summary_table, file.path(out_dir, "body_condition_stopover_summary.csv"), row.names = FALSE)

sink(file.path(out_dir, "body_condition_model_summary.txt"))
cat("PRIMARY ANALYSIS: body condition = residual body mass relative to wing chord\n")
cat("====================================================================\n")
print(summary(fit_primary_simple))
cat("\n")
print(summary(fit_primary_full))
cat("\n\nSENSITIVITY ANALYSIS: body condition = residual body mass relative to tarsus length\n")
cat("===============================================================================\n")
print(summary(fit_tarsus_simple))
cat("\n")
print(summary(fit_tarsus_full))
cat("\n\nSample sizes\n")
cat("------------\n")
cat("Primary wing-based analysis n = ", nrow(wing_analysis), "\n", sep = "")
cat("Sensitivity tarsus-based analysis n = ", nrow(tarsus_analysis), "\n", sep = "")
sink()

# -------------------------------------------------------------------
# Export figures and show immediately
# -------------------------------------------------------------------
primary_plot <- make_plot(
  plot_df = wing_analysis,
  x_label = "Body condition at capture (wing-corrected residual mass)",
  out_file = file.path(out_dir, "s1_primary_wing_initial_stopover.png"),
  show_plot = TRUE
)

sensitivity_plot <- make_plot(
  plot_df = tarsus_analysis,
  x_label = "Body condition at capture (tarsus-corrected residual mass)",
  out_file = file.path(out_dir, "s1_sensitivity_tarsus_initial_stopover.png"),
  show_plot = TRUE
)

cat("Done.\n")
cat("Recommended supplementary figure: s1_primary_wing_initial_stopover.png\n")
cat("Sensitivity figure: s1_sensitivity_tarsus_initial_stopover.png\n")
cat("Summary table: body_condition_stopover_summary.csv\n")
cat("Model summary: body_condition_model_summary.txt\n")

# -------------------------------------------------------------------
# Open the main figure automatically
# -------------------------------------------------------------------
primary_png <- file.path(out_dir, "s1_primary_wing_initial_stopover.png")
if (file.exists(primary_png)) {
  cat("Opening main figure:", primary_png, "\n")
  try(browseURL(primary_png), silent = TRUE)
} else {
  cat("Could not find the main figure to open automatically.\n")
}
