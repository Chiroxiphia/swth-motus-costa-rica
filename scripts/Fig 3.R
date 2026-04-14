# =========================================================
# FIGURE 3 — Minimum spring stopover windows in Costa Rica
# Final manuscript-style version
# Uses the same analytical birds as Figure 2
# Expected subset: 28 birds (16 Finca Cántaros, 12 Rancho Quemado)
# =========================================================

# ---------------------------
# 0) Packages
# ---------------------------
library(tidyverse)
library(lubridate)
library(ggplot2)
library(openxlsx)
library(stringr)
library(janitor)
library(grid)

# ---------------------------
# 1) Select files manually
# ---------------------------
cat("Select the BANDING Excel file...\n")
band_file <- file.choose()

cat("Select the DETECTIONS Excel file...\n")
det_file <- file.choose()

# ---------------------------
# 2) Helper functions
# ---------------------------
clean_tag <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\.0$", "") %>%
    str_replace_all("[^0-9]", "")
}

pick_first <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

normalize_text <- function(x) {
  x %>%
    as.character() %>%
    iconv(from = "", to = "ASCII//TRANSLIT") %>%
    str_to_lower() %>%
    str_squish()
}

harmonize_site <- function(x) {
  x_raw  <- as.character(x)
  x_norm <- normalize_text(x_raw)
  
  case_when(
    str_detect(x_norm, "casa\\s*bont") ~ "Finca Cántaros",
    str_detect(x_norm, "finca") & str_detect(x_norm, "cat|cant|bont|botan") ~ "Finca Cántaros",
    str_detect(x_norm, "cataros|cantaros") ~ "Finca Cántaros",
    str_detect(x_norm, "osa|rancho|quemado") ~ "Rancho Quemado",
    TRUE ~ x_raw
  )
}

parse_motus_datetime <- function(x) {
  x_chr <- as.character(x)
  
  # 1) Standard datetime string
  out <- suppressWarnings(lubridate::ymd_hms(x_chr, tz = "UTC", quiet = TRUE))
  
  # 2) Unix timestamp fallback
  idx_na <- is.na(out) & !is.na(x_chr) & x_chr != ""
  if (any(idx_na)) {
    x_num <- suppressWarnings(as.numeric(x_chr[idx_na]))
    out[idx_na] <- suppressWarnings(as.POSIXct(x_num, origin = "1970-01-01", tz = "UTC"))
  }
  
  # 3) Generic POSIXct fallback
  idx_na2 <- is.na(out) & !is.na(x_chr) & x_chr != ""
  if (any(idx_na2)) {
    out[idx_na2] <- suppressWarnings(as.POSIXct(x_chr[idx_na2], tz = "UTC"))
  }
  
  out
}

# ---------------------------
# 3) Read files
# ---------------------------
band_raw <- openxlsx::read.xlsx(
  band_file,
  startRow = 2,
  colNames = TRUE
)

det_raw <- openxlsx::read.xlsx(
  det_file,
  colNames = TRUE
)

names(band_raw) <- janitor::make_clean_names(names(band_raw))
names(det_raw)  <- janitor::make_clean_names(names(det_raw))

cat("\nBanding columns:\n")
print(names(band_raw))

cat("\nDetection columns:\n")
print(names(det_raw))

# ---------------------------
# 4) Detect useful columns
# ---------------------------
band_tag_col <- pick_first(
  names(band_raw),
  c("tag", "tag_id", "tagdeployid", "tag_deploy_id")
)

band_site_col <- pick_first(
  names(band_raw),
  c("tagging_site", "site", "capture_site", "station_site", "field_site", "banding_site")
)

det_tag_col <- pick_first(
  names(det_raw),
  c("tagdeployid", "tag_deploy_id", "tag", "tag_id")
)

det_date_col <- pick_first(
  names(det_raw),
  c("ts_corrected", "tscorrected", "ts", "timestamp", "datetime", "date")
)

det_country_col <- pick_first(
  names(det_raw),
  c("country", "recv_country", "deploy_country")
)

det_filter_col <- pick_first(
  names(det_raw),
  c("motus_filter", "motusfilter")
)

det_lat_col <- pick_first(
  names(det_raw),
  c("recv_deploy_lat", "recvdeploylat", "latitude", "lat", "recv_lat")
)

det_lon_col <- pick_first(
  names(det_raw),
  c("recv_deploy_lon", "recvdeploylon", "longitude", "lon", "recv_lon")
)

cat("\nDetected columns:\n")
cat("band_tag_col   =", band_tag_col, "\n")
cat("band_site_col  =", band_site_col, "\n")
cat("det_tag_col    =", det_tag_col, "\n")
cat("det_date_col   =", det_date_col, "\n")
cat("det_country_col=", det_country_col, "\n")
cat("det_filter_col =", det_filter_col, "\n")
cat("det_lat_col    =", det_lat_col, "\n")
cat("det_lon_col    =", det_lon_col, "\n")

if (is.na(band_tag_col)) stop("Could not identify the tag column in the banding file.")
if (is.na(band_site_col)) stop("Could not identify the tagging site column in the banding file.")
if (is.na(det_tag_col)) stop("Could not identify the tag column in the detections file.")
if (is.na(det_date_col)) stop("Could not identify the date column in the detections file.")

# ---------------------------
# 5) Build site map from banding file
# ---------------------------
site_map <- band_raw %>%
  transmute(
    tag = clean_tag(.data[[band_tag_col]]),
    tagging_site_raw = as.character(.data[[band_site_col]])
  ) %>%
  mutate(
    tagging_site = harmonize_site(tagging_site_raw)
  ) %>%
  filter(!is.na(tag), tag != "", !is.na(tagging_site), tagging_site != "") %>%
  distinct(tag, tagging_site)

cat("\nRaw tagging-site values from banding:\n")
print(sort(unique(na.omit(as.character(band_raw[[band_site_col]])))))

cat("\nHarmonized tagging-site counts from banding file:\n")
print(table(site_map$tagging_site, useNA = "ifany"))

# ---------------------------
# 6) Prepare detections
# ---------------------------
det <- det_raw %>%
  transmute(
    tag = clean_tag(.data[[det_tag_col]]),
    det_datetime = parse_motus_datetime(.data[[det_date_col]]),
    country = if (!is.na(det_country_col)) as.character(.data[[det_country_col]]) else NA_character_,
    motus_filter = if (!is.na(det_filter_col)) suppressWarnings(as.numeric(.data[[det_filter_col]])) else NA_real_,
    lat = if (!is.na(det_lat_col)) suppressWarnings(as.numeric(.data[[det_lat_col]])) else NA_real_,
    lon = if (!is.na(det_lon_col)) suppressWarnings(as.numeric(.data[[det_lon_col]])) else NA_real_
  ) %>%
  mutate(det_date = as.Date(det_datetime)) %>%
  filter(!is.na(tag), tag != "", !is.na(det_date))

cat("\nNumber of detections after date parsing:\n")
print(nrow(det))

cat("\nPreview parsed datetimes:\n")
print(det %>% select(tag, det_datetime, det_date) %>% head(10))

# ---------------------------
# 7) Keep spring 2024
# ---------------------------
det_spring <- det %>%
  filter(
    det_date >= as.Date("2024-03-01"),
    det_date <= as.Date("2024-05-31")
  )

cat("\nNumber of detections in spring window:\n")
print(nrow(det_spring))

# ---------------------------
# 8) Keep motusfilter == 1 (or NA if absent)
# ---------------------------
det_spring_filt <- det_spring %>%
  filter(is.na(motus_filter) | motus_filter == 1)

cat("\nNumber of detections after motusfilter step:\n")
print(nrow(det_spring_filt))

# ---------------------------
# 9) Restrict to Costa Rica
# ---------------------------
cr_bbox <- list(
  lon_min = -86.5,
  lon_max = -82.0,
  lat_min = 8.0,
  lat_max = 11.5
)

has_country_info <- !all(is.na(det_spring_filt$country))
has_coords <- !all(is.na(det_spring_filt$lat)) & !all(is.na(det_spring_filt$lon))

if (has_country_info) {
  det_cr <- det_spring_filt %>%
    filter(str_to_lower(country) %in% c("costa rica", "costa_rica", "cr"))
  cat("\nCosta Rica filter used: COUNTRY column\n")
} else if (has_coords) {
  det_cr <- det_spring_filt %>%
    filter(
      lon >= cr_bbox$lon_min,
      lon <= cr_bbox$lon_max,
      lat >= cr_bbox$lat_min,
      lat <= cr_bbox$lat_max
    )
  cat("\nCosta Rica filter used: LAT/LON bounding box\n")
} else {
  det_cr <- det_spring_filt
  cat("\nWarning: no usable country or coordinate columns found; using all spring detections.\n")
}

cat("\nNumber of detections kept after Costa Rica restriction:\n")
print(nrow(det_cr))

# ---------------------------
# 10) Summarize detections in Costa Rica
# ---------------------------
stopover_summary <- det_cr %>%
  group_by(tag) %>%
  summarise(
    first_det_cr = min(det_date, na.rm = TRUE),
    last_det_cr  = max(det_date, na.rm = TRUE),
    min_stopover_days = as.numeric(last_det_cr - first_det_cr) + 1,
    .groups = "drop"
  ) %>%
  left_join(site_map, by = "tag")

cat("\nSummary counts before >=10 day filter:\n")
print(stopover_summary %>% count(tagging_site))

cat("\nPreview stopover_summary:\n")
print(
  stopover_summary %>%
    select(tag, tagging_site, first_det_cr, last_det_cr, min_stopover_days) %>%
    arrange(tagging_site, first_det_cr, tag)
)

# ---------------------------
# 11) Build Figure 2 analytical subset
# Figure 3 must use these same birds
# ---------------------------
fig2_dat <- stopover_summary %>%
  filter(min_stopover_days >= 10) %>%
  filter(!is.na(tagging_site)) %>%
  filter(tagging_site %in% c("Finca Cántaros", "Rancho Quemado")) %>%
  arrange(tagging_site, first_det_cr, last_det_cr, tag)

cat("\nFigure 2 / Figure 3 shared counts after >=10 day filter:\n")
print(fig2_dat %>% count(tagging_site))

cat("\nShared analytical birds:\n")
print(fig2_dat %>% select(tag, tagging_site, first_det_cr, last_det_cr, min_stopover_days))

write.csv(fig2_dat, "Figure2_and_Figure3_shared_stopover_data_spring2024.csv", row.names = FALSE)

# ---------------------------
# 12) Build Figure 3 data
# ---------------------------
fig3_dat <- fig2_dat %>%
  mutate(
    tagging_site = factor(tagging_site, levels = c("Finca Cántaros", "Rancho Quemado"))
  ) %>%
  arrange(tagging_site, first_det_cr, last_det_cr, tag) %>%
  group_by(tagging_site) %>%
  mutate(tag_f = factor(tag, levels = rev(unique(tag)))) %>%
  ungroup()

cat("\nFigure 3 counts:\n")
print(fig3_dat %>% count(tagging_site))

cat("\nFigure 3 preview:\n")
print(fig3_dat %>% select(tag, tagging_site, first_det_cr, last_det_cr, min_stopover_days))

write.csv(fig3_dat, "Figure3_data_minimum_stopover_windows_spring2024.csv", row.names = FALSE)

if (nrow(fig3_dat) == 0) {
  stop("fig3_dat is empty after filtering. Check diagnostics above.")
}

# ---------------------------
# 13) Colors
# ---------------------------
site_colors <- c(
  "Finca Cántaros" = "#1b9e77",
  "Rancho Quemado" = "#2c7fb8"
)

# ---------------------------
# 14) Build plot
# ---------------------------
p_fig3 <- ggplot(
  fig3_dat,
  aes(
    x = first_det_cr,
    xend = last_det_cr,
    y = tag_f,
    yend = tag_f,
    color = tagging_site
  )
) +
  geom_segment(linewidth = 1.2, lineend = "round") +
  facet_wrap(~ tagging_site, ncol = 1, scales = "free_y") +
  scale_color_manual(values = site_colors, guide = "none") +
  scale_x_date(
    limits = c(as.Date("2024-03-01"), as.Date("2024-05-31")),
    breaks = as.Date(c("2024-03-01", "2024-04-01", "2024-05-01")),
    date_labels = "%b",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Date",
    y = "Tag ID"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.6),
    strip.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 7),
    axis.line = element_line(linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.9),
    panel.spacing = unit(0.18, "lines")
  )

print(p_fig3)

# ---------------------------
# 15) Save outputs
# ---------------------------
ggsave(
  "Figure3_minimum_stopover_windows_spring2024.png",
  plot = p_fig3,
  width = 6.8,
  height = 7.6,
  dpi = 600
)

ggsave(
  "Figure3_minimum_stopover_windows_spring2024.pdf",
  plot = p_fig3,
  width = 6.8,
  height = 7.6
)

ggsave(
  "Figure3_minimum_stopover_windows_spring2024.tiff",
  plot = p_fig3,
  width = 6.8,
  height = 7.6,
  dpi = 600,
  compression = "lzw"
)