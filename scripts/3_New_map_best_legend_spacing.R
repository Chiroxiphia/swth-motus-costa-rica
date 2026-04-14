# ============================================================
# FIGURE_MAP_HIGHLIGHTED_TRACKS_SPRING2024_POLISHED.R
# Clean manuscript figure: spring 2024 only
# - Broad continental extent with Great Lakes visible
# - Short legend at right
# - Legend moved away from highlighted lines
# - Opens exported figure automatically on Mac/Windows/Linux
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(lubridate)
  library(sf)
  library(readxl)
  library(rnaturalearth)
  library(geosphere)
  library(scales)
})

Sys.setenv(TZ = "UTC")

# ---------------------------
# 1) SETTINGS
# ---------------------------

highlight_ids <- c(52993, 53001, 53002)
tagProjID_keep <- 280

highlight_labels <- c(
  "52993 - Rapid progression",
  "53001 - Long detectability",
  "53002 - Highest latitude"
)
names(highlight_labels) <- as.character(highlight_ids)

date_start <- "2024-03-01"
date_end   <- "2024-08-31"

max_step_km <- 2500

# Keep broad extent so Great Lakes remain visible
xlim <- c(-130, -58)
ylim <- c(5, 60)

plot_title <- "Migratory tracks of Swainsons Thrush inferred from Motus detections"
plot_subtitle <- "Spring 2024 detections only"

# Automatic opening after export
open_output <- TRUE
open_file <- "png"   # options: "png" or "pdf"

# ---------------------------
# 2) HELPERS
# ---------------------------

to_num_clean <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\s+", "", x)
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

parse_dt_any <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  if (inherits(x, "Date")) return(as.POSIXct(x, tz = "UTC"))

  x0 <- as.character(x)

  dt <- suppressWarnings(as.POSIXct(x0, tz = "UTC"))
  if (all(is.na(dt))) dt <- suppressWarnings(ymd_hms(x0, tz = "UTC", quiet = TRUE))
  if (all(is.na(dt))) dt <- suppressWarnings(ymd_hm(x0, tz = "UTC", quiet = TRUE))
  if (all(is.na(dt))) dt <- suppressWarnings(ymd(x0, tz = "UTC", quiet = TRUE))

  dt
}

# ---------------------------
# 3) LOAD FILE
# ---------------------------

message("Select the Motus detections Excel file...")
input_path <- file.choose()

# Output: save next to the selected Excel file
out_dir <- file.path(dirname(input_path), "graficos")
out_png <- file.path(out_dir, "Figure_migratory_tracks_highlighted_spring2024_polished.png")
out_pdf <- file.path(out_dir, "Figure_migratory_tracks_highlighted_spring2024_polished.pdf")

d <- read_excel(input_path, guess_max = 100000) |>
  as.data.frame()

# ---------------------------
# 4) REQUIRED COLUMNS
# ---------------------------

needed <- c(
  "ts",
  "motusTagID",
  "tagDeployID",
  "runID",
  "runLen",
  "motusFilter",
  "tagProjID",
  "recvDeployLat",
  "recvDeployLon",
  "recvDeployName"
)

missing <- setdiff(needed, names(d))
if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

# ---------------------------
# 5) CLEAN + FILTER
# ---------------------------

df <- d %>%
  mutate(
    det_time = parse_dt_any(ts),
    recvDeployLat = to_num_clean(recvDeployLat),
    recvDeployLon = to_num_clean(recvDeployLon),
    tagDeployID   = to_num_clean(tagDeployID),
    motusTagID    = to_num_clean(motusTagID),
    tagProjID     = to_num_clean(tagProjID),
    runLen        = to_num_clean(runLen),
    motusFilter   = to_num_clean(motusFilter),
    recvDeployName = if_else(
      is.na(recvDeployName) | recvDeployName == "",
      paste(recvDeployLat, recvDeployLon, sep = ":"),
      as.character(recvDeployName)
    )
  ) %>%
  filter(
    !is.na(det_time),
    !is.na(recvDeployLat),
    !is.na(recvDeployLon),
    !is.na(tagDeployID),
    between(recvDeployLat, -90, 90),
    between(recvDeployLon, -180, 180),
    tagProjID == tagProjID_keep,
    is.na(motusFilter) | motusFilter == 1
  ) %>%
  filter(
    det_time >= as_datetime(date_start, tz = "UTC"),
    det_time <  as_datetime(date_end, tz = "UTC") + days(1)
  )

if (nrow(df) == 0) {
  stop("No detections remained after filtering.")
}

# ---------------------------
# 6) REDUCE REPEATED POINTS
# ---------------------------

df_path <- df %>%
  group_by(
    tagDeployID,
    motusTagID,
    runID,
    recvDeployName,
    recvDeployLat,
    recvDeployLon
  ) %>%
  summarize(
    det_time = mean(det_time),
    max_runLen = suppressWarnings(max(runLen, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(tagDeployID, det_time)

if (nrow(df_path) == 0) {
  stop("No path data available after summarizing detections.")
}

# ---------------------------
# 7) QC: REMOVE EXTREME JUMPS
# ---------------------------

df_path_qc <- df_path %>%
  group_by(tagDeployID) %>%
  arrange(det_time, .by_group = TRUE) %>%
  mutate(
    prev_lon = lag(recvDeployLon),
    prev_lat = lag(recvDeployLat),
    step_km = geosphere::distHaversine(
      matrix(c(prev_lon, prev_lat), ncol = 2),
      matrix(c(recvDeployLon, recvDeployLat), ncol = 2)
    ) / 1000
  ) %>%
  filter(is.na(step_km) | step_km <= max_step_km) %>%
  ungroup()

if (nrow(df_path_qc) == 0) {
  stop("No data left after QC filtering.")
}

# ---------------------------
# 8) BUILD POINTS + TRACKS
# ---------------------------

pts_sf <- st_as_sf(
  df_path_qc,
  coords = c("recvDeployLon", "recvDeployLat"),
  crs = 4326,
  remove = FALSE
)

tracks_sf <- df_path_qc %>%
  arrange(tagDeployID, det_time) %>%
  group_by(tagDeployID) %>%
  filter(n() >= 2) %>%
  summarize(
    geometry = st_sfc(
      st_linestring(as.matrix(cbind(recvDeployLon, recvDeployLat))),
      crs = 4326
    ),
    .groups = "drop"
  ) %>%
  st_as_sf()

tracks_sf <- tracks_sf %>%
  mutate(
    track_class = if_else(tagDeployID %in% highlight_ids, "Highlighted", "Other"),
    tag_label = as.character(tagDeployID)
  )

pts_sf <- pts_sf %>%
  mutate(
    track_class = if_else(tagDeployID %in% highlight_ids, "Highlighted", "Other"),
    tag_label = as.character(tagDeployID)
  )

tracks_sf$tag_label <- factor(
  as.character(tracks_sf$tag_label),
  levels = as.character(highlight_ids)
)

pts_sf$tag_label <- factor(
  as.character(pts_sf$tag_label),
  levels = as.character(highlight_ids)
)

# ---------------------------
# 9) BASEMAP
# ---------------------------

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# ---------------------------
# 10) PLOT
# ---------------------------

p <- ggplot() +
  geom_sf(data = world, fill = "gray96", color = "gray78", linewidth = 0.25) +
  geom_sf(
    data = tracks_sf %>% filter(track_class == "Other"),
    color = "gray70",
    linewidth = 0.55,
    alpha = 0.35,
    show.legend = FALSE
  ) +
  geom_sf(
    data = tracks_sf %>% filter(track_class == "Highlighted"),
    aes(color = tag_label),
    linewidth = 1.7,
    alpha = 0.98,
    show.legend = TRUE
  ) +
  geom_sf(
    data = pts_sf %>% filter(track_class == "Highlighted"),
    aes(color = tag_label),
    size = 2.6,
    alpha = 0.98,
    show.legend = FALSE
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(
    title = plot_title,
    subtitle = plot_subtitle,
    x = NULL,
    y = NULL,
    color = "Highlighted individuals"
  ) +
  scale_color_discrete(labels = highlight_labels) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major = element_line(color = "gray88", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray45", linewidth = 0.6, fill = NA),

    plot.title = element_text(face = "plain", size = 18),
    plot.subtitle = element_text(size = 12),
    plot.margin = margin(t = 12, r = 16, b = 10, l = 10),

    legend.position = c(0.84, 0.52),
    legend.justification = c(0, 0.5),
    legend.direction = "vertical",
    legend.box.background = element_blank(),
    legend.background = element_rect(fill = alpha("white", 0.90), color = NA),
    legend.key = element_rect(fill = NA, color = NA),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9.5),
    legend.spacing.y = unit(0.22, "cm")
  )

print(p)

# ---------------------------
# 11) SAVE
# ---------------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(out_png, plot = p, width = 12, height = 8, dpi = 600, bg = "white")
ggsave(out_pdf, plot = p, width = 12, height = 8, device = cairo_pdf, bg = "white")

message("Saved PNG: ", normalizePath(out_png))
message("Saved PDF: ", normalizePath(out_pdf))

# ---------------------------
# 12) OPEN EXPORTED FILE
# ---------------------------

if (open_output) {
  file_to_open <- if (tolower(open_file) == "pdf") out_pdf else out_png
  file_to_open <- normalizePath(file_to_open, mustWork = FALSE)

  if (file.exists(file_to_open)) {
    sysname <- Sys.info()[["sysname"]]

    if (sysname == "Darwin") {
      system2("open", shQuote(file_to_open), wait = FALSE)
    } else if (sysname == "Windows") {
      shell.exec(file_to_open)
    } else if (sysname == "Linux") {
      system2("xdg-open", shQuote(file_to_open), wait = FALSE)
    } else {
      message("Figure saved, but automatic opening is not configured for this OS: ", file_to_open)
    }
  } else {
    warning("Could not open file because it was not found: ", file_to_open)
  }
}
