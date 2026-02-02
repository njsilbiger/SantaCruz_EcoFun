# Compile temperature and light data - detided
# Tena Dhayalan
# 2026-02-02


library(data.table)
library(lubridate)
library(here)

# --------- file paths ----------
coops_file       <- here("Data","CO-OPS_Tidal_Height.csv")
tidal_file       <- here("Data","tidal_height_2.csv")
full_data_dir    <- here("Data", "Temperature", "Spring_2024", "full data")
output_dir       <- file.path(full_data_dir, "detided")
tz          <- "America/Los_Angeles"  # PDT / user's timezone
match_tol_secs   <- 0
# names 
coops_date_col   <- "Date"          # CO-OPS date column name
coops_time_col   <- "Time_LST_LDT"  # CO-OPS time column name
coops_height_col <- "Height_m"      # CO-OPS height column name
plot_datetime_col <- c("Date-Time (PDT)")

# ---------- function -----------

# ---- CO-OPS ----
coops <- fread(coops_file)
coops[, datetime := parse_date_time(
  paste(Date, Time_LST_LDT),
  orders = c("Ymd HMS","Ymd HM","mdY HMS","mdY HM","Y-m-d HMS","Y-m-d HM","Y/m/d HMS","Y/m/d HM"),
  tz = tz
)]
coops <- coops[!is.na(datetime), .(datetime, coops_height = as.numeric(Height_m))]
setorder(coops, datetime)
setkey(coops, datetime)

# ---- tidal heights (constant per plot) ----
tidal <- fread(tidal_file)
# column names
tidal_plot_col   <- "Plot"         
tidal_height_col <- "tidal_height" 

tidal[, plot_id := as.character(get(tidal_plot_col))]
tidal[, plot_tidal_height := as.numeric(get(tidal_height_col))]
tidal <- unique(tidal[!is.na(plot_id) & !is.na(plot_tidal_height), .(plot_id, plot_tidal_height)], by="plot_id")
tidal_lookup <- setNames(tidal$plot_tidal_height, tidal$plot_id)

# ---- loop files ----
files <- list.files(full_data_dir, pattern="\\.csv$", full.names=TRUE)
files <- files[!grepl("/detided/", files)]

for (f in files) {
  dt <- fread(f)
  
  # plot id from the Plot column
  pid <- as.character(dt$Plot[1])
  th <- unname(tidal_lookup[pid])
  
  if (is.na(th)) {
    warning("No tidal height for Plot='", pid, "' in tidal_height.csv; skipping ", basename(f))
    next
  }
  
  # parse datetime
  dt[, datetime := parse_date_time(`Date-Time (PDT)`,
                                   orders = c("Ymd HMS","Ymd HM",
                                              "mdY HMS","mdY HM",
                                              "Y-m-d HMS","Y-m-d HM",
                                              "m/d/Y HMS","m/d/Y HM"),
                                   tz = tz)]
  dt <- dt[!is.na(datetime)]
  if (nrow(dt) == 0) {
    warning("All datetimes failed to parse for ", basename(f))
    next
  }
  setorder(dt, datetime)
  
  # nearest CO-OPS height for each timestamp, tolerance set to 0
  dt_times <- data.table(datetime = dt$datetime)
  setkey(dt_times, datetime)
  coops_nearest <- coops[dt_times, roll="nearest"]
  
  dt[, coops_height := coops_nearest$coops_height]
  dt[, coops_dt_diff := abs(as.numeric(difftime(datetime, coops_nearest$datetime, units="secs")))]
  if (is.finite(match_tol_secs)) dt[coops_dt_diff > match_tol_secs, coops_height := NA_real_]
  
  # filter
  detided <- dt[!is.na(coops_height) & (coops_height > th)]
  
  # remove helper columns before writing
  detided[, c("datetime", "coops_height", "coops_dt_diff") := NULL]
  
  out_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(f)), "_detided.csv")
  )
  fwrite(detided, out_file)
  
  message(basename(f), " | Plot=", pid, " | th=", th,
          " | kept=", nrow(detided), " / ", nrow(dt), " | wrote ", out_file)
}

message("Done. Detided files are in: ", output_dir)