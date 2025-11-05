##############################################
# Parallel timelag calculation for 3 gases
# using future.apply (robust, stable)
##############################################

library(RFlux)
library(forecast)
library(future)
library(future.apply)
library(progressr)

# -------------------
# PARALLEL & PROGRESS SETTINGS

plan(multisession, workers = parallel::detectCores() - 1)
options(future.globals.maxSize = 8e9)  # allow large data

handlers(handler_txtprogressbar())

# -------------------
# MAIN SETTINGS - ADJUST FOR YUOUR OWN PROJECT!

year <- "2025" # REPLACE WITH YOUR OWN!
input_folder  <- "I:/Sinikula/raw/" # REPLACE WITH YOUR OWN!
output_folder <- "I:/Sinikula/timelag_opt/" # REPLACE WITH YOUR OWN!

gas_settings <- list(
  list(name = "N2O", scalar_column = 4),
  list(name = "CH4", scalar_column = 7),
  list(name = "CO2", scalar_column = 9)
)

wind_x_col <- 10
wind_y_col <- 11
wind_z_col <- 12
temp_col   <- 13

if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
  cat("Output folder created:", output_folder, "\n")
} else {
  cat("Output folder already exists:", output_folder, "\n")
}

# -------------------
# FUNCTIONS

mad_despike <- function(data, column_index, threshold = 5) {
  column <- data[[column_index]]
  median_val <- median(column, na.rm = TRUE)
  mad_val <- median(abs(column - median_val), na.rm = TRUE)
  limit <- threshold * mad_val
  column[abs(column - median_val) > limit] <- NA
  data[[column_index]] <- column
  return(data)
}

double_coordinate_rotation <- function(data, x_col, y_col, z_col) {
  u_mean <- mean(data[[x_col]], na.rm = TRUE)
  v_mean <- mean(data[[y_col]], na.rm = TRUE)
  w_mean <- mean(data[[z_col]], na.rm = TRUE)
  theta <- atan2(v_mean, u_mean)
  u_rot <- data[[x_col]] * cos(theta) + data[[y_col]] * sin(theta)
  v_rot <- -data[[x_col]] * sin(theta) + data[[y_col]] * cos(theta)
  data[[x_col]] <- u_rot
  data[[y_col]] <- v_rot
  phi <- atan2(w_mean, mean(data[[x_col]], na.rm = TRUE))
  u_rot2 <- data[[x_col]] * cos(phi) + data[[z_col]] * sin(phi)
  w_rot2 <- -data[[x_col]] * sin(phi) + data[[z_col]] * cos(phi)
  data[[x_col]] <- u_rot2
  data[[z_col]] <- w_rot2
  return(data)
}

# -------------------
# LOAD FILE LIST

file_list <- list.files(input_folder, pattern = "\\.dat$", full.names = TRUE)
cat("Found", length(file_list), "files to process.\n")

# -------------------
# MAIN PARALLEL LOOP WITH PROGRESS BAR

results_list <- with_progress({
  
  p <- progressor(along = file_list)
  
  future_lapply(
    file_list,
    future.seed = TRUE,   # FIX RANDOMNESS WARNINGS
    function(file) {
      
      filename <- tools::file_path_sans_ext(basename(file))
      p(sprintf("Processing %s", filename))
      
      # Try-catch protects worker crashes
      tryCatch({
        
        data <- read.table(file, header = TRUE, sep = "\t")
        
        # Basic despike once
        data <- mad_despike(data, wind_z_col)
        data <- mad_despike(data, wind_x_col)
        data <- mad_despike(data, wind_y_col)
        data <- mad_despike(data, temp_col)
        
        valid_rows <- is.finite(data[[temp_col]]) &
          is.finite(data[[wind_z_col]])
        data <- data[valid_rows, ]
        
        # Not enough data to estimate timelag
        if (nrow(data) < 100) {
          message(filename, " skipped (too few rows)")
          return(NULL)
        }
        
        data <- double_coordinate_rotation(data, wind_x_col, wind_y_col, wind_z_col)
        
        gas_outputs <- list()
        
        # ----------- Gas loop -----------
        for (gas in gas_settings) {
          gas_name      <- gas$name
          scalar_column <- gas$scalar_column
          
          data_g <- mad_despike(data, scalar_column)
          vals <- data_g[[scalar_column]]
          scalar_var <- var(vals, na.rm = TRUE)
          
          # Skip gas if no variation
          if (is.na(scalar_var) || scalar_var == 0) {
            message(filename, " (", gas_name, ") has no variance, skipping.")
            next
          }
          
          # Save plot
          png_filename <- file.path(
            output_folder,
            paste0(filename, "_", gas_name, "_tlag.png")
          )
          
          png(png_filename, width = 1200, height = 900)
          
          # try timelag calculation
          tlag_output <- try({
            tlag_detection(
              scalar_var = data_g[[scalar_column]],
              tsonic_var = data_g[[temp_col]],
              w_var      = data_g[[wind_z_col]],
              mfreq = 10,
              wdt   = 5,
              lws   = 0,
              uws   = 10,
              Rboot = 29,
              plot.it = TRUE
            )
          }, silent = TRUE)
          
          dev.off()
          
          # error handling inside the gas loop
          if (inherits(tlag_output, "try-error")) {
            message(filename, " (", gas_name, ") tlag_detection failed")
            next
          }
          
          df <- as.data.frame(tlag_output)
          df <- cbind(Filename = filename, df)
          gas_outputs[[gas_name]] <- df
        }
        
        return(gas_outputs)
        
      }, error = function(e) {
        message("Worker crashed on: ", filename, " → ", e$message)
        return(NULL)
      })
    })
})

# -------------------
# COLLECT RESULTS PER GAS

results_by_gas <- list()
for (g in gas_settings) {
  results_by_gas[[g$name]] <- list()
}

for (res in results_list) {
  if (is.null(res)) next
  for (g in names(res)) {
    results_by_gas[[g]][[res[[g]]$Filename[1]]] <- res[[g]]
  }
}

# -------------------
# SAVE RESULTS PER GAS

for (g in gas_settings) {
  gas_name <- g$name
  message("\nSaving results for ", gas_name, " ...")
  
  combined <- do.call(rbind, results_by_gas[[gas_name]])
  
  # Add HDI (95% highest density interval)
  # pwb_uci - pwb_lci 
  # ----------------------------------------------
  if (all(c("pwb_uci","pwb_lci") %in% names(combined))) {
    combined$HDI <- combined$pwb_uci - combined$pwb_lci
  } else {
    combined$HDI <- NA
  }
  
  # ----------------------------------------------
  # Add datetime parsed from Filename
  # ----------------------------------------------
  parts <- strsplit(combined$Filename, "_")
  
  # Extract components (site, year, month, day, HHMM)
  combined$year   <- as.integer(sapply(parts, `[`, 2))
  combined$month  <- as.integer(sapply(parts, `[`, 3))
  combined$day    <- as.integer(sapply(parts, `[`, 4))
  hm              <- sapply(parts, `[`, 5)
  
  combined$hour   <- as.integer(substr(hm, 1, 2))
  combined$minute <- as.integer(substr(hm, 3, 4))
  
  # Create POSIXct datetime (assumes local time)
  combined$datetime <- as.POSIXct(
    sprintf("%04d-%02d-%02d %02d:%02d:00",
            combined$year,
            combined$month,
            combined$day,
            combined$hour,
            combined$minute),
    format = "%Y-%m-%d %H:%M:%S"
  )
  
  # ----------------------------------------------

  
  # Optional: remove temporary columns if you like
   combined$year <- combined$month <- combined$day <- combined$hour <- combined$minute <- NULL
  
  
  out_csv <- file.path(
    output_folder,
    paste0(gas_name, "_tlag_results_", year, ".csv")
  )
  
  write.csv(combined, out_csv, row.names = FALSE)
  
  message(" → ", out_csv)
}

message("\nALL DONE ✅\n")
