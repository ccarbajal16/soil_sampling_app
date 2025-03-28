# Soil Sampling Points Selection App
# A Shiny application for selecting soil sampling points using CLHS and spatial coverage methods

# Load required packages
library(shiny)
library(shinydashboard)
library(terra)
library(sf)
library(raster)
library(clhs)
library(spcosa)
library(sp)
library(rJava)
library(tidyverse)
library(tmap)
library(DT)
library(shinyWidgets)
library(leaflet)

# UI
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "Soil Sampling Tool"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Input", tabName = "data", icon = icon("upload")),
      menuItem("CLHS Method", tabName = "clhs", icon = icon("layer-group")),
      menuItem("Spatial Coverage", tabName = "spatcov", icon = icon("map")),
      menuItem("Results", tabName = "results", icon = icon("download"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Data Input Tab
      tabItem(tabName = "data",
        fluidRow(
          box(
            title = "Input Configuration", width = 12, status = "primary",
            collapsible = TRUE,
            column(6,
              radioButtons("inputType", "Choose input method:",
                          choices = c("Upload raster files" = "upload"),
                          selected = "upload"),
              conditionalPanel(
                condition = "input.inputType == 'upload'",
                fileInput("rasterFiles", "Upload raster files (GeoTIFF format)",
                         multiple = TRUE, accept = c(".tif", ".tiff")),
                helpText("Upload multiple GeoTIFF files. All files must have the same extent and projection.")
              )
            ),
            column(6,
              radioButtons("areaType", "Select area of interest:",
                          choices = c("Use raster extent" = "raster",
                                     "Upload boundary" = "upload"),
                          selected = "raster"),
              conditionalPanel(
                condition = "input.areaType == 'upload'",
                fileInput("boundaryFile", "Upload boundary file (.gpkg, .geojson)", 
                         accept = c(".gpkg", ".geojson"))
              )
            )
          )
        ),
        fluidRow(
          box(
            title = "Data Summary", width = 12, status = "success",
            actionButton("checkData", "Check Data", icon = icon("check"), 
                        class = "btn-success"),
            hr(),
            htmlOutput("dataSummary"),
            verbatimTextOutput("extensionCheck")
          )
        )
      ),
      
      # CLHS Method Tab
      tabItem(tabName = "clhs",
        fluidRow(
          box(
            title = "CLHS Parameters", width = 6, status = "primary",
            numericInput("clhsSize", "Number of sample points:", 60, min = 1, max = 1000),
            numericInput("clhsIter", "Number of iterations:", 2000, min = 100, max = 10000),
            checkboxInput("useCost", "Use cost surface for sampling", TRUE),
            conditionalPanel(
              condition = "input.useCost == true",
              selectInput("costLayer", "Select cost layer:", choices = NULL)
            ),
            numericInput("costThreshold", "Cost threshold for accessibility:", 1.5, min = 0),
            actionButton("runClhs", "Run CLHS Sampling", icon = icon("play"), 
                        class = "btn-success")
          ),
          box(
            title = "CLHS Method Info", width = 6, status = "info",
            HTML("<p>The Conditioned Latin Hypercube Sampling (CLHS) method:</p>
                 <ul>
                 <li>Ensures good representation of the full range of each variable</li>
                 <li>Maintains the correlation structure between variables</li>
                 <li>Can incorporate a cost surface to prioritize accessible locations</li>
                 <li>Creates a statistically representative sample of multivariate distributions</li>
                 </ul>")
          )
        ),
        fluidRow(
          box(
            title = "CLHS Results", width = 12, status = "success",
            tabsetPanel(
              tabPanel("Raster Map", plotOutput("clhsRasterMap", height = 500)),
              tabPanel("Point Data", DTOutput("clhsTable"))
            )
          )
        )
      ),
      
      # Spatial Coverage Tab
      tabItem(tabName = "spatcov",
        fluidRow(
          column(6,
            box(
              title = "Spatial Coverage Parameters", width = 12, status = "primary",
              numericInput("nStrata", "Number of strata:", 100, min = 1, max = 1000),
              numericInput("nTry", "Number of optimization tries:", 10, min = 1, max = 100),
              numericInput("nGridCells", "Number of grid cells:", 12000, min = 1000),
              radioButtons("samplingType", "Sampling method:",
                          choices = c("Centroid" = "centroid",
                                     "Random" = "random",
                                     "Infill" = "infill"),
                          selected = "centroid")
            )
          ),
          column(6,
            box(
              title = "Run Method", width = 12, status = "primary",
              conditionalPanel(
                condition = "input.samplingType == 'random'",
                numericInput("nPointsPerStratum", "Points per stratum:", 1, min = 1, max = 10)
              ),
              conditionalPanel(
                condition = "input.samplingType == \"infill\"",
                fileInput("priorPoints", "Upload prior points file (.gpkg, .shp)", 
                         accept = c(".gpkg", ".shp", ".geojson")),
                checkboxInput("useLocalPrior", "Use existing prior points", TRUE),
                conditionalPanel(
                  condition = "input.useLocalPrior == true",
                  selectInput("localPriorPoints", "Select prior points:", 
                             choices = list.files(path = "data", pattern = "prior.*\\.(gpkg|shp)$", 
                                                full.names = TRUE))
                )
              ),
              actionButton("runSpatCov", "Run Spatial Coverage Sampling", 
                          icon = icon("play"), class = "btn-success")
            )
          )
        ),
        fluidRow(
          box(
            title = "Spatial Coverage Method Info", width = 12, status = "info", collapsible = TRUE, collapsed = TRUE,
            HTML("<p>The Spatial Coverage sampling method:</p>
                 <ul>
                 <li>Divides the area into compact strata of equal size</li>
                 <li>Ensures good spatial distribution of sampling points</li>
                 <li>Can take existing sampling points into account (infill mode)</li>
                 <li>Optimizes spatial coverage of the entire area</li>
                 </ul>")
          )
        ),
        fluidRow(
          box(
            title = "Spatial Coverage Results", width = 12, status = "success",
            tabsetPanel(
              tabPanel("Map", leafletOutput("spatcovMap", height = 500)),
              tabPanel("Point Data", DTOutput("spatcovTable"))
            )
          )
        )
      ),
      # Results Tab
      tabItem(tabName = "results",
        fluidRow(
          column(6,
            box(
              title = "Export Results", width = 12, status = "primary",
              radioButtons("exportMethod", "Select export format:",
                          choices = c("CSV" = "csv",
                                     "GeoPackage" = "gpkg"),
                          selected = "csv"),
              radioButtons("exportData", "Select data to export:",
                          choices = c("CLHS points" = "clhs",
                                     "Spatial Coverage points" = "spatcov",
                                     "Both methods" = "both"),
                          selected = "both"),
              textInput("exportPrefix", "File name prefix:", "sampling_points"),
              downloadButton("downloadData", "Download Data", class = "btn-success"),
              hr()
            )
          ),
          column(6,
            box(
              title = "Comparison of Methods", width = 12, status = "success",
              tabsetPanel(
                tabPanel("Map", leafletOutput("comparisonMap", height = 400)),
                tabPanel("Statistics", plotOutput("comparisonPlot", height = 400))
              )
            )
          )
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Reactive values to store data
  values <- reactiveValues(
    rasters = NULL,
    boundary = NULL,
    clhsPoints = NULL,
    spatcovPoints = NULL,
    strataObj = NULL,
    dataChecked = FALSE
  )
  
  # Update cost layer choices when rasters are loaded
  observe({
    req(values$rasters)
    raster_names <- names(values$rasters)
    updateSelectInput(session, "costLayer", choices = raster_names)
  })
  
  # Check data button handler
  observeEvent(input$checkData, {
    # Load raster data based on input method
    if (input$inputType == "upload" && !is.null(input$rasterFiles)) {
      temp_files <- input$rasterFiles$datapath
      values$rasters <- tryCatch({
        rast(temp_files)
      }, error = function(e) {
        showNotification(paste("Error loading uploaded rasters:", e$message), type = "error")
        return(NULL)
      })
    }
    
    # Load boundary based on input method
    if (input$areaType == "upload" && !is.null(input$boundaryFile)) {
      values$boundary <- tryCatch({
        st_read(input$boundaryFile$datapath, quiet = TRUE)
      }, error = function(e) {
        showNotification(paste("Error loading uploaded boundary:", e$message), type = "error")
        return(NULL)
      })
    } else if (input$areaType == "raster" && !is.null(values$rasters)) {
      # Create boundary from raster extent
      ext <- ext(values$rasters)
      values$boundary <- tryCatch({
        st_sf(geometry = st_as_sfc(st_bbox(c(xmin = ext[1], xmax = ext[2], 
                                           ymin = ext[3], ymax = ext[4])), 
                                 crs = crs(values$rasters, proj = TRUE)))
      }, error = function(e) {
        showNotification(paste("Error creating boundary from raster extent:", e$message), type = "error")
        return(NULL)
      })
    }
    
    # Check if all required data is loaded
    if (is.null(values$rasters) || is.null(values$boundary)) {
      values$dataChecked <- FALSE
      showNotification("Missing data! Please check your inputs.", type = "error")
    } else {
      values$dataChecked <- TRUE
      showNotification("Data loaded successfully!", type = "message")
    }
  })
  # Data summary output
  output$dataSummary <- renderUI({
    req(values$rasters, values$boundary)
    raster_names <- names(values$rasters)
    raster_dims <- dim(values$rasters)
    boundary_area <- st_area(values$boundary)
    HTML(paste0(
      "<h4>Loaded Data Summary:</h4>",
      "<p><b>Raster layers:</b> ", paste(raster_names, collapse = ", "), "</p>",
      "<p><b>Dimensions:</b> ", raster_dims[1], " x ", raster_dims[2], " (", 
      format(raster_dims[1] * raster_dims[2], big.mark = ","), " cells)</p>",
      "<p><b>Spatial resolution:</b> ", res(values$rasters)[1], " x ", res(values$rasters)[2], "</p>",
      "<p><b>Area of interest size:</b> ", format(boundary_area, units = "ha"), "</p>"
    ))
  })
  # Extension check output
  output$extensionCheck <- renderPrint({
    req(values$rasters)
    # Check if all rasters have the same extent
    r_ext <- ext(values$rasters)
    cat("Raster extent check:\n")
    cat("Xmin:", r_ext[1], "\n")
    cat("Xmax:", r_ext[2], "\n")
    cat("Ymin:", r_ext[3], "\n")
    cat("Ymax:", r_ext[4], "\n")
    cat("\nRaster dimensions check:\n")
    cat("Dimensions:", paste(dim(values$rasters), collapse = " x "), "\n")
    cat("Resolution:", paste(res(values$rasters), collapse = " x "), "\n")
    cat("\nCRS:", crs(values$rasters, proj = TRUE), "\n")
    # If boundary is available, check if it's compatible with rasters
    if (!is.null(values$boundary)) {
      b_crs <- st_crs(values$boundary)$epsg
      r_crs <- terra::crs(values$rasters, describe=TRUE)$epsg
      cat("\nBoundary CRS:", b_crs, "\n")
      # Ensure both CRS values are not NULL before comparison
      if (!is.null(b_crs) && !is.null(r_crs)) {
        cat("CRS match:", ifelse(b_crs == r_crs, "Yes ✓", "No ✗"), "\n")
        if (b_crs != r_crs) {
          cat("Note: The boundary will be reprojected to match raster CRS when processing.\n")
        }
      } else {
        cat("CRS match: Cannot determine (one or both CRS values are NULL)\n")
      }
    }
  })
  # Run CLHS sampling
  observeEvent(input$runClhs, {
    req(values$dataChecked, values$rasters, values$boundary)
    
    withProgress(message = 'Running CLHS sampling...', value = 0, {
      incProgress(0.2, detail = "Preparing data")
      
      # Create regular grid of points
      set.seed(5)
      sample_points <- tryCatch({
        spatSample(values$rasters, size = 1000, method = "regular", xy = TRUE)
      }, error = function(e) {
        showNotification(paste("Error in sampling points:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(sample_points)) return(NULL)
      
      # Filter points within boundary
      incProgress(0.2, detail = "Filtering points")
      points_sf <- st_as_sf(sample_points, coords = c("x", "y"), 
                           crs = crs(values$rasters, proj = TRUE))
      points_in_boundary <- st_intersects(points_sf, values$boundary, sparse = FALSE)[,1]
      sample_points <- sample_points[points_in_boundary, ]
      
      # Filter out NA values
      if (input$useCost && !is.null(input$costLayer) && input$costLayer %in% names(sample_points)) {
        sample_points <- sample_points[!is.na(sample_points[[input$costLayer]]), ]
      } else {
        # Filter out rows with any NA values
        sample_points <- na.omit(sample_points)
      }
      
      # Create spatial points data frame
      incProgress(0.2, detail = "Creating spatial points")
      s_spdf <- SpatialPointsDataFrame(
        coords = sample_points[, 1:2],
        data = sample_points[, 3:ncol(sample_points)]
      )
      
      # Run CLHS
      incProgress(0.2, detail = "Running CLHS algorithm")
      clhs_result <- tryCatch({
        if (input$useCost && !is.null(input$costLayer) && input$costLayer %in% names(s_spdf@data)) {
          clhs(s_spdf, size = input$clhsSize, progress = FALSE, 
               iter = input$clhsIter, cost = input$costLayer, simple = FALSE)
        } else {
          clhs(s_spdf, size = input$clhsSize, progress = FALSE, 
               iter = input$clhsIter, simple = FALSE)
        }
      }, error = function(e) {
        showNotification(paste("Error in CLHS algorithm:", e$message), type = "error")
        return(NULL)
      })
      
      incProgress(0.2, detail = "Processing results")
      
      if (!is.null(clhs_result)) {
        # Extract points
        subset_idx <- clhs_result$index_samples
        selected_points <- sample_points[subset_idx, ]
        
        # Add accessibility classification
        if (input$useCost && !is.null(input$costLayer) && input$costLayer %in% names(selected_points)) {
          selected_points$accessibility <- ifelse(
            selected_points[[input$costLayer]] > input$costThreshold, 
            "less accessible", "more accessible"
          )
        } else {
          selected_points$accessibility <- "not evaluated"
        }
        
        # Add ID column
        selected_points$ID <- seq_len(nrow(selected_points))
        
        # Store result
        values$clhsPoints <- selected_points
        values$clhs_result <- clhs_result
        values$s_spdf <- s_spdf
        values$subset_idx <- subset_idx
        
        showNotification("CLHS sampling completed successfully!", type = "message")
      }
    })
  })
  # CLHS raster map
  output$clhsRasterMap <- renderPlot({
    req(values$rasters, values$clhsPoints)
    
    # Select the DEM layer, or the first layer if DEM is not available
    dem_layer_name <- ifelse("dem" %in% names(values$rasters), "dem", names(values$rasters)[1])
    dem_layer <- values$rasters[[dem_layer_name]]
    
    # Plot the DEM with points overlaid
    terra::plot(dem_layer, main = "Sampling Points on Terrain")
    
    # Add points with different colors based on accessibility
    if ("accessibility" %in% names(values$clhsPoints)) {
      # Split points by accessibility
      more_accessible <- values$clhsPoints[values$clhsPoints$accessibility == "more accessible", ]
      less_accessible <- values$clhsPoints[values$clhsPoints$accessibility == "less accessible", ]
      
      # Add more accessible points in green
      if (nrow(more_accessible) > 0) {
        terra::points(more_accessible[, c("x", "y")], col = "green", pch = 16, cex = 1.5)
      }
      
      # Add less accessible points in red
      if (nrow(less_accessible) > 0) {
        terra::points(less_accessible[, c("x", "y")], col = "red", pch = 16, cex = 1.5)
      }
      
      # Add legend with clearer position and styling
      legend("topright", 
             legend = c("More accessible", "Less accessible"), 
             col = c("green", "red"), 
             pch = 16, 
             cex = 1.2,
             pt.cex = 1.5,
             bg = "white",
             box.lty = 1,
             title = "Accessibility")
    } else {
      # Add all points in blue if no accessibility info
      terra::points(values$clhsPoints[, c("x", "y")], col = "blue", pch = 16, cex = 1.5)
      
      # Add legend for blue points only
      legend("topright", 
             legend = "Sampling points", 
             col = "blue", 
             pch = 16, 
             cex = 1.2,
             pt.cex = 1.5,
             bg = "white",
             box.lty = 1)
    }
  })
  # Run Spatial Coverage sampling
  observeEvent(input$runSpatCov, {
    req(values$dataChecked, values$boundary)
    
    withProgress(message = 'Running Spatial Coverage sampling...', value = 0, {
      incProgress(0.2, detail = "Preparing data")
      
      # Convert boundary to spatial polygon for spcosa
      poly_sp <- tryCatch({
        as(values$boundary, "Spatial")
      }, error = function(e) {
        showNotification(paste("Error converting boundary:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(poly_sp)) return(NULL)
      
      # Create strata
      incProgress(0.2, detail = "Creating strata")
      
      strata <- tryCatch({
        if (input$samplingType == "infill") {
          # Load prior points
          prior_points <- NULL
          
          if (input$useLocalPrior && !is.null(input$localPriorPoints)) {
            prior_points <- tryCatch({
              st_read(input$localPriorPoints, quiet = TRUE)
            }, error = function(e) {
              showNotification(paste("Error loading prior points:", e$message), type = "error")
              return(NULL)
            })
          } else if (!is.null(input$priorPoints) && !is.null(input$priorPoints$datapath)) {
            prior_points <- tryCatch({
              st_read(input$priorPoints$datapath, quiet = TRUE)
            }, error = function(e) {
              showNotification(paste("Error loading uploaded prior points:", e$message), type = "error")
              return(NULL)
            })
          }
          
          if (is.null(prior_points)) {
            showNotification("Prior points required for infill sampling!", type = "error")
            return(NULL)
          }
          
          # Convert to sp
          prior_points_sp <- as(prior_points, "Spatial")
          
          # Create strata with prior points
          stratify(poly_sp, prior_points_sp,
                  nStrata = input$nStrata, 
                  nTry = input$nTry,
                  nGridCells = input$nGridCells)
        } else {
          # Create strata without prior points
          stratify(poly_sp,
                  nStrata = input$nStrata, 
                  nTry = input$nTry,
                  nGridCells = input$nGridCells)
        }
      }, error = function(e) {
        showNotification(paste("Error creating strata:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(strata)) return(NULL)
      
      values$strataObj <- strata
      
      # Sample points
      incProgress(0.4, detail = "Sampling points")
      
      sample_points <- tryCatch({
        if (input$samplingType == "random") {
          spsample(strata, n = input$nPointsPerStratum)
        } else {
          spsample(strata)
        }
      }, error = function(e) {
        showNotification(paste("Error sampling points:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(sample_points)) return(NULL)
      
      # Convert to data frame and then to sf
      incProgress(0.2, detail = "Processing results")
      points_df <- as(sample_points, "data.frame")
      points_sf <- st_as_sf(points_df, coords = c("x1", "x2"))
      st_crs(points_sf) <- st_crs(values$boundary)
      
      # Add ID column
      points_sf$ID <- seq_len(nrow(points_sf))
      
      # Store result
      values$spatcovPoints <- points_sf
      
      showNotification("Spatial Coverage sampling completed successfully!", type = "message")
    })
  })
  # CLHS map
  output$clhsMap <- renderTmap({
    req(values$clhsPoints, values$boundary)
    
    # Convert CLHS points to sf
    clhs_sf <- st_as_sf(values$clhsPoints, coords = c("x", "y"), 
                        crs = crs(values$rasters, proj = TRUE))
    
    # Create tmap object with boundary and points
    tm <- tm_shape(values$boundary) +
          tm_borders(col = "black", lwd = 2, alpha = 0.7) +
          tm_shape(clhs_sf)
    
    # Add points with styling based on accessibility if available
    if ("accessibility" %in% names(clhs_sf)) {
      tm <- tm + tm_dots(col = "accessibility", 
                        palette = c("green", "red"),
                        size = 0.1,
                        title = "Accessibility",
                        popup.vars = c("ID", "accessibility"),
                        id = "ID")
    } else {
      tm <- tm + tm_dots(col = "blue", 
                        size = 0.1,
                        popup.vars = "ID",
                        id = "ID")
    }
    
    # Add basemap
    tm + tm_basemap("Esri.WorldImagery") +
      tm_layout(legend.outside = TRUE,
                legend.position = c("right", "bottom"))
  })
  # CLHS table
  output$clhsTable <- renderDT({
    req(values$clhsPoints)
    # Create data table with selected columns
    point_data <- as.data.frame(values$clhsPoints)
    # Select columns to display
    display_cols <- c("ID", "x", "y", "accessibility")
    # Add at most 5 more environmental variables
    env_vars <- setdiff(names(point_data), c("ID", "x", "y", "accessibility"))
    display_cols <- c(display_cols, env_vars[1:min(5, length(env_vars))])
    datatable(point_data[, display_cols, drop = FALSE],
              options = list(pageLength = 10, scrollX = TRUE))
  })
  # CLHS variable distribution plot
  output$clhsPlot <- renderPlot({
    req(values$clhsPoints, values$clhs_result)
    # Check if objective function values exist and are finite
    obj_values <- values$clhs_result$objective_function
    if (is.null(obj_values) || length(obj_values) == 0 || all(!is.finite(obj_values))) {
      # Return an empty plot with a message if no valid data
      plot(0, 0, type = "n", axes = FALSE, ann = FALSE)
      text(0, 0, "No valid objective function data available", cex = 1.5)
      return()
    }
    # Filter out any non-finite values
    obj_values <- obj_values[is.finite(obj_values)]
    if (length(obj_values) == 0) {
      # Return an empty plot with a message if all values were non-finite
      plot(0, 0, type = "n", axes = FALSE, ann = FALSE)
      text(0, 0, "No valid objective function data available", cex = 1.5)
      return()
    }
    # Plot objective function convergence with valid data
    plot(obj_values, type = "l", col = "blue", lwd = 2, 
         xlab = "Iteration", ylab = "Objective function value",
         main = "CLHS Optimization Progress")
  })
  # Also improve the comparison plot to handle potential non-finite values
  output$comparisonPlot <- renderPlot({
    req(values$rasters)
    # Prepare data for plotting
    plot_data <- data.frame(Method = character(), Value = numeric(), Variable = character())
    # Extract values for CLHS points if available
    if (!is.null(values$clhsPoints) && nrow(values$clhsPoints) > 0) {
      clhs_sf <- st_as_sf(values$clhsPoints, coords = c("x", "y"), 
                         crs = crs(values$rasters, proj = TRUE))
      raster_names <- names(values$rasters)
      for (var_name in raster_names[1:min(3, length(raster_names))]) {
        if (var_name %in% names(values$clhsPoints)) {
          # Filter out non-finite values
          values_to_add <- values$clhsPoints[[var_name]]
          values_to_add <- values_to_add[is.finite(values_to_add)]
          if (length(values_to_add) > 0) {
            plot_data <- rbind(plot_data, data.frame(
              Method = "CLHS",
              Value = values_to_add,
              Variable = var_name
            ))
          }
        }
      }
    }
    # Extract values for Spatial Coverage points if available
    if (!is.null(values$spatcovPoints) && nrow(values$spatcovPoints) > 0) {
      # Extract values at point locations
      points_coords <- st_coordinates(values$spatcovPoints)
      raster_names <- names(values$rasters)
      for (var_name in raster_names[1:min(3, length(raster_names))]) {
        point_values <- terra::extract(values$rasters[[var_name]], points_coords)
        if (!is.null(point_values) && nrow(point_values) > 0 && var_name %in% names(point_values)) {
          # Filter out non-finite values
          values_to_add <- point_values[[var_name]]
          values_to_add <- values_to_add[is.finite(values_to_add)]
          if (length(values_to_add) > 0) {
            plot_data <- rbind(plot_data, data.frame(
              Method = "Spatial Coverage",
              Value = values_to_add,
              Variable = var_name
            ))
          }
        }
      }
    }
    # Create comparison plot if data is available
    if (nrow(plot_data) > 0) {
      ggplot(plot_data, aes(x = Method, y = Value, fill = Method)) +
        geom_boxplot() +
        facet_wrap(~ Variable, scales = "free_y") +
        theme_bw() +
        labs(title = "Comparison of Variable Distributions",
             x = "Sampling Method", y = "Value") +
        theme(legend.position = "bottom")
    } else {
      # Empty plot with message
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Run both methods to see comparison") +
        theme_void()
    }
  })
  # Spatial coverage map
  output$spatcovMap <- renderLeaflet({
    req(values$spatcovPoints)
    
    # Convert to 4326 for leaflet
    points_4326 <- st_transform(values$spatcovPoints, 4326)
    
    # Convert boundary to 4326 as well
    boundary_4326 <- st_transform(values$boundary, 4326)
    
    # Create leaflet map with boundary and points
    map <- leaflet() %>%
      addProviderTiles(providers$Esri.WorldImagery) %>%
      # Add the boundary polygon
      addPolygons(data = boundary_4326, 
                 weight = 4,
                 color = "#000000",
                 fillOpacity = 0,
                 label = "Study area")
    
    # Add the points
    map <- map %>%
      addCircleMarkers(data = points_4326, radius = 5, 
                      color = "blue",
                      fillOpacity = 1,
                      popup = ~paste("ID:", ID),
                      label = ~as.character(ID))
    
    # Add a legend for the points
    map %>%
      addLegend(
        position = "bottomright",
        title = "Sampling Points",
        colors = "blue",
        labels = "Sampling locations",
        opacity = 0.7
      )
  })
  
  # Spatial coverage table
  output$spatcovTable <- renderDT({
    req(values$spatcovPoints)
    # Create data table
    point_data <- cbind(
      ID = values$spatcovPoints$ID,
      st_coordinates(values$spatcovPoints)
    )
    datatable(point_data,
              options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # Strata plot - Keep this function but it won't be displayed in the UI
  output$strataPlot <- renderPlot({
    req(values$strataObj, values$spatcovPoints)
    
    # Create a clean plot area
    par(mar = c(4, 4, 2, 2))
    
    # Plot strata
    plot(values$strataObj, main = paste(input$nStrata, "compact strata with sampling points"))
    
    # Add points based on the selected sampling type
    if (input$samplingType == "centroid") {
      # For centroid sampling, highlight the centroids
      points_coords <- st_coordinates(values$spatcovPoints)
      points(points_coords, col = "red", pch = 16, cex = 1.5)
      # Add legend for points
      legend("topright", 
             legend = "Centroid sampling points", 
             col = "red", 
             pch = 16, 
             cex = 0.8,
             pt.cex = 1.5,
             bg = "white")
    } else if (input$samplingType == "random") {
      # For random sampling
      points_coords <- st_coordinates(values$spatcovPoints)
      points(points_coords, col = "blue", pch = 16, cex = 1.5)
      # Add legend for random points
      legend("topright", 
             legend = "Random sampling points", 
             col = "blue", 
             pch = 16, 
             cex = 0.8,
             pt.cex = 1.5,
             bg = "white")
    } else if (input$samplingType == "infill") {
      # For infill sampling, show both existing and new points
      points_coords <- st_coordinates(values$spatcovPoints)
      points(points_coords, col = "green", pch = 16, cex = 1.5)
      # Show prior points if available
      if (exists("prior_points_sp", envir = environment(server))) {
        prior_coords <- coordinates(prior_points_sp)
        points(prior_coords, col = "purple", pch = 17, cex = 1.5)
        # Legend with both point types
        legend("topright", 
               legend = c("New infill points", "Existing points"), 
               col = c("green", "purple"), 
               pch = c(16, 17),
               cex = 0.8,
               pt.cex = 1.5,
               bg = "white")
      } else {
        # Legend with only new points
        legend("topright", 
               legend = "Infill sampling points", 
               col = "green", 
               pch = 16,
               cex = 0.8,
               pt.cex = 1.5,
               bg = "white")
      }
    }
  })
  # Comparison map
  output$comparisonMap <- renderLeaflet({
    # Create base map with boundary
    map <- leaflet() %>%
      addProviderTiles(providers$Esri.WorldImagery)
    
    # Add boundary if available
    if (!is.null(values$boundary)) {
      boundary_4326 <- st_transform(values$boundary, 4326)
      map <- map %>%
        addPolygons(data = boundary_4326, 
                   weight = 2,
                   color = "#000000",
                   fillOpacity = 0,
                   label = "Study area")
    }
    
    # Add CLHS points if available
    if (!is.null(values$clhsPoints)) {
      clhs_sf <- st_as_sf(values$clhsPoints, coords = c("x", "y"), 
                         crs = crs(values$rasters, proj = TRUE))
      clhs_sf <- st_transform(clhs_sf, 4326)
      
      map <- map %>%
        addCircleMarkers(data = clhs_sf, radius = 5, 
                        color = "red",
                        popup = ~paste("CLHS Point ID:", ID),
                        label = ~paste("CLHS", ID),
                        group = "CLHS Points")
    }
    
    # Add Spatial Coverage points if available
    if (!is.null(values$spatcovPoints)) {
      points_4326 <- st_transform(values$spatcovPoints, 4326)
      
      map <- map %>%
        addCircleMarkers(data = points_4326, radius = 5, 
                        color = "blue",
                        popup = ~paste("Spatial Coverage Point ID:", ID),
                        label = ~paste("SC", ID),
                        group = "Spatial Coverage Points")
    }
    
    # Add layers control for points only
    groups <- c()
    if (!is.null(values$clhsPoints)) groups <- c(groups, "CLHS Points")
    if (!is.null(values$spatcovPoints)) groups <- c(groups, "Spatial Coverage Points")
    
    if (length(groups) > 0) {
      map <- map %>%
        addLayersControl(
          overlayGroups = groups,
          options = layersControlOptions(collapsed = FALSE)
        )
    }
    
    return(map)
  })
  # Download handler for data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$exportPrefix, "_", Sys.Date(), ".", input$exportMethod)
    },
    content = function(file) {
      export_data <- NULL
      
      # Determine which data to export
      if (input$exportData == "clhs" && !is.null(values$clhsPoints) && nrow(values$clhsPoints) > 0) {
        export_data <- st_as_sf(values$clhsPoints, coords = c("x", "y"), 
                               crs = crs(values$rasters, proj = TRUE))
      } else if (input$exportData == "spatcov" && !is.null(values$spatcovPoints) && nrow(values$spatcovPoints) > 0) {
        export_data <- values$spatcovPoints
      } else if (input$exportData == "both") {
        # Combine both datasets if available
        if (!is.null(values$clhsPoints) && nrow(values$clhsPoints) > 0 && 
            !is.null(values$spatcovPoints) && nrow(values$spatcovPoints) > 0) {
          
          clhs_sf <- st_as_sf(values$clhsPoints, coords = c("x", "y"), 
                             crs = crs(values$rasters, proj = TRUE))
          clhs_sf$method <- "CLHS"
          
          spatcov_sf <- values$spatcovPoints
          spatcov_sf$method <- "SpatialCoverage"
          
          # Match CRS if needed
          if (st_crs(clhs_sf) != st_crs(spatcov_sf)) {
            spatcov_sf <- st_transform(spatcov_sf, st_crs(clhs_sf))
          }
          
          # Combine (keeping only common columns)
          common_cols <- intersect(names(clhs_sf), names(spatcov_sf))
          if (length(common_cols) > 0) {
            export_data <- rbind(
              clhs_sf[, common_cols],
              spatcov_sf[, common_cols]
            )
          } else {
            # If no common columns, just bind with different column sets
            clhs_sf$type <- "CLHS"
            spatcov_sf$type <- "SpatialCoverage"
            export_data <- rbind(clhs_sf, spatcov_sf)
          }
        } else if (!is.null(values$clhsPoints) && nrow(values$clhsPoints) > 0) {
          export_data <- st_as_sf(values$clhsPoints, coords = c("x", "y"), 
                                 crs = crs(values$rasters, proj = TRUE))
          export_data$method <- "CLHS"
        } else if (!is.null(values$spatcovPoints) && nrow(values$spatcovPoints) > 0) {
          export_data <- values$spatcovPoints
          export_data$method <- "SpatialCoverage"
        }
      }
      
      if (is.null(export_data)) {
        showNotification("No data available to export!", type = "error")
        return()
      }
      
      # Export based on selected format
      if (input$exportMethod == "csv") {
        # For CSV, extract coordinates
        export_df <- cbind(
          as.data.frame(export_data),
          st_coordinates(export_data)
        )
        export_df$geometry <- NULL
        write.csv(export_df, file, row.names = FALSE)
      } else if (input$exportMethod == "gpkg") {
        st_write(export_data, file, quiet = TRUE)
      }
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)