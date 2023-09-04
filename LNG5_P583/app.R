library(shinydashboard)
library(htmlwidgets)
library(shiny)
library(leaflet)
library(leaflet.extras)
library(ggplot2)
library(zeallot)
library(shinyBS)
library(shinyWidgets)
library(RColorBrewer)

coords.list <- setNames(lapply(paste0("data/coords", c(1:8)), function(f) readRDS(file = f)), nm = paste0(c(1:8)))
expr.data.list <- setNames(lapply(paste0("data/expr.data", c(1:8)), function(f) readRDS(file = f)), nm = paste0(c(1:8)))
expr.data.list <- lapply(expr.data.list, as.matrix)
cell.data.list <- setNames(lapply(paste0("data/cell.data", c(1:8)), function(f) readRDS(file = f)), nm = paste0(c(1:8)))
gene.choices <- readRDS("data/gene.choices")


palette.select <- list("Spectral" = rev(brewer.pal(n = 9, name = "Spectral")),
                    "viridis" = "viridis", "magma" = "magma", "Blues" = "Blues", "Greens" = "Greens",
                    "Oranges" = "Oranges", "OrRd" = "OrRd", "Purples" = "Purples", "Reds" = "Reds",
                    "Jet" = c("darkblue", "cyan", "yellow", "red", "darkred"), 
                    "GrRdBlk" = c("lightgray", "mistyrose", "red", "darkred", "black"),
                    "GrRd" = c("lightgray", "mistyrose", "red", "darkred"))



ui <- dashboardPage(
  
  header = dashboardHeader(tags$li(class = "dropdown",
                                   tags$style(".main-header {max-height: 60px}"),
                                   tags$style(".main-header .logo {height: 60px}")),
                           title = "DiscovAir \ndata explorer"),
  
  sidebar = dashboardSidebar(
    tags$style(".left-side, .main-sidebar {padding-top: 60px}"),
    width = 250,
    column(width = 12,
           selectizeInput(
             inputId = "gene",
             label = "gene",
             choices = NULL),
           selectInput(
             inputId = "cell",
             label = "cell",
             choices = rownames(cell.data.list[[1]]),
             selected = "AT1"),
           sliderInput(
             inputId = "alpha",
             label = "Opacity", value = 1,
             min = 0, max = 1, step = 0.01
           ),
           sliderInput(
             inputId = "maxcutoff",
             label = "Max cutoff", value = 1,
             min = 0.5, max = 1, step = 0.01
           ),
           fluidRow(
             column(width = 4, checkboxInput(
               inputId = "revpal", 
               label = "reverse palette", 
               value = FALSE)),
             column(width = 4, checkboxInput(
               inputId = "scalealpha", 
               label = "opacity gradient", 
               value = FALSE))
           ),
           actionBttn(inputId = "go2", style = "bordered", size = "s", label = "Value histogram"),
           selectInput(
             inputId = "pal",
             label = "palette",
             choices = names(palette.select),
             selected = "Jet")
    )
    
  ),
  
  body = dashboardBody(

    fluidRow(
      tabBox(
        title = "",
        id = "tabset1", width = "80%",
        tabPanel(
          value = 1, 
          title = "section 1", 
          column(width = 12, leafletOutput("map1", width = "100%", height = "1000px"))),
        tabPanel(
          value = 2, 
          title = "section 2", 
          column(width = 12, leafletOutput("map2", width = "100%", height = "1000px"))),
        tabPanel(
          value = 3, 
          title = "section 3", 
          column(width = 12, leafletOutput("map3", width = "100%", height = "1000px"))),
        tabPanel(
          value = 4, 
          title = "section 4", 
          column(width = 12, leafletOutput("map4", width = "100%", height = "1000px"))),
        tabPanel(
          value = 5, 
          title = "section 5", 
          column(width = 12, leafletOutput("map5", width = "100%", height = "1000px"))),
        tabPanel(
          value = 6, 
          title = "section 6", 
          column(width = 12, leafletOutput("map6", width = "100%", height = "1000px"))),
        tabPanel(
          value = 7, 
          title = "section 7", 
          column(width = 12, leafletOutput("map7", width = "100%", height = "1000px"))),
        tabPanel(
          value = 8, 
          title = "section 8", 
          column(width = 12, leafletOutput("map8", width = "100%", height = "1000px")))
      ),
      bsModal(
        id = "modalExample2", 
        title = "histogram", "go2", 
        size = "large", 
        plotOutput("histogram"))
    )
  )
)

server <- function(input, output, session) {
  
  updateSelectizeInput(session, 'gene', selected = "SCGB3A1", choices = gene.choices, server = TRUE)
  
  # Create reactive value (responds on input)
  rv <- reactiveValues(lastBtn = "gene")
  rv2 <- reactiveValues(curgene = "SCGB3A1")
  
  
  # Set last button to "cell" when a cell type is selected
  observeEvent(input$cell, {
    if (input$cell > 0 ) {
      rv$lastBtn = "cell"
    }
  })
  
  
  # Set last button to "gene" when a gene is selected
  observeEvent(input$gene, {
    if (input$gene > 0 ) {
      rv$lastBtn = "gene"
    }
  })
  
  
  # Generate data
  get_data <- function() {
    # Selected tab
    i <- as.integer(input$tabset1)
    # Get coordinates for selected tab
    data <- coords.list[[paste0(i)]]
    # Collect data based on input type; gene, factor or cell type
    if (rv$lastBtn == "gene") {
      variable <- input$gene
      if (input$gene %in% gene.choices) {
        rv2$curgene <- variable
      } else {
        variable <- rv2$curgene
      }
      data[, paste0(variable, "_raw")] <- expr.data.list[[paste0(i)]][variable, ]
      data[, variable] <- expr.data.list[[paste0(i)]][variable, ]
    } else if (rv$lastBtn == "cell") {
      variable <- input$cell
      data[, paste0(variable, "_raw")] <- cell.data.list[[paste0(i)]][variable, ]
      data[, variable] <- cell.data.list[[paste0(i)]][variable, ]
    }
    
    # Clip data if max cut-off if specified, i.e. if max cut-off is lower than 1
    if (input$maxcutoff < 1) {
      data[, variable][data[, variable] > quantile(data[, variable], probs = input$maxcutoff)] <- quantile(data[, variable], probs = input$maxcutoff)
    }
    # Create a column with opacity values for each spot scaled by the value to be visualized
    # The maximum opcaity value is set by alpha 
    # If no scaling is applied, all spots will have the same opacity specified by alpha
    if (input$scalealpha) {
      data[, "alpha"] <- scales::rescale(data[, variable], to = c(0, input$alpha))
    } else {
      data[, "alpha"] <- input$alpha
    }
    # reverse color palette if specified
    if (input$revpal %in% c(T, F)) {
      pal <- colorNumeric(
        palette = palette.select[[input$pal]],
        alpha = TRUE, 
        reverse = input$revpal,
        domain = data[, variable])
      pal_rev <- colorNumeric(
        palette = palette.select[[input$pal]],
        alpha = TRUE, reverse = !input$revpal,
        domain = data[, variable])
    }
    return(list(data, variable, pal, pal_rev))
  }
  
  
  # render a histogram pop-up for the selected values to be visualzied
  output$histogram <- renderPlot({
    c(data, variable, pal, pal_rev) %<-% get_data()
    ggplot(data, aes_string(paste0("`", variable, "_raw`"))) + 
      geom_histogram(fill = "lightgray", color = "black", bins = 30) +
      theme(panel.background = element_blank()) +
      labs(x = variable) +
      geom_vline(aes(xintercept = input$maxcutoff*max(data[, paste0(variable, "_raw")]), color = "cutoff_threshold"), linetype = "longdash") +
      scale_color_manual(values = c("cutoff_threshold" = "black"))
  })

  # Create leaflets for two tabs and add tiles
  # Output id is map1 and map2
  lapply(c(1:8), function(i) {
    output[[paste0("map", i)]] <- renderLeaflet({
      m <- leaflet(
                    options = leafletOptions(crs = leafletCRS(crsClass = "L.CRS.Simple"), preferCanvas = TRUE)) %>% 
        addTiles(urlTemplate = paste0("map", i, "/tiles/{z}/{x}/{y}.png"), 
                 options = tileOptions(continuousWorld = FALSE, 
                               tms = TRUE, 
                               tileSize = "256", 
                               minZoom = "0", 
                               maxZoom = "4")) %>%
        setView(lat = 0, lng = 0, zoom = 0) %>%
        setMaxBounds(lng1 = 0, lat1 = 0, lng2 = 596, lat2 = -601) %>% # TODO: find a way to automatically set bounds around tiled image
        addStyleEditor() %>%
        addDrawToolbar(targetGroup = "draw",
                       polylineOptions = FALSE,
                       circleOptions = FALSE,
                       markerOptions = FALSE,
                       editOptions = editToolbarOptions(selectedPathOptions = selectedPathOptions())) %>%
        fitBounds(lng1 = 0, lat1 = 0, lng2 = 596, lat2 = -601)
      return(m)
    })
  })
  

  # Map circles on map 1
  observeEvent(ifelse(req(input$tabset1 == "1"), list(input$map1_zoom, 
                                                      input$gene, 
                                                      input$cell,
                                                      input$scalealpha,
                                                      input$alpha,
                                                      input$maxcutoff,
                                                      input$revpal,
                                                      input$pal), list()), {
    i <- 1
    c(data, variable, pal, pal_rev) %<-% get_data()
    proxy <- leafletProxy(paste0("map", i))
    radi = (1.7*2.5)*2^(input[[paste0("map", i, "_zoom")]] - 1)
    proxy %>% 
      clearMarkers() %>% 
      clearControls() %>%
      addCircleMarkers(radius = radi, 
                       lat = -coords.list[[paste0(i)]]$pixel_y, 
                       lng = coords.list[[paste0(i)]]$pixel_x, 
                       fillColor = pal(data[, variable]), 
                       group = "Spots", 
                       fillOpacity = data[, "alpha"], 
                       stroke = FALSE) %>%
      addLegend(pal = pal_rev, 
                className = "legendbox",
                position = "topright", 
                title = variable, 
                values = data[, variable], 
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  })

  # Map circles on map 2
  observeEvent(ifelse(req(input$tabset1 == "2"), 
                      list(input$map2_zoom, 
                           input$cell,
                           input$gene, 
                           input$scalealpha,
                           input$alpha,
                           input$maxcutoff,
                           input$revpal,
                           input$pal), 
                      list()), {
    i <- 2
    c(data, variable, pal, pal_rev) %<-% get_data()
    proxy <- leafletProxy(paste0("map", i))
    radi = (1.7*2.5)*2^(input[[paste0("map", i, "_zoom")]] - 1)
    proxy %>% 
      clearMarkers() %>% 
      clearControls() %>%
      addCircleMarkers(radius = radi, 
                       lat = -coords.list[[paste0(i)]]$pixel_y, 
                       lng = coords.list[[paste0(i)]]$pixel_x, 
                       fillColor = pal(data[, variable]), 
                       group = "Spots", 
                       fillOpacity = data[, "alpha"], 
                       stroke = FALSE) %>%
      addLegend(pal = pal_rev, 
                className = "legendbox",
                position = "topright", 
                title = variable, 
                values = data[, variable], 
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  })
  
  # Map circles on map 3
  observeEvent(ifelse(req(input$tabset1 == "3"), 
                      list(input$map3_zoom, 
                           input$gene, 
                           input$cell,
                           input$scalealpha,
                           input$alpha,
                           input$maxcutoff,
                           input$revpal,
                           input$pal), 
                      list()), {
    i <- 3
    c(data, variable, pal, pal_rev) %<-% get_data()
    proxy <- leafletProxy(paste0("map", i))
    radi = (1.7*2.5)*2^(input[[paste0("map", i, "_zoom")]] - 1)
    proxy %>% 
      clearMarkers() %>% 
      clearControls() %>%
      addCircleMarkers(radius = radi, 
                       lat = -coords.list[[paste0(i)]]$pixel_y, 
                       lng = coords.list[[paste0(i)]]$pixel_x, 
                       fillColor = pal(data[, variable]), 
                       group = "Spots", 
                       fillOpacity = data[, "alpha"], 
                       stroke = FALSE) %>%
      addLegend(pal = pal_rev, 
                className = "legendbox",
                position = "topright", 
                title = variable, 
                values = data[, variable], 
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  })
  
  # Map circles on map 4
  observeEvent(ifelse(req(input$tabset1 == "4"), 
                      list(input$map4_zoom, 
                           input$gene, 
                           input$cell,
                           input$scalealpha,
                           input$alpha,
                           input$maxcutoff,
                           input$revpal,
                           input$pal), 
                      list()), {
                        i <- 4
                        c(data, variable, pal, pal_rev) %<-% get_data()
                        proxy <- leafletProxy(paste0("map", i))
                        radi = (1.7*2.5)*2^(input[[paste0("map", i, "_zoom")]] - 1)
                        proxy %>% 
                          clearMarkers() %>% 
                          clearControls() %>%
                          addCircleMarkers(radius = radi, 
                                           lat = -coords.list[[paste0(i)]]$pixel_y, 
                                           lng = coords.list[[paste0(i)]]$pixel_x, 
                                           fillColor = pal(data[, variable]), 
                                           group = "Spots", 
                                           fillOpacity = data[, "alpha"], 
                                           stroke = FALSE) %>%
                          addLegend(pal = pal_rev, 
                                    className = "legendbox",
                                    position = "topright", 
                                    title = variable, 
                                    values = data[, variable], 
                                    labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
                      })
  
  # Map circles on map 5
  observeEvent(ifelse(req(input$tabset1 == "5"), 
                      list(input$map5_zoom, 
                           #input$var, 
                           input$gene, 
                           input$cell,
                           input$scalealpha,
                           input$alpha,
                           input$maxcutoff,
                           input$revpal,
                           input$pal), 
                      list()), {
    i <- 5
    c(data, variable, pal, pal_rev) %<-% get_data()
    proxy <- leafletProxy(paste0("map", i))
    radi = (1.7*2.5)*2^(input[[paste0("map", i, "_zoom")]] - 1)
    proxy %>% 
      clearMarkers() %>% 
      clearControls() %>%
      addCircleMarkers(radius = radi, 
                       lat = -coords.list[[paste0(i)]]$pixel_y, 
                       lng = coords.list[[paste0(i)]]$pixel_x, 
                       fillColor = pal(data[, variable]), 
                       group = "Spots", 
                       fillOpacity = data[, "alpha"], 
                       stroke = FALSE) %>%
      addLegend(pal = pal_rev, 
                className = "legendbox",
                position = "topright", 
                title = variable, 
                values = data[, variable], 
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  })
  
  # Map circles on map 6
  observeEvent(ifelse(req(input$tabset1 == "6"), 
                      list(input$map6_zoom, 
                           input$gene, 
                           input$cell,
                           input$scalealpha,
                           input$alpha,
                           input$maxcutoff,
                           input$revpal,
                           input$pal), 
                      list()), {
    i <- 6
    c(data, variable, pal, pal_rev) %<-% get_data()
    proxy <- leafletProxy(paste0("map", i))
    radi = (1.7*2.5)*2^(input[[paste0("map", i, "_zoom")]] - 1)
    proxy %>% 
      clearMarkers() %>% 
      clearControls() %>%
      addCircleMarkers(radius = radi, 
                       lat = -coords.list[[paste0(i)]]$pixel_y, 
                       lng = coords.list[[paste0(i)]]$pixel_x, 
                       fillColor = pal(data[, variable]), 
                       group = "Spots", 
                       fillOpacity = data[, "alpha"], 
                       stroke = FALSE) %>%
      addLegend(pal = pal_rev, 
                className = "legendbox",
                position = "topright", 
                title = variable, 
                values = data[, variable], 
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  })
  
  # Map circles on map 7
  observeEvent(ifelse(req(input$tabset1 == "7"), 
                      list(input$map7_zoom, 
                           input$gene, 
                           input$cell,
                           input$scalealpha,
                           input$alpha,
                           input$maxcutoff,
                           input$revpal,
                           input$pal), 
                      list()), {
    i <- 7
    c(data, variable, pal, pal_rev) %<-% get_data()
    proxy <- leafletProxy(paste0("map", i))
    radi = (1.7*2.5)*2^(input[[paste0("map", i, "_zoom")]] - 1)
    proxy %>% 
      clearMarkers() %>% 
      clearControls() %>%
      addCircleMarkers(radius = radi, 
                       lat = -coords.list[[paste0(i)]]$pixel_y, 
                       lng = coords.list[[paste0(i)]]$pixel_x, 
                       fillColor = pal(data[, variable]), 
                       group = "Spots", 
                       fillOpacity = data[, "alpha"], 
                       stroke = FALSE) %>%
      addLegend(pal = pal_rev, 
                className = "legendbox",
                position = "topright", 
                title = variable, 
                values = data[, variable], 
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  })
  
  # Map circles on map 8
  observeEvent(ifelse(req(input$tabset1 == "8"), 
                      list(input$map8_zoom, 
                           input$gene, 
                           input$cell,
                           input$scalealpha,
                           input$alpha,
                           input$maxcutoff,
                           input$revpal,
                           input$pal), 
                      list()), {
    i <- 8
    c(data, variable, pal, pal_rev) %<-% get_data()
    proxy <- leafletProxy(paste0("map", i))
    radi = (1.7*2.5)*2^(input[[paste0("map", i, "_zoom")]] - 1)
    proxy %>% 
      clearMarkers() %>% 
      clearControls() %>%
      addCircleMarkers(radius = radi, 
                       lat = -coords.list[[paste0(i)]]$pixel_y, 
                       lng = coords.list[[paste0(i)]]$pixel_x, 
                       fillColor = pal(data[, variable]), 
                       group = "Spots", 
                       fillOpacity = data[, "alpha"], 
                       stroke = FALSE) %>%
      addLegend(pal = pal_rev, 
                className = "legendbox",
                position = "topright", 
                title = variable, 
                values = data[, variable], 
                labFormat = labelFormat(transform = function(x) sort(x, decreasing = TRUE)))
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
