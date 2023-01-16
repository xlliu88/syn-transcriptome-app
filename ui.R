## in this version: 
## 1. file input works
## 2. datatable and count table works
## 3. datatable sortable 
## 4. customizable gene names
## 5. customizable colors
## 6. uncoupled auto plotsize and manual plotsize adjust
## 7. default customized color to grays
## 8. image download with ggsave, all major format supported

## to do:
## 1. gene input re-arrangement - family/file input changable
## 2. make tables downloadable
## 3. re-align download buttons
## 4. show statistics on plot


## build user interface ------
ui <- fluidPage(
 sidebarLayout(
     sidebarPanel(
        
        ## Gene Input/selection
        fluidRow(
            column(6, 
                tabsetPanel(id = "input_method",
                     tabPanel("Search",
                             textInput(inputId = "searchterm",
                                       label = "Your Search:",
                                       placeholder = "type something here; support regex"),
                             actionButton(inputId = "searchButton",
                                          label = "Search")
                      ),
                
                    tabPanel("Gene_Family",
                                     selectInput(inputId = "genefamily",
                                                  label = "Select a Gene Family",
                                                  choices = c("", levels(as.factor(GENE_FAMILY$Family))),
                                                  selected = ""),
                
                                     selectInput(inputId = "subfamily",
                                                 label = "Select a subfamily",
                                                 choices = c("", levels(as.factor(GENE_FAMILY$Sub_Family))),
                                                 selected = "")
                      ),
      
                    tabPanel("File",
                                    fileInput(inputId = "file",
                                              label = "Upload your file",
                                              accept = c(".txt",".csv", "tsv"),
                                              placeholder = "accept '.csv', '.txt'", 
                                              multiple = F)
                      ),
                    tabPanel("Filter",
                             br(),
                             actionLink("build_filter", label = "Filter Genes"),
                             br(),
                             actionLink("filter_hist", label = "Filting History")
                    ))),
                        
            column(6, 
                   #div(style = "background-color:gray;", 
                   textAreaInput(inputId = "searchresult",
                                   label = "Search Result",
                                   value = "",
                                   placeholder = "your search result will show here",
                                   height = "100%"),
                   actionButton(inputId = "addresultButton",
                                  label = "Add Genes")
             )#)

        ),
       
        textAreaInput(inputId = "genelist",
                      label = "Selected Genes", 
                      value = "",
                      placeholder = gl.placeholder,
                      cols = 1,
                      resize = "both"),
        actionButton(inputId = "clear",
                     label = "Clear"),
        
        ## dataset selection
        fluidRow(
          column(4, 
                #radioButtons(inputId = "dataset",
                radioButtons(inputId = "dataset",
                               label = "Select Dataset",
                               choices = c("EIS", "LCM"),
                               selected = "EIS",
                               inline = F)),
          column(4,
                checkboxGroupInput(inputId = "genotype",
                                   label = "Select Genotypes", 
                                   choices = c("WT", "clv_triple"),
                                   selected = c("WT", "clv_triple"))),
          column(4, 
                checkboxGroupInput(inputId = "treatment",
                                   label = "Select Treatment", 
                                   choices = c("Ctrl","HsCLE2", "BCN"),
                                   selected = c("Ctrl", "HsCLE2", "BCN")))
         ),
       
        ## plot option selection
        fluidRow(
               column(6,
                        h5(strong("Plot Options")),
                        radioButtons(inputId = "plottype", 
                                     label = "ploy type",
                                     choices = c("bar", "point", "heatmap"),
                                     selected = "bar",
                                     inline = TRUE),
                      
                        radioButtons(inputId = "groupby",
                                     label = "group by: ",
                                     choices = c("Genotype", "Treatment"),
                                     selected = "Genotype",
                                     inline = TRUE),
                       
                        radioButtons(inputId = "showlegend",
                                     label = "legend position",
                                     choices = c("right", "top", "bottom", "none"),
                                     inline = TRUE),

                        checkboxInput(inputId = "uselog",
                                     label = "use log2 scale",
                                     value = TRUE),
                       
                        checkboxInput(inputId = "showstat",
                                     label = "show statistics",
                                     value = TRUE),
                        checkboxInput(inputId = "facet",
                                      label = "Facet by Gene",
                                      value = TRUE),
                        numericInput(inputId = "ncols",
                                       label = div(style = "font-size:12px;", "columns:"),
                                       value = NA,
                                       width = "100px")

                      ),
               
                      column(6,
                            sliderInput(inputId = "fontsize",
                                        label = "Font Size",
                                        min = 1,
                                        max = 30, value = 12),
                            sliderInput(inputId = "image_w",
                                        label = "Image Width",
                                        min = 1, max = 100, value = 50),
                            sliderInput(inputId = "image_h",
                                        label = "Image Hight",
                                        min = 1, max = 100, value = 50),
                            checkboxInput(inputId = "usemycolor",
                                      label = "Customize colors",
                                      value = F),
                            uiOutput(outputId = "colorpicker")
                      )
        
        ),
        hr(),
        div(style = " text-align: left;",
                numericInput(inputId = "cutoff", 
                 label = div(style = "font-size:12px;", "Count cutoff (Average)"),
                 value = 10, 
                 width = "100px")),
        textInput(inputId = "plottitle",
                  label = "Title: "),
        textInput(inputId = "xlab",
                  label = "x-axis title"),
        textInput(inputId = "ylab",
                  label = "y-axis title")
     ),
     
     mainPanel(
             tabsetPanel(
                     tabPanel("Plot",                     
                              h5(""),
                              fluidRow(
                                  column(2, 
                                       actionButton(inputId = "update_plot",
                                                       label = "Update Plot")),
                                  # column(2,
                                  #        h5("place holder")),
                                  column(10,
                                         downloadButton(outputId = "download_img",
                                                              label = "Download"),
                                         selectInput(inputId = "img_format",
                                                      label = "Format",
                                                      choices = DEVICES,
                                                      selected = "pdf",
                                                      multiple = F,
                                                     width = "80px")
                                         )

                              ),
                              hr(),
                              uiOutput("plot"),
                              hr(),
                              htmlOutput(outputId = "plotInfo")
                             ),   

                     tabPanel("Results",
                              h5(""),
                              fluidRow(
                                      column(4, 
                                              selectInput(inputId = "comparisons",
                                                          label = "Select Comparison(s)",
                                                          choices = c("show all", unname(unlist(P_LABELS))),
                                                          multiple = T)
                                                          ),
                                      column(8, 
                                             h5(strong("Download Data")),
                                             div(style="display:inline-block",
                                                 downloadButton(outputId = "download_curr",
                                                                label = "Current Table"), 
                                                 style="float:left"),
                                             div(style="display:inline-block",
                                                 downloadButton(outputId = "download_all",
                                                                label = "All Tables"), 
                                                 style="float:left")
                                             )
                                      ),
                              hr(),
                              uiOutput("results"),
                              ),
                     
                     tabPanel("Counts",
                              h5(""),
                              actionButton(inputId = "download_counts",
                                           label = "Download Table"),
                              hr(),
                              h5("Normalized Count Table"),
                              uiOutput("counts"))
              )),
         
             position = "left"
   ) #sidebarlayout
) #fluidpage
