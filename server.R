## 20210526: rollback to EIS and LCM seperate plots

## to run the server -----------------
server <- function(input, output, session){
  ## main code
  # GENES_FOR_PLOT <- tibble::tibble()
  
  ## get parameters
  re.genelist <- reactive(input2Id(input$genelist))
  re.title <- reactive(input$plottitle)
  re.xlab <- reactive(input$xlab)
  re.ylab <- reactive(input$ylab)
  re.gt <- reactive(input$genotype)
  re.trt <- reactive(input$treatment)
  re.groupby <- reactive(input$groupby)
  re.cutoff <- reactive(input$cutoff)
  re.ncol <- reactive(input$ncols)
  re.pty <- reactive(input$plottype)
  re.uselog <- reactive(input$uselog)
  re.showl <- reactive(input$showlegend) 
  re.font <- reactive(input$fontsize)
  re.family <- reactive(input$genefamily)
  re.subf <- reactive(input$subfamily)
  re.facet <- reactive({ switch(input$facet + 1, NULL, "Gene") })
  re.legend <- reactive(input$showlegend)
  re.usemycolor <- reactive(input$usemycolor)
  re.format <- reactive(input$img_format)
  re.imgsize <- reactive(list(w = input$img_save_w, 
                              h = input$img_save_h))
  

  re.colors <- reactive({
     fcts <- levels(re.dds()[[re.groupby()]])
     message(sprintf("!!! selected factors: %s", str_c(fcts, collapse = ", ")))
     res <- switch(input$usemycolor + 1, 
            #NULL,
            rev(generateColors(length(fcts), theme = "bw")),
            sapply(fcts, function(x) input[[sprintf("%s_color", x)]]) )
     message(sprintf("!!! colors: %s", str_c(res, collapse = ", ")))
     return(res)
    })
  
  # assign selected dds object
  re.dds <- reactive({
    switch(input$dataset,
            EIS = dds.EIS,
            LCM = dds.LCM) %>%
    ddsSlice(., factorName = "Genotype", factorLevel = input$genotype) %>% 
    ddsSlice(., factorName = "Treatment", factorLevel = input$treatment)
    })
 
   # filter out result tables
  re.res_tbls <- reactive({
    res_tbls <- switch(input$dataset,
                       EIS = EIS.RES,
                       LCM = LCM.RES)
    tbl_names <- names(res_tbls)

    res_tbls <- lapply(res_tbls, function(x) {x %>%
        filter(Gene_id %in% re.genelist())})
    names(res_tbls) <- tbl_names

    return(res_tbls)
  })
  
  # return selected table names, in the form of "wec", "wbc", etc...
  re.tbl_sel <- reactive({
      p_labels <- unlist(P_LABELS)
      sel <- input$comparisons
      
      if("show all" %in% sel) {
        tbl_sel <- names(p_labels)
      } else {
        tbl_sel <- names(p_labels[p_labels %in% sel])
      }
      return(tbl_sel)
  })
   
  # calculate plot size based on gene number and plot type
  # unit in px
  plotsize <- eventReactive(input$update_plot, {
      ngene <- length(input2Id(input$genelist))
      ncol <- re.ncol()
      ngt <- length(input$genotype)
      ntrt <- length(input$treatment)

      if(ngene > 0) {
          if(re.pty() %in% c("point", "bar")) {
                 if(!is.null(re.facet())) {
                     plot_col <- ifelse(is.na(ncol), ncolOptimize(ngene), ncol)
                     plot_col <- ifelse(plot_col > ngene, ngene, plot_col)
                     plot_row <- floor(ngene/plot_col)
                     plot_w <- 120 * (plot_col + 1.5)
                     plot_h <- 120 * (plot_row + 1.5)
                 } else {
                     plot_w <- 20 * ntrt * ngene + 200
                     plot_h <- 300 * ngt
                 }
            } else if(re.pty() == "heatmap") {
               plot_w <- 20 * 3 * ngt * ntrt * 3 + 100
               plot_h <- 16 * ngene * 1.3 + 20
            }
        } else {
            plot_w <- plot_h <- 300
        }

        return(list(w = plot_w,
                    h = plot_h))
  })
  
  # mannul adjustment of plot size
  plotsize_adj <- reactive({
     return(list(w = as.integer(plotsize()$w * input$image_w/50),
                 h = as.integer(plotsize()$h * input$image_h/50)
          ))
  })
  
  ###################################################################
  ## observe inputs #################################################
  ###################################################################
   
  ## genes search box
  # to initiate search
  # return a STRING in a format that,
  #    each row has one entry that contains
  #    gene ID and gene name, separated by \t
  observeEvent(input$searchButton,
    {
      term <- input$searchterm
      res <- data.frame()
      
      if(!grepl("[A-Za-a]", term)) res <- NULL
      else if(nchar(term) < 2) res <- NULL
      # else if(grepl("^GO:\\d{5,9}", term)) res <- genesInGo(term)
      else res <- geneSearch(term)
      
      if(is.null(res)) res_show <- "Your search does not yield any result."
      else {
           res_show <- str_c(res[[1]], res[[2]], sep = "\t")
           res_show <- str_c(res_show, collapse = "\n")
          
           if(nrow(res) > 200) {
             msg <- "Over 200 genes found.\nplease refine your search\n\n"
             res_show <- sprintf("%s%s", msg, res_show)
           }

      }
      
      updateTextAreaInput(session,
                          inputId = "searchresult",
                          label = sprintf("Genes Found (%s)", nrow(res)),
                          value = res_show)
  })
  
  # update gene list based on gene family
  # update subfamily selection options
  observe({
      req(input$genefamily)
      sub <- GENE_FAMILY$Sub_Family[GENE_FAMILY$Family == re.family()]
      #if(any(is.na(sub))) sub <- ""i
      
      subfamilies <- levels(as.factor(sub))
      nsub <- length(subfamilies)
      add_opt <- case_when(
                          nsub == 1 ~ "",
                          nsub == 2 ~ "Both",
                          nsub > 2 ~ "All"
                          )
      
      label <- sprintf("Select a SubFamily for %s", re.family())
      updateSelectInput(session, 
                        inputId = "subfamily",
                        label = label,
                        choices = c(subfamilies, add_opt),
                        selected = NULL)
  })
  
  # update family genes to gene list textbox
  observe({
    req(input$subfamily)
    if(re.subf() %in% c("Both", "All")) {
        gf_df <- filter(GENE_FAMILY, Family == re.family())
    } else {
        gf_df <- filter(GENE_FAMILY, Family == re.family(), Sub_Family == re.subf())
    }
    gf_df <- droplevels(gf_df)
    genes <- toupper(gf_df$Genomic_Locus_Tag)
    # if(any(is.na(genes))) genes <- ""
    
    res <- geneSearch(genes, asis = T, search_fields = "TAIR")
    res_show <- str_c(res[[1]], res[[2]], sep = "\t")
    res_show <- str_c(res_show, collapse = "\n")
    updateTextAreaInput(session, 
                      inputId = "searchresult",
                      label = sprintf("Genes Found (%d)", nrow(res)), 
                      value = res_show)  
  })
  
  # processing file upload
  observeEvent(input$file, {
    req(input$file)
    ext <- tools::file_ext(input$file$name)
    if(!ext %in% c("csv", "txt", "tsv")) {
       res_show <- "Invalid file type; Please upload a .csv or .txt file"
       ngene <- 0
    } else {
       text <- readLines(input$file$datapath, warn = F)
       text <- str_c(text, collapse = ", ")
       genes <- trimws(unique(input2Id(text)))
       res <- geneSearch(genes, search_fields = "TAIR", asis = T)
       ngene <- nrow(res)
       res_show <- str_c(res[[1]], res[[2]], sep = "\t")
       res_show <- str_c(res_show, collapse = "\n")
       if(ngene == 0) res_show <- "No gene found in the file."
    }
 
    updateTextAreaInput(session,
                        inputId = "searchresult",
                        label = sprintf("Genes Found in File (%d)", ngene),
                        value = res_show)
  })
  
  ## to add search result to input box
  observeEvent(input$addresultButton,  {   
          search_raw <- input$searchresult
          search_genes <- search_raw %>% trimws %>% str_split(., "\n") %>% unlist()
                          #unlist(str_split(trimws(search_raw), "\n"))
          search_genes <- trimws(search_genes[grepl("AT[0-5MC]G[0-9]{5}", search_genes, ignore.case = T)])
          search_genes <- ifelse(grepl("\t", search_genes), search_genes, str_c(search_genes, "\t"))
          if(is.null(search_genes)){ #} | is.na(search_genes) | search_genes == ""){
              genes_add <- ""
          } else if (length(search_genes) == 0) {
              genes_add <- ""
          } else if(is.na(search_genes[1])) {
             genes_add <- ""
          } else if(search_genes[1] == "") {
            genes_add <- ""
          } else {
              df_genes <-  unlist(str_split(search_genes, "\t"))%>%
                matrix(ncol = 2, byrow = T) %>%
                trimws %>%
                as_tibble()
              colnames(df_genes) <- c("Gene_id", "Gene_name")
              
              GENES_FOR_PLOT <<- bind_rows(GENES_FOR_PLOT, df_genes) %>%
                mutate(idx = seq_along(1:length(Gene_id))) %>%
                arrange(dplyr::desc(idx)) %>%
                filter(!duplicated(Gene_id)) %>%
                select("Gene_id", "Gene_name")
              
              genes_add <- df_genes$Gene_id

          }

          if(!input$genelist == "") {
            old_genes <- unlist(str_split(input$genelist, ","))
            genes <- c(old_genes, genes_add) %>%
              trimws()%>%
              unique()
          } else {
            genes <- genes_add
          }
          
          updateTextAreaInput(session,
                              inputId = "genelist",
                              label = sprintf("Selected Genes (%d)", length(genes)), 
                              value = str_c(genes, collapse = ", ")
                              )
        })
      
  # to clear input box
  observeEvent(input$clear,  {
    
            GENES_FOR_PLOT <<- tibble()
            updateTextAreaInput(session,
                                inputId = "genelist",
                                label = "Selected Genes (0)",
                                value = "")
    })
 
  # update treatment options based on dataset selection
  observe({
      dset <- input$dataset
      
      getfctlvls <- function(x) {
        sprintf("dds.%s", dset) %>% get() %>% .[[x]] %>% levels()
        # levels(get(sprintf("dds.%s", dset))[[x]])
      }
      trts <- getfctlvls("Treatment")
      gts <- getfctlvls("Genotype")
      
      updateCheckboxGroupInput(session,
                               inputId = "genotype",
                               label = "Genotype",
                               choices = gts,
                               selected = gts)
      updateCheckboxGroupInput(session,
                               inputId = "treatment",
                               label = "Treatment",
                               choices = trts,
                               selected = trts)
  })
 
                        
 # to generate colors select panel
 output$colorpicker <- renderUI({
   if(re.usemycolor()) {
     
        fcts <- switch(input$groupby,
                       Genotype = input$genotype,
                       Treatment = input$treatment)

        cols <- generateColors(length(fcts), theme = "bw") %>% rev()
        h5(strong("Choose Color"))
        lapply(1:length(fcts), function(x) colourInput(inputId = sprintf('%s_color', fcts[x]),
                                              label = fcts[x],
                                              value = cols[x],
                                              showColour = "both", #"background",
                                              palette =  "square", #"limited",
                                              allowTransparent = T,
                                              returnName = F))
                                                                                           
   }
 })   
 
 ## plot
 plotInput <- reactive({
   
    plotCounts2(re.dds(), 
                 re.genelist(), 
                 intgroup = re.groupby(), 
                 pty = re.pty(), 
                 use.log = re.uselog(),  
                 countthreshold = re.cutoff(),
                 ncol = re.ncol(),
                 facet = re.facet(),
                 font = re.font(), # * FONT_CONVERT,
                 main = re.title(), 
                 xlab = re.xlab(), 
                 ylab = re.ylab(),
                 cols = re.colors(),
                 legend = re.legend())
 })
 
 # update the plote
 observeEvent(input$update_plot, {
   output$plotx <- renderPlot({
     isolate(print(plotInput()))
   })
 })
 
 
# observeEvent(input$update_plot,{
#       updateActionButton(session,
#                          inputId = "update_plot",
#                          label = sprintf("Update Plot (%d)", input$update_plot))
#  }) 
 
  # render plot area
  output$plot <- renderUI({
     plotOutput("plotx",  width = plotsize_adj()$w, height = plotsize_adj()$h)
  })
  
  ## download plot
  output$download_img <- downloadHandler(
    filename = function() { 
      label <- now() %>% 
        as.character %>% 
        str_replace_all(., ":", "") %>% 
        str_replace(., " ", "_")
      
      fn <- sprintf("Rplot_%s.%s", 
                    label,
                    re.format())
      return(fn)
      },
    
    content = function(file) {
      message("downloading...", file)
      ggsave(file,
             plot = plotInput(),
             device = re.format(),
             width = plotsize_adj()$w/DPI,
             height = plotsize_adj()$h/DPI,
             units = PLOT_UNIT,
             scale = 0.5, #
             1/FONT_CONVERT,
             dpi = DPI)
      }
  )
  
  output$plotInfo <- renderUI({
    # "some text here"
    HTML(paste("plot size: ", 
               sprintf("%d x %d (px)", plotsize_adj()$w, plotsize_adj()$h),
               sprintf("%.2f x %.2f (in)", plotsize_adj()$w/DPI, plotsize_adj()$h/DPI),
               sep = "<br/>"))
  })
  
   
  ## result table output
  observeEvent(input$comparisons, {
     # d <- re.res2show()
     lapply(re.tbl_sel(), 
             function(x) {
                      data_table <- re.res_tbls()[[x]] %>%
                        datatable(rownames = F,
                                  options = list(autoWidth = TRUE)) 
                     output[[x]] <- renderDT(data_table) 
      })
  }) 
  
  ## update select dropdown menu
  observe({
      res_names <- names(re.res_tbls())
      p_labels <- unname(unlist(P_LABELS[res_names]))
      
      updateSelectInput(session,
                        inputId = "comparisons",
                        label = "Select Comparison(s)",
                        choices = c("show all", p_labels)
                        )
  })
  
  # render result table outputs
  output$results <- renderUI({
    
      ntabs <- length(re.tbl_sel())
      p_labels <- P_LABELS[re.tbl_sel()]
      msg_0 <- "No comparison selected"
      res_tabs <- list()
      
      if(ntabs == 0) { # when no selection
            res_tabs[[1]] <- tabPanel(" ",
                                      h5(msg_0),
                                      hr())
      } else {
            for (t in 1:ntabs) {
                comps <- p_labels[[t]]
                cap <- sprintf("%s Result", comps)
                res_tabs[[t]] <- tabPanel(comps,
                                          h5(strong(cap)),
                                          DTOutput(re.tbl_sel()[t]),
                                          hr())
            }
      }
      
      do.call(tabsetPanel, res_tabs) 

  })
  
 ## prepare count table
 observe({
 
   dds <- re.dds()
     
   if(is.null(dds)) {
      output[["cnts"]] <- renderDT({data.frame(nrow = 0, ncol = 0)})
   } else {

       dds <- dds[row.names(dds) %in% re.genelist(), ]
       count_table <- counts(dds, normalized = T)
       count_table <- round(count_table, 2)
       output[["cnts"]] <- renderDT(count_table)
      }

 })
 
 ## render count table UI
 output$counts <- renderUI({
     DTOutput("cnts")
   })

}
