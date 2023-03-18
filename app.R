library(shiny)
library(shinythemes)
library(ape)
library(vegan)
library(plyr)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)
library(magrittr)
library(ggplot2)
library(ggpubr)
#library(ampvis)
library(BiocManager)
library(Biobase)
options(repos = BiocManager::repositories())

options(shiny.maxRequestSize=150*1024^2)
shinyApp(
  ui = fluidPage(
    navbarPage(
      theme = shinytheme("flatly"),
      "FGanalyst",
      tabPanel("Import Data",
               sidebarPanel(
                 fileInput("file1", "FunGuilds table:"),
                 fileInput("file2", "Mapping File:"),
                 
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("FunGuilds OTU table",
                            DT::dataTableOutput("otutable")
                   ),
                   tabPanel("Mapping File",
                            DT::dataTableOutput("mapping")
                   ),
                   tabPanel("FunGuilds TAX table",
                            DT::dataTableOutput("taxtable"))
                 )
               )
      ),
      tabPanel("Ploting",
               sidebarLayout(
                 sidebarPanel(
                   uiOutput("metachoices"),
                   uiOutput("guildchoices")
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel("FunGuilds-barplot",
                              sidebarPanel(
                                radioButtons("ftype", "Select the file type", c("pdf", "png", "jpeg")),
                                downloadButton('dwd','Download graph')),
                              mainPanel(plotOutput("fgps"),height = "20%")
                     ),
                     tabPanel("Richness plot",
                              plotOutput("rn"),height = "20%",
                              radioButtons("rtype", "Select the file type", c("pdf", "png", "jpeg")),
                              downloadButton('dwdrn','Download graph')
                     ),
                     tabPanel("Ordination plot",
                              plotOutput("od"),height = "20%",
                              radioButtons("otype", "Select the file type", c("pdf", "png", "jpeg")),
                              downloadButton('dwdod','Download graph')
                     ),    
                     tabPanel("Heatmap plot",
                              plotOutput("hm"),height = "20%",
                              radioButtons("htype", "Select the file type", c("pdf", "png", "jpeg")),
                              downloadButton('dwdhm','Download graph')
                     ))
                   
                 ))) ,tabPanel("User guide", 
                               sidebarPanel(textOutput("text"),  width = "100px", height="300px"),
                               mainPanel(uiOutput("Table")
                                         ,uiOutput('Mapping')))
    )),
  server = function(input, output) {
    
    url <- a("Download", href="https://drive.google.com/drive/folders/1shKSl1vuNCosuFu_hF5yUW7yj40Mj6TY?usp=share_link")
    output$Table <- renderUI({
      tagList("Example Files", url)})

    
    file_name <- reactive({
      inFile <- input$file1
      
      if (is.null(inFile))
        return(NULL)
      
      return (stringi::stri_extract_first(str = inFile$name, regex = ".*(?=\\.)"))
    })
    
    output$text <- renderText('User required to up load FunGuilds table file and Mapping file !! FunGuilds table file cannot have # or number in the begining of it header otherwise it cannot be read and cause errors"')
    #user metadata input
    dfmeta <- reactive({
      FGmeta <- read.table(input$file2$datapath, header=T,row.names=1, sep="\t")
    })
    
    output$metachoices <- renderUI({
      if (is.null(input$file2)){return()}
      dfmeta <- dfmeta()
      selectInput("meta", "Choose metadata", as.list(colnames(dfmeta)))
    })
    
    #user guilds input
    dffg <- reactive({
      FG <- read.table(input$file1$datapath,header=T,sep="\t",row.names=1)
      x <- select(FG, Trophic.Mode, Guild, Growth.Morphology)
      
    })
    
    output$guildchoices <- renderUI({
      if (is.null(input$file1)){return()}
      dffg <- dffg()
      selectInput("guild", "Choose guilds", as.list(colnames(dffg)))
    })
    
    
    #FunGuilds OTU table output
    output$otutable <- DT::renderDataTable({
      if (is.null(input$file1)){return()}
      
      FG <- read.table(input$file1$datapath,header=T,sep="\t",row.names=1)
      FGotus <- select(FG, -(Taxonomy:Citation.Source))
      FGotumat <- as(as.matrix(FGotus), "matrix")
      FGOTU = otu_table(FGotumat, taxa_are_rows = TRUE)
      DT::datatable(FGOTU, options = list(scrollX = T))
    })
    
    #Mapping File output
    output$mapping <- DT::renderDataTable({
      if (is.null(input$file2)){return()}
      
      readmapping <- read.table(input$file2$datapath, header=TRUE,row.names=1, sep="\t")
      sampleData <- sample_data(readmapping)
      DT::datatable(sampleData, options = list(scrollX = T))
    })
    
    
    #FunGuilds TAX table output
    output$taxtable <- DT::renderDataTable({
      if (is.null(input$file1)){return()}
      FG <- read.table(input$file1$datapath,header=T,sep="\t",row.names=1)
      FGtaxmat <- select(FG, Trophic.Mode, Guild, Growth.Morphology)
      FGtaxmat <- as(as.matrix(FGtaxmat),"matrix")
      FGTAX = tax_table(FGtaxmat)
      DT::datatable(FGTAX, options = list(scrollX = T))
    })
    
    
    #Construct FunGuilds-Phyloseq object
    output$fgps <-renderPlot({
      if (is.null(input$file1)){return()}
      
      FG <- read.table(input$file1$datapath,header=T,sep="\t",row.names=1)
      FGotus <- select(FG, -(Taxonomy:Citation.Source))
      FGotumat <- as(as.matrix(FGotus), "matrix")
      FGOTU <<- otu_table(FGotumat, taxa_are_rows = TRUE)
      
      readmapping <<- read.table(input$file2$datapath, header=TRUE,row.names=1, sep="\t")
      #readmappingfg <- select(readmapping, input$meta)
      sampleData <- sample_data(readmapping)
      
      FGtaxmat <- select(FG, Trophic.Mode, Guild, Growth.Morphology)
      FGtaxmat <- as(as.matrix(FGtaxmat),"matrix")
      FGTAX <<- tax_table(FGtaxmat)
      
      fgps <- phyloseq(FGOTU,FGTAX,sampleData)
      
      fgps.prune <<- prune_taxa(taxa_sums(fgps) > 1, fgps)
      
      fgps.prune.no.na <- subset_taxa(fgps.prune, Trophic.Mode!="-")
      

      
      
      
      #psTopNOTUs = names(sort(taxa_sums(fgps.prune.no.na), TRUE)[1:100])
      #ps100 = prune_taxa(psTopNOTUs, fgps.prune.no.na)
      #ps100.t = transform_sample_counts(ps100, function(x) x / sum(x) )
      #colSums(otu_table(ps100.t)[, 1:5] )
      
      #ps100.tm <- merge_samples(ps100.t, input$meta)
      #sample_data(ps100.tm)$Site <- factor(sample_names(ps100.tm))
      #ps100.tmg <- tax_glom(ps100.tm, input$guild)
      #ps100.tmg = transform_sample_counts(ps100.tmg, function(x) 100 * x/sum(x))
      
      #plot_bar(ps100.tmg, input$meta, fill = input$guild, title = input$meta) + theme(plot.title = element_text(hjust = 0.5, size = 22)) +
      #  theme(legend.title = element_text(size = 18), legend.text = element_text(size = 18), axis.text = element_text(size = 18), 
      #        axis.title = element_text(size = 18))
      fgps_plot <<- ggplot(data = psmelt(fgps.prune.no.na), mapping = aes_string(x = input$meta ,y = "Abundance", fill = input$guild )) + geom_bar(stat="identity", position="fill") + 
        scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", 
                                   "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
      fgps_plot
    }, height = 800)
    output$dwd <- downloadHandler(
      filename = function(){
        paste(file_name(),"barplot", input$ftype, sep = ".")
      },
      content = function(file){
        if(input$ftype=="pdf"){pdf(file)}
        if(input$ftype=="png"){png(file)}
        if(input$ftype=="jpeg"){jpeg(file)}
        plot(fgps_plot)
        dev.off()
      })
    output$rn <- renderPlot({plot_richness(fgps.prune , x = input$meta, color = input$meta, measures=c("Chao1") ) + geom_boxplot() + theme_bw() + ggtitle("Richness plot") + 
        theme(plot.title = element_text(hjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))})
    output$dwdrn <- downloadHandler(
      filename = function(){
        paste(file_name(),"Richness plot", input$rtype, sep = ".")
      },
      content = function(file){
        if(input$rtype=="pdf"){pdf(file)}
        if(input$rtype=="png"){png(file)}
        if(input$rtype=="jpeg"){jpeg(file)}
        plot({plot_richness(fgps.prune , x = input$meta, color = input$meta, measures=c("Chao1") ) + geom_boxplot() + theme_bw() + ggtitle("Richness plot") + 
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))})
        dev.off()
      })
    output$od <- renderPlot({
      fgps.prune.ord <- ordinate(fgps.prune, "DCA", "unifrac")
      od_plot <<- plot_ordination(fgps.prune, fgps.prune.ord, type = input$meta, color = input$meta) + theme_bw()
      od_plot
    })
    output$dwdod <- downloadHandler(
      filename = function(){
        paste(file_name(),"Ordination plot", input$otype, sep = ".")
      },
      content = function(file){
        if(input$otype=="pdf"){pdf(file)}
        if(input$otype=="png"){png(file)}
        if(input$otype=="jpeg"){jpeg(file)}
        plot(od_plot)
        dev.off()
      })
    output$hm <- renderPlot({
      sampleData <- sample_data(readmapping)
      
      physeq = phyloseq(FGOTU, FGTAX)
      
      random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
      
      physeq1 <<- merge_phyloseq(physeq,sampleData , random_tree)
      
      
      plot_heatmap(physeq1,
                   col_metadata = input$meta,
                   row_metadata = input$guild)
    })
    output$dwdhm <- downloadHandler(
      filename = function(){
        paste(file_name(),"Heatmap plot", input$htype, sep = ".")
      },
      content = function(file){
        if(input$htype=="pdf"){pdf(file)}
        if(input$htype=="png"){png(file)}
        if(input$htype=="jpeg"){jpeg(file)}
        plot(plot_heatmap(physeq1))
        dev.off()
      })
    
  }
)