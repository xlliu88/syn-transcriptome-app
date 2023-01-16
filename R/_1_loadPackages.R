
## import libraries
message("loading packages.....")
suppressMessages(library(RColorBrewer))
suppressMessages(library(stringr))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(pheatmap))
suppressMessages(library(VennDiagram))
suppressMessages(library(eulerr))
#suppressMessages(library(shiny))
suppressMessages(library(DT))

suppressMessages(library(LSD))           # a colorful plot tool

suppressMessages(library(DESeq2))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(biomaRt))
suppressMessages(library(org.At.tair.db))

#suppressMessages(library(tidyr))
#suppressMessages(library(tidyverse))
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(readr))
suppressMessages(library(colourpicker))
suppressMessages(library(lubridate))
suppressMessages(library(svglite))

select <- dplyr::select
desc <- dplyr::desc
dbi.select <- AnnotationDbi::select
tibble <- tibble::tibble()
