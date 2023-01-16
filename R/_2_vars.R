## to set global variables
message("Setting global variables...")

gl.placeholder <- "accept list of genes seperated by\n\tcomma (,), space, or tab.\n"

### constants
 GENES_FOR_PLOT <- tibble::tibble()
DDS <- NULL
RES_TBLS <- list()
DEVICES <- c("pdf", "svg", "eps", "ps", "png",  "jpeg", "tiff")

### plot parameters
DPI <- 72
FONT_CONVERT <- 72/DPI
PLOT_UNIT <- "in"

## pairwise comparison name mapping
P_LABELS <- list(wec="WT_CLE2.vs.WT_Ctrl",
                wbc="WT_BCN.vs.WT_Ctrl",
                
                vec="clv_CLE2.vs.clv_Ctrl",
                vbc="clv_BCN.vs.clv_Ctrl",
                
                vwc="clv_Ctrl.vs.WT_Ctrl",
                vwe="clv_CLE2.vs.WT_CLE2",
                vwb="clv_BCN.vs.WT_BCN",
                
                clvCLE = "interaction_clvCLE2",
                clvBCN = "interaction_clvBCN")
					  
### count tables, condition list, and annotaiton file

### file path and names	
CNT_FILES <- list(path="data/salmon_counts",
              EIS.cnts = "EIS_Counts_Per_Sample.txt",
      				EIS.conditions = "EIS_samples.txt",
      				LCM.cnts = "LCM_Counts_Per_Sample.txt",
      				LCM.conditions = "LCM_samples.txt")

A_FILES <- list(path="data",
          			  gff3 = "Arabidopsis_thaliana.TAIR10.42.gff3",
          			  annotation = "Arabidopsis_thaliana.TAIR10.42.annotation",
          			  genefamily = "gene_families_20201221_update.txt")
			  
## EIS results files, without sample EIS_7-5 (HsCLE2 treated clv triple mutant, rep 1)
EIS_RES_FILES <- list(path = "data/deseq2output",
                        
                        wec = "EIS.WT_HsCLE2.vs.WT_Ctrl.csv", 
                        wbc = "EIS.WT_BCN.vs.WT_Ctrl.csv",
                                         
                        vec = "EIS.ccr_HsCLE2.vs.ccr_Ctrl.csv", 
                        vbc = "EIS.ccr_BCN.vs.ccr_Ctrl.csv",
                         
                        vwc = "EIS.ccr_Ctrl.vs.WT_Ctrl.csv", 
                        vwe = "EIS.ccr_HsCLE2.vs.WT_HsCLE2.csv",  
                        vwb = "EIS.ccr_BCN.vs.WT_BCN.csv", 
                                       
                        clvCLE = "EIS.interaction_ccr.HsCLE2.csv",
                        clvBCN = "EIS.interaction_ccr.BCN.csv")
                    
## LCM results files
LCM_RES_FILES <- list(path="data/deseq2output",
            					wbc = "LCM.WT_BCN.vs.WT_Ctrl.csv",
            					vbc = "LCM.ccr_BCN.vs.ccr_Ctrl.csv",
            					vwc = "LCM.ccr_Ctrl.vs.WT_Ctrl.csv",
            					vwb = "LCM.ccr_BCN.vs.WT_BCN.csv",
            					clvBCN = "LCM.interaction_ccr.BCN.csv")

