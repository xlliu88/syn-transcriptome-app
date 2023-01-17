message("import utility module...")
## convert input strings to TAIR ID
## input should be a character vector
input2Id <- function(input) {
  if(input == "") return("")
  input <- geneNamesFromText(input)
  ids <- sapply(input, function(x) get(str_c(getIdType(x), "2Tair"))(x)) %>%
          unlist()
  no_ids <- input[which(is.na(ids))]
  if(length(no_ids) > 0) {
      ids <- ids[!is.na(ids)]
      msg <- sprintf("  Following input don't match to any gene: \n\t%s", 
                     str_c(no_ids, collapse = ", "))
      message(msg)
  }
  
  return(unique(ids))
}

geneNamesFromText <- function(text) {
  delim <- c(",", " ", "\t", "\n")
  genelist <- text
  for (d in delim) {
    if(!any(grepl(d, genelist))) next
    genelist <- sapply(genelist, function(x) trimws(unlist(str_split(x, d))))
    genelist <- as.vector(unname(unlist(genelist)))
  }
  genelist <- genelist[!genelist == " "]
  return(genelist)
}

## take a string and test the id type
## return one of:
##    "symbol", "tair", "seqref", "entrezid"
##    "void" if didn't find any match
getIdType <- function(id) {
  type_pattern <- c(tair = "^AT[1-5CM]G[0-9]{5}",
                       refseq = "^N[MPR]_[0-9]{5,12}",
                       entrezid = "^[0-9]{5,12}")
  
  idx <- which(sapply(type_pattern, function(x) grepl(x, id, ignore.case = T)))
  
  if(length(idx) > 0 ) return(names(type_pattern)[idx][1])
  
  symbols <- keys(org.At.tair.db, keytype = "SYMBOL")
  if(id %in% symbols) return("symbol")
  
  return("void")
  
}

tair2Tair <- function(tair) {
  return(toupper(tair))
}

void2Tair <- function(void) {
  return(NA)
}

## convert gene symbols to TAIR ID
symbol2Tair <- function(symbol) {
  id <- AnnotationDbi::select(org.At.tair.db, 
                              keys = symbol, 
                              columns = c("SYMBOL", "TAIR"),
                              keytype = "SYMBOL")
  id <- unique(id$TAIR)
  return(id)
}

## convert gene refseq_id to TAIR ID
refseq2Tair <- function(refseq) {
  id <- AnnotationDbi::select(org.At.tair.db, 
                              keys = refseq, 
                              columns = c("REFSEQ", "TAIR"),
                              keytype = "REFSEQ")
  id <- unique(id$TAIR)
  return(id)
}

## convert gene ENTREZID to TAIR ID
entrezid2Tair <- function(entrezid) {
  id <- AnnotationDbi::select(org.At.tair.db, 
                              keys = entrezid, 
                              columns = c("ENTREZID", "TAIR"),
                              keytype = "ENTREZID")
  id <- unique(id$TAIR)
  return(id)
}

importResults <- function(result_files, pattern = NA) {
  # result <- list()
  p_labels <- P_LABELS[names(P_LABELS) %in% names(result_files)]
  if(!is.na(pattern[1])) {
    p_labels <- p_labels[sapply(pattern, function(x) which(grepl(pattern, names(p_labels))))]
  }
  
  result <- lapply(names(p_labels), function(x) {
    file_name <- result_files[[x]]
    sep <- ifelse(grepl(".txt$", file_name), "\t", ",")
    res <- read_delim(file.path(result_files$path, result_files[[x]]),
                      col_types = cols(),
                      delim = sep, 
                      quote = "") %>%      
            arrange(Gene_id) %>%
            select(Gene_id, Gene_name, baseMean, log2FC = log2FoldChange, everything()) %>%
            #rename(log2FC = log2FoldChange) %>%
            mutate(baseMean = round(baseMean, 2),
                   log2FC = round(log2FC, 2),
                   lfcSE = round(lfcSE, 2),
                   stat = round(stat, 2),
                   pvalue = ifelse(is.na(pvalue), 1, pvalue),
                   pvalue = signif(pvalue, 1),
                   padj = ifelse(is.na(padj), 1, padj),
                   padj = signif(padj, 1))

    return(res)
  })
  names(result) <- names(p_labels)
  return(result)
}

## takes tair ID or gene name; 
## returns a dataframe with two columns: TairID SYMBOL
geneSearch <- function(term, db = org.At.tair.db, 
                       search_fields = c("TAIR", "SYMBOL", "GENENAME"), 
                       asis = F, ignore_case = T) {

    columns <- unique(c("TAIR", "SYMBOL", search_fields))
    ##if(grepl("^GO:\\d{5,9}", term)) term <- str_replace(term, "GO:", "")
    fun <- function(x) str_subset(keys(db, x), regex(term, ignore_case))
    key_res <- switch(asis + 1,
                      lapply(search_fields, fun),
                      list(term)
                     )
    
    fun2 <- function(x) AnnotationDbi::select(db, key_res[[x]], columns, search_fields[[x]])
    res2 <- lapply(1:length(key_res), fun2)
    
    res3 <- do.call("bind_rows", res2) %>%
      filter(!duplicated(TAIR)) %>%
      select(TAIR, SYMBOL) %>%
      droplevels()
   res3$SYMBOL[is.na(res3$SYMBOL)] <- ""
   
  return(res3)
}

## input can be a string or a vector
## format = c("any", "name", "tair")
geneParse <- function(input, search_fields) {
   
  res <- lapply(input, function(x) geneSearch(x, search_fields))
  res2 <- do.call(bind_rows, res) %>%
    filter(!duplicated(TAIR))
  
  return(res2)
}

genesInGo <- function(GOID) {
## takes an GO.ID
## return Genes annotated with given GO.ID
    GOID <- GOID %>%
        str_replace(., "GO:", "") %>%
        trimws()
    genes <- AnnotationDbi::select(org.At.tair.db, keys = GOID, keytype = "GO", columns = "TAIR", )
    return (genes)
}
## ui functions
## to create a ui for filtering
## under constructions
uiFilter <- function(filter_idx) {
  ids <- str_c(c("op", "filter_from", "filter_by", "criteria"), filter_idx)
  tags$div(
    tags$div(class = "divleft", selectInput(ids[1], label = "Op", choices = c("AND", "OR"))),
    tags$div(class = "divleft", selectInput(ids[2], label = "Filter_from", choices = c("Sample", "Compare"))),
    tags$div(class = "divleft", selectInput(ids[3], label = "Filter_by", choices = c("AND", "OR"))),
    tags$div(class = "divleft", textInput(ids[4], label = "Criteria", choices = c("AND", "OR"))))
  
}
