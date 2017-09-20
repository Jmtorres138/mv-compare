
"%&%" <- function(a,b) paste0(a,b)

library("data.table")
library("dplyr")
library("Homo.sapiens")

server.dir <- "/well/got2d/jason/"
work.dir <- server.dir %&% "projects/mv-compare/"
input.dir <- work.dir %&% "input_files/"
sig.results.file <- input.dir %&% "metaxcan_full_results_0.01_only.txt"
sig.df <- fread(sig.results.file)



# Create reference of Gene Ontology (GO) terms for each significant gene

get_go_info <- function(genename){
  goid.vec <- tryCatch({(select(Homo.sapiens,keys=genename,keytype="SYMBOL",columns="GO"))$GO},
                       warning=function(war){return(NA)},
                       error=function(err){
                         return(NA)
                       })

  if (is.na(goid.vec)==FALSE){
    go.df <- select(GO.db,keys=goid.vec,columns=columns(GO.db)[1:4])
    gene_name <- rep(genename,dim(go.df)[1])
    go.df <- cbind(gene_name,go.df)
  } else{
    go.df <- data.frame("gene_name"=NA,"DEFINITION"=NA,"GOID"=NA,"ONTOLOGY"=NA,"TERM"=NA)
  }
  return(go.df)
}

get_go_ids <- function(ensid){
  goid.vec <- tryCatch({(select(Homo.sapiens,keys=ensid,keytype="ENSEMBL",columns="GO"))$GO},
                       warning=function(war){return(NA)},
                       error=function(err){
                         return(NA)
                       })
  return(goid.vec)
}

build_sig_go_df <- function(){
  out.df <- c()
  pb <- txtProgressBar(min=0,max=length(unique(sig.df$gene)),initial=1,style=3)
  for (i in 1:length(unique(sig.df$gene_name))){
    setTxtProgressBar(pb,i)
    gene <- unique(sig.df$gene_name)[i]
    print(gene)
    df <- suppressMessages(get_go_info(gene))
    out.df <- rbind(out.df,df)
  }
  return(out.df)
}

sig.go.df <- build_sig_go_df()
sig.go.df$gene_name <- as.character(sig.go.df$gene_name)
saveRDS(sig.go.df,work.dir%&%"sig.go.df.RDS")


# GO Enrichment Analysis

#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")
library("topGO")
library("GO.db")

go_enrich_metab <- function(metab,ont="BP"){
  # Performs GO enrichment on the set of genes significantly associated with a specified metablite
  sig.genes <- unique(filter(sig.df,metabolite==metab)$gene)
  ensgenes <- unique(keys(Homo.sapiens,keytype="ENSEMBL"))
  geneList <- as.factor(as.integer(ensgenes %in% sig.genes))
  names(geneList) <- ensgenes
  # Build topGOdata object
  GOdata <- new("topGOdata", ontology = ont,
                allGenes = geneList,annotationFun = annFUN.org,
                mapping="org.Hs.eg.db",ID="ensembl")
  result.fisher <- runTest(GOdata,algorithm="classic",statistic="fisher")
  pval.fisher <- sort(score(result.fisher))
  pval.adjust <- p.adjust(pval.fisher,method="BH")
  sig.fisher <- pval.fisher[pval.fisher < 0.01]
  out.df <- data.frame("metabolite"=rep(metab,length(sig.fisher)),
                       "ontology"=rep(ont,length(sig.fisher)),
                       "GO"=names(sig.fisher),
                       "pval"=sig.fisher,
                       "pval.adj"=pval.adjust[1:length(sig.fisher)],
                       stringsAsFactors = FALSE)
  return(out.df)
}

build_go_results_df <- function(){
  metab.vec <- unique(sig.df$metabolite)
  out.df <- c()
  pb <- txtProgressBar(min=0,max=length(metab.vec),initial=1,style=3)
  for (i in 1:length(metab.vec)){
    setTxtProgressBar(pb,i)
    metab <- metab.vec[i]
    print(metab)
    bp.df <- suppressMessages(go_enrich_metab(metab,ont="BP"))
    mf.df <- suppressMessages(go_enrich_metab(metab,ont="MF"))
    cc.df <- suppressMessages(go_enrich_metab(metab,ont="CC"))
    out.df <- rbind(out.df,bp.df)
    out.df <- rbind(out.df,mf.df)
    out.df <- rbind(out.df,cc.df)
  }
  return(out.df)
}

go.results.df <- build_go_results_df()
write.table(x=go.results.df,file=work.dir%&%"metab.sig-go.results.txt",sep="\t",quote=FALSE,row.names=F)

# Useful function genesInTerm(GOdata, sel.terms) $ get the gene annotations in a set GO terms
# Term Stat gets useful GO statistics: termStat(GOdata, sel.terms)
