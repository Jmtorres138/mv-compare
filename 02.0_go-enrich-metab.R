
"%&%" <- function(a,b) paste0(a,b)

library("data.table")
library("dplyr")
library("Homo.sapiens")
#source("http://bioconductor.org/biocLite.R")
#biocLite("topGO")
library("topGO")
library("GO.db")

server.dir <- "/well/got2d/jason/"
work.dir <- server.dir %&% "projects/mv-compare/"
input.dir <- work.dir %&% "input_files/"
output.dir <- work.dir %&% "output_files/"
sig.results.file <- input.dir %&% "metaxcan_full_results_0.01_only.txt"
sig.df <- fread(sig.results.file)

args = commandArgs(trailingOnly=TRUE)
metab <- args[1]
ont <- args[2]


#Assess Gene's contribution to gene set enrichment 


go_enrich_metab_nodrop <- function(metab,gene.vec,ont="BP",pval.thresh=0.01){
  # gene.vec is vector of ensid ids for significant genes or set of interest 
  ensgenes <- unique(keys(Homo.sapiens,keytype="ENSEMBL"))
  geneList <- as.factor(as.integer(ensgenes %in% gene.vec))
  names(geneList) <- ensgenes 
  # Build topGOdata object 
  GOdata <- new("topGOdata", ontology = ont, 
                allGenes = geneList,annotationFun = annFUN.org, 
                mapping="org.Hs.eg.db",ID="ensembl")
  result.fisher <- runTest(GOdata,algorithm="classic",statistic="fisher")
  pval.fisher <- sort(score(result.fisher))
  pval.adjust <- p.adjust(pval.fisher,method="BH")
  sig.fisher <- pval.fisher[pval.fisher < pval.thresh]
  out.df <- data.frame("metabolite"=rep(metab,length(sig.fisher)),
                       "ontology"=rep(ont,length(sig.fisher)),
                       "GO"=names(sig.fisher),
                       "pval"=sig.fisher,
                       "pval.adj"=pval.adjust[1:length(sig.fisher)],
                       stringsAsFactors = FALSE,
                       "dropped.gene"=rep("full",length(sig.fisher)),
                       "impact"=rep("baseline",length(sig.fisher)))
  return(out.df)
}

go_enrich_metab_drop <- function(metab,gene.vec,dropped.gene,nodrop.df,ont="BP"){
  # gene.vec is vector of ensid ids for significant genes or set of interest 
  sub.vec <- sig.genes[sig.genes!=dropped.gene]
  ensgenes <- unique(keys(Homo.sapiens,keytype="ENSEMBL"))
  geneList <- as.factor(as.integer(ensgenes %in% sub.vec))
  names(geneList) <- ensgenes 
  # Build topGOdata object 
  GOdata <- new("topGOdata", ontology = ont, 
                allGenes = geneList,annotationFun = annFUN.org, 
                mapping="org.Hs.eg.db",ID="ensembl")
  result.fisher <- runTest(GOdata,algorithm="classic",statistic="fisher")
  pval.fisher <- sort(score(result.fisher))
  pval.adjust <- p.adjust(pval.fisher,method="BH")
  keep.fisher <- pval.fisher[names(pval.fisher) %in% nodrop.df$GO]
  keep.fisher <- keep.fisher[order(match(names(keep.fisher),nodrop.df$GO))]
  keep.adjust <- pval.adjust[names(pval.adjust) %in% nodrop.df$GO]
  keep.adjust <- keep.adjust[order(match(names(keep.adjust),nodrop.df$GO))]
  impact <- keep.fisher < nodrop.df$pval
  impact <- ifelse(impact==TRUE,"more.sig","less.sig")
  out.df <- data.frame("metabolite"=rep(metab,length(keep.fisher)),
                       "ontology"=rep(ont,length(keep.fisher)),
                       "GO"=names(keep.fisher),
                       "pval"=keep.fisher,
                       "pval.adj"=keep.adjust,
                       stringsAsFactors = FALSE,
                       "dropped.gene"=rep(dropped.gene,length(keep.fisher)),
                       "impact"=impact)
  return(out.df)
}

build_go_file <- function(metab,ont){
  out.df <- c()
  sig.genes <- unique(filter(sig.df,metabolite==metab)$gene)
  nodrop.df <- go_enrich_metab_nodrop(metab,sig.genes,ont,pval.thresh=0.01)
  out.df <- rbind(out.df,nodrop.df)
  pb <- txtProgressBar(min=0,max=length(sig.genes),style=3)
  for (i in 1:length(sig.genes)){
    setTxtProgressBar(pb,i)
    dgene <- sig.genes[i]
    build.df <- go_enrich_metab_drop(metab,sig.genes,dgene,nodrop.df,ont) 
    out.df <- rbind(out.df,build.df)
  }
  savename <- output.dir %&% metab %&% "_GO_" %&% ont %&% ".txt"
  write.table(out.df,file=savename,sep="\t",row.names=F,quote=F)
  return(out.df)
}

df <- build_go_file(metab,ont)
