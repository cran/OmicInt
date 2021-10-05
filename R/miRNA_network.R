#' @title miRNA_network

#' @description miRNA_network function allows to assess how many genes are regulated by the same miRNA. Note if you supply too many genes the function will take longer to run.
#'
#' @param genes Requires a  gene list (HGNC gene symbol); class list of strings
#' @return  a heatmap plot for found interactions and a list of miRNA and regulated genes. The list output value can be used for downstream analyses. Classes returned - a plot and a list
#' @importFrom RCurl getURL
#' @importFrom pheatmap pheatmap
#' @import methods
#' @import utils
#' @examples
#' \dontrun{
#' # basic usage of miRNA_network
#' return_df<-miRNA_network(c("PIP4K2A","MOB1A","PHACTR2","MDM2","YWHAG" ,"RAB31"  ))
#' head(return_df)
#' }
#' @export
miRNA_network<-function(genes){

  #access data


  #class data preparation

  classes_url <- RCurl::getURL("https://gitlab.com/Algorithm379/databases/-/raw/main/HS_protein_classes_curated.csv")
  classes <- utils::read.csv(text = classes_url)
  #prepare data frame
  gene_annotations<-as.data.frame(genes)
  gene_annotations$"Class"<-ifelse(gene_annotations$"genes"%in%classes$"Gene",classes$"Class","NA")

  #download the data from curated databases
  location_url <- RCurl::getURL("https://gitlab.com/Algorithm379/databases/-/raw/main/Subcellular.locationmerged_protein_data.csv")
  location_df <- utils::read.csv(text = location_url)
  gene_annotations$"Location"<-ifelse(gene_annotations$"genes"%in%location_df$"Symbol",location_df$"Subcellular.location","NA")
  gene_annotations$"Location"<-ifelse( is.na(gene_annotations$"Location"),"NA", gene_annotations$"Location")

  rownames(gene_annotations)<-gene_annotations$"genes"
  gene_annotations<-gene_annotations[,-1]

  #miRNA data preparation
  miRNA_url <- RCurl::getURL("https://gitlab.com/Algorithm379/databases/-/raw/main/miRNA_df_validated.csv")
  miRNA <- utils::read.csv(text = miRNA_url)

  #prepare data frame
  miRNA_list<-list()

  #pre-filter miRNA data
  index_list<-c()
  for(gene in genes){

    index<-which(gene==miRNA$"Target_Symbol")
    index_list<-c(index,index_list)
  }
  miRNA_filtered<-miRNA[index_list,]
  #only miRNA data containing genes are reported
  for(miRNA_var in miRNA_filtered$"mature_miRNA"){

    #what genes are regulated by miRNA
    gene_list<-unique(miRNA_filtered$"Target_Symbol"[which(miRNA_var==miRNA_filtered$"mature_miRNA")])

    #what genes are contained in the data from the genes regulated by miRNA
    gene_list<-gene_list[which(gene_list%in%genes)]

    miRNA_list[[miRNA_var]]<-gene_list

  }

  #build a matrix
  matrix<-matrix(0, nrow=length(genes),ncol=length(miRNA_list))
  rownames(matrix)<-genes
  colnames(matrix)<-names(miRNA_list)

  #fill the matrix

  for(gene_var in rownames(matrix)){


    for(miRNA_var in names(miRNA_list)){

      if(gene_var%in%miRNA_list[[miRNA_var]]){
      matrix[gene_var, miRNA_var]<-1}
    }

  }



  pheatmap::pheatmap(matrix, annotation_row=gene_annotations, main = "Interactor heatmap for all genes and regulatory miRNA", legend = FALSE)

  return(miRNA_list)
}
