#' @title cluster_links
#'
#' @description Function to select an optimal number of clusters and a model to be fitted during the EM phase of clustering for Gaussian Mixture Models. The function provides summaries and helps to visualise  gene clusters based on generated data using score_genes function. Weighed gene expression is clustered based on a specific disease score which can be either the association or specificity for a disease, i.e., if the gene has known links to disease phenotypes or how specific it is when describing a pathology. The function also provides scatter plots and dimension reduction plots to analyse the clusters and features in the experimental data.
#'
#' @param data data frame containing processed expression file from score_genes with LFCscore; Class - data frame
#' @param max_range number of clusters to consider during model selection; default 20 clusters. Class - integer
#' @param type type of score to consider which can be either "association" or "specificity"; default "association". Class - string
#' @param clusters number of clusters to test not based on the best BIC output, user also needs to supply modelNames; class - integer
#' @param modelNames can only be supplied when clusters are also specified, this option will model based on the user parameters; class- string
#'
#' @return A data frame object that contains a summary of clusters; class - data frame
#' @importFrom plotly ggplotly
#' @import mclust
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @ImpportFrom ggplot2 ggplot
#' @ImportFrom ggplot2 ggtitle
#' @import utils
#' @import methods
#' @examples
#'  \dontrun{
#' path_to_test_data<- system.file("extdata", "test_data.tabular", package="OmicInt")
#' # basic usage of cluster_links
#' df<-utils::read.table(path_to_test_data)
#' df<-cluster_links(df)
#' head(df)
#' }
#' @export
cluster_links<-function(data, max_range=20, type="association", clusters=NULL, modelNames=NULL){

  #prepare data frame based on a score selected
  if(type=="association"){
  df<-as.data.frame(data$"Association_score")
  df<-cbind(df,data$"log2FoldChange")
  colnames(df)<-c("Association_score","log2FoldChange")
  rownames(df)<-data$"Symbol"}

  if(type=="specificity"){
    df<-as.data.frame(data$"Specificity_score")
    df<-cbind(df,data$"log2FoldChange")
    colnames(df)<-c("Specificity_score","log2FoldChange")
    rownames(df)<-data$"Symbol"}

  #calculate Bayesian information criterion and plot different GMM models
  BIC <- mclust::mclustBIC(df, G=1:max_range)
  plot(BIC,with="BIC")
  #report the best cluster value
  print(summary(BIC))

  #To select other BIC values based on the report
  if(is.null(clusters)){
    model <- mclust::Mclust(df, x = BIC)}
  if(!is.null(clusters)&&!is.null(modelNames)){
    model <- mclust::Mclust(df,G=clusters, modelNames=modelNames)}

  #prepare a model for reporting and plotting by extracting relevant information
  model_report<-as.data.frame(model$"data")
  model_report$"Cluster"<-as.factor(model$"classification")
  model_report$"Symbol"<-names(model$"classification")
  #Plot GMM

  #to avoid namescape conflicts
  if(type=="association"){
  Association_score<-model_report$"Association_score"
  Cluster<-model_report$"Cluster"
  LFC<-model_report$"log2FoldChange"
  Symbol<-model_report$"Symbol"

  plot<-ggplot2::ggplot(model_report, ggplot2::aes(x=LFC, y=Association_score, color=Cluster, key=Symbol))+ggplot2::geom_point(alpha=0.5,size=2)+ggplot2::ggtitle(label="Cluster distribution for a gene set")
  plot<-plotly::ggplotly(plot)
  methods::show(plot)
  }
  if(type=="specificty"){
    Specificity_score<-model_report$"Specificity_score"
    Cluster<-model_report$"Cluster"
    LFC<-model_report$"log2FoldChange"
    Symbol<-model_report$"Symbol"

    plot<-ggplot2::ggplot(model_report, ggplot2::aes(x=LFC, y=Specificity_score, color=Cluster, key=Symbol))+ggplot2::geom_point(alpha=0.5,size=2)+ggplot2::ggtitle(label="Cluster distribution for a gene set")
    plot<-plotly::ggplotly(plot)
    methods::show(plot)
  }

  #Dimension reduction based clustering visualisation

  model_dir <- mclust::MclustDR(model)
  print(summary(model_dir))
  plot(model_dir, what = "scatterplot", main="Distribution for feature dimension reduction")


  return(model_report)


}
