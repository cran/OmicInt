#' @title cluster_genes
#'
#' @description Function helps to select an optimal number of clusters and a model to be fitted during the EM phase of clustering for Gaussian Mixture Models. The function provides summaries and helps to visualise  gene clusters based on generated data using score_genes function. Weighed gene expression is clustered based on the interactome complexity, i.e., the number of known interactors according to STRING DB, with a cutoff of 700 for the score threshold. The function also provides scatter plotting and dimension reduction plots to analyse the clusters and features in the experimental data.
#'
#' @param data data frame containing processed expression file from score_genes with LFCscore; class - data frame
#' @param max_range number of clusters to consider during model selection; default 20 clusters; class - integer
#' @param clusters number of clusters to test not based on the best BIC output, user also needs to supply modelNames; class - integer
#' @param modelNames can only be supplied when clusters are also specified, this option will model based on the user parameters; class - string
#'
#' @return A data frame object that contains a summary of clusters as well as clustering and summary plots
#' @importFrom plotly ggplotly
#' @import mclust
#' @importFrom RCurl getURL
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @ImpportFrom ggplot2 ggtitle
#' @ImportFrom ggplot2 ggplot
#' @import utils
#' @import methods
#' @examples
#'  \dontrun{
#' path_to_test_data<- system.file("extdata", "test_scores.tabular", package="OmicInt")
#' # basic usage of cluster_genes
#' df<-utils::read.table(path_to_test_data)
#' df<-cluster_genes(df)
#' head(df)
#' }
#' @export
cluster_genes<-function(data, max_range=20, clusters=NULL, modelNames=NULL){

  #prepare data frame
  df<-as.data.frame(data$"Interactors")
  df<-cbind(df,data$"LFCscore")
  colnames(df)<-c("Interactors","LFCscore")
  rownames(df)<-data$"Symbol"

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
  Interactors<-model_report$"Interactors"
  Cluster<-model_report$"Cluster"
  LFCscore<-model_report$"LFCscore"
  Symbol<-model_report$"Symbol"

  plot<-ggplot2::ggplot(model_report, ggplot2::aes(x=LFCscore, y=Interactors, color=Cluster, key=Symbol))+ggplot2::geom_point(alpha=0.5,size=2)+ggplot2::ggtitle(label="Cluster distribution for a gene set")
  plot<-plotly::ggplotly(plot)
  methods::show(plot)

  #Dimension reduction based clustering visualisation

  model_dir <- mclust::MclustDR(model)
  print(summary(model_dir))
  plot(model_dir, what = "scatterplot", main="Distribution for feature dimension reduction")


  return(model_report)


}
