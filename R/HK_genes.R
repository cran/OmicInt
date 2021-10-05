#' @title HK_genes

#' @description HK_genes function provides a way to visualise how housekeeping genes changed throughout the conditions under the investigation. Depending on the number of conditions separate plots will be generated. Function requires a path variable to a normalised count data file.
#'
#' @param data Requires a path variable to a data file of normalised scores in CSV format (comma separated); class - string
#' @param meta Requires a path variable to a data file of metadata  in CSV format (comma separated); class - string
#' @return  multiple plots; class - plots
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by_at
#' @importFrom dplyr ungroup
#' @importFrom stringr str_detect
#' @importFrom dplyr summarize
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_point
#' @ImpportFrom ggplot2 geom_violin
#' @ImportFrom ggplot2 geom_line
#' @ImportFrom ggplot2 ggplot
#' @Import viridis
#' @import methods
#' @import utils
#' @examples
#' \dontrun{
#' path_to_test_data<- system.file("extdata", "normalised_counts.csv", package="OmicInt")
#' path_to_meta_data<- system.file("extdata", "meta_data.csv", package="OmicInt")
#' # basic usage of HK_genes
#' HK_genes(path_to_test_data,path_to_meta_data)}
#' @export
HK_genes<-function(data, meta){

  #HK list for homo sapiens
  HK_list<-c("ACTB","GAPDH","PGK1","PPIA","RPLP0","ARBP","B2M","YWHAZ","SDHA","TFRC","GUSB","HMBS","HPRT1","TBP")


  #access data
  data<-utils::read.csv(data, header=TRUE)
  meta<-utils::read.csv(meta, header=TRUE)

  #transform and merge data
  data_transform<-reshape2::melt(data,"Symbol")
  data_transform<-merge(data_transform,meta,by.x="variable",by.y="Sample_ID")


  #assess the number of conditions supplied
  conditions<-colnames(meta)[str_detect(colnames(meta),"Condition")]

  #plot data

  for(condition in conditions){
  for(gene in HK_list){

    if(!(gene%in%data_transform$"Symbol")){next}
    df<-data_transform[which(gene==data_transform$"Symbol"),]

    df_mean <- dplyr::group_by_at(df, .vars=condition)
    df_mean <-dplyr::ungroup(dplyr::summarize(df_mean, average = mean(value)))

    #to avoid plotting conflicts
    Condition_1<-df[,condition]
    value<-df$"value"
    Condition_1_mean<-as.data.frame(df_mean)[,condition]
    average<-df_mean$"average"

    p<- ggplot2::ggplot(df, ggplot2::aes(Condition_1,value, fill=Condition_1)) +ggplot2::geom_violin(trim=FALSE, alpha=0.5)+ggplot2::geom_point(data = df_mean,  mapping = ggplot2::aes(x = Condition_1_mean, y = average),color="red") +ggplot2::geom_line(data = df_mean, mapping = aes(x = Condition_1_mean, y = average,group=1))+viridis::scale_fill_viridis(discrete = TRUE)+ggtitle(label=gene)

    methods::show(p)

  }}




}
