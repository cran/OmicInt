#' @title plot_3D_distribution
#'
#' @description Function allows to explore 3D distribution between the number of interactors, LFCscore and p.adj values. Function takes a data frame provided by score_genes function.
#'
#' @param data a data frame containing processed expression file from score_genes with LFCscore; class - data frame
#' @param type default value is "association", the user can select how to color data points depending on association or specificity score (e.g., selecting "specificity"); class - string
#' @return function returns an interactive plot; class - plot
#' @importFrom  plotly plot_ly
#' @importFrom  plotly layout
#' @import utils
#' @examples
#' \dontrun{
#' path_to_test_data<- system.file("extdata", "test_data.tabular", package="OmicInt")
#' # basic usage of plot_3D_distribution
#' df<-utils::read.table(path_to_test_data)
#' plot_3D_distribution(df)}
#' @export
plot_3D_distribution<-function(data, type="association"){



  #plot structures

  if(type=="association"){
  plot<-plotly::plot_ly(data=data,x=data$"LFCscore", y=data$"Interactors", z=data$"pvalue", type="scatter3d", mode="markers", size=0.7,color=data$"Association_score", text=~data$"Symbol")

  plotly::layout(plot,scene = list(xaxis = list(title = "LFCscore"), yaxis = list(title = "Interactors"), zaxis = list(title ="p.adj")))}


  if(type=="specificity"){
    plot<-plotly::plot_ly(data=data,x=data$"LFCscore", y=data$"Interactors", z=data$"pvalue", type="scatter3d", mode="markers", size=0.7,color=data$"Specificity_score", text=~data$"Symbol")

    plotly::layout(plot,scene = list(xaxis = list(title = "LFCscore"), yaxis = list(title = "Interactors"), zaxis = list(title ="p.adj")))}



}
