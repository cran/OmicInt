#' @title pattern_plots

#' @description pattern_plots function uses a subsetted pattern data from the function pattern_search. The function plots distribution plots as well as a selected set of genes and how they changed patterns.
#'
#' @param data Requires a  data frame of normalised scores subsetted from pattern_search function. Class - string
#' @param meta Requires a  path variable to a data frame of metadata  in CSV format. Class - string
#' @param Condition Requires a condition  name to select if there are multiple conditions in meta data file, default "Condition_1". Conditions need to match between pattern_search and pattern_plot functions.Class - string
#' @param low the lowest value for the expression value; class - integer or float
#' @param high the highest value for the expression value; class -integer or float
#' @return  function plots multiple plots, class - plots
#' @import RColorBrewer
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom plotly plot_ly
#' @importFrom tidyr spread
#' @importFrom gtools permutations
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_point
#' @ImpportFrom ggplot2 geom_violin
#' @ImportFrom ggplot2 geom_line
#' @ImportFrom ggplot2 ggplot
#' @ImportFrom stats aggregate
#' @Import viridis
#' @import methods
#' @import utils
#' @examples
#' \dontrun{
#' path_to_test_data<- system.file("extdata","subsetted_data.csv", package="OmicInt")
#' path_to_meta_data<- system.file("extdata", "meta_data.csv", package="OmicInt")
#' # basic usage of pattern_search
#' pattern_plots(path_to_test_data,path_to_meta_data, 20, 10000)}
#' @export
pattern_plots<-function(data, meta, low=NA,high=NA, Condition="Condition_1"){


  if((is.na(low))||(is.na(high))){stop("You need to supply low an d high parameters for plotting")}

  #access data and select data

  meta<-utils::read.csv(meta, header=TRUE)
  meta<-meta[,c("Sample_ID",Condition)]
  colnames(meta)<-c("Sample_ID","Condition")

  data_temp<-data
  rownames(data_temp)<-data$"Symbol"
  data_temp<-data_temp[meta$"Sample_ID"]

  data_temp$"Mean"<-apply(data_temp,1,mean)
  data_temp<-data_temp[data_temp$"Mean">=low,]
  data_temp<-data_temp[data_temp$"Mean"<=high,]

  symbol_vals<-rownames(data_temp)

  data<-data[which(data$"Symbol"%in%symbol_vals),]



  #transform and merge data
  data_transform<-reshape2::melt(data,"Symbol")
  data_transform<-merge(data_transform,meta,by.x="variable",by.y="Sample_ID")

  temp_df_mean <-stats::aggregate(.~ Condition+Symbol, data_transform, mean, na.rm = TRUE)

  #plot data

      df<-data_transform

      df_mean <- dplyr::group_by(df, Condition)
      df_mean <-dplyr::ungroup(dplyr::summarize(df_mean, average = mean(value)))

      #to avoid plotting conflicts
      Condition<-df[,"Condition"]
      value<-df$"value"
      Condition_mean<-as.data.frame(df_mean)[,"Condition"]
      average<-df_mean$"average"

      p<- ggplot2::ggplot(df, ggplot2::aes(Condition,value, fill=Condition)) +ggplot2::geom_violin(trim=FALSE, alpha=0.5)+ggplot2::geom_point(data = df_mean,  mapping = ggplot2::aes(x = Condition_mean, y = average),color="red") +ggplot2::geom_line(data = df_mean, mapping = aes(x = Condition_mean, y = average,group=1))+viridis::scale_fill_viridis(discrete = TRUE)

      methods::show(p)

#prepare colours
      qual_col_pals <- RColorBrewer::brewer.pal.info[which(RColorBrewer::brewer.pal.info$"category"%in%c('qual')),] #max number of colours 335, setting for qual gives 74
      col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$"maxcolors", rownames(qual_col_pals)))

      colors <- col_vector[1:nlevels(as.factor(data$"Symbol"))]

    plot_ly(data = temp_df_mean, x = ~Condition, y = ~value, color = ~Symbol, colors = colors, type="bar")


}
