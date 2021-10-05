#' @title pattern_search

#' @description pattern_search function searches for gene patterns that were upregulated or downregulated throughout the conditions when comparing to the geometric mean across all conditions. The geometric mean serves as a base value to compare across multiple conditions if more complex patterns exist and also allows for a universal baseline. Function takes path variables to  data frames for normalised gene counts and meta data file (CSV format) as well as an additional variable that describes the name of a column that contains the condition under the investigation.
#'
#' @param data Requires a path variable to a data frame of normalised scores in CSV format; class - string
#' @param meta Requires a path variable to a data frame of metadata  in CSV format; class - string
#' @param Condition Requires a condition  name to select if there are multiple conditions in meta data file, default "Condition_1"; class - string
#' @return  a list variable which contains a pattern list with pattern names and associated genes; class - list
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
#' @importFrom tidyr spread
#' @importFrom gtools permutations
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 geom_point
#' @ImpportFrom ggplot2 geom_violin
#' @ImportFrom ggplot2 geom_line
#' @ImportFrom ggplot2 ggplot
#' @ImportFrom viridis scale_fill_viridis
#' @ImportFrom stats aggregate
#' @ImportFrom tidyselect all_of
#' @import methods
#' @import utils
#' @examples
#' \dontrun{
#' path_to_test_data<- system.file("extdata", "normalised_counts.csv", package="OmicInt")
#' path_to_meta_data<- system.file("extdata", "meta_data.csv", package="OmicInt")
#' # basic usage of pattern_search
#' pattern_search(path_to_test_data,path_to_meta_data)}
#' @export
pattern_search<-function(data, meta, Condition="Condition_1"){

  #geometric mean function
  geometric_mean<-function(array){
    n<-length(array)
    gm_val<-1
    for(val in array){

      gm_val<-gm_val*val
    }
    gm<-gm_val^(1/n)
    return(gm)
  }


  #access data
  data<-utils::read.csv(data, header=TRUE,row.names = 1)
  meta<-utils::read.csv(meta, header=TRUE)

  #select and rename to avoid conflicts when variables are passed dynamically
  meta<-meta[,c("Sample_ID",Condition)]
  colnames(meta)<-c("Sample_ID","Condition")

  #transform and merge data
  data_transform<-reshape2::melt(data,"Symbol")
  data_transform<-merge(data_transform,meta,by.x="variable",by.y="Sample_ID")

  df_mean <- as.data.frame(dplyr::group_by_at(data_transform,tidyselect::all_of( c("Condition", "Symbol"))))
  df_mean<-df_mean[order(df_mean$"Symbol"),c(2,3,4)]

  df_mean <-stats::aggregate(.~ Condition+Symbol, df_mean, mean, na.rm = TRUE)

  df_mean<-tidyr::spread(df_mean, 1,3)

  #prepare pattern groups
  con_subclass<-levels(as.factor(meta$"Condition"))
  pattern<-gtools::permutations(n = 2, r =length(con_subclass), v = c("up","down"), repeats.allowed = TRUE)

  #prepare pattern list
  names<-apply(pattern, 1, paste, collapse="_")
  pattern_groups<- vector(mode = "list", length = length(names))

  names(pattern_groups)<-names


  #plot data
  for(index in seq(1,nrow(df_mean))){

    #prepare a temporary data frame for the evaluation
    df_temp<-df_mean[index,]
    gene_symbol_temp<-df_temp$"Symbol"
    rownames(df_temp)<-df_temp$"Symbol"
    df_temp<-df_temp[,-1]
    mean_temp<-as.double(apply(df_temp,1,geometric_mean))

    #set changes based on the mean across all conditions
    temp_cond<-c()
    for(cond in colnames(df_temp)){

      if(df_temp[1,cond]>=mean_temp){
        temp_cond<-c(temp_cond,"up")
      }else if(df_temp[1,cond]<mean_temp){
        temp_cond<-c(temp_cond,"down")
      }

    }
    temp_cond<-paste(temp_cond, collapse="_")

    if(is.null(pattern_groups[temp_cond])){

      pattern_groups[[temp_cond]]<-c(gene_symbol_temp)
    }
    else{
      pattern_groups[[temp_cond]]<-c(pattern_groups[[temp_cond]],gene_symbol_temp)
    }




  }


  #print summary

  df<-as.data.frame(summary(pattern_groups)[,1])
  colnames(df)<-"Gene count"
  print("Condition subclasses")
  print(con_subclass)
  print(df)

  #return pattern list
  return(pattern_groups)

}
