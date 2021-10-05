#' @title score_genes
#'
#' @description Function collects data from STRINGDB and disease association databases to  scale as well as prepare additional score integration. Function returns a data frame with calculated scores for downstream analyses.
#'
#' @param data requires a path variable to CSV file containing gene names as row names and a column with LFC values in csv format (comma separated). Columns must contain values: 'Symbol', "log2FoldChange", and 'pvalue'. Class - data frame
#' @param alpha default value returns "association" which gives a score from 0 to 1 based on how strongly the gene is associated with a disease or pathological phenotype; other options are "specificity" - to give values based on how specific the gene is for a given disease and "geometric" - to give a geometric score of both association and specificity. Class - string
#' @param beta default  FALSE; if TRUE, please supply data with column beta that contains information on gene associations from single cell studies. Class - string
#' @param gamma default FALSE; if TRUE, please supply data with column gamma that contains information on gene associations from proteome studies. Class - string
#' @return a data frame with  calculated score values for the downstream analyses; class - data frame
#'
#' @importFrom  RCurl getURL
#' @import utils
#' @examples
#' \dontrun{
#' path_to_test_data<- system.file("extdata", "data.csv", package="OmicInt")
#' #basic usage of score_genes function
#' df<-score_genes(path_to_test_data)
#' head(df)}
#' @export
score_genes<-function(data, alpha="association",beta=FALSE, gamma=FALSE){

#download the data from curated databases
string_url <- RCurl::getURL("https://gitlab.com/Algorithm379/databases/-/raw/main/Interactor_number_STRING_HS.csv")
string_db <- utils::read.csv(text = string_url)

disgenet_url <- RCurl::getURL("https://gitlab.com/Algorithm379/databases/-/raw/main/Gene_disease_scores_HS.csv")
disgenet_db <- utils::read.csv(text = disgenet_url)

data<-utils::read.csv(data, header = TRUE)

#Imput check

entered_names<-colnames(data)

if(!("Symbol"%in%entered_names)){stop("There is no 'Symbol' column in your table, please check and resubmit")}

if(!("log2FoldChange"%in%entered_names)){stop("There is no 'log2FoldChange' column in your table, please check and resubmit")}

if(!("pvalue"%in%entered_names)){stop("There is no 'pvalue' column in your table, please check and resubmit")}

if((beta==TRUE)&&!("beta"%in%entered_names)){stop("There is no 'beta' column in your table, but you selected option TRUE, please check and resubmit")}

if((gamma==TRUE)&&!("gamma"%in%entered_names)){stop("There is no 'gamma' column in your table, but you selected option TRUE, please check and resubmit")}

#Prepare columns

data$"Interactors"<-ifelse(data$"Symbol"%in%string_db$"Gene_name",string_db$"interactor_number",0)

data$"Association_score"<-ifelse(data$"Symbol"%in%disgenet_db$"Gene_name",disgenet_db$"harmonic_mean_association",0)

data$"Specificity_score"<-ifelse(data$"Symbol"%in%disgenet_db$"Gene_name",disgenet_db$"harmonic_mean_specificity",0)

#Prepare scores based on selected parameters

if(beta==FALSE && gamma==FALSE){

if((alpha=="association")){
       data$"LFCscore"<-data$"log2FoldChange"*(1+data$"Association_score")
}
if(alpha=="specificity"){
      data$"LFCscore"<-data$"log2FoldChange"*(1+data$"Specificity_score")
}
if(alpha=="geometric"){
  data$"LFCscore"<-data$"log2FoldChange"*(1+(data$"Association_score"*data$"Specificity_score")^(1/2))
}

}

if(beta==TRUE && gamma==FALSE){

  if((alpha=="association")){
    data$"LFCscore"<-data$"log2FoldChange"*(1+data$"Association_score"+data$"beta")
  }
  if(alpha=="specificity"){
    data$"LFCscore"<-data$"log2FoldChange"*(1+data$"Specificity_score"+data$"beta")
  }
  if(alpha=="geometric"){
    data$"LFCscore"<-data$"log2FoldChange"*(1+(data$"Association_score"*data$"Specificity_score")^(1/2)+data$"beta")
  }

}

if(beta==FALSE && gamma==TRUE){

  if((alpha=="association")){
    data$"LFCscore"<-data$"log2FoldChange"*(1+data$"Association_score"+data$"gamma")
  }
  if(alpha=="specificity"){
    data$"LFCscore"<-data$"log2FoldChange"*(1+data$"Specificity_score"+data$"gamma")
  }
  if(alpha=="geometric"){
    data$"LFCscore"<-data$"log2FoldChange"*(1+(data$"Association_score"*data$"Specificity_score")^(1/2)+data$"gamma")
  }

}



if(beta==TRUE && gamma==TRUE){

  if((alpha=="association")){
    data$"LFCscore"<-data$"log2FoldChange"*(1+data$"Association_score"+data$"gamma"+data$"beta")
  }
  if(alpha=="specificity"){
    data$"LFCscore"<-data$"log2FoldChange"*(1+data$"Specificity_score"+data$"gamma"+data$"beta")
  }
  if(alpha=="geometric"){
    data$"LFCscore"<-data$"log2FoldChange"*(1+(data$"Association_score"*data$"Specificity_score")^(1/2)+data$"gamma"+data$"beta")
  }

}

return (data)}
