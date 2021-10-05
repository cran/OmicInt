#' @title cluster_heatmap
#'
#' @description cluster_heatmap uses information mined from STRING database to map experimental, referenced, and inferred interactions to see if there are any interactors in the set of significantly changed genes. This heatmap provides clustered visualisation of all genes and the genes that have shared interactions.
#'
#' @param data requires a data frame containing gene names as row names and a column with LFC values. Class - data frame
#' @return heatmap; class - plot
#'
#' @Import STRINGdb
#' @importFrom stats heatmap
#' @import utils
#' @examples
#'  \dontrun{
#' path_to_test_data<- system.file("extdata", "test_data.tabular", package="OmicInt")
#' # basic usage of cluster_heatmap
#' df<-utils::read.table(path_to_test_data)
#' cluster_heatmap(df)
#' }
#' @export
cluster_heatmap<-function(data){

  #download the data from curated databases
  string_db<-STRINGdb::STRINGdb$new(version="11.0",species=9606,score_threshold=700)

  #add protein interactor  number

  #this function maps data between string database and the gene list

  db_interactors<-string_db$map(data,"Symbol",removeUnmappedRows = FALSE)

  db_interactors$"Interactor_list"<-lapply(db_interactors$"STRING_id", function(x){string_db$get_neighbors(x)})

  #create a matrix

  matrix<-matrix(0, length(db_interactors$"Symbol"),length(db_interactors$"Symbol"))
  colnames(matrix)<-db_interactors$"Symbol"
  rownames(matrix)<-db_interactors$"Symbol"

  #search for matches
  interactor_list<-list()

  for(index in seq(1:nrow(db_interactors))){

    id_list<-unlist(db_interactors[index,"Interactor_list"])
    symbol<-db_interactors[index,"Symbol"]
    #check if ID is shared with other genes in the list of differentially expressed genes and STRING Database interactor list
    if(any(db_interactors$"STRING_id"%in%id_list)){

      matches<-db_interactors$"STRING_id"[db_interactors$"STRING_id"%in%id_list]
      for(match in matches){
        match_symbol<-db_interactors[which(db_interactors$"STRING_id"==match),"Symbol"]
        interactor_list[[symbol]]<-match_symbol

        matrix[symbol, match_symbol]<-1 #match is counted once

      }

    }  }

  # plot all genes
  stats::heatmap(matrix, scale = "none", col = c("lightyellow","red"), main = "Interactor heatmap for all genes")

  #plot only interacting

  interactor_matrix<-matrix[names(interactor_list),unique(unlist(interactor_list))]

  stats::heatmap(interactor_matrix, scale = "none", col = c("lightyellow","red"), main = "Interactor heatmap")

}


