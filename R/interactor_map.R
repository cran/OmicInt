#' @title interactor_map
#'
#' @description interactor_map uses information mined from STRING database to map  experimental, predicted,  or referenced interactions  to see if there are  any interactors in the set of significantly changed genes and how they are linked. The function requires a data frame prepared by score_genes. The output is a plot depicting interaction map.
#'
#' @param data requires a data frame containing gene names as row names and a column with LFC values; class - data frame
#' @return interaction map/plot; class - plot
#'
#' @importFrom igraph graph_from_adjacency_matrix degree delete.vertices
#' @Import STRINGdb
#' @import RColorBrewer
#' @importFrom  dendextend find_k
#' @importFrom  dendextend color_labels
#' @importFrom  dendextend color_branches
#' @importFrom  dendextend colored_bars
#' @import cluster
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme
#' @ImpportFrom ggplot2 element_text
#' @ImportFrom ggplot2 geom_col
#' @ImportFrom ggplot2 ggplot
#' @import RColorBrewer
#' @importFrom graphics legend
#' @import utils
#' @importFrom stats  hclust
#' @importFrom stats as.dendrogram
#' @importFrom  RCurl getURL
#' @examples
#'  \dontrun{
#' path_to_test_data<- system.file("extdata", "test_data.tabular", package="OmicInt")
#' # basic usage of interactor_map
#' df<-utils::read.table(path_to_test_data)
#' interactor_map(df)
#' }
#' @export
interactor_map<-function(data){


  #helper functions
  min_max<-function(x){
    (x-min(x))/(max(x)-min(x))
  }


  #download the data from curated databases
  string_db<-STRINGdb::STRINGdb$new(version="11.0",species=9606,score_threshold=700)

  classes_url <- RCurl::getURL("https://gitlab.com/Algorithm379/databases/-/raw/main/HS_protein_classes_curated.csv")
  classes <- utils::read.csv(text = classes_url)

  #prepare data frame
  data$"Class"<-ifelse(data$"Symbol"%in%classes$"Gene",classes$"Class","NA")


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


  #prepare network list
  undirected_network <- igraph::graph_from_adjacency_matrix(matrix, mode="undirected", weighted=NULL,diag = FALSE)

  #remove unconnected nodes
  non_connected<- which(igraph::degree(undirected_network)==0)
  undirected_network= igraph::delete.vertices(undirected_network, non_connected)




  qual_col_pals <- RColorBrewer::brewer.pal.info[which(RColorBrewer::brewer.pal.info$"category"%in%c('qual')),] #max number of colours 335, setting for qual gives 74
  col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$"maxcolors", rownames(qual_col_pals)))

  colors <- col_vector[1:nlevels(as.factor(data$"Class"))]

  plot(undirected_network,edge.arrow.size=0.2, edge.curved=0.3,vertex.size=6,edge.color="red", vertex.label.cex=0.6, vertex.color=colors)
  graphics::legend("topleft",legend=levels(as.factor(data$"Class")), col = colors,  pch=15,cex=0.7)



}



