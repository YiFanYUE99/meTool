#' Title a function to make correlation network between genes and metabolites
#'
#' @param gene a csv file containing gene abundance without group// example in data/gene.csv
#' @param metabolites a csv file containing peak concentration and group//example in data/metabolitesg.csv
#' @param genenode the color of genes nodes
#' @param metabnode the color of metabolites nodes
#' @param poscol if correlation no.>0, the color of the edge
#' @param negcol if correlation no.<0, the color of the edge
#' @import tidygraph
#' @import igraph
#' @import ggraph
#' @import psych
#' @return the correlationnetwork
#' @export
#'
#' @examples
#'      #setwd("D:/R_work/meTool")
#'      #a<-cornet("data/gene.csv","data/metabolites.csv",geneedge="green",metabedge="pink",genenode="blue",metabnode="yellow",poscol="pink",negcol="green")
#'

cornet<-function(gene, metabolites,genenode,metabnode,poscol,negcol){
  genes<-read.table(gene, header = TRUE,row.names = 1,check.names = FALSE,sep = ",")
  metabolites<-read.table(metabolites,header = TRUE,row.names = 1,sep = ",",check.names = FALSE)
  all<-cbind(genes,metabolites)
  all<-na.omit(all)
  all1<-all[,1:(dim(all)[2]-1)]
  a<-cor(all1,method = "pearson")#列之间计算相关性
  #remove unrelated data
  for (i in 1:dim(a)[1]){
    for (j in 1:dim(a)[2]){
      if (a[i,j]<0.6 & a[i,j]>-0.6){
        a[i,j]<-NA
      }
    }
  }
  #make nodes table
  nodes<-matrix(nrow = (dim(all)[2]-1),ncol = 4)
  nodes[,1]<-rownames(a)
  nodes[,2]<-c(rep("genes",dim(genes)[2]),rep("metabolites",dim(metabolites)[2]-1))
  for (i in 1:dim(all1)[2]){
    nodes[i,3]<-mean(all1[,i])
  }
  colnames(nodes)<-c("node","type","abundance","color")
  nodes<-as.data.frame(nodes)
  nodes[,3]<-as.numeric(nodes[,3])#更改第三列数据结构从chr变numeric
  nodes[,4]<-c(rep("geneedge",dim(genes)[2]),rep("metabedge",dim(metabolites)[2]-1))

  #make edges table
  edges<-matrix(nrow = dim(a)[1]*dim(a)[2],ncol = 4)
  k=1
  for (i in 1:dim(a)[1])
  {
    for (j in 1:i)
    {
      if (!is.na(a[i,j]) & a[i,j]!=1)
      {
        if (a[i,j]>0){
          edges[k,1]<-rownames(a)[i]
          edges[k,2]<-colnames(a)[j]
          edges[k,3]<-"positive"
          edges[k,4]<-a[i,j]
          k<-k+1
        }else{
          edges[k,1]<-rownames(a)[i]
          edges[k,2]<-colnames(a)[j]
          edges[k,3]<-"negative"
          edges[k,4]<-a[i,j]
          k<-k+1
        }
      }
    }
  }

  edges<-as.data.frame(edges)
  colnames(edges)<-c("from", "to", "relation","coeffi")
  edges[,4]<-as.numeric(edges[,4])
  edges<-na.omit(edges)#delete redundant lines

  # make correlation network
  G = graph_from_data_frame(d=edges, v = nodes) %>% as_tbl_graph()


  #visualize met-genes
  pic<-G %>%
    ggraph(layout = 'circle')+#see ??layout; stress circle are good ones
    #parameter of edges
    geom_edge_link(aes(edge_width = abs(edges$coeffi),
                       color=relation,
                       edge_alpha= ifelse(abs(edges$coeffi)>0.9,0.9,0.4)),
                   linetype = "solid")+#the edge color depends on edge relation
    scale_edge_width(breaks = seq(0.2,1,0.2),
                     label=seq(0.2,1,0.2),
                     range=c(0.2,2),
                     guide = guide_legend(title = "|r|",
                                          title.theme = element_text(
                                            size=10,
                                            face="italic",
                                            family = "mono",
                                            colour = "gray1"),
                                          label.theme = element_text(
                                            size = 10,
                                            color = "grey1",
                                            family = "mono",
                                            face = "italic")))+
    scale_edge_alpha(guide= "none")+
    scale_edge_color_manual(values = c('positive' = poscol, 'negative' = negcol),
                            guide = guide_legend(title = "Correlation",
                                                 title.theme = element_text(
                                                   size = 10,
                                                   face = "italic",
                                                   family = "mono",
                                                   color = "grey1"),
                                                 label.theme = element_text(
                                                   size = 10,
                                                   color = "grey1",
                                                   family = "mono",
                                                   face ="italic")
                            ))+
    #设置点的参数
    geom_node_point(aes(color = as.factor(type) , size = abundance),
                    alpha=1)+#aes(color=,size=);points color by ; points size on ;need to be continuous data
    scale_color_manual(values = c('genes'=genenode, 'metabolites'=metabnode),
                       guide = guide_legend(title = "type",
                                            title.theme = element_text(
                                              size = 10,
                                              face="italic",
                                              family = "mono",
                                              colour = "grey1"),
                                            label.theme = element_text(
                                              size = 10,
                                              colour = "grey1",
                                              family = "mono",
                                              face="italic")))+
    scale_size(
      range = c(2,10),
      guide=guide_legend(title = "Relative Abundance",
                         title.theme = element_text(
                           size = 10,
                           face="italic",
                           family = "mono",
                           color = "gray1"),
                         label.theme = element_text(
                           size = 10,family = "mono",colour = "grey1",face = "italic")))+
    geom_node_text(aes(x=x,y=y,
                       #label= ifelse(nodes$type=="genes",as.character(nodes$node),"")
                       label=nodes$node,
                       family = "serif"
    ),
    size=5,
    color= "grey2",
    show.legend= F)+
    theme(
      text = element_text(family = "mono", size = 12),
      legend.title.align = 0,#legend左对齐
      legend.title = element_text(family = "mono", size = 8),
      legend.text.align =0,#legend text左对齐
      legend.text = element_text(family = "mono", size = 8),
      legend.spacing.x = unit(0.1,"cm"),#图例各个元素距离
      legend.spacing.y = unit(0.1,"cm"),
      panel.background = element_blank(),#delete background.
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )


  return(pic)
}
