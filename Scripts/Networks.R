#!/usr/local/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
file<-args[1]
ppis<-read.table(file,sep="\t",header=T)
library(igraph)
g<-graph_from_data_frame(ppis[,1:2],directed=F)
g2<-simplify(g)
nl<-length(names(V(g2)))
if(nl<=100){ l=3 } else if(nl>100 & nl<=200){ l=2 } else if(nl>200 & nl<=400){ l=1.5 } else if(nl>400 & nl<=650){ l=1 } else if(nl>650 & nl<=1000){ l=0.8 } else if(nl>1000 & nl<=2000){ l=0.5 } else { l=0.2 }
if(l==3){ew=0.4} else if(l==2){ew=0.2} else if(l==1.5){ew=0.15} else if(l==1){ew=0.1} else if(l==0.8){ew=0.05} else if(l==0.5){ew=0.001} else {ew=0.0001}
V(g2)$size<-1+(degree(g2)-1)*0.3
V(g2)$color<-"lightgreen"
V(g2)$frame.color<-"lightgreen"
out_ntwk<-gsub("_PPIs_gs.txt","_PPIs_Network.pdf",file)
pdf(out_ntwk,height=5.5,width=5.5)
plot(g2,vertex.label.cex=0.15,edge.width=ew,layout=layout_nicely(g2))
dev.off()
deg<-degree(g2,mode="all")
try({deg2<-data.frame(names(deg),deg)},silent=TRUE)
out_deg<-gsub("_PPIs_gs.txt","_Network_DD.txt",file)
write.table(deg2,out_deg,sep="\t",quote=F,col.names=F,row.names=F)

