### Merge SV breakpoints from split reads and discordant read pairs.
### Usage: Rscript Merge.SR.PEM.signals.R args.list.Rdata SV-type chromosome-1 [chromosome-2]

args=commandArgs(T)
options(scipen=999)

suppressMessages(library(data.table))
suppressMessages(library(igraph))

merge.SR.PEM.signals=function(args.list,chrs,pattern)
{
libinfo = args.list$lib.info
insertsize = libinfo[[1]] * 3 + libinfo[[3]] - libinfo[[2]] * 2
clipped.file = paste(args.list$OutDir,"/CLIP/",pattern,".",paste(chrs,collapse = "."),".res.txt",sep = "")
pem.file = paste(args.list$OutDir,"/DISC/",pattern,".",paste(chrs,collapse = "."),".res.txt",sep = "")
if(length(readLines(clipped.file,n = 1))>0 & length(readLines(pem.file,n = 1))>0)
{
clipped.data = fread(clipped.file,sep = '\t',showProgress = F)
c.label = 1:nrow(clipped.data)
pem.data = fread(pem.file,sep = '\t',showProgress = F)
p.label = 1:nrow(pem.data)
c.file = paste(args.list$OutDir,"/CLIP/tmp.",pattern,".",paste(chrs,collapse = "."),".bedpe",sep = '')
p.file = paste(args.list$OutDir,"/DISC/tmp.",pattern,".",paste(chrs,collapse = "."),".bedpe",sep = '')
o.file = paste(args.list$OutDir,"/output.",pattern,".",paste(chrs,collapse = "."),".bedpe",sep = '')
fwrite(data.table(clipped.data$V1,clipped.data$V2,clipped.data$V2+1,clipped.data$V3,clipped.data$V4,clipped.data$V4+1,clipped.data$V10,c.label),c.file,sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)
fwrite(data.table(pem.data$V1,pem.data$V2,pem.data$V2 + 1,pem.data$V3,pem.data$V4,pem.data$V4 + 1,pem.data$V9,p.label),p.file,sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)
overlap.cmd=system(paste(args.list$bedtools,"pairtopair -slop",insertsize,"-a",p.file,"-b",c.file,">",o.file))
if(length(readLines(o.file))>0)
{
overlap=unique(fread(o.file,sep = '\t',showProgress = F))
}else
{
overlap=data.table(matrix(,nrow = 0,ncol = 18))
}
tmp.obj=file.remove(c(c.file,p.file,o.file))

if(nrow(overlap)==0)
{
Clipped.only.SV = clipped.data
PEM.only.SV = pem.data
overlap.SV = NULL
}

}else if(length(readLines(clipped.file,n=1))>0 & length(readLines(pem.file,n=1))==0)
{
Clipped.only.SV = fread(clipped.file,sep = '\t',showProgress = F)
PEM.only.SV = data.table(matrix(nrow = 0,ncol = 9))
overlap.SV = data.table(matrix(,nrow = 0,ncol = 18))
overlap = data.table(matrix(nrow = 0,ncol = 18))
}else if(length(readLines(clipped.file,n=1))==0 & length(readLines(pem.file,n=1))>0)
{
Clipped.only.SV = data.table(matrix(nrow = 0,ncol = 10))
PEM.only.SV = fread(pem.file,sep = '\t',showProgress = F)
overlap.SV = data.table(matrix(nrow = 0,ncol = 15))
overlap = data.table(matrix(nrow = 0,ncol = 18))
}else
{
Clipped.only.SV = data.table(matrix(nrow = 0,ncol = 10))
PEM.only.SV = data.table(matrix(nrow = 0,ncol = 9))
overlap.SV = data.table(matrix(nrow = 0,ncol = 15))
overlap = data.table(matrix(nrow = 0,ncol = 18))
}


if(nrow(overlap)>0)
{
overlap1 = overlap[!is.na(match(overlap$V16,names(which(table(overlap$V16)==1)))),]
overlap2 = overlap[is.na(match(overlap$V16,names(which(table(overlap$V16)==1)))),]
edges = unique(cbind(paste("A",overlap2$V8,sep = "_"),paste("B",overlap2$V16,sep = "_")))
g = graph_from_edgelist(edges,directed = F)
clu = components(g)
g.clu = groups(clu)
edges = gsub("A_|B_","",edges)
A.list = sapply(g.clu,function(x){gsub("A_","",x[grep("A",x)])})
A.cls = rep(1:length(A.list),sapply(A.list,length))
A.vec = as.numeric(unlist(A.list))
edge.cls = A.cls[match(as.numeric(edges[,1]),A.vec)]
edge.count = apply(edges,1,function(x){
x = as.numeric(x)
y = clipped.data$V5[x[2]] + pem.data$V5[x[1]]
})
# edge.dist=abs(as.numeric(overlap2[,3])-as.numeric(overlap2[,11]))+abs(as.numeric(overlap2[,6])-as.numeric(overlap2[,14]))

max.row=function(mat)
{
mat = apply(mat,2,as.numeric)
mat = data.table(mat)
colnames(mat) = c("V1","V2","V3")

mat = mat[order(mat$V2,mat$V3,decreasing = T),]
mat = mat[match(unique(mat$V2),mat$V2),]
}



if(nrow(overlap2)>0)
{
combn.pairs = max.row(cbind(edges,edge.count))
combn.label = paste(combn.pairs$V1,combn.pairs$V2,sep = "_")
overlap2.label = gsub(" ","",paste(overlap2$V8,overlap2$V16,sep = "_"))
overlap2 = overlap2[match(combn.label,overlap2.label),]
}

overlap.unique = rbind(overlap1,overlap2)
Clipped.only.SV = clipped.data[-unique(overlap$V16),]
PEM.only.SV = pem.data[-unique(overlap$V8),]
pos = as.matrix(overlap.unique[,c(9,10,12,13)])
strand.clipped = clipped.data$V7[overlap.unique$V16]
strand.drp = pem.data$V8[overlap.unique$V8]
direction = clipped.data$V8[overlap.unique$V16]
read.names.clipped = gsub("\\/1|\\/2","",clipped.data$V9[overlap.unique$V16])
patterns = overlap.unique$V7
read.names.drp = pem.data$V6[overlap.unique$V8]
read.names = data.table(t(sapply(strsplit(paste(read.names.clipped,read.names.drp,sep=";"),";"),function(x){
x = unique(x)
y = c(length(x),paste(x,collapse = ";"))
})))
read.names$V1 = as.numeric(read.names$V1)


colnames(read.names) = c("MRC","Reads")
mapq.clipped = sapply(strsplit(as.character(clipped.data$V6[overlap.unique$V16]),";"),function(x){sum(as.numeric(x))})
mq20.clipped = sapply(strsplit(as.character(clipped.data$V6[overlap.unique$V16]),";"),function(x){sum(as.numeric(x)>=20)})

mapq.drp = sapply(strsplit(as.character(pem.data$V7[overlap.unique$V8]),";"),function(x){sum(as.numeric(x))})
mq20.drp = sapply(strsplit(as.character(pem.data$V7[overlap.unique$V8]),";"),function(x){sum(as.numeric(x)>=20)})

mrc.clipped = clipped.data$V5[overlap.unique$V16]
mrc.drp = pem.data$V5[overlap.unique$V8]
pos[mq20.clipped<=0 & mq20.clipped<mq20.drp,] = as.matrix(overlap.unique[mq20.clipped<=0 & mq20.clipped<mq20.drp,c(1,2,4,5)])

overlap.SV = data.table(pos,
                        patterns,
						strand.clipped,
						strand.drp,
						direction,
						read.names,
						mrc.clipped,
						mrc.drp,
						mapq.clipped,
						mapq.drp,
						mq20=mq20.clipped+mq20.drp)

colnames(overlap.SV) = paste("V",1:ncol(overlap.SV),sep = "")
overlap.SV$V2 = as.numeric(overlap.SV$V2)
overlap.SV$V4 = as.numeric(overlap.SV$V4)

}

if(nrow(Clipped.only.SV)>0)
{
clipped.only.mapq.list = lapply(strsplit(as.character(Clipped.only.SV$V6),";"),function(x){as.numeric(x)})
Clipped.only.out = cbind(Clipped.only.SV[,c(1:4,10,7)],
                         NA,
						 Clipped.only.SV[,c(8,5,9,5)],
						 0,
						 sapply(clipped.only.mapq.list,sum),
						 0,
						 sapply(clipped.only.mapq.list,function(x){sum(x>=20)}))

colnames(Clipped.only.out) = paste("V",1:ncol(Clipped.only.out),sep="")

}else
{
Clipped.only.out = data.table(matrix(nrow = 0,ncol = 15))
}


if(nrow(PEM.only.SV)>0)
{
pem.only.mapq.list = sapply(strsplit(as.character(PEM.only.SV$V7),";"),function(x){as.numeric(x)})
PEM.only.out = cbind(PEM.only.SV[,c(1:4,9)],
                     NA,
					 PEM.only.SV[,8],
					 NA,
					 PEM.only.SV[,5:6],
					 0,
					 PEM.only.SV[,5],
					 0,
					 sapply(pem.only.mapq.list,sum),
					 sapply(pem.only.mapq.list,function(x){sum(x>=20)}))

colnames(PEM.only.out)=paste("V",1:ncol(PEM.only.out),sep="")

}else
{
PEM.only.out=data.table(matrix(nrow=0,ncol=15))
}

combn.SVs=rbind(overlap.SV,
                Clipped.only.out,
				PEM.only.out)
}



c2format=function(svs)
{
if(nrow(svs)>0)
{
c.strand = t(sapply(strsplit(as.character(svs$V6),";"),function(x){
y=c(sum(x=="+-"),
    sum(x=="-+"),
	sum(x=="++"),
	sum(x=="--"))
}))

c.strand[is.na(c.strand)] = 0
c.strand=paste("CLIP.strand=",paste(c.strand[,1],c.strand[,2],c.strand[,3],c.strand[,4],sep = ","),sep = "")

p.strand = t(sapply(strsplit(as.character(svs$V7),";"),function(x){
y=c(sum(x=="+-"),
    sum(x=="-+"),
	sum(x=="++"),
	sum(x=="--"))
}))

p.strand[is.na(p.strand)] = 0
p.strand=paste("DISC.strand=",paste(p.strand[,1],p.strand[,2],p.strand[,3],p.strand[,4],sep = ","),sep = "")

c.direction = t(sapply(strsplit(as.character(svs$V8),";"),function(x){
y=c(sum(x=="RL"),
    sum(x=="LR"),
	sum(x=="RR"),
	sum(x=="LL"))
}))

c.direction[is.na(c.direction)] = 0
c.direction = paste("CLIP.direction=",paste(c.direction[,1],c.direction[,2],c.direction[,3],c.direction[,4],sep = ","),sep = "")


c.mapq = svs$V13/svs$V11
c.mapq[is.na(c.mapq)] = 0
c.mapq = paste("CLIP.MAPQ=",round(c.mapq),sep = "")
p.mapq = svs$V14/svs$V12
p.mapq[is.na(p.mapq)] = 0
p.mapq = paste("DISC.MAPQ=",round(p.mapq),sep = "")

MRC = paste("MRC=",gsub(" ","",svs$V9),sep = "")
c.MRC = paste("CLIP.MRC=",gsub(" ","",svs$V11),sep = "")
p.MRC = paste("DISC.MRC=",gsub(" ","",svs$V12),sep = "")
mq20 = paste("MQ20=",svs$V15,sep = "")
fields = paste(MRC,c.MRC,c.strand,c.direction,c.mapq,p.MRC,p.strand,p.mapq,mq20,sep = ":")
y = cbind(svs[,1:5],MRC,fields)

}else
{
y = data.table(matrix(,nrow = 0,ncol = 10))
}
y
}

args.list.file = args[1]
load(args.list.file)
pattern = args[2]
if(pattern=="TRA")
{
chrs = args[3:4]
}else
{
chrs = args[3]
}

combn.SVs = c2format(merge.SR.PEM.signals(args.list,chrs,pattern))
out.filename = paste(args.list$OutDir,"/Output/",pattern,".",paste(chrs,collapse = "."),".txt",sep = "")
fwrite(combn.SVs,out.filename,sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)
