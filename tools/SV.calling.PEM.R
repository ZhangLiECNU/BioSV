### SV calling using discordant read pairs (PEM information).
### The two-step clustering method is used to cluster discordant 
### read pairs, and predicts the candidate SV breakpoints. 
### Usage: Rscript SV.calling.PEM.R args.list.Rdata SV-type chromosome-1 [chromosome-2]

args=commandArgs(T)
#load(args[1])
options(scipen=999)
suppressMessages(library(data.table))
suppressMessages(library(igraph))

ClusterByPos=function(vec,min.dis,multi.samples = FALSE)
{
names(vec) = 1:length(vec)
sort.vec = sort(vec)
diff.sort.vec = diff(sort.vec)
split.index = which(diff.sort.vec > min.dis)
start.index = c(1,split.index + 1)
end.index = c(split.index,length(sort.vec))
cls.list = lapply(apply(cbind(start.index,end.index),1,function(x){y = list(as.numeric(names(sort.vec)[x[1]:x[2]])) }),unlist)
cls.list
}

ClusterFirstStep=function(drps,insertion,strand=TRUE)
{
if(strand)
{
chroms = paste(drps$V1,drps$V3,drps$V7,drps$V8,sep="_")
}else
{
chroms = paste(drps$V1,drps$V3,sep="_")
}
uniq.chroms = unique(chroms)
hclust.cls = list()
for(i in 1:length(uniq.chroms))
{
time0 = proc.time()[3]
tmp.drps = drps[chroms==uniq.chroms[i],]
index1 = ClusterByPos(tmp.drps$V2,insertion)
index2 = ClusterByPos(tmp.drps$V4,insertion)
cls1 = rep(1:length(index1),sapply(index1,length))
cls2 = rep(1:length(index2),sapply(index2,length))
cls1 = cls1[match(1:length(cls1),unlist(index1))]
cls2 = cls2[match(1:length(cls2),unlist(index2))]
cls = paste(cls1,cls2,sep = "_")
filtered.cls = names(which(table(cls) > 500))
tmp.drps1 = tmp.drps[is.na(match(cls,filtered.cls)),]
tmp.drps1 = cbind(tmp.drps1,cluster = cls[is.na(match(cls,filtered.cls))])

tmp.drps2 = tmp.drps[!is.na(match(cls,filtered.cls)),]
tmp.drps2 = cbind(tmp.drps2,cluster = cls[!is.na(match(cls,filtered.cls))])
if(nrow(tmp.drps2) > 0)
{
ods = unlist(tapply(tmp.drps2$V6,tmp.drps2$cluster,function(t){order(t,decreasing = T)}))
rows = unlist(tapply(1:nrow(tmp.drps2),tmp.drps2$cluster,list))[ods<=500]
tmp.drps2 = tmp.drps2[rows,]
res = rbind(tmp.drps1,tmp.drps2)
}else
{
res = tmp.drps1
}


hclust.cls[i] = list(res)
}
names(hclust.cls) = uniq.chroms
hclust.cls
}

coord2cliques=function(coord,min.ratio = .05,max.ratio = 0.5,insertsize)
{
pairs = which(as.matrix(dist(coord,method="maximum"))<=insertsize,arr.ind = T)
pairs = matrix(pairs[pairs[,1] < pairs[,2],],ncol = 2)
if(length(pairs) > 0)
{
pairs = matrix(apply(pairs,2,as.character),ncol = 2)
g = max_cliques(graph_from_edgelist(pairs,directed = F))
if(length(g) > 1)
{
size = sapply(g,length)
g = g[order(size)]
overlap = t(sapply(g[-length(g)],function(x){c(length(intersect(x,g[[length(g)]])),length(x))}))
overlap = matrix(overlap,ncol = 2)
overlap1 = overlap[,1]/overlap[,2]
second.idx = ifelse(length(which(overlap1<=min.ratio)) > 0,max(which(overlap1<=min.ratio)),NA)
second.idx = ifelse(is.na(second.idx) & size[length(size)]==size[length(size)-1] & overlap1[length(overlap1)]<max.ratio,length(size)-1,second.idx)

max.overlap.idx = ifelse(length(which(overlap1>=max.ratio)) > 0,max(which(overlap1>=max.ratio)),NA)
max.clique = rownames(coord)[unique(as.numeric(names(unlist(g[c(max.overlap.idx,length(g))]))))]
second.max.clique = rownames(coord)[unique(as.numeric(names(unlist(g[second.idx]))))]
y = list(max.clique,second.max.clique)
}else
{
max.clique = rownames(coord)[g[[1]]]
y = list(max.clique)
}
}else
{
y = NULL
}
y
}

reads2cliques=function(mat.forcluster,insertsize,max.reads = 500)
{
time = proc.time()[3]
cls.test = names(which(table(mat.forcluster$cluster)<=max.reads))
mat.forcluster = mat.forcluster[!is.na(match(mat.forcluster$cluster,cls.test)),]
pos = paste(mat.forcluster$V2,mat.forcluster$V4,sep="\t")
names(pos) = mat.forcluster$V5
top2cliques = tapply(pos,mat.forcluster$cluster,function(x){
time = proc.time()[3]
coord = t(sapply(strsplit(x,"\t"),as.numeric))
y = coord2cliques(coord,0.05,0.5,insertsize)
})
time1 = proc.time()[3] - time
top2cliques = top2cliques[sapply(top2cliques,length) > 0]
top2cliques
}

disbyfactor=function(numbers,factor,strand)
{
if(strand=="+")
{
y = tapply(as.numeric(numbers),factor,max)
}else
{
y = tapply(as.numeric(numbers),factor,min)
}
y
}

sumbyfactor=function(numbers,factor)
{
y = tapply(as.numeric(numbers),factor,sum)
}

catbyfactor=function(strs,factor)
{
y = tapply(strs,factor,function(x){paste(x,collapse = ";")})
}

uniquebyfactor=function(numbers,factor)
{
y = tapply(numbers,factor,unique)
}

countbyfactor=function(numbers,factor)
{
y = tapply(numbers,factor,length)
}

ClusterSecondStep=function(mat,insertsize)
{
dups = mat$V5[duplicated(mat$V5)]
mat = mat[is.na(match(mat$V5,dups)),]
len1 = tapply(mat$V2,mat$cluster,function(x){diff(range(x))})
len2 = tapply(mat$V4,mat$cluster,function(x){diff(range(x))})
mat.cliques = mat[!is.na(match(mat$cluster,names(which(len1<=insertsize & len2<=insertsize)))),]
mat.forcluster = mat[is.na(match(mat$cluster,names(which(len1<=insertsize & len2<=insertsize)))),]
if(nrow(mat.forcluster)>0)
{
reads.clusters = reads2cliques(mat.forcluster,insertsize)
if(length(reads.clusters)>0)
{
reads.vec = unlist(reads.clusters)
test.for.list = sapply(reads.clusters,is.list)
reads.clusters[!test.for.list] = sapply(reads.clusters[!test.for.list],function(x){list(list(x))})
reads.count = unlist(sapply(reads.clusters,function(x){y=sapply(x,length)}))
reads.count = reads.count[reads.count > 0]
cls = rep(1:length(reads.count),reads.count)
mat.forcluster = mat.forcluster[match(reads.vec,mat.forcluster$V5),]
mat.forcluster$cluster = cls
}else
{
mat.forcluster = NULL
}

}else
{
mat.forcluster = NULL
}
pooled.cliques = rbind(mat.cliques,mat.forcluster)
if(nrow(pooled.cliques) > 0)
{
strand = c(pooled.cliques$V7[1],pooled.cliques$V8[1])
drp.svs = data.table(uniquebyfactor(pooled.cliques$V1,pooled.cliques$cluster),
                     disbyfactor(pooled.cliques$V2,pooled.cliques$cluster,strand[1]),
                     uniquebyfactor(pooled.cliques$V3,pooled.cliques$cluster),
					 disbyfactor(pooled.cliques$V4,pooled.cliques$cluster,strand[2]),
                     countbyfactor(pooled.cliques$V5,pooled.cliques$cluster),
					 catbyfactor(pooled.cliques$V5,pooled.cliques$cluster),
					 catbyfactor(pooled.cliques$V6,pooled.cliques$cluster),
                     catbyfactor(paste(pooled.cliques$V7,pooled.cliques$V8,sep=""),pooled.cliques$cluster))
}else{
drp.svs = data.table()}
colnames(drp.svs) = paste("V",1:ncol(drp.svs),sep = "")
drp.svs = drp.svs
}

merge.clusters=function(clusters)
{
chr1 = unique(clusters$V1)
pos1 = round(tapply(paste(clusters$V2,as.numeric(clusters$V5),sep = "_"),clusters$cluster,function(t){
pos.tmp = sapply(strsplit(t,"_"),function(strs){as.numeric(strs[1])})
rc.tmp = sapply(strsplit(t,"_"),function(strs){as.numeric(strs[2])})
pos.out = round(sum(pos.tmp*rc.tmp)/sum(rc.tmp))
}))
chr2 = unique(clusters$V3)
pos2 = round(tapply(paste(clusters$V4,as.numeric(clusters$V5),sep = "_"),clusters$cluster,function(t){
pos.tmp = sapply(strsplit(t,"_"),function(strs){as.numeric(strs[1])})
rc.tmp = sapply(strsplit(t,"_"),function(strs){as.numeric(strs[2])})
pos.out = round(sum(pos.tmp*rc.tmp)/sum(rc.tmp))
}))
rc = tapply(as.numeric(clusters$V5),clusters$cluster,sum)
rnm = tapply(clusters$V6,clusters$cluster,function(t){paste(t,collapse = ";")})
mapq = tapply(clusters$V7,clusters$cluster,function(x){paste(x,collapse = ";")})
strand = tapply(clusters$V8,clusters$cluster,function(t){paste(t,collapse = ";")})
y = as.data.table(cbind(chr1,pos1,chr2,pos2,rc,rnm,mapq,strand,clusters$cluster[1]))
}


BalancedSV=function(clusters,insertsize)
{
strands = sapply(strsplit(names(clusters),"_"),function(x){x})[3:4,]
paired.idx = unique(t(apply(which(apply(strands,2,function(x){apply(strands,2,function(t){sum(t==x)})})==0,arr.ind = T),1,sort)))

if(length(clusters)>length(paired.idx) & length(paired.idx)>0)
{
unable.to.merge = eval(as.call(c(rbind,clusters[-c(paired.idx)])))
colnames(unable.to.merge) = paste("V",1:ncol(unable.to.merge),sep="")
}else if (length(clusters)>length(paired.idx) & length(paired.idx)==0)
{
unable.to.merge = eval(as.call(c(rbind,clusters)))
}else
{
unable.to.merge = NULL
}

if(nrow(paired.idx)>0)
{
res = apply(paired.idx,1,function(idx){
tmp.clusters = ClusterFirstStep(eval(as.call(c(rbind,clusters[idx]))),insertsize,FALSE)[[1]]
dup.clusters = tmp.clusters$cluster[duplicated(tmp.clusters$cluster)]
tmp.clusters.unpaired = tmp.clusters[is.na(match(tmp.clusters$cluster,dup.clusters)),]
tmp.clusters.paired = tmp.clusters[!is.na(match(tmp.clusters$cluster,dup.clusters)),]
if(nrow(tmp.clusters.paired)>0)
{
paired.clusters = tapply(1:nrow(tmp.clusters.paired),tmp.clusters.paired$cluster,function(row.idx){
dist.bkps = as.matrix(dist(as.matrix(tmp.clusters.paired[row.idx,c(2,4)]),method = "maximum"))
dist.bkps = as.data.table(cbind(which(dist.bkps>=0,arr.ind=T),dist.bkps[which(dist.bkps>=0)]))
dist.bkps = dist.bkps[dist.bkps$row<dist.bkps$col,]
dist.bkps = dist.bkps[dist.bkps$V3==min(dist.bkps$V3),]
merged.paired.clusters = apply(as.matrix(dist.bkps),1,function(t){
merge.clusters(tmp.clusters.paired[row.idx[t[1:2]],])
})
unmerged.paired.clusters = tmp.clusters.paired[row.idx[-c(as.matrix(dist.bkps)[,-3])],]
uniq.clusters = eval(as.call(c(rbind,merged.paired.clusters)))
colnames(uniq.clusters) = paste("V",1:ncol(uniq.clusters),sep = "")
colnames(unmerged.paired.clusters) = paste("V",1:ncol(unmerged.paired.clusters),sep = "")
final.uniq.clusters = rbind(uniq.clusters,unmerged.paired.clusters)
final.uniq.clusters
})
paired.clusters = eval(as.call(c(rbind,paired.clusters)))
}else
{
paired.clusters = NULL
}
colnames(tmp.clusters.unpaired) = paste("V",1:ncol(tmp.clusters.unpaired),sep = "")
merged.SVs = rbind(tmp.clusters.unpaired,paired.clusters)
})

res = eval(as.call(c(rbind,res)))[,-9]
colnames(res) = paste("V",1:ncol(res),sep = "")
}else
{
res = NULL
}
res = rbind(res,unable.to.merge)
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
infile = paste(args.list$OutDir,"/DISC/",pattern,".",paste(chrs,collapse="."),".txt",sep = "")
if(file.exists(infile))
{
drp.breakpoints.tmp = unique(fread(infile,sep = '\t',showProgress = F))
libinfo = args.list$lib.info
insertsize = libinfo[[1]] * 3 + libinfo[[3]] - libinfo[[2]] * 2
drp.clusters = ClusterFirstStep(drp.breakpoints.tmp,insertsize)

tmp.SVs = lapply(drp.clusters,function(x){
y = ClusterSecondStep(x,insertsize)
})


if(length(tmp.SVs)>1)
{
###
# tmp.SVs=ClusterFirstStep(eval(as.call(c(rbind,tmp.SVs))),insertsize,FALSE)
tmp.SVs = BalancedSV(tmp.SVs,insertsize)

}else
{
tmp.SVs = tmp.SVs[[1]]
}


PEM.SV = as.data.table(cbind(tmp.SVs,pattern))
colnames(PEM.SV) = paste("V",1:ncol(PEM.SV),sep = "")

}else
{
PEM.SV = data.table(matrix(,nrow = 0,ncol = 9))
}
outfile = paste(args.list$OutDir,"/DISC/",pattern,".",paste(chrs,collapse = "."),".res.txt",sep = "")
fwrite(PEM.SV,outfile,sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)



