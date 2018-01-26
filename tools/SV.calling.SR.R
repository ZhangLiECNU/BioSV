### SV calling using split reads.
### The two-step clustering method is used to cluster split reads, 
### and predicts the candidate SV breakpoints. 
### Usage: Rscript SV.calling.SR.R args.list.Rdata SV-type chromosome-1 [chromosome-2]

args=commandArgs(T)
options(scipen=999)
suppressMessages(library(data.table))
suppressMessages(library(igraph))
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
	
	maxbyfactor=function(numbers,factor)
	{
	y = tapply(as.numeric(numbers),factor,max)
	}
	
	mostbyfactor=function(numbers,factor)
	{
	y = tapply(as.numeric(numbers),factor,function(t){
	max(as.numeric(names(which(table(t)==max(table(t))))))
	})
	}
	
	leastbyfactor=function(numbers,factor)
	{
	y = tapply(as.numeric(numbers),factor,function(t){
	min(as.numeric(names(which(table(t)==max(table(t))))))
	})
	}
	
	minbyfactor=function(numbers,factor)
	{
	y = tapply(as.numeric(numbers),factor,min)
	}
	
	sumbyfactor=function(numbers,factor)
	{
	y = tapply(as.numeric(numbers),factor,sum)
	}
	
	catbyfactor=function(strs,factor)
	{
	y = tapply(strs,factor,function(x){paste(x,collapse=";")})
	}
	
	uniquebyfactor=function(numbers,factor)
	{
	y = tapply(numbers,factor,unique)
	}
	
	countbyfactor=function(numbers,factor)
	{
	y = tapply(numbers,factor,length)
	}

	merge.bkp=function(sv,label)
	{
	MAPQ = tapply(apply(sv[,c(6,8)],1,function(x){min(as.numeric(x))}),label,function(x){paste(x,collapse = ";")})
	strand = tapply(paste(sv$V7,sv$V9,sep = ""),label,function(x){paste(x,collapse = ";")})
	direction = tapply(paste(sv$V10,sv$V11,sep = ""),label,function(x){paste(x,collapse = ";")})
	read.name = tapply(sv$V5,label,function(x){paste(x,collapse = ";")})
	read.count = as.numeric(table(label))
	res = cbind(sv[match(levels(label),as.matrix(label)),c(1:4)],read.count,MAPQ,strand,direction,read.name)
	rownames(res) = NULL
	res
	}

	ClusterByPos=function(vec,min.dis,multi.samples = FALSE)
	{
	names(vec) = 1:length(vec)
	sort.vec = sort(vec)
	diff.sort.vec = diff(sort.vec)
	split.index = which(diff.sort.vec > min.dis)
	start.index = c(1,split.index + 1)
	end.index = c(split.index,length(sort.vec))
	cls.list = lapply(apply(cbind(start.index,end.index),1,function(x){y=list(as.numeric(names(sort.vec)[x[1]:x[2]])) }),unlist)
	cls.list
	}
	
	coord2cliques=function(coord,min.ratio = .05,max.ratio = 0.5,insertsize)
	{
	pairs = which(as.matrix(dist(coord,method = "maximum"))<=insertsize,arr.ind=T)
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
	second.idx = ifelse(length(which(overlap1<=min.ratio))>0,max(which(overlap1<=min.ratio)),NA)
	second.idx = ifelse(is.na(second.idx) & size[length(size)]==size[length(size) - 1] & overlap1[length(overlap1)] < max.ratio,length(size) - 1,second.idx)
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
	cls.test = names(which(table(mat.forcluster$V10)<=max.reads))
	mat.forcluster = mat.forcluster[!is.na(match(mat.forcluster$V10,cls.test)),]
	pos = paste(mat.forcluster$V2,mat.forcluster$V4,sep = "\t")
	names(pos) = mat.forcluster$V9
	top2cliques = tapply(pos,mat.forcluster$V10,function(x){
		time = proc.time()[3]
		coord = t(sapply(strsplit(x,"\t"),as.numeric))
		y = coord2cliques(coord,0.05,0.5,insertsize)
		})
	top2cliques = top2cliques[sapply(top2cliques,length) > 0]
	top2cliques
	}

	ClusterFirstStep=function(clippeds,insertion)
	{
	dims = dim(clippeds)
	chroms = paste(clippeds$V1,clippeds$V3,sep = "_")
	uniq.chroms = unique(chroms)

	hclust.cls = sapply(uniq.chroms,function(x){
	tmp.clippeds = clippeds[chroms==x,]
	index1 = ClusterByPos(tmp.clippeds$V2,insertion)
	index2 = ClusterByPos(tmp.clippeds$V4,insertion)
	cls1 = rep(1:length(index1),sapply(index1,length))
	cls2 = rep(1:length(index2),sapply(index2,length))
	cls1 = cls1[match(1:length(cls1),unlist(index1))]
	cls2 = cls2[match(1:length(cls2),unlist(index2))]
	cls = paste(cls1,cls2,sep="_")
	res = list(cbind(tmp.clippeds,cluster=cls))
	res
	})

	hclust.cls
	}
	ClusterSecondStep=function(mat,min.dist=15)
	{
	#dups=mat[duplicated(mat[,5]),5]
	#mat=mat[is.na(match(mat[,5],dups)),]
	len1=tapply(as.numeric(mat$V2),mat$V10,function(x){diff(range(x))})
	len2=tapply(as.numeric(mat$V4),mat$V10,function(x){diff(range(x))})
	mat.cliques=mat[!is.na(match(mat$V10,names(which(len1<=min.dist & len2<=min.dist)))),]
	mat.forcluster=mat[is.na(match(mat$V10,names(which(len1<=min.dist & len2<=min.dist)))),]
	if(nrow(mat.forcluster)>0)
	{
	reads.clusters=reads2cliques(mat.forcluster,min.dist)
	if(length(reads.clusters)>0)
	{
	reads.vec=unlist(reads.clusters)
	reads.count=unlist(sapply(reads.clusters,function(x){sapply(x,length)}))
	reads.count=reads.count[reads.count>0]

	cls=rep(1:length(reads.count),reads.count)
	mat.forcluster=mat.forcluster[match(reads.vec,mat.forcluster$V9),]
	mat.forcluster$V10=cls
	}else
	{
	mat.forcluster=NULL
	}

	}else
	{
	mat.forcluster=NULL
	}
	pooled.cliques=rbind(mat.cliques,mat.forcluster)
	clipped.svs=data.table(uniquebyfactor(pooled.cliques$V1,pooled.cliques$V10),mostbyfactor(pooled.cliques$V2,pooled.cliques$V10),
				  uniquebyfactor(pooled.cliques$V3,pooled.cliques$V10),leastbyfactor(pooled.cliques$V4,pooled.cliques$V10),
				  sumbyfactor(pooled.cliques$V5,pooled.cliques$V10),catbyfactor(pooled.cliques$V6,pooled.cliques$V10),
				  catbyfactor(pooled.cliques$V7,pooled.cliques$V10),catbyfactor(pooled.cliques$V8,pooled.cliques$V10),
				  catbyfactor(pooled.cliques$V9,pooled.cliques$V10))
	colnames(clipped.svs)=paste("V",1:ncol(clipped.svs),sep="")
	clipped.svs
	}


sc.mismatch = 20
infile = paste(args.list$OutDir,"/CLIP/",pattern,".",paste(chrs,collapse = "."),".txt",sep = "")
if(file.exists(infile))
{
SR.breakpoints.tmp = fread(infile,sep = '\t',showProgress = F)
label = as.factor(paste(SR.breakpoints.tmp$V1,
                        SR.breakpoints.tmp$V2,
						SR.breakpoints.tmp$V3,
						SR.breakpoints.tmp$V4,
						sep = "_"))

SV.merge = cbind(SR.breakpoints.tmp[,1:4],
               V5=1,
			   V6=apply(as.matrix(SR.breakpoints.tmp[,c(6,8)]),1,min),
			   V7=paste(SR.breakpoints.tmp$V7,SR.breakpoints.tmp$V9,sep = ""),
			   V8=paste(SR.breakpoints.tmp$V10,SR.breakpoints.tmp$V11,sep = ""),
			   V9=SR.breakpoints.tmp$V5)

clipped.cluster = ClusterFirstStep(SV.merge,sc.mismatch)
clipped.cluster = lapply(clipped.cluster,function(x){
		cls = x$cluster
		cls.filtered = names(which(table(cls)>500))
        x1 = x[is.na(match(cls,cls.filtered)),]
		x2 = x[!is.na(match(cls,cls.filtered)),]
		if(nrow(x2) > 0)
			{
				ods = unlist(tapply(x2$V6,x2$cluster,function(t){order(t,decreasing = T)}))
				rows = unlist(tapply(1:nrow(x2),x2$cluster,list))[ods<=500]
				x2 = x2[rows,]
				x = rbind(x1,x2)
			}else
			{
				x=x1
			}
			x
})

colnames(clipped.cluster[[1]]) = paste("V",1:ncol(clipped.cluster[[1]]),sep = "")

SV.cluster = lapply(clipped.cluster,function(x){
	y = ClusterSecondStep(x,sc.mismatch)
	})

clipped.SV = cbind(SV.cluster[[1]],pattern)

}else
{
clipped.SV = data.table(matrix(ncol = 10,nrow = 0))
}
outfile = paste(args.list$OutDir,"/CLIP/",pattern,".",paste(chrs,collapse = "."),".res.txt",sep = "")
fwrite(clipped.SV,outfile,sep = '\t',row.names = F,col.names = F,quote = F)


