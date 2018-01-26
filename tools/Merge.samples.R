### Merge SVs from multiple samples and generate joint SV breakpoints.
### Usage: Rscript Merge.samples.R SVfile-1,SVfile-2,...,SVfile-n sample-1,sample-2,...,sample-n outfile

args=commandArgs(T)
options(scipen=999)
suppressMessages(library(data.table))

ClusterByPos=function(vec,min.dis,multi.samples=FALSE)
{
names(vec)=1:length(vec)
sort.vec=sort(vec)
diff.sort.vec=diff(sort.vec)
split.index=which(diff.sort.vec>min.dis)
start.index=c(1,split.index+1)
end.index=c(split.index,length(sort.vec))
cls.list=lapply(apply(cbind(start.index,end.index),1,function(x){y=list(as.numeric(names(sort.vec)[x[1]:x[2]])) }),unlist)
cls.list
}

SepReadsByGW=function(drps,insertion,strand=TRUE)
{
if(strand)
{
chroms=paste(drps$V1,drps$V3,drps$V7,drps$V8,sep="_")
}else
{
chroms=paste(drps$V1,drps$V3,sep="_")
}
uniq.chroms=unique(chroms)
hclust.cls=list()
for(i in 1:length(uniq.chroms))
{
time0=proc.time()[3]
tmp.drps=drps[chroms==uniq.chroms[i],]
index1=ClusterByPos(tmp.drps$V2,insertion)
index2=ClusterByPos(tmp.drps$V4,insertion)
cls1=rep(1:length(index1),sapply(index1,length))
cls2=rep(1:length(index2),sapply(index2,length))
cls1=cls1[match(1:length(cls1),unlist(index1))]
cls2=cls2[match(1:length(cls2),unlist(index2))]
cls=paste(cls1,cls2,sep="_")
filtered.cls=names(which(table(cls)>500))
tmp.drps=tmp.drps[is.na(match(cls,filtered.cls)),]
res=cbind(tmp.drps,cluster=cls[is.na(match(cls,filtered.cls))])
hclust.cls[i]=list(res)
}
names(hclust.cls)=uniq.chroms
hclust.cls
}

undet.svformat=function(clipped.out,pem.out,rows)
{
MRC=rep(0,length(rows))
CLIP.MRC=rep(0,length(rows))
CLIP.strand=rep("0,0,0,0",length(rows))
CLIP.direction=rep("0,0,0,0",length(rows))
CLIP.MAPQ=rep(0,length(rows))
DISC.MRC=rep(0,length(rows))
DISC.strand=rep("0,0,0,0",length(rows))
DISC.MAPQ=rep(0,length(rows))
MQ20=rep(0,length(rows))

if(length(clipped.out)>0)
{
MRC[match(as.numeric(names(table(clipped.out[,8]))),rows)]=MRC[match(as.numeric(names(table(clipped.out[,8]))),rows)]+table(clipped.out[,8])
CLIP.MRC[match(as.numeric(names(table(clipped.out[,8]))),rows)]=CLIP.MRC[match(as.numeric(names(table(clipped.out[,8]))),rows)]+table(clipped.out[,8])
CLIP.strand.direction.stat=t(sapply(tapply(clipped.out[,16],clipped.out[,8],function(x){
tmp=t(sapply(strsplit(x,"_"),function(t){t}))
strand=c(sum(tmp[,1]=="+-"),sum(tmp[,1]=="-+"),sum(tmp[,1]=="++"),sum(tmp[,1]=="--"))
direction=c(sum(tmp[,2]=="RL"),sum(tmp[,2]=="LR"),sum(tmp[,2]=="RR"),sum(tmp[,2]=="LL"))
y=c(strand,direction)
}),function(m){m}))

CLIP.strand[match(as.numeric(rownames(CLIP.strand.direction.stat)),rows)]=paste(CLIP.strand.direction.stat[,1],CLIP.strand.direction.stat[,2],CLIP.strand.direction.stat[,3],CLIP.strand.direction.stat[,4],sep=",")
CLIP.direction[match(as.numeric(rownames(CLIP.strand.direction.stat)),rows)]=paste(CLIP.strand.direction.stat[,5],CLIP.strand.direction.stat[,6],CLIP.strand.direction.stat[,7],CLIP.strand.direction.stat[,8],sep=",")
CLIP.MAPQ.tmp=tapply(as.numeric(clipped.out[,15]),clipped.out[,8],function(x){
round(median(x))
})
CLIP.MAPQ[match(as.numeric(names(CLIP.MAPQ.tmp)),rows)]=CLIP.MAPQ.tmp
CLIP.MQ20=tapply(as.numeric(clipped.out[,15]),clipped.out[,8],function(x){
sum(x>=20)
})
MQ20[match(as.numeric(names(CLIP.MQ20)),rows)]=CLIP.MQ20
}
if(length(pem.out)>0)
{
MRC[match(as.numeric(names(table(pem.out[,8]))),rows)]=MRC[match(as.numeric(names(table(pem.out[,8]))),rows)]+table(pem.out[,8])
DISC.MRC[match(as.numeric(names(table(pem.out[,8]))),rows)]=DISC.MRC[match(as.numeric(names(table(pem.out[,8]))),rows)]+table(pem.out[,8])
DISC.strand.tmp=tapply(pem.out[,16],pem.out[,8],function(x){
y=paste(sum(x=="+-"),sum(x=="-+"),sum(x=="++"),sum(x=="--"),sep=",")
})
DISC.strand[match(as.numeric(names(DISC.strand.tmp)),rows)]=DISC.strand.tmp
DISC.MAPQ.tmp=tapply(as.numeric(pem.out[,15]),pem.out[,8],function(x){round(median(x))})
DISC.MAPQ[match(as.numeric(names(DISC.MAPQ.tmp)),rows)]=DISC.MAPQ.tmp

DISC.MQ20=tapply(as.numeric(pem.out[,15]),pem.out[,8],function(x){sum(x>=20)})
MQ20[match(as.numeric(names(DISC.MQ20)),rows)]=DISC.MQ20

}
sv.format=paste("MRC=",MRC,":CLIP.MRC=",CLIP.MRC,":CLIP.strand=",CLIP.strand,":CLIP.direction=",CLIP.direction,":CLIP.MAPQ=",CLIP.MAPQ,":DISC.MRC=",DISC.MRC,
":DISC.strand=",DISC.strand,":DISC.MAPQ=",DISC.MAPQ,":MQ20=",MQ20,sep="")
}



files=unlist(strsplit(args[1],"\\,"))
samples=sapply(strsplit(unlist(strsplit(args[2],"\\,")),"\\/"),function(x){x[length(x)]})
files.exists=sapply(files,function(x){length(readLines(x,n=1))>0})
outfile=args[3]


if(sum(files.exists)>0)
{
svs=lapply(files[files.exists],function(x){
y=fread(x,sep='\t',showProgress=F)
})

insertsize=sapply(strsplit(files,"\\/"),function(x){rdata=paste(c(x[-((length(x)-1):length(x))],"args.list.Rdata"),collapse="/");load(rdata);y=get(ls())$lib.info[[1]]})*3
bedtools=sapply(strsplit(files[1],"\\/"),function(x){rdata=paste(c(x[-((length(x)-1):length(x))],"args.list.Rdata"),collapse="/");load(rdata);y=get(ls())$bedtools})
clipped.reads=gsub("\\/Output\\/","/CLIP/",files)
pem.reads=gsub("\\/Output\\/","/DISC/",files)
temp.outfiles=gsub("Output/","",files)
samples.table=rep(samples[files.exists],sapply(svs,nrow))

svs.table=as.data.table(cbind(eval(as.call(c(rbind,svs))),V8=samples.table))
colnames(svs.table)=paste("V",1:ncol(svs.table),sep="")
svs.table.sort=svs.table[order(svs.table$V2,svs.table$V4),]
merged.svs.multi.samples=SepReadsByGW(svs.table,min(insertsize),F)[[1]]
MRC=as.numeric(gsub("MRC=","",merged.svs.multi.samples$V6))
merged.svs.multi.samples=merged.svs.multi.samples[order(merged.svs.multi.samples$cluster,MRC,decreasing=T),]
bkps=merged.svs.multi.samples[match(unique(merged.svs.multi.samples$cluster),merged.svs.multi.samples$cluster),1:5]
sv.format=merged.svs.multi.samples$V7
names(sv.format)=merged.svs.multi.samples$V8
sample.fields=t(sapply(tapply(sv.format,merged.svs.multi.samples$cluster,function(x){
sample=names(x)
fmt=x
y=x[match(samples,sample)]
 }),function(x){x}))

 sample.fields=sample.fields[unique(merged.svs.multi.samples$cluster),]
 colnames(sample.fields)=samples
 res=cbind(bkps,sample.fields)
 labels=paste(merged.svs.multi.samples$V8,merged.svs.multi.samples$cluster,sep="\t")
 if(length(which(table(labels)>1))>0)
 {
 repeat.bkps=merged.svs.multi.samples[!is.na(match(labels,names(which(table(labels)>1)))),]
 repeat.bkps=repeat.bkps[order(repeat.bkps$cluster,repeat.bkps$V8,as.numeric(gsub("MRC=","",repeat.bkps$V6)),decreasing=T),]
 res1=lapply(unique(repeat.bkps$cluster),function(x){
	tmp=repeat.bkps[repeat.bkps$cluster==x,]
	y=tmp[match(samples,tmp$V8)+1,]
 })
 res1=eval(as.call(c(rbind,res1)))
 res1=res1[!is.na(res1$V1),]
 sv.format.repeat=res1$V7
 names(sv.format.repeat)=res1$V8
 sample.fields.repeat=t(sapply(tapply(sv.format.repeat,1:nrow(res1),function(x){
sample=names(x)
fmt=x
y=x[match(samples,sample)]
 }),function(x){x}))
 colnames(sample.fields.repeat)=samples
 res2=cbind(res1[,1:5],sample.fields.repeat)
 }else
 {
 res2=NULL
 }
 res=rbind(res,res2)
 colnames(res)=paste("V",1:ncol(res),sep="")
 for(i in 6:ncol(res))
 {
	rows=which(is.na(res[[paste("V",i,sep="")]]))
	if(length(rows)>0){
	bkp.na=res[rows,]
	bkp.bedpe=cbind(bkp.na[,1:2],bkp.na$V2+1,bkp.na[,3:4],bkp.na$V4+1,bkp.na$V5,rows)
	fwrite(bkp.bedpe,temp.outfiles[i-5],sep='\t',row.names=F,col.names=F,quote=F,showProgress=F)
	tmp.clipped.out=gsub("txt$","bedpe",clipped.reads[i-5])
	tmp.pem.out=gsub("txt$","bedpe",pem.reads[i-5])
	cmd1=system(paste("awk '{OFS=\"\t\";mapq=$6;if($8<mapq){mapq=$8};print $1,$2,$2+1,$3,$4,$4+1,mapq,$7$9\"_\"$10$11}'",clipped.reads[i-5],">",tmp.clipped.out))
	cmd1.1.out=t(sapply(strsplit(system(paste(bedtools,"pairtopair -slop",insertsize[i-5],"-a",temp.outfiles[i-5],"-b",tmp.clipped.out),intern=T),"\t"),function(x){x}))
	cmd2=system(paste("awk '{OFS=\"\t\";print $1,$2,$2+1,$3,$4,$4+1,$6,$7$8}'",pem.reads[i-5],">",tmp.pem.out))
	cmd2.1.out=t(sapply(strsplit(system(paste(bedtools,"pairtopair -slop",insertsize[i-5],"-a",temp.outfiles[i-5],"-b",tmp.pem.out),intern=T),"\t"),function(x){x}))
	undet.sv.format=undet.svformat(cmd1.1.out,cmd2.1.out,rows)
	res[[paste("V",i,sep="")]][rows]=undet.sv.format
	tmp.obj=file.remove(c(tmp.clipped.out,tmp.pem.out,temp.outfiles[i-5]))
	}
}
}else
{
res=data.table(matrix(nrow=0,ncol=length(samples)+6))
}
 fwrite(res,outfile,sep='\t',row.names=F,col.names=F,quote=F,showProgress=F)
 
 
 
 

