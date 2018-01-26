### Extract SV reads, including discordant read pairs and split reads, and proper reads from bam file.
### Usage: Rscript Extract.SV.Reads.R args.list.Rdata chromosome 

args=commandArgs(T)
### args.list chr1
options(scipen = 999)
args.list.file = args[1]
chr = args[2]
load(args.list.file)
samtools = args.list$samtools
bedtools = args.list$bedtools
bamfile = args.list$BamFile
args.list$MaxRC = 500
suppressMessages(library(data.table))

SR.reformat = function(split.reads,chr,args.list)
{
out1 = paste(args.list$OutDir,"/",chr,"/split.reads.tmp.bed",sep = "")
out2 = paste(args.list$OutDir,"/",chr,"/SR.paired.tmp.bed",sep = "")

cmd1 = paste("awk '{print $1,$2,$3,$4,$5,$6 >>",paste("\"",out1,"\"",sep = ""),"; if($7~/SA:Z/){print $0 >>",paste("\"",out2,"\"",sep = ""),"} }'",split.reads)

tmp.obj = system(cmd1)
SR.data = fread(out1,sep = " ",showProgress = F)
colnames(SR.data) = c("name","flag","chr","pos","mapq","cigar")
proper.test = bitwAnd(SR.data$flag,0x2)
SR.data$proper = proper.test
read1.test = bitwAnd(SR.data$flag,0x40)
read1.test[read1.test>0] = 1
read1.test[read1.test==0] = 2
SR.data$name = paste(SR.data$name,read1.test,sep = "/")
strand.test = bitwAnd(SR.data$flag,0x10)
strand = ifelse(strand.test > 0,"-","+")
SR.data$strand = strand
rm.S.H = gsub("[0-9]+S$|[0-9]+H$","",gsub("^[0-9]+S|^[0-9]+H","",SR.data$cigar))
map.len = sapply(strsplit(rm.S.H,"[A-Z]"),function(x){sum(as.numeric(x))})
SR.data$end = SR.data$pos + map.len - 1
SR.paired.tmp = fread(out2,sep = ' ',showProgress = F)
read1.test.tmp = bitwAnd(SR.paired.tmp$V2,0x40)
read1.test.tmp[read1.test.tmp > 0] = 1
read1.test.tmp[read1.test.tmp==0] = 2
SR.paired.tmp$V1 = paste(SR.paired.tmp$V1,read1.test.tmp,sep = "/")
SR.paired.data = SR.data[!is.na(match(SR.data$name,SR.paired.tmp$V1)),c("chr","pos","end","name","mapq","strand","cigar")]

split.reads.paired.out = paste(args.list$OutDir,"/",chr,"/SR.paired.bed",sep = "")

fwrite(SR.paired.data,split.reads.paired.out,sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)

SR.data = SR.data[SR.data$proper==0,c("chr","pos","end","name","mapq","strand","cigar")]
SR.data$name = gsub("\\/1$|\\/2$","",SR.data$name)
SR.data = SR.data[grep("S",SR.data$cigar),]
}


DRP.reformat = function(improper.reads,SR.improper,chr,args.list)
{
IPP = fread(improper.reads,sep = " ",showProgress = F)
colnames(IPP) = c("name","flag","chr","pos","mapq","cigar","mate")
IPP = IPP[IPP$cigar!="*",]
strand.test = bitwAnd(IPP$flag,0x10)
strand = ifelse(strand.test > 0,"-","+")
IPP$strand = strand
strand.test = bitwAnd(IPP$flag,0x10)
strand = ifelse(strand.test > 0,"-","+")
IPP$strand = strand
rm.S.H = gsub("[0-9]+S$|[0-9]+H$","",gsub("^[0-9]+S|^[0-9]+H","",IPP$cigar))
map.len = sapply(strsplit(rm.S.H,"[A-Z]"),function(x){sum(as.numeric(x))})
IPP$end = IPP$pos + map.len - 1

IPP = IPP[,c("chr","pos","end","name","mapq","strand","cigar")]
combn.ipp = rbind(IPP,SR.improper)
dup.name = combn.ipp$name[duplicated(combn.ipp$name)]
combn.ipp.dup = combn.ipp[!is.na(match(combn.ipp$name,dup.name)),]
combn.ipp.dup = combn.ipp.dup[order(combn.ipp.dup$name),]
mapq = apply(cbind(combn.ipp.dup$mapq[seq(1,nrow(combn.ipp.dup),by = 2)],
                   combn.ipp.dup$mapq[seq(2,nrow(combn.ipp.dup),by = 2)]),1,min)
				   
paired.drp = cbind(combn.ipp.dup[seq(1,nrow(combn.ipp.dup),by = 2),c("chr","pos","end")],
                   combn.ipp.dup[seq(2,nrow(combn.ipp.dup),by = 2),c("chr","pos","end","name")],
				   mapq,
				   combn.ipp.dup$strand[seq(1,nrow(combn.ipp.dup),by = 2)],
				   combn.ipp.dup$strand[seq(2,nrow(combn.ipp.dup),by = 2)])
				 
fwrite(paired.drp,paste(args.list$OutDir,"/",chr,"/IntraChrom.bedpe",sep = ""),sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)
improper.unpaired = combn.ipp[is.na(match(combn.ipp$name,dup.name)),c("chr","pos","end","name","mapq","strand")]
fwrite(improper.unpaired,paste(args.list$OutDir,"/",chr,"/InterChrom.bedpe",sep = ""),sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)

}
	ClusterByPos = function(vec,min.dis,multi.samples=FALSE)
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

	drp2cluster = function(drps,insertion,strand = TRUE)
	{
	dims = dim(drps)
	if(strand)
	{
	chroms = paste(drps[,1],drps[,3],drps[,7],drps[,8],sep = "_")
	}else
	{
	chroms = paste(drps[,1],drps[,3],sep = "_")
	}
	uniq.chroms = unique(chroms)

	hclust.cls = sapply(uniq.chroms,function(x){
	time0=proc.time()[3]

	tmp.drps = matrix(drps[chroms==x,],ncol = dims[2])
	pos = matrix(apply(matrix(tmp.drps[,c(2,4)],ncol = 2),2,as.numeric),ncol = 2)

	index1 = ClusterByPos(pos[,1],insertion)
	index2 = ClusterByPos(pos[,2],insertion)
	cls1 = rep(1:length(index1),sapply(index1,length))
	cls2 = rep(1:length(index2),sapply(index2,length))
	cls1 = cls1[match(1:length(cls1),unlist(index1))]
	cls2 = cls2[match(1:length(cls2),unlist(index2))]
	cls = paste(cls1,cls2,sep = "_")
	res = list(cbind(tmp.drps,cls))
	res
	})
	hclust.cls
	}


Short.del = function(SR.data,short.deletions,lib.info,args.list)
{
SD.data=fread(short.deletions,sep = ' ',showProgress = F)
colnames(SD.data) = c("name","flag","chr","pos","mapq","cigar")
strand.test = bitwAnd(SD.data$flag,0x10)
strand = ifelse(strand.test > 0,"-","+")
SD.data$strand = strand
map.len = sapply(strsplit(SD.data$cigar,"[A-Z]"),function(x){sum(as.numeric(x))})
SD.data$end = SD.data$pos + map.len - 1
SR.data$name = gsub("\\/1$|\\/2$","",SR.data$name)

dup.reads = c(SD.data$name[!is.na(match(SD.data$name,SR.data$name))],SD.data$name[duplicated(SD.data$name)])
if(length(dup.reads) > 0)
{
SD.data = rbind(SR.data[,c("chr","pos","end","name","mapq","strand","cigar")],SD.data[,c("chr","pos","end","name","mapq","strand","cigar")])

SD.data.paired = SD.data[!is.na(match(SD.data$name,dup.reads)),]
SD.data.paired = SD.data.paired[order(SD.data.paired$name),]
SD.data.paired.r1 = SD.data.paired[seq(1,nrow(SD.data.paired),by = 2),]
colnames(SD.data.paired.r1) = paste(colnames(SD.data.paired.r1),"_1",sep = "")
SD.data.paired.r2 = SD.data.paired[seq(2,nrow(SD.data.paired),by = 2),]
colnames(SD.data.paired.r2) = paste(colnames(SD.data.paired.r2),"_2",sep = "")
mapq = apply(cbind(SD.data.paired.r1[,"mapq_1"],SD.data.paired.r2[,"mapq_2"]),1,min)
combn.SD.data.paired = cbind(SD.data.paired.r1[,c("chr_1","pos_1","end_1")],
                             SD.data.paired.r2[,c("chr_2","pos_2","end_2","name_2")],
							 mapq,
							 SD.data.paired.r1$strand_1,
							 SD.data.paired.r2$strand_2)
drps = as.matrix(combn.SD.data.paired[,c( "chr_1" ,"end_1",  "chr_2",  "pos_2",  "name_2")])
insertsize = lib.info[[1]]*3
sd.cluster = drp2cluster(drps,insertsize,F)[[1]]
hc.cluster = names(which(table(sd.cluster[,ncol(sd.cluster)]) > 1  & table(sd.cluster[,ncol(sd.cluster)]) < args.list$MaxRC))
sd.cluster.filtered = sd.cluster[!is.na(match(sd.cluster[,ncol(sd.cluster)],hc.cluster)),]
combn.SD.data.paired = combn.SD.data.paired[match(sd.cluster.filtered[,5],combn.SD.data.paired$name),]
fwrite(combn.SD.data.paired,paste(args.list$OutDir,"/",chr,"/IntraChrom.bedpe",sep = ""),sep = '\t',row.names = F,col.names = F,quote = F,append = T,showProgress = F)
}

}

PP.unique = function(proper.reads,lib.info,chr)
{
read.len = lib.info[[2]]
proper.reads.out = gsub("tmp.proper.txt$","proper.txt",proper.reads)
cmd = system(paste("uniq -c ",proper.reads,"|awk '{OFS=\"\t\";print ",paste("\"",chr,"\"",sep = ""),",$2,$2+",read.len,"-1,","$1}' >",proper.reads.out,sep = ""))
file.remove(proper.reads)
}

cat(paste("##INFO:",date(),": Extract SV reads from the ",chr," of ",bamfile,"\n",sep = ""))

split.reads = paste(args.list$OutDir,"/",chr,"/split.reads.txt",sep = "")
short.deletions = paste(args.list$OutDir,"/",chr,"/short.deletions.txt",sep = "")
improper.reads = paste(args.list$OutDir,"/",chr,"/improper.txt",sep = "")
proper.reads = paste(args.list$OutDir,"/",chr,"/tmp.proper.txt",sep = "")
test.files.exist = sapply(c(split.reads,short.deletions,improper.reads),function(x){if(file.exists(x)){file.remove(x)}})
lib.info = args.list$lib.info
insertsize = lib.info[[3]] + 3*lib.info[[1]]
if(args.list$BiosvMode=="call")
{
scripts = paste(samtools,"view", bamfile ,chr ,"| awk -F \"\t\" '{if($6~/S|H/){print $1,$2,$3,$4,$5,$6,$16 >>",paste("\"",split.reads,"\"",sep = "") ,"}else{ flen=sqrt($9**2); if(and($2,0x2)){if(flen>",insertsize,"){print $1,$2,$3,$4,$5,$6 >>",paste("\"",short.deletions,"\"",sep=""), "}else{print $4 >>",paste("\"",proper.reads,"\"",sep=""),"}}else{print $1,$2,$3,$4,$5,$6,$7 >>",paste("\"",improper.reads,"\"",sep = ""),"}  }} '")
tmp.obj = system(scripts)

## return the soft-clipped reads
SR.data = SR.reformat(split.reads,chr,args.list)
## return the improper reads, the pair of which mapped to the same chromosome, but unpaired.
DRP.unpaired.data = DRP.reformat(improper.reads,SR.data,chr,args.list)
## return proper pairs, but unpaired
SD.unpaired.data = Short.del(SR.data,short.deletions,lib.info,args.list)
}else
{
scripts = paste(samtools,"view -q 20 -f 0x2", bamfile ,chr ,"| awk -F \"\t\" '{ flen=sqrt($9**2); if(flen<=",insertsize,"){print $4 >>",paste("\"",proper.reads,"\"",sep=""),"}}'")
tmp.obj = system(scripts)
}
tmp.obj = PP.unique(proper.reads,lib.info,chr)








