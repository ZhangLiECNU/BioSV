### Arrange the extracted SV reads by chromosome
### Usage: Rscript Arrange.SV.Reads.R args.list.Rdata

args=commandArgs(T)
args.list.file=args[1]
options(scipen=999)
load(args.list.file)
suppressMessages(library(data.table))

samtools = args.list$samtools
bedtools = args.list$bedtools
bam = args.list$BamFile
getTestChroms=function(samtools,bamfile)
{
chrs = system(paste(samtools,"view -H",bamfile),intern = T)
chrs = gsub("^SN:","",sapply(strsplit(chrs[grep("@SQ",chrs)],"\t"),function(x){x[2]}))
test.chroms = c(1:22,"X","Y",paste("chr",c(1:22,"X","Y"),sep = ""))
chrs = intersect(chrs,test.chroms)
}


Split.dup=function(data.table)
{
name = data.table$name
dup.name = name[duplicated(name)]
dup.data.table = data.table[!is.na(match(data.table$name,dup.name)),]
dup.data.table = dup.data.table[order(dup.data.table$name),]
uniq.data.table = data.table[is.na(match(data.table$name,dup.name)),]
res = list(uniq.data.table,dup.data.table)
}

bed2bedpe=function(bed)
{
read1 = bed[seq(1,nrow(bed),by = 2),]
read2 = bed[seq(2,nrow(bed),by = 2),]
mapq = apply(cbind(read1$mapq,read2$mapq),1,min)
res = cbind(read1[,1:3],read2[,1:4],mapq,read1[,6],read2[,6])
}

SVreadByChrom=function(svread.dir,Discordant.Reads.bedpe,clipped.SH.aln,test.chroms,args.list)
{
tmp.sv.reads.bed = paste(svread.dir,"/tmp.sv.reads.bed",sep = "")
tmp.obj = system(paste("awk '{OFS=\"\t\";print $1,$2,$3,1,$8;print $4,$5,$6,1,$8 }'",Discordant.Reads.bedpe,"| cat -",clipped.SH.aln," |awk '{if($5>=20){print $1,$2,$3}}'" ,">",tmp.sv.reads.bed))
tmp.obj = sapply(test.chroms,function(chr){
chr1 = paste("\"",chr,"\"",sep = "")
out = paste(svread.dir,"/",chr,".bed",sep = "")
tmp = system(paste("awk '{if($1==",chr1,"){print $1,$2,$3}}'",tmp.sv.reads.bed,">",out))
datainR = fread(out,sep = ' ',showProgress = F)
colnames(datainR) = c("chr","start","end")
datainR = cbind(datainR[order(datainR$start,datainR$end),],1)
fwrite(datainR,out,sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)
})
}


clip.dir = paste(args.list$OutDir,"/CLIP",sep = "")
disc.dir = paste(args.list$OutDir,"/DISC",sep = "")
svread.dir = paste(args.list$OutDir,"/SvReads",sep = "")
sv.output.dir = paste(args.list$OutDir,"/Output",sep = "")
tmp.obj = sapply(c(clip.dir,disc.dir,svread.dir,sv.output.dir),function(x){
if(!dir.exists(x))
{
dir.create(x)
}
})

test.chroms = getTestChroms(samtools,bam)
InterChrom = paste(args.list$OutDir,"/",test.chroms,"/InterChrom.bedpe",sep = "")
InterChrom.sorted = paste(args.list$OutDir,"/",test.chroms,"/InterChrom.sorted.bed",sep = "")

IntraChrom = paste(args.list$OutDir,"/",test.chroms,"/IntraChrom.bedpe",sep = "")
SR.paired = paste(args.list$OutDir,"/",test.chroms,"/SR.paired.bed",sep = "")

clipped.SH.aln = paste(clip.dir,"/clipped.SH.aln",sep = "")
if(args.list$rm.high.cov)
{
cmd.SR.merge = system(paste("cat",paste(SR.paired,collapse = " "),"|",bedtools,"intersect -v -a - -b",args.list$exclude_bed,">",clipped.SH.aln))
}else
{
cmd.SR.merge = system(paste("cat",paste(SR.paired,collapse = " "),">",clipped.SH.aln))
}
tmp.obj = sapply(SR.paired,function(f){file.remove(f)}) ### SR.paired.bed removed



InterChrom.data=sapply(1:length(InterChrom),function(i){
cmd1 = paste("sort -k4",InterChrom[i],">",InterChrom.sorted[i])
tmp.obj = system(cmd1)
file.remove(InterChrom[i]) ## InterChrom.bedpe removed
})

Discordant.Reads.raw.bedpe = paste(disc.dir,"/Discordant.Reads.raw.bedpe",sep="")

for(i in 1:(length(InterChrom.sorted)-1))
{
	for(j in (i+1):length(InterChrom.sorted))
	{
	cmd1 = paste("join --nocheck-order -j 4 -o 1.1,1.2,1.3,2.1,2.2,2.3,1.4,1.5,2.5,1.6,2.6 ",InterChrom.sorted[i],InterChrom.sorted[j],"|awk '{OFS=\"\t\";mapq=$8;if($9<$8){mapq=$9};print $1,$2,$3,$4,$5,$6,$7,mapq,$10,$11}'",">>",Discordant.Reads.raw.bedpe)
	tmp.obj = system(cmd1)
	}
}
tmp.obj=sapply(InterChrom.sorted,function(f){file.remove(f)}) ## InterChrom.sorted.bedpe removed

Discordant.Reads.bedpe = paste(disc.dir,"/Discordant.Reads.bedpe",sep = "")
if(args.list$rm.high.cov)
{
cmd.DRP.merge = system(paste("cat",paste(c(IntraChrom,Discordant.Reads.raw.bedpe),collapse = " "),"|",bedtools,"pairtobed -type neither -a - -b",args.list$exclude_bed,">",Discordant.Reads.bedpe))
}else
{
cmd.DRP.merge = system(paste("cat",paste(c(IntraChrom,Discordant.Reads.raw.bedpe),collapse = " "),">",Discordant.Reads.bedpe))
}
tmp.obj = sapply(c(IntraChrom,Discordant.Reads.raw.bedpe),function(f){file.remove(f)}) ## InterChrom.bedpe and Discordant.Reads.raw.bedpe removed
tmp.obj = SVreadByChrom(svread.dir,Discordant.Reads.bedpe,clipped.SH.aln,test.chroms,args.list)
tmp.files = unlist(sapply(c("improper.txt","short.deletions.txt","split.reads.tmp.bed","split.reads.txt","SR.paired.tmp.bed","tmp.proper.txt"),
function(f){paste(args.list$OutDir,"/",test.chroms,"/",f,sep = "")}))
tmp.obj = sapply(tmp.files,function(f){if(file.exists(f)){file.remove(f)}})





