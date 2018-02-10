args=commandArgs(T)

options(scipen = 999)
{ ### load R packages
suppressMessages(library(data.table))
suppressMessages(library(igraph))
suppressMessages(library(doParallel))
suppressMessages(library(tools))
}

{ ### parse arguments
commandlines=function(args)
{
biosv.mode = args[1]
bamfile = unlist(strsplit(args[grep("^-b$|^--bam$",args)+1],"\\,"))
bamfile = sapply(bamfile,file_path_as_absolute)
out.dir = args[grep("^-o$|^--output$",args)+1]
if(length(out.dir)==0)
{
out.dir = "BioSV-output"
}
if(!dir.exists(out.dir))
{
dir.create(out.dir)
}else
{
file.remove(out.dir)
dir.create(out.dir)
}
out.dir = file_path_as_absolute(out.dir)
threads = args[grep("^-t$|^--thread$",args)+1]
if(length(threads)==0)
{
threads = 1
}else
{
threads = as.numeric(threads)
}


min.MRC = args[grep("^^--minMRC$",args)+1]
if(length(min.MRC)==0)
{
min.MRC = 4
}else
{
min.MRC = as.numeric(min.MRC)
}


remove.high.cov = grep("^-e$",args)
if(length(remove.high.cov)==0)
{
remove.high.cov = FALSE
}else
{
remove.high.cov = TRUE
}

prepared.SVs = grep("^--bedpe$",args)
if(length(prepared.SVs)>0)
{
prepared.SVs = args[grep("^--bedpe$",args)+1]
}else
{
prepared.SVs = NULL
}

subdir = paste(out.dir,"/",1:length(bamfile),".",gsub("\\.bam$","",sapply(strsplit(bamfile,"\\/"),function(x){x[length(x)]})),sep="")
mkdirs = sapply(subdir,function(x){
if(!dir.exists(x))
{
dir.create(x)
}
})

out= list( BiosvMode = biosv.mode,
           BamFile = bamfile,
           OutDir=out.dir,
           Threads = threads,
		   MinMRC=min.MRC,
		   rm.high.cov=remove.high.cov,
		   prepared.SVs=prepared.SVs)
out
}

print.help=function()
{
cat("            \n")
cat("  BioSV: Breakpoint-based identification of Structure Variation\n")
cat("            \n")
cat("  Usage: BioSV.main.R [call/genotype] [options]\n")
cat("  Options:\n")
cat("  -b <file>, --bam=<file>, <required>\n")
cat("        Input bam files (sorted bam files by bwa mem, seperated by comma)\n")
cat("  --bedpe <file>, <optional>\n")
cat("        SVs prepared to be genotyped when the mode is set to genotype\n")
cat("  -o <directory>, --output=<directory>\n")
cat("        Directory for output files [BioSV-output]\n")
cat("  -t <int>, --thread=<int>\n")
cat("        Number of threads\n")
cat("  -e <flag>\n")
cat("        Exclude high coverage regions\n")
cat("  --minMRC <int>\n")
cat("        Minimal mutant read counts [4]\n")
cat("  -h, --help\n")
cat("        Print this help message and exit\n")
}

option.cmd = args[grep("^-",args)]
unidentified.opt = setdiff(option.cmd,c("-b","--bam","-o","-e","--output","-t","--threads","--minMRC","-h","--help","--bedpe"))
if(length(unidentified.opt)>0)
{
cat(paste("##INFO:","Unidentified options:",paste(unidentified.opt,collapse=", "),"\n"))
quit(save="no")
}

if(length(grep("^-h$|^--help$",args))>0)
{
print.help()
quit(save = "no")
}


if(length(grep("^call$|^genotype$",args[1]))==0)
{
print.help()
cat("  Error: the mode should be specified as call or genotype\n")
quit(save = "no")
}


if(length(grep("^-b$|^--bam$",args))==0)
{
print.help()
cat("  Error: bam file should been specified by -b/--bam\n")
quit(save = "no")
}else
{
tmp.bamfile = unlist(strsplit(args[grep("^-b$|^--bam$",args)+1],"\\,"))
test.for.bam = file.exists(tmp.bamfile)
if(sum(!test.for.bam)>0)
{
cat(paste("  Error: no such bam files:",paste(tmp.bamfile[!test.for.bam],collapse=", "),"\n"))
print.help()
quit(save = "no")
}
}

}

LibInfo=function(args.list,samtools,bedtools)
{
bam = args.list$BamFile
chrs = system(paste(samtools,"view -H",bam),intern=T)
chrs = matrix(t(gsub("^SN:|^LN:","",sapply(strsplit(chrs[grep("@SQ",chrs)],"\t"),function(x){x[2:3]}))),ncol = 2)
chrs = data.frame(chr = chrs[,1], len = as.numeric(chrs[,2]))
test.chroms = intersect(c(1:22,"X","Y",paste("chr",c(1:22,"X","Y"),sep = "")),chrs$chr)
chrs = chrs[match(test.chroms,chrs$chr),]
regions = apply(chrs,1,function(x){
pos = paste(x[1],":",round(runif(1,1,as.numeric(x[2]))),sep = "")
})
lines = round(1000000/nrow(chrs)) + 1
reads.info = list()
for(i in 1:length(regions))
{
reads.info[i] = list(t(sapply(strsplit(system(paste(samtools," view -f 2 -F 0x800 -f 0x40 ",bam,regions[i] ,"| awk '{OFS=\"\t\";print $6,$9}' | head -n", lines),intern=T),"\t"),function(x){x})))
}

reads.info = eval(as.call(c(rbind,reads.info)))
reads.info = reads.info[setdiff(1:nrow(reads.info),grep("I|D|S",reads.info[,1])),]
read.len = as.numeric(gsub("[A-Z]","",reads.info)[1])

insertsize = abs(as.numeric(reads.info[,2]))
y=list(sigma = floor(mad(insertsize)*1.4826),
       read.len = read.len,
	   med = median(insertsize))

}

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

getTestChroms=function(samtools,bamfile)
{
chrs = system(paste(samtools,"view -H",bamfile),intern = T)
chrs = gsub("^SN:","",sapply(strsplit(chrs[grep("@SQ",chrs)],"\t"),function(x){x[2]}))
test.chroms = c(1:22,"X","Y",paste("chr",c(1:22,"X","Y"),sep = ""))
chrs = intersect(chrs,test.chroms)
}

BayesGenotype=function(count.mat,args.list)
{
prior = c(0.1,0.4)
pattern = count.mat[,1]
count.mat = apply(count.mat[,-1],2,as.numeric)
mrc = count.mat[,1]
count.mat = cbind(count.mat[,1],apply(count.mat[,2:3],1,max))
homo.pt = count.mat[,1] > count.mat[,2]
count.mat[homo.pt,] = count.mat[homo.pt,2:1]
log.lik.homo = dbinom(count.mat[,1],apply(count.mat,1,sum),prior[1],log = T)
log.lik.het = dbinom(count.mat[,1],apply(count.mat,1,sum),prior[2],log = T)
log.lik.homo[pattern=="INV" |pattern=="TRA"]=dbinom(count.mat[pattern=="INV" |pattern=="TRA",1],apply(count.mat,1,sum)[pattern=="INV" |pattern=="TRA"],prior[1]/2,log = T)
log.lik.het[pattern=="INV" |pattern=="TRA"]=dbinom(count.mat[pattern=="INV" |pattern=="TRA",1],apply(count.mat,1,sum)[pattern=="INV" |pattern=="TRA"],prior[2]/2,log = T)
genotype.index = apply(cbind(log.lik.homo,log.lik.het),1,which.max)
genotype = rep("0/1",nrow(count.mat))
genotype[genotype.index==1 & homo.pt]="1/1"
genotype[genotype.index==1 & !homo.pt & mrc<args.list$MinMRC]="0/0"
genotype
}

sortBed=function(bedfile)
{
beddata = fread(bedfile,sep = '\t',showProgress = F)
ods = order(beddata$V2,beddata$V3)
sorted.beddata = beddata[ods,]
fwrite(sorted.beddata,bedfile,sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)
}
###
PPcount=function(args.list,merged.svs.file)
{
samtools = args.list$samtools
bedtools = args.list$bedtools
bam = args.list$BamFile
test.chroms = getTestChroms(samtools,bam)
outfile = paste(args.list$OutDir,"/PP.count.bed",sep = "")
if(file.exists(outfile))
{
tmp.obj = file.remove(outfile)
}
PPs.files = paste(args.list$OutDir,"/",test.chroms,"/proper.txt",sep = "")
PP.count = sapply(test.chroms,function(chr){
chr1 = paste("\"",chr,"\"",sep = "")
tmp.svs = paste(args.list$OutDir,"/",chr,"/tmp.svs.bed",sep = "")
cmd1 = system(paste("awk '{if($1==",chr1,"){print $0}}'", merged.svs.file,">",tmp.svs))
tmp.obj = sortBed(tmp.svs)
tmp.pp = paste(args.list$OutDir,"/",chr,"/proper.txt",sep = "")
tmp.out = paste(args.list$OutDir,"/",chr,"/tmp.count.bed",sep = "")
cmd2 = system(paste(bedtools,"intersect -sorted  -wo -f 1 -a",tmp.svs,"-b",tmp.pp,"|",bedtools,"groupby -i - -g 4 -c 8 -o sum  >",tmp.out))
count = fread(tmp.out,sep = '\t',header = F,showProgress = F)
tmp.obj = file.remove(c(tmp.svs,tmp.out))
count = list(count)
})
PP.count = as.matrix(eval(as.call(c(rbind,PP.count))))
}


CalDepth=function(BioSV.output,PPcount.mat,mrc.mat,VAF = 0.05,args.list.allbam,pattern = "DEL")
{
PPcount.mat = apply(PPcount.mat,2,function(x){apply(matrix(x,ncol = 2,byrow = F),1,min)})
vaf.total = mrc.mat / (mrc.mat + PPcount.mat)
test.for.caldepth = which(BioSV.output$SVTYPE==pattern & apply(vaf.total>=VAF,1,sum) > 0 & (apply(mrc.mat >= args.list.allbam$MinMRC,1,sum)==0 | BioSV.output$FILTER=="LowQuality" ))
if(length(test.for.caldepth) > 0)
{
lowq.deletions = BioSV.output[test.for.caldepth,]
interval.left = cbind(lowq.deletions[[1]],lowq.deletions[[2]] - 50,lowq.deletions[[2]] - 1,paste("L",1:nrow(lowq.deletions),sep = "_"))
interval.right = cbind(lowq.deletions[[4]],lowq.deletions[[6]] + 1,lowq.deletions[[6]] + 50,paste("R",1:nrow(lowq.deletions),sep = "_"))
interval.center = cbind(lowq.deletions[[1]],round(lowq.deletions[[3]] / 2 + lowq.deletions[[6]] / 2 - 25),round(lowq.deletions[[3]] / 2 + lowq.deletions[[6]] / 2 + 25),paste("C",1:nrow(lowq.deletions),sep = "_"))
interval.total = rbind(interval.left,interval.right,interval.center)
interval.total = interval.total[order(interval.total[,1],as.numeric(interval.total[,2])),]
interval.total[as.numeric(interval.total[,2]) <= 0,2] = 1
test.chroms = sort(unique(interval.total[,1]))
interval.total = data.table(V1 = interval.total[,1],
                            V2 = as.numeric(interval.total[,2]),
							V3 = as.numeric(interval.total[,3]),
							V4 = interval.total[,4])

n.samples=length(args.list.allbam$BamFile)

log2ratio.mat = matrix(0, nrow = length(test.for.caldepth), ncol = n.samples)
for(i in 1:n.samples)
{
args.list = args.list.allbam
bam = args.list$BamFile[i]
bam.name = unlist(strsplit(bam,"\\/"))
bam.name = bam.name[length(bam.name)]
outdir = paste(args.list$OutDir,"/",i,".",gsub("\\.bam$","",bam.name),sep = "")
args.list$BamFile = bam
args.list$OutDir = outdir

read.depth = lapply(test.chroms,function(chr){
sv.file = paste(args.list$OutDir,"/",chr,"/",pattern,".",chr,".breakpoints.bed",sep = "")
tmp = interval.total[interval.total$V1==chr,]
write.table(tmp,sv.file,sep = '\t',row.names = F,col.names = F,quote = F)
pp.file = paste(args.list$OutDir,"/",chr,"/proper.txt",sep = "")
drp.file = paste(args.list$OutDir,"/SvReads/",chr,".bed",sep = "")
depth.file = paste(args.list$OutDir,"/",chr,"/",pattern,".",chr,".depth.bed",sep = "")
cmd.read.depth1 = system(paste(bedtools,"intersect  -sorted -wo -a",sv.file,"-b",pp.file,"|",bedtools,"groupby -g 1,2,3,4 -c 8,9 -o collapse,collapse >",depth.file))
cmd.read.depth2 = system(paste(bedtools,"intersect  -sorted -wo -a",sv.file,"-b",drp.file,"|",bedtools,"groupby -g 1,2,3,4 -c 8,9 -o collapse,collapse >>",depth.file))
read.depth = fread(depth.file,showProgress = F)
base.count = apply(read.depth[,5:6],1,function(x){
tmp = matrix(sapply(strsplit(as.character(x),","),as.numeric),ncol = 2)
bc = sum(tmp[,1]*tmp[,2])
})
read.depth$base.count = base.count
read.depth = read.depth[,-c(5:6)]
tmp.obj = file.remove(c(sv.file,depth.file))
read.depth = read.depth
})
read.depth = eval(as.call(c(rbind,read.depth)))
read.depth$avg.base = as.numeric(read.depth$base.count)/(read.depth$V3-read.depth$V2 + 1)
bkp.label = read.depth$V4
combn.rd = tapply(as.numeric(read.depth$avg.base),bkp.label,sum)
combn.label = t(sapply(strsplit(names(combn.rd),"_"),function(x){x}))

avg.depth = matrix(0,nrow=length(test.for.caldepth),ncol = 3)
avg.depth[as.numeric(combn.label[combn.label[,1]=="L",2]),1] = combn.rd[combn.label[,1]=="L"]
avg.depth[as.numeric(combn.label[combn.label[,1]=="C",2]),2] = combn.rd[combn.label[,1]=="C"]
avg.depth[as.numeric(combn.label[combn.label[,1]=="R",2]),3] = combn.rd[combn.label[,1]=="R"]

avg.depth = avg.depth + .01

log2ratio = apply(cbind(log2(avg.depth[,2] / avg.depth[,1]),log2(avg.depth[,2] / avg.depth[,3])),1,max)

log2ratio.mat[,i] = log2ratio

}

if(pattern=="DEL")
{
lowq2pass = test.for.caldepth[apply(log2ratio.mat < -1,1,sum) > 0]
}else if(pattern=="DUP")
{
lowq2pass = test.for.caldepth[apply(log2ratio.mat > 0.6 , 1 , sum) > 0]
}

BioSV.output$FILTER[lowq2pass] = "PASS"
for(i in 1:n.samples)
{
BioSV.output[[i + 8]][test.for.caldepth]=paste(BioSV.output[[i+8]][test.for.caldepth],":log2ratio=",log2ratio.mat[,i],sep = "")
}
}
BioSV.output

}




{ #### loading tools used for SV calling
exclude_bed = paste(tools_dir,"/lumpy.exclude.bed",sep = "")
samtools = paste(tools_dir,"/samtools",sep = "")
bedtools = paste(tools_dir,"/bedtools",sep = "")
extract.SVreads.rscript = paste(tools_dir,"/Extract.SV.Reads.R",sep = "")
arrange.SVreads.rscript = paste(tools_dir,"/Arrange.SV.Reads.R",sep = "")
Sr2Bkp.rscript = paste(tools_dir,"/Split.Reads.to.Bkp.R",sep = "")
Drp2Bkp.rscript = paste(tools_dir,"/Drp.Reads.to.Bkp.R",sep = "")
PEM.calling.rscript = paste(tools_dir,"/SV.calling.PEM.R",sep = "") 
SR.calling.rscript = paste(tools_dir,"/SV.calling.SR.R",sep = "")
Merge.SR.PEM.signal.rscript = paste(tools_dir,"/Merge.SR.PEM.signals.R",sep = "")
Merge.samples.rscript = paste(tools_dir,"/Merge.samples.R",sep = "")
}






{ ### SV calling
cat(paste("##INFO:",date(),": BioSV is running\n",sep = ""))
args.list.allbam = commandlines(args)
args.list.allbam$rscript = rscript
n.samples = length(args.list.allbam$BamFile)
cat(paste("##INFO:",date(),": A total of ",n.samples," bam ", ifelse(n.samples==1,"file","files"), " detected\n",sep = ""))
args.list.allbam$samtools = samtools
args.list.allbam$bedtools = bedtools

all.svlist = list()
chroms.svtype.list = list()
for(i in 1:n.samples)
{
args.list = args.list.allbam
bam = args.list$BamFile[i]
cat(paste(paste("##INFO:",date(),":",sep=""),"SV calling for",bam,"\n"))
bam.name = unlist(strsplit(bam,"\\/"))
bam.name = bam.name[length(bam.name)]
outdir = paste(args.list$OutDir,"/",i,".",gsub("\\.bam$","",bam.name),sep = "")
args.list$BamFile = bam
args.list$OutDir = outdir
args.list$exclude_bed = exclude_bed
test.chroms = getTestChroms(samtools,bam)
args.list$test.chroms = test.chroms
lib.info = LibInfo(args.list,samtools,bedtools)
args.list$lib.info = lib.info
args.list.outfile = paste(outdir,"/args.list.Rdata",sep = "")
save(args.list,file = paste(outdir,"/args.list.Rdata",sep = ""))

tmp.obj = sapply(test.chroms,function(chr){
if(!dir.exists(paste(args.list$OutDir,"/",chr,sep = ""))){
dir.create(paste(args.list$OutDir,"/",chr,sep = ""))
}
})

cl = makeCluster(min(args.list$Threads,length(test.chroms)))
registerDoParallel(cl)  
clusterExport(cl = cl, varlist = c("extract.SVreads.rscript","rscript","args.list.outfile"))
tmp.obj = parLapply(cl,sample(test.chroms,length(test.chroms),replace = F),function(chr){
cmd = paste(rscript,extract.SVreads.rscript,args.list.outfile,chr)
tmp.obj = system(cmd)
})
stopCluster(cl = cl)
if(args.list.allbam$BiosvMode=="call")
{
cat(paste("##INFO:",date(),": Arrange SV reads by chromosome\n",sep = ""))
cmd = system(paste(rscript,arrange.SVreads.rscript,args.list.outfile))

cat(paste("##INFO:",date(),": Determine SV breakpoints of SV reads\n",sep = ""))

script.vec = c(Sr2Bkp.rscript,Drp2Bkp.rscript)
cl = makeCluster(2)
registerDoParallel(cl)  
clusterExport(cl = cl, varlist = c("script.vec","rscript","args.list.outfile"))
tmp.obj = parLapply(cl,1:2,function(i){
cmd = paste(rscript,script.vec[i],args.list.outfile)
tmp.obj = system(cmd)
})
stopCluster(cl = cl)

cat(paste("##INFO:",date(),": SV calling using split reads and discordant read pairs, respectively\n",sep = ""))
command.sv.calling = c(paste(rscript,SR.calling.rscript,args.list.outfile,"DEL",test.chroms),
                       paste(rscript,PEM.calling.rscript,args.list.outfile,"DEL",test.chroms),
                       paste(rscript,SR.calling.rscript,args.list.outfile,"DUP",test.chroms),
                       paste(rscript,PEM.calling.rscript,args.list.outfile,"DUP",test.chroms),
                       paste(rscript,SR.calling.rscript,args.list.outfile,"INV",test.chroms),
                       paste(rscript,PEM.calling.rscript,args.list.outfile,"INV",test.chroms))

chroms.combn = t(combn(sort(test.chroms),2))
command.sv.calling = c(command.sv.calling,
                       paste(rscript,SR.calling.rscript,args.list.outfile,"TRA",chroms.combn[,1],chroms.combn[,2]),
                       paste(rscript,PEM.calling.rscript,args.list.outfile,"TRA",chroms.combn[,1],chroms.combn[,2]))

cl = makeCluster(min(args.list$Threads,length(command.sv.calling)))
registerDoParallel(cl)  
clusterExport(cl = cl, varlist = c("command.sv.calling"))
tmp.obj = parLapply(cl,1:length(command.sv.calling),function(i){
tmp.obj=system(command.sv.calling[i])
})
stopCluster(cl = cl)


cat(paste("##INFO:",date(),": Merge SV breakpoints from split reads and discordant read pairs\n",sep=""))
command.signal.merge = c(paste(rscript,Merge.SR.PEM.signal.rscript,args.list.outfile,"DEL",test.chroms),
                         paste(rscript,Merge.SR.PEM.signal.rscript,args.list.outfile,"DUP",test.chroms),
                         paste(rscript,Merge.SR.PEM.signal.rscript,args.list.outfile,"INV",test.chroms),
                         paste(rscript,Merge.SR.PEM.signal.rscript,args.list.outfile,"TRA",chroms.combn[,1],chroms.combn[,2]))

cl = makeCluster(args.list$Threads)
registerDoParallel(cl)  
clusterExport(cl = cl, varlist = c("command.signal.merge"))
tmp.obj = parLapply(cl,1:length(command.signal.merge),function(i){
tmp.obj = system(command.signal.merge[i])
})
stopCluster(cl = cl)

sv.files = c(paste(args.list$OutDir,"/Output/DEL.",test.chroms,".txt",sep = ""),
             paste(args.list$OutDir,"/Output/DUP.",test.chroms,".txt",sep = ""),
             paste(args.list$OutDir,"/Output/INV.",test.chroms,".txt",sep = ""),
			 paste(args.list$OutDir,"/Output/TRA.",chroms.combn[,1],".",chroms.combn[,2],".txt",sep = ""))
chroms.svtype = c(paste("DEL",test.chroms,sep = "_"),
                  paste("DUP",test.chroms,sep = "_"),
				  paste("INV",test.chroms,sep = "_"),
                  paste("TRA",chroms.combn[,1],chroms.combn[,2],sep = "_"))
file.e.test = file.exists(sv.files)
chroms.svtype.list[i] = list(data.table(V1 = chroms.svtype[file.e.test],V2 = sv.files[file.e.test]))
}else if(args.list.allbam$BiosvMode=="genotype")
{
cat(paste("##INFO:",date(),": Proper reads of ",bam," are prepared","\n",sep = ""))
if(i==n.samples)
{
prepared.SVs = fread(args.list.allbam$prepared.SVs,sep = '\t',header = T)
}
}

}
}

{###merge multiple samples
if(args.list.allbam$BiosvMode=="call")
{
if(n.samples > 1)
{
pooled.chroms.svtype = unique(c(unlist(sapply(chroms.svtype.list,function(x){x$V1}))))

files.mat=matrix(,nrow = length(pooled.chroms.svtype),ncol = n.samples)

for(i in 1:n.samples)
{
files.mat[match(chroms.svtype.list[[i]]$V1,pooled.chroms.svtype),i] = chroms.svtype.list[[i]]$V2
}

if(!dir.exists(paste(args.list.allbam$OutDir,"/multi-samples",sep = ""))){
dir.create(paste(args.list.allbam$OutDir,"/multi-samples",sep = ""))}

outfiles = paste(args.list.allbam$OutDir,"/multi-samples/",pooled.chroms.svtype,".txt",sep = "")

all.samples = paste(paste(1:length(args.list.allbam$BamFile),".",sapply(strsplit(args.list.allbam$BamFile,"\\/"),function(t){t[length(t)]}),sep = ""),collapse = ",")

sample.merge.commandlines = apply(files.mat,1,function(x){
files = paste(x,collapse = ",")
cmd = paste(rscript,Merge.samples.rscript,files,all.samples)
})

sample.merge.commandlines = paste(sample.merge.commandlines,outfiles)

cl = makeCluster(args.list.allbam$Threads)
registerDoParallel(cl)  

clusterExport(cl = cl, varlist = c("sample.merge.commandlines"))

tmp.obj = parLapply(cl,1:length(sample.merge.commandlines),function(i){
tmp.obj = system(sample.merge.commandlines[i])
})
stopCluster(cl = cl)

merged.multi.samples = lapply(outfiles,function(tmp.file){
if(length(readLines(tmp.file,n = 1))>0)
{
out.data=fread(tmp.file,sep='\t',showProgress = F)
}else
{
out.data = NULL
}
out.data
})

merged.multi.samples = eval(as.call(c(rbind,merged.multi.samples)))
}else
{

merged.multi.samples = lapply(chroms.svtype.list[[1]]$V2,function(tmp.file){
if(length(readLines(tmp.file,n = 1))>0)
{
out.data = fread(tmp.file,sep = '\t',showProgress = F)
}else
{
out.data = data.table(matrix(nrow = 0,ncol = 6 + n.samples))
}
out.data
})
merged.multi.samples = eval(as.call(c(rbind,merged.multi.samples)))
merged.multi.samples = merged.multi.samples[,-6]
colnames(merged.multi.samples) = paste("V",1:ncol(merged.multi.samples),sep = "")

}
}else if(args.list.allbam$BiosvMode=="genotype")
{
merged.multi.samples = prepared.SVs[,c(1,2,4,5,7)]
colnames(merged.multi.samples) = paste("V",1:ncol(merged.multi.samples),sep = "")
for(i in 1:n.samples)
{
merged.multi.samples = cbind(merged.multi.samples,paste("MRC=",prepared.SVs[[i + 7]],sep = ""))
}
colnames(merged.multi.samples) = paste("V",1:ncol(merged.multi.samples),sep = "")

}

}



{
if(args.list.allbam$BiosvMode=="call")
{
test.for.del.width = ((merged.multi.samples$V4 - merged.multi.samples$V2)>=50 & merged.multi.samples$V5=="DEL") | merged.multi.samples$V5!="DEL"
merged.multi.samples = merged.multi.samples[test.for.del.width,]
sample.fields = matrix(as.matrix(merged.multi.samples[,6:ncol(merged.multi.samples)]),ncol = n.samples)

mq20.mat = matrix(
apply(sample.fields,2,function(x){
y = as.numeric(sapply(strsplit(x,":"),function(t){y=gsub("MQ20=","",t[grep("MQ20=",t)]);ifelse(length(y)==1,y,0)}))
}),
ncol = n.samples)

mrc.mat=matrix(
apply(sample.fields,2,function(x){
y = as.numeric(sapply(strsplit(x,":"),function(t){y=gsub("^MRC=","",t[grep("^MRC=",t)]);ifelse(length(y)==1,y,0)}))
}),
ncol=n.samples)

merged.multi.samples = merged.multi.samples[apply(mq20.mat > 0,1,sum) > 0,]
mrc.mat = matrix(mrc.mat[apply(mq20.mat > 0,1,sum) > 0,],ncol = n.samples)
mq20.mat = matrix(mq20.mat[apply(mq20.mat > 0,1,sum) > 0,],ncol = n.samples)
if(args.list.allbam$rm.high.cov)
{
cat(paste("##INFO:",date(),": Remove regions with high coverage...\n",sep = ""))
merged.svs = cbind(merged.multi.samples[,1],merged.multi.samples[,2],merged.multi.samples$V2 + 1,merged.multi.samples[,3],merged.multi.samples[,4],merged.multi.samples$V4 + 1,1:nrow(merged.multi.samples))
temp.bedpe = paste(args.list.allbam$OutDir,"/temp.svs.bedpe",sep = "")
write.table(merged.svs,temp.bedpe,sep = '\t',row.names = F,col.names = F,quote = F)
svs.excluded = t(sapply(strsplit(system(paste(bedtools,"pairtobed -type neither -a",temp.bedpe,"-b",exclude_bed),intern = T),"\t"),function(x){x}))
if(nrow(svs.excluded) > 0)
{
merged.multi.samples = merged.multi.samples[as.numeric(svs.excluded[,7]),]
mrc.mat = mrc.mat[as.numeric(svs.excluded[,7]),]
mq20.mat = mq20.mat[as.numeric(svs.excluded[,7]),]
}
}
}else
{
mrc.mat = as.matrix(merged.multi.samples[,6:ncol(merged.multi.samples)])
mrc.mat = apply(gsub("MRC=","",mrc.mat),2,as.numeric)
mq20.mat = mrc.mat

}

merged.svs.bkp.left = data.table(V1 = merged.multi.samples$V1,
                                 V2 = merged.multi.samples$V2,
								 V3=merged.multi.samples$V2,
								 V4=paste("L",1:nrow(merged.multi.samples),sep="_"))

merged.svs.bkp.left$V3[merged.multi.samples$V5=="DEL"] = merged.svs.bkp.left$V3[merged.multi.samples$V5=="DEL"] + 20
merged.svs.bkp.left$V2[merged.multi.samples$V5=="DUP"] = merged.svs.bkp.left$V2[merged.multi.samples$V5=="DUP"] - 20
merged.svs.bkp.left$V2[merged.multi.samples$V5=="INV"] = merged.svs.bkp.left$V2[merged.multi.samples$V5=="INV"] - 10
merged.svs.bkp.left$V3[merged.multi.samples$V5=="INV"] = merged.svs.bkp.left$V3[merged.multi.samples$V5=="INV"] + 10
merged.svs.bkp.left$V2[merged.multi.samples$V5=="TRA"] = merged.svs.bkp.left$V2[merged.multi.samples$V5=="TRA"] - 10
merged.svs.bkp.left$V3[merged.multi.samples$V5=="TRA"] = merged.svs.bkp.left$V3[merged.multi.samples$V5=="TRA"] + 10

merged.svs.bkp.right=data.table(V1=merged.multi.samples$V3,
                                V2=merged.multi.samples$V4,
								V3=merged.multi.samples$V4,
								V4=paste("R",1:nrow(merged.multi.samples),sep="_"))

merged.svs.bkp.right$V2[merged.multi.samples$V5=="DEL"] = merged.svs.bkp.right$V2[merged.multi.samples$V5=="DEL"] - 20
merged.svs.bkp.right$V3[merged.multi.samples$V5=="DUP"] = merged.svs.bkp.right$V3[merged.multi.samples$V5=="DUP"] + 20
merged.svs.bkp.right$V2[merged.multi.samples$V5=="INV"] = merged.svs.bkp.right$V2[merged.multi.samples$V5=="INV"] - 10
merged.svs.bkp.right$V3[merged.multi.samples$V5=="INV"] = merged.svs.bkp.right$V3[merged.multi.samples$V5=="INV"] + 10
merged.svs.bkp.right$V2[merged.multi.samples$V5=="TRA"] = merged.svs.bkp.right$V2[merged.multi.samples$V5=="TRA"] - 10
merged.svs.bkp.right$V3[merged.multi.samples$V5=="TRA"] = merged.svs.bkp.right$V3[merged.multi.samples$V5=="TRA"] + 10

merged.svs.bkp=rbind(merged.svs.bkp.left,merged.svs.bkp.right)



merged.svs.file = paste(args.list.allbam$OutDir,"/merged.svs.bkp.bed",sep = "")
merged.svs.bkp$V2[merged.svs.bkp$V2<0] = 0
fwrite(merged.svs.bkp,merged.svs.file,sep = '\t',row.names = F,col.names = F,quote = F,showProgress = F)

cat(paste("##INFO:",date(),": Calculate proper read count\n",sep = ""))

PPcount.mat=matrix(0,nrow = nrow(merged.svs.bkp),ncol = n.samples)
for(i in 1:n.samples)
{
args.list = args.list.allbam
bam = args.list$BamFile[i]
bam.name = unlist(strsplit(bam,"\\/"))
bam.name = bam.name[length(bam.name)]
outdir = paste(args.list$OutDir,"/",i,".",gsub("\\.bam$","",bam.name),sep = "")
args.list$BamFile = bam
args.list$OutDir = outdir
PP.count = PPcount(args.list,merged.svs.file)
PPcount.mat[match(PP.count[,1],merged.svs.bkp$V4),i] = as.numeric(PP.count[,2])
}

for(i in 1:n.samples)
{
args.list = args.list.allbam
bam = args.list$BamFile[i]
cat(paste("##INFO:",date(),": Genotype SVs for ",bam,"\n",sep = ""))
bam.name = unlist(strsplit(bam,"\\/"))
bam.name = bam.name[length(bam.name)]
outdir = paste(args.list$OutDir,"/",i,".",gsub("\\.bam$","",bam.name),sep = "")
args.list$BamFile = bam
args.list$OutDir = outdir
mrc.i = mrc.mat[,i]
count.mat = as.matrix(cbind(merged.multi.samples[,5],mrc.i,matrix(PPcount.mat[,i],ncol = 2,byrow = F)))
genotype = BayesGenotype(count.mat,args.list)
genotype.format = paste("GT=",genotype,sep = "")
merged.multi.samples[[paste("V",i+5,sep="")]] = paste(genotype.format,
                                                      paste("PP1=",as.numeric(count.mat[,3]),sep = ""),
													  paste("PP2=",as.numeric(count.mat[,4]),sep = ""),
													  merged.multi.samples[[paste("V",5+i,sep = "")]],
													  sep = ":")
}


filter.test = ifelse(apply(mq20.mat>=args.list.allbam$MinMRC,1,sum) > 0,"PASS","LowQuality")

BioSV.output = cbind(merged.multi.samples[,1:2],
                     merged.multi.samples$V2+1,
					 merged.multi.samples[,3:4],
					 merged.multi.samples$V4+1,
					 merged.multi.samples[,5],
					 filter.test,
					 merged.multi.samples[,6:ncol(merged.multi.samples)])

bamfiles = sapply(strsplit(args.list.allbam$BamFile,"\\/"),function(x){x[length(x)]})

colnames(BioSV.output) = c("#Chrom-1","Start-1","End-1","Chrom-2","Start-2","End-2","SVTYPE","FILTER",paste("FORMAT-",bamfiles,sep=""))

rownames(BioSV.output)=NULL

if(args.list.allbam$BiosvMode=="call")
{
BioSV.output = CalDepth(BioSV.output,PPcount.mat,mrc.mat,VAF = 0.05,args.list.allbam)
BioSV.output.withRD=CalDepth(BioSV.output,PPcount.mat,mrc.mat,VAF = 0.05,args.list.allbam,"DUP")

}else
{
BioSV.output.withRD = BioSV.output
}
headers = c("##fileformat=BEDPE","##source=BioSV","##SVTYPE=<ID=DEL,DUP,INV,TRA>","##FILTER=<ID=PASS,Description=\"More than minimal mutant read count and higher than minimal MAPQ\">",
"##FILTER=<ID=LowQuality,Description=\"Otherwise\">", "##FORMAT=<ID=GT,genotype>","##FORMAT=<ID=PP1,Proper read pair count for left breakpoint>","##FORMAT=<ID=PP2,Proper read pair count for right breakpoint>",
"##FORMAT=<ID=MRC,mutant read count>","##FORMAT=<ID=CLIP.MRC,clipped mutant read count>","##FORMAT=<ID=CLIP.strand, Description=\"strand for clipped reads(+-,-+,++,--)\">",
"##FORMAT=<ID=CLIP.direction,Description=\"clipped reads for the inference of breakpoint direction (RL,LR,RR,LL)\">","##FORMAT=<ID=DISC.MRC,Description=\"discordant read pair count\">",
"##FORMAT=<ID=DISC.strand,Description=\"strands for discordant reads pairs (+-,-+,++,--), which can also indicate the direction of breakpoints (RL,LR,RR,LL)\">","##FORMAT=<ID=DISC.MAPQ,Description=\"discordant read pairs MAPQ\"")
headers = c(headers,paste(c("#Chrom-1","Start-1","End-1","Chrom-2","Start-2","End-2","SVTYPE","FILTER",args.list.allbam$BamFile),collapse = "\t"))
write.table(headers,paste(args.list.allbam$OutDir,"/BioSV.output.bedpe",sep = ""),sep = '\t',row.names = F,col.names = F,quote = F)
write.table(BioSV.output.withRD,paste(args.list.allbam$OutDir,"/BioSV.output.bedpe",sep = ""),sep = '\t',row.names = F,col.names = F,quote = F,append = T)
cat(paste("##INFO:",date(),": BioSV finished its work\n",sep = ""))
save.image(paste(args.list.allbam$OutDir,"/BioSV.output.Rdata",sep = ""))

}


