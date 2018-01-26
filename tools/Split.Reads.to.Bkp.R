### Determine the breakpoints of split reads based on the alignment information,
### such as chromosome, position, strand, and SV type.
### Usage: Rscript Split.Reads.to.Bkp.R args.list.Rdata

args.list.file=commandArgs(T)
    options(scipen = 999)
    load(args.list.file)
	## load libraries
	suppressMessages(library(data.table))
	suppressMessages(library(igraph))
	suppressMessages(library(doParallel))

	samtools = args.list$samtools
	bedtools = args.list$bedtools
	aln2SVlist = function(clipped.SH.aln)
	{
	SH.pairs = fread(clipped.SH.aln,sep = '\t',showProgress = F)
	chroms = sort(unique(c(SH.pairs$V1)))
	soft.aln = SH.pairs[grep("S",SH.pairs$V7),]
	hard.aln = SH.pairs[grep("H",SH.pairs$V7),]
	dup.soft = soft.aln$V4[duplicated(soft.aln$V4)]
	dup.hard = hard.aln$V4[duplicated(hard.aln$V4)]
	uniq.soft = soft.aln$V4[is.na(match(soft.aln$V4,dup.soft))]
	uniq.hard = hard.aln$V4[is.na(match(hard.aln$V4,dup.hard))]
	SH.pairs.read.names = uniq.soft[!is.na(match(uniq.soft,uniq.hard))]

	soft.aln = soft.aln[match(SH.pairs.read.names,soft.aln$V4),]
	hard.aln = hard.aln[match(SH.pairs.read.names,hard.aln$V4),]

	soft.cigar.str = sapply(strsplit(soft.aln$V7,"[0-9]"),function(x){x[x!=""]})
	hard.cigar.str = sapply(strsplit(hard.aln$V7,"[0-9]"),function(x){x[x!=""]})
	soft.cigar.base = strsplit(soft.aln$V7,"[A-Z]")
	hard.cigar.base = strsplit(hard.aln$V7,"[A-Z]")
	soft.cigar.nchar = nchar(soft.aln$V7)
	soft.end.test = substr(soft.aln$V7,soft.cigar.nchar,soft.cigar.nchar)=="S"
	hard.cigar.nchar = nchar(hard.aln$V7)
	hard.end.test = substr(hard.aln$V7,hard.cigar.nchar,hard.cigar.nchar)=="H"

	soft.S.count = sapply(soft.cigar.str,function(x){sum(x=="S")})

	hard.H.count = sapply(hard.cigar.str,function(x){sum(x=="H")})

	direction.hard = rep("L",nrow(hard.aln))
	direction.soft = direction.hard

	direction.soft[soft.S.count==1] = ifelse(soft.end.test[soft.S.count==1],"R","L")
	direction.hard[hard.H.count==1] = ifelse(hard.end.test[hard.H.count==1],"R","L")
	direction.soft[soft.S.count>1] = ifelse(sapply(soft.cigar.base[soft.S.count>1],function(x){
		x = as.numeric(x)
		y = x[1]<x[length(x)]
		}),"R","L")

	direction.hard[hard.H.count>1] = ifelse(sapply(hard.cigar.base[hard.H.count>1],function(x){
		x = as.numeric(x)
		y = x[1]<x[length(x)]
		}),"R","L")
	soft.breakpoint = soft.aln$V2
	soft.breakpoint[direction.soft=="R"] = soft.aln$V3[direction.soft=="R"]
	hard.breakpoint = hard.aln$V2
	hard.breakpoint[direction.hard=="R"] = hard.aln$V3[direction.hard=="R"]


	sv.breakpoints = cbind(soft.aln[,1],
	                       soft.breakpoint,
						   hard.aln[,1],
						   hard.breakpoint,
						   soft.aln[,4:6],
						   hard.aln[,5:6],
						   direction.soft,
						   direction.hard)
	colnames(sv.breakpoints) = paste("V",1:ncol(sv.breakpoints),sep = "")
	chrom.test = sv.breakpoints$V1==sv.breakpoints$V3
	pos.test = sv.breakpoints$V2<sv.breakpoints$V4
	strand.test = sv.breakpoints$V7==sv.breakpoints$V9
	direction.test = sv.breakpoints$V10==sv.breakpoints$V11

	TRA = sv.breakpoints[!chrom.test ,]
	INV = sv.breakpoints[chrom.test & !strand.test,]
	DEL_DUP.split = sv.breakpoints[chrom.test & strand.test,]
	DEL_DUP.split[DEL_DUP.split$V10=="L",] = DEL_DUP.split[DEL_DUP.split$V10=="L",c(3,4,1,2,5,8,9,6,7,11,10)]
	DEL = DEL_DUP.split[DEL_DUP.split$V2 < DEL_DUP.split$V4,]
	DUP = DEL_DUP.split[DEL_DUP.split$V2 > DEL_DUP.split$V4,]

	INV[INV$V2 > INV$V4,] = INV[INV$V2 > INV$V4,c(3,4,1,2,5,8,9,6,7,11,10)]
	chr1.index = match(TRA$V1,chroms)
	chr2.index = match(TRA$V3,chroms)
	TRA[chr1.index > chr2.index,] = TRA[chr1.index > chr2.index,c(3,4,1,2,5,8,9,6,7,11,10)]
	SV.list = list( DEL = DEL,
	                DUP = DUP,
					INV = INV,
					TRA = TRA)
	}
	

	merge.bkp=function(sv,label)
	{
	MAPQ = tapply(apply(sv[,c(6,8)],1,function(x){min(as.numeric(x))}),label,function(x){paste(x,collapse = ";")})
	strand = tapply(paste(sv$V7,sv$V9,sep = ""),label,function(x){paste(x,collapse = ";")})
	direction = tapply(paste(sv$V10,sv$V11,sep = ""),label,function(x){paste(x,collapse = ";")})
	read.name = tapply(sv$V5,label,function(x){paste(x,collapse = ";")})
	read.count = table(label)
	res = cbind(sv[match(levels(label),as.matrix(label)),c(1:4)],read.count,MAPQ,strand,direction,read.name)
	rownames(res) = NULL
	res
	}
	

	bam = args.list$BamFile
	output_dir = args.list$OutDir
	threads = args.list$Threads
	sc.mismatch = args.list$lib.info[[2]] - 60
	output_dir.clipped = paste(output_dir,"/CLIP",sep = "")
	clipped.SH.aln = paste(output_dir.clipped,"/clipped.SH.aln",sep = "")
	SV.merge = aln2SVlist(clipped.SH.aln)
	for(svtype in names(SV.merge))
	{
		if(svtype!="TRA")
		{
		chroms = unique(SV.merge[[svtype]]$V1)
		tmp.obj = sapply(chroms,function(x){
		outfile = paste(output_dir.clipped,"/",svtype,".",x,".txt",sep = "")
		fwrite(SV.merge[[svtype]][SV.merge[[svtype]]$V1==x,],outfile,sep = "\t",row.names = F,col.names = F,quote = F)
		})
		}else
		{
		tmp.chroms = unique(SV.merge[[svtype]][,c(1,3)])
		tmp.obj = apply(tmp.chroms,1,function(x){
		outfile = paste(output_dir.clipped,"/",svtype,".",x[1],".",x[2],".txt",sep = "")
		fwrite(SV.merge[[svtype]][SV.merge[[svtype]]$V1==x[1] & SV.merge[[svtype]]$V3==x[2] ,],outfile,sep = "\t",row.names = F,col.names = F,quote = F)
		})
		}
	}


