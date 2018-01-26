### Determine the breakpoints of discordant read pairs based on the alignment information,
### such as chromosome, position, strand, and SV type.
### Usage: Rscript Drp.Reads.to.Bkp.R args.list.Rdata

args.list.file=commandArgs(T)
options(scipen=999)

	suppressMessages(library(data.table))
	suppressMessages(library(igraph))
	suppressMessages(library(doParallel))

	load(args.list.file)
	samtools=args.list$samtools
	bedtools=args.list$bedtools
	drp2breakpoint=function(drps)
	{
	breakpoints.left = drps[,c(1,3)]
	breakpoints.left$V3[drps$V9=="-"] = drps$V2[drps$V9=="-"]
	breakpoints.right = drps[,c(4,5)]
	breakpoints.right$V5[drps$V10=="+"] = drps$V6[drps$V10=="+"]
	res = cbind(breakpoints.left,breakpoints.right,drps[,7:10])
	colnames(res) = paste("V",1:ncol(res),sep = "")
	res
	}

	output_dir.PEM = paste(args.list$OutDir,"/DISC",sep = "")
	insert.info = args.list$lib.info
	insertsize=  floor(insert.info[[1]]*3)
	SV.drps = as.data.frame(fread(paste(output_dir.PEM,"/Discordant.Reads.bedpe",sep=""),sep = '\t',showProgress= F ))
	
	chroms = sort(unique(c(SV.drps$V1,SV.drps$V4)))
	chrom.test = SV.drps$V1==SV.drps$V4
	strand.test = SV.drps$V9==SV.drps$V10
	DEL_DUP.drp = SV.drps[chrom.test & !strand.test,]
	DEL_DUP.drp[DEL_DUP.drp$V9=="-",] = DEL_DUP.drp[DEL_DUP.drp$V9=="-",c(4:6,1:3,7:8,10:9)]
	test.DEL = as.numeric(DEL_DUP.drp$V3) <= as.numeric(DEL_DUP.drp$V5)
	test.DUP = as.numeric(DEL_DUP.drp$V6) <= as.numeric(DEL_DUP.drp$V2)
	DEL.drp = DEL_DUP.drp[test.DEL,]
	DUP.drp = DEL_DUP.drp[test.DUP,]
	
	INV.drp = SV.drps[chrom.test & strand.test,]
	INV.drp = rbind(INV.drp[as.numeric(INV.drp$V3)<=as.numeric(INV.drp$V5),],
	                INV.drp[as.numeric(INV.drp$V6)<=as.numeric(INV.drp$V2),c(4:6,1:3,7:8,10:9)])
	TRA.drp=SV.drps[!chrom.test,]
	test.chrom.order=match(TRA.drp$V1,chroms) > match(TRA.drp$V4,chroms)
	TRA.drp[test.chrom.order,] = TRA.drp[test.chrom.order,c(4:6,1:3,7:8,10:9)]
	INV.drp.bkp = drp2breakpoint(INV.drp)
	INV.drp.bkp[INV.drp.bkp$V2 > INV.drp.bkp$V4,c(2,4)] = INV.drp.bkp[INV.drp.bkp$V2 > INV.drp.bkp$V4,c(4,2)]
	drp.breakpoints=list(DEL = drp2breakpoint(DEL.drp),
	                     DUP = drp2breakpoint(DUP.drp),
						 INV = INV.drp.bkp,
						 TRA = drp2breakpoint(TRA.drp))
	for(svtype in names(drp.breakpoints))
	{
		if(svtype!="TRA")
		{
		tmp.obj = sapply(chroms,function(chr){
		outfile = paste(output_dir.PEM,"/",svtype,".",chr,".txt",sep = "")
		fwrite(drp.breakpoints[[svtype]][drp.breakpoints[[svtype]]$V1==chr,],outfile,sep = "\t",row.names = F,col.names = F,quote = F)
		})
		}else
		{
		tmp.chroms = unique(drp.breakpoints[[svtype]][,c(1,3)])
		tmp.obj = apply(tmp.chroms,1,function(x){
		outfile = paste(output_dir.PEM,"/",svtype,".",x[1],".",x[2],".txt",sep = "")
		fwrite(drp.breakpoints[[svtype]][drp.breakpoints[[svtype]]$V1==x[1] & drp.breakpoints[[svtype]]$V3==x[2] ,],outfile,sep="\t",row.names=F,col.names=F,quote=F)
		})
		}
	}


