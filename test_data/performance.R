overlap=t(sapply(strsplit(system("./tools/bedtools pairtopair -slop 100 -a ./test_data/simulate/simu.bedpe -b ./test_data/BioSV_output/predicted.hc.bedpe",intern=T),"\t"),function(x){x}))
overlap=overlap[overlap[,11]==overlap[,23],]
overlap[overlap=="1/0"]="0/1"
simu=as.matrix(read.table("./test_data/simulate/simu.bedpe",sep='\t'))
hc.pred=as.matrix(read.table("./test_data/BioSV_output/predicted.hc.bedpe",sep='\t'))
sens=table(unique(overlap[,c(7,11)])[,2])/table(simu[,11])
prec=table(unique(overlap[,c(7+12,11+12)])[,2])/table(hc.pred[,11])
acc=sapply(c("DEL","DUP","INV","TRA"),function(x){m=overlap[overlap[,11]==x,c(12,24)];sum(m[,1]==m[,2])/nrow(m)})
perf=rbind(sensitivity=sens,precision=prec,accuracy=acc)
write.table(perf,"./test_data/performance.txt",sep='\t',quote=F)


