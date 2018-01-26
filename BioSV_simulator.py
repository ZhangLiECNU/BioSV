# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 15:16:15 2016

@author: DingWB
"""

import pysam
import os,sys
from optparse import OptionParser
import pandas as pd
import random
from scipy import stats

"""
# =============================================================================

           S T A R T   D E F I N I T I O N S 

# =============================================================================
"""
# =============================================================================
def prepare_optparser():
    usage = "usage: python %prog -g <genome> -o <out_dir> [options]"
    description = "BioSV_simulator:Simulate structural variation \
    from diploid genomes and output paired-end fastq files"
    optparser = OptionParser(version="%prog -- 1.0",description=description,\
    usage=usage,add_help_option=False)
    optparser.add_option("-g","--Genome",dest="Reference",type="string",\
    help="Input reference genome file",metavar="<file>")
    optparser.add_option("-o","--out dir",dest="out_dir",type="string",\
    help="output directory",metavar="<path>")
    optparser.add_option("-p","--prefix",dest="prefix",type="string",\
    default="simulate",help="prefix of out file",metavar="<string>")
    optparser.add_option("-T",dest="Tumor",action="store_true",default=False,\
    help="generate tumor file [False]")
    optparser.add_option("-b","--bed",dest="bed",type="string",\
    help="Using supplied bed file instead of randomly generated bed file",\
    metavar="<file>")
    optparser.add_option("-n","--num",dest="n",type="string",\
    default="2500,0,2500,2500,2500",help="the approximate number of \
    structral variation events, seperated by comma in the order of \
    ('DEL','INS','INV','TRA','DUP') [2500,0,2500,2500,2500]")
    optparser.add_option("-m","--mu",dest="mu",type="int",default=2,\
    help="mu of poisson distribution of DUP copy number [2]",metavar="int")
    optparser.add_option("-c","--chrom",dest="chrom",type="string",\
    default=os.path.abspath(sys.path[0]+"/data/chromosomes.txt"),\
    help="Only simulate structral variation in the chromosomes of this \
    file(per chromosome per line),default is %s"%os.path.abspath(sys.path[0]+\
    "/data/chromosomes.txt"))
    optparser.add_option("--sv_len",dest="sv_len",type="string",\
    default="100,10000",help="minimum and maximum length of variation, seperated by comma [100,10000]")
    optparser.add_option("--wgsim","--wgsim",dest="wgsim",type="string",\
    default=os.path.abspath(sys.path[0]+"/tools/wgsim"),help="wgsim path [%s]"\
    %(os.path.abspath(sys.path[0]+"/tools/wgsim")),metavar="file")
    optparser.add_option("-e","--error rate",dest="error",type="float",\
    default=0.02,help="base error rate [0.020]",metavar="float")
    optparser.add_option("-d","--distance",dest="distance",type="int",\
    default=350,help="outer distance between the two ends [350]",metavar="int")
    optparser.add_option("-s","--stand deviation",dest="std",type="int",\
    default=50,help="standard deviation [50]",metavar="int")
    optparser.add_option("-N",dest="N",type="int",default=150000000,\
    help="number of read pairs,default is coverage of 10X [150000000]",metavar="int")
    optparser.add_option("-l",dest="read_length",type="int",default=100,\
    help="length of the reads [100]",metavar="int")
    optparser.add_option("-r",dest="rate",type="float",default=0.001,\
    help="rate of mutations [0.0010]",metavar="float")
    optparser.add_option("-R",dest="fraction",type="float",default=0.15,\
    help="fraction of indels [0.15]",metavar="float")
    optparser.add_option("-X",dest="probability",type="float",default=0.30,\
    help="probability an indel is extended [0.30]",metavar="float")
    optparser.add_option("-S",dest="seed",type="int",default=-1,\
    help="seed for random generator [-1]",metavar="int")
    optparser.add_option("-H",dest="haplotype",action="store_true",\
    default=False,help="haplotype mode [False]",metavar="bool")
    optparser.add_option("-h","--help",action="help",\
    help="Show this help message and exit")
    return optparser
# =============================================================================
class RunError(Exception):
    pass
# =============================================================================
def pyshell(cmd):
    sys.stdout.write(cmd+'\n')
    e=os.system(cmd)
    if e!=0:
        raise RunError("Wrong![CMD]\n%s\nRun failed"%cmd)
# =============================================================================
def generateBed(genome,chrom_file,n,sv_len):
    with open(chrom_file,'r') as f:
        chroms=f.read().strip().split('\n')
    D={}
    f=pysam.FastaFile(genome)
    for chrom in chroms:
        start=1
        while f.fetch(chrom,start-1,start)=='N':
            if f.fetch(chrom,start+2000-1,start+2000)=='N':
                start+=2000
            start+=1
        end=f.get_reference_length(chrom)
        while f.fetch(chrom,end-1,end)=='N':
            if f.fetch(chrom,end-2000-1,end-2000)=='N':
                end-=2000
            end-=1
        D[chrom]=[start,end]
    patterns=['DEL']*n[0]+['INS']*n[1]+['INV']*n[2]+['TRA']*n[3]+['DUP']*n[4]
    L=[]
    for pattern in patterns:
        chrom=random.choice(chroms)
        length=random.randint(sv_len[0],sv_len[1])
        genotype=random.choice(['0/1','1/0','1/1'])
        start=random.randint(D[chrom][0],D[chrom][1]-length)
        end=start+length-1
        sg=random.choice(['somatic','germline'])
        if pattern=="DEL":
            L.append([chrom,start,end,"*","*","*","*",pattern,genotype,sg])
        elif pattern=='INS':
            chrom2=random.choice(chroms)
            start2=random.randint(D[chrom2][0],D[chrom2][1]-length)
            end2=start2+length-1
            strand=random.choice(['+','-'])
            L.append([chrom,start,start,chrom2,start2,end2,strand,pattern,genotype,sg])
        elif pattern=='INV':
            L.append([chrom,start,end,"*","*","*","*",pattern,genotype,sg])
        elif pattern=='TRA':
            chrom2=random.choice([ch for ch in chroms if ch!=chrom])
            start2=random.randint(D[chrom2][0],D[chrom2][1]-length)
            end2=start2+length-1
            strand=random.choice(['+','-'])
            L.append([chrom,start,end,chrom2,start2,end2,strand,pattern,genotype,sg])
            L.append([chrom2,start2,end2,chrom,start,end,strand,pattern,genotype,sg])
        elif pattern=="DUP":
            L.append([chrom,end+1,end+1,chrom,start,end,'+',pattern,genotype,sg])
    f.close()
    return pd.DataFrame(L,columns=['#chrom1','start1','end1','chrom2','start2','end2','strand','pattern','genotype','sg'])
# =============================================================================
def rev_comp(seq):
    D = {
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c',
    'k' : 'm',
    'm' : 'k',
    'r' : 'y',
    'y' : 'r',
    's' : 's',
    'w' : 'w',
    'b' : 'v',
    'v' : 'b',
    'h' : 'd',
    'd' : 'h',
    'n' : 'n',
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
    'K' : 'M',
    'M' : 'K',
    'R' : 'Y',
    'Y' : 'R',
    'S' : 'S',
    'W' : 'W',
    'B' : 'V',
    'V' : 'B',
    'H' : 'D',
    'D' : 'H',
    'N' : 'N'}
    r=[D[i] for i in seq]
    return ''.join(r[::-1])
# =============================================================================
def bed2bedpe(bed,bedpe):
    L=[]
    with open(bed,'r') as f:
        lines=f.readlines()
        for line in lines[1:]:
            chrom1,start1,end1,chrom2,start2,end2,strand,pattern,genotype,sg=line.strip().split('\t')
            start1=int(start1)
            end1=int(end1)
            if pattern=='TRA':
                if strand=='+':
                    L.append([chrom1,start1-1,chrom2,int(start2),'+',strand,pattern,genotype,sg])
                    L.append([chrom2,int(end2),chrom1,end1+1,strand,'+',pattern,genotype,sg])
                else:
                    L.append([chrom1,start1-1,chrom2,int(end2),'+',strand,pattern,genotype,sg])
                    L.append([chrom2,int(start2),chrom1,end1+1,strand,'+',pattern,genotype,sg])
            elif pattern=='INS':
                if strand =='+':
                    L.append([chrom1,start1-1,chrom2,int(start2),'+',strand,pattern,genotype,sg])
                    L.append([chrom2,int(end2),chrom1,end1,strand,'+',pattern,genotype,sg])
                else:
                    L.append([chrom1,start1-1,chrom2,int(end2),'+',strand,pattern,genotype,sg])
                    L.append([chrom2,int(start2),chrom1,end1,strand,'+',pattern,genotype,sg])
            elif pattern=='DEL':
                L.append([chrom1,start1-1,chrom1,end1+1,'+','+',pattern,genotype,sg])
            elif pattern=='INV':
                L.append([chrom1,start1-1,chrom1,end1,'+','-',pattern,genotype,sg])
#                L.append([chrom1,start1,chrom1,end1+1,'-','+',pattern,genotype])
            else:#DUP
                L.append([chrom1,end1,chrom1,int(start2),'+','+',pattern,genotype,sg])
    df=pd.DataFrame(L,columns=['chrI','breakpointI','chrII','breakpointII','strandI','strandII','pattern','genotype','sg'])
    df['startI']=df.breakpointI.apply(int)-1
    df['startII']=df.breakpointII.apply(int)-1
    df[['chrI','startI','breakpointI','chrII','startII','breakpointII','strandI','strandII','pattern'\
    ,'genotype','sg']].rename(columns={'chrI':'#chrI'}).to_csv(bedpe,header=True,index=False,sep='\t')
# =============================================================================
def writeModifiedGenome(gnome,prefix,array,names,mu,TN,haplotype=False):
    f=pysam.FastaFile(gnome)
    f1=open(prefix+'_'+names[0]+'.fa','w')
    if haplotype:
        f2=open(prefix+'_'+names[1]+'.fa','w')
    f0=open(prefix+'_'+TN+'.bed','w')
    f0.write("\t".join(['#chrom1','start1','end1','chrom2','start2','end2','strand','pattern','genotype','sg'])+'\n')
    chrom=None
    flag=False
    while True:
        try:
            line=array.pop(0)
        except:
            break
        if line[0]!=chrom:
            if flag:
                f1.write(f.fetch(chrom,start1,)+'\n')
                if haplotype:
                    f2.write(f.fetch(chrom,start2,)+'\n')
            flag=True
            chrom=line[0]
            start1=0
            start2=0
            f1.write(">%s\n"%(line[0]+"_"+names[0]))
            if haplotype:
                f2.write(">%s\n"%(line[0]+"_"+names[1]))
#            print(line[0])
        if line[1]<=start1 or line[1]<=start2:
            continue
        if line[1]<10:
            continue
        if 'N' in f.fetch(chrom,int(line[1])-10,int(line[2])+10):
            continue
        if line[4]!='*' and int(line[4])<10:
            continue
        if line[4]!='*' and 'N' in f.fetch(line[3],int(line[4])-10,int(line[5])+10):
            continue
        chrom,pattern,genotype,sg=line[0],line[7],line[8].split('/'),line[9]
        if sg=='somatic' and TN=='N':
            continue
        f0.write('\t'.join([str(i) for i in line])+'\n')
        f1.write(f.fetch(chrom,start1,line[1]-1))
        if haplotype:
            f2.write(f.fetch(chrom,start2,line[1]-1))
        if pattern=='DEL':
            start1=line[2] if genotype[0]=='1' else line[1]-1
            start2=line[2] if genotype[1]=='1' else line[1]-1
        elif pattern=='INS':
            seq=f.fetch(line[3],line[4]-1,line[5]) if line[6]=='+' \
            else rev_comp(f.fetch(line[3],line[4]-1,line[5]))
            if genotype[0]=='1':
                f1.write(seq)
            if haplotype and genotype[1]=='1':
                f2.write(seq)
            start1=start2=line[2]-1
        elif pattern=='DUP':
            seq=f.fetch(line[3],int(line[4])-1,int(line[5])) if line[6]=='+' \
            else rev_comp(f.fetch(line[3],int(line[4])-1,int(line[5])))
            copy_num=stats.poisson.rvs(mu,1)#poisson distribution
            if genotype[0]=='1':
                for i in range(copy_num):
                    f1.write(seq)
            if haplotype and genotype[1]=='1':
                for i in range(copy_num):
                    f2.write(seq)
            start1=start2=line[2]-1
        elif pattern=='INV':
            seq=rev_comp(f.fetch(chrom,line[1]-1,line[2])) if genotype[0]=='1' \
            else f.fetch(chrom,line[1]-1,line[2])
            f1.write(seq)
            if haplotype:
                seq=rev_comp(f.fetch(chrom,line[1]-1,line[2])) if genotype[1]=='1' \
                else f.fetch(chrom,line[1]-1,line[2])
                f2.write(seq)
            start1=start2=line[2]
        elif pattern=="TRA":
            seq=f.fetch(line[3],int(line[4])-1,int(line[5])) if line[6]=='+' \
            else rev_comp(f.fetch(line[3],int(line[4])-1,int(line[5])))
            seq1=seq if genotype[0]=='1' else f.fetch(chrom,line[1]-1,line[2])
            f1.write(seq1)
            if haplotype:
                seq2=seq if genotype[1]=='1' else f.fetch(chrom,line[1]-1,line[2])
                f2.write(seq2)
            start1=start2=line[2]
    f1.write(f.fetch(chrom,start1,)+'\n')
    if haplotype:
        f2.write(f.fetch(chrom,start2,)+'\n')
#            print(chrom)
    f.close()
    f0.close()
    f1.close()
    bed2bedpe(prefix+'_'+TN+'.bed',prefix+'_'+TN+'.bedpe')
    if haplotype:
        os.system("cat %s %s > %s"%(prefix+'_'+names[0]+'.fa',prefix+'_'+names[1]+'.fa',prefix+'_'+TN+'.fa'))
        os.remove(prefix+'_'+names[0]+'.fa')
        os.remove(prefix+'_'+names[1]+'.fa')
        f2.close()
# =============================================================================
"""
# =============================================================================

           E N D    O F    D E F I N I T I O N S 

# =============================================================================
"""
if __name__ == "__main__":
    optparser=prepare_optparser()
    (options,args) = optparser.parse_args()
    if not (options.Reference and options.out_dir):
        optparser.print_help()
        sys.exit("-G -O must be assigned")
    if not os.path.exists(options.out_dir):
        os.mkdir(options.out_dir)
    options.n=[int(i) for i in options.n.strip().split(',')]
    options.sv_len=[int(i) for i in options.sv_len.strip().split(',')]
    options.Reference=os.path.abspath(options.Reference)
    if not options.bed:
        df=generateBed(options.Reference,options.chrom,options.n,options.sv_len)
    else:
        df=pd.read_table(os.path.abspath(options.bed),header=None)
        df.columns=['#chrom1','start1','end1','chrom2','start2','end2','strand','pattern','genotype','sg']
    df=df.sort(['#chrom1','start1','end1','pattern'])
    sys.stdout.write("Generate simulated fasta file\n")
    prefix=os.path.abspath(options.out_dir)+"/"+options.prefix
    writeModifiedGenome(options.Reference,prefix,df.values.tolist(),['NA','NB'],options.mu,'N',haplotype=True)
    if options.Tumor:
        writeModifiedGenome(options.Reference,prefix,df.values.tolist(),['TA','TB'],options.mu,'T',haplotype=True)
    haplotype="-h" if options.haplotype else ""
    pyshell(options.wgsim+" -e %s -d %s -s %s -N %s -1 %s -2 %s -r %s -R %s -X %s -S %s %s %s %s %s > %s"\
    %(options.error,options.distance,options.std,options.N,options.read_length,options.read_length,options.rate,options.fraction,options.probability,\
    options.seed,haplotype,prefix+"_N.fa",prefix+"_N1.fq",prefix+"_N2.fq",prefix+"_Nsnp.txt"))
    
    if options.Tumor:
        pyshell(options.wgsim+" -e %s -d %s -s %s -N %s -1 %s -2 %s -r %s -R %s -X %s -S %s %s %s %s %s > %s"\
        %(options.error,options.distance,options.std,options.N,options.read_length,options.read_length,options.rate,options.fraction,options.probability,\
        options.seed,haplotype,prefix+"_T.fa",prefix+"_T1.fq",prefix+"_T2.fq",prefix+"_Tsnp.txt"))
