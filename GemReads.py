#!/usr/bin/python -W ignore

# Copyright (C) 2010, Kerensa McElroy.
# kerensa@unsw.edu.au

# This file is part of the sequence simulator GemSIM. 
# It is used to calculate simulated next-gen sequencing
# reads, based on a set of reference genomes and an 
# error model.

# GemSIM is free software; it may be redistributed and 
# modified under the terms of the GNU General Public 
# License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option)
# any later version.

# GemSIM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, without even the implied 
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more 
# details.

# You should have recieved a copy of the GNU General Public
# License along with GemSIM. If not, see 
# http://www.gnu.org/licenses/. 


import sys
import getopt
import random
import numpy
import bisect
import cPickle
import gzip
import os
import logging
import logging.handlers
import Experiment as E


# Make a global logging object.
rdlog=logging.getLogger("RDlog")

# Set logging level, and write everything to a file
rdlog.setLevel(logging.DEBUG)
LOG_FILENAME='./rd.log'
h=logging.FileHandler(LOG_FILENAME,'w')
f=logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
h.setFormatter(f)
rdlog.addHandler(h)

def getRef(refFile):
    """Returns a genome reference, and the length of that reference."""
    refDict={}
    comDict={}
    choiceList=[]
    hdList=[]
    ref=''
    num=0
    try:
        f=open(refFile)
    except IOError:
        rdlog.error('Cannot find reference file ' +refFile+'. Please check pathname.')
        sys.exit('Cannot find reference file '+refFile+'. Please check pathname.')
    i=f.readline()
    head=i[1:51].rstrip()
    i=f.readline().rstrip()    
    while i:
        if i[0]!='>':
            ref+=i.rstrip()
            i=f.readline()
        else:
            if head in hdList:
                num+=1
                head=head+str(num)
            ref=ref.upper()
            for l in 'RYMKSWHBVD':
                ref=ref.replace(l,'N')
            name=refFile.split('/')[-1]
            refDict[head]=(ref,name)
            comDict[head]=Comp(ref)
            choiceList.append((head,len(ref)))
            hdList.append(head)
            head=i[1:51].rstrip()
            i=f.readline()    
            ref=''
    ref=ref.upper()
    for l in 'RYMKSWHBVD':
        ref=ref.replace(l,'N')
    name=refFile.split('/')[-1]
    refDict[head]=(ref,name)
    comDict[head]=Comp(ref)
    choiceList.append((head,len(ref))) 
    rdlog.debug('Reference file successfully parsed.')
    return refDict,comDict, choiceList

def getMet(directory,genDir, abund):
    """Returns the genomes in a directory, with weights."""
    refDict={}
    comDict={}
    genDict={}
    choiceList=[]
    try:
        a=open(abund)
    except IOError:
        rdlog.error('Cannot find species abundance file ' +abund+'. Please check pathname.')
        sys.exit('Cannot find species abundance file ' +abund+'. Please check pathname.')
    ab={}
    sum=0.0
    for i in a:
        parts=i.rstrip().split('\t')
        val=float(parts[1])
        ab[parts[0]]=val
        sum+=val
    key=ab.keys()
    for k in key:
        val=ab[k]
        val=val/sum
        ab[k]=val
    path1='./'+directory
    path2='./'+genDir
    for k in key:
        file=os.path.join(path1,k)
        rD, cD,choice=getRef(file)

        # resolve names with "." in basename
        # will be caught later if errors
        file_parts = k.split('.')
        if len(file_parts) > 1:
            stem = ".".join(file_parts[0:len(file_parts)-1])
        else:
            stem=k.split('.')[0]

        print stem
        try:
            gen=os.path.join(path2,stem+'.txt')
            gens=gParse(gen)
            gens=bisect_gens(gens)
        except:
            gens=genRef('')
            print stem
        for r,v in choice:
            val=v*ab[k]
            refDict[r]=rD[r]
            comDict[r]=cD[r]
            choiceList.append((r,val))
        genDict[k]=gens
    return refDict,genDict,comDict,choiceList

def readGen1(ref,cRef,refLen,readLen,genos,inter,mx1,insD1,delD1,gQ,bQ,iQ,qual,circular,hd):
    """Generates a random read of desired length from a reference."""
    ind=random.randint(0,(refLen-1))
    dir=random.randint(1,2)
    readPlus=int(readLen*1.5)
    if circular:
        if dir==1:
            end=ind+readPlus
            if end<=refLen:
                read=ref[ind:end]
            else:
                read=ref[ind:]+ref[:end-refLen]
            if genos!='':
                read=mutate(read,ind,genos,refLen,1,readPlus,hd)
            read, quals=mkErrors(read,readLen,mx1,insD1,delD1,gQ,bQ,iQ,qual)
        elif dir==2:
            end=ind-readPlus+1
            if end>=0:
                read=cRef[end:ind+1]
            else:
                read=cRef[end:]+cRef[:ind+1]
            if genos!='':
                read=mutate(read,end,genos,refLen,2,readPlus,hd)
            read=read[::-1]
            read,quals=mkErrors(read,readLen,mx1,insD1,delD1,gQ,bQ,iQ,qual)
    else:
        if dir==1:
            frag=ind+inter
            end=ind+readPlus
            while frag>=refLen or end>=refLen:
                ind=random.randint(0,(refLen-inter))
                frag=ind+inter
                end=ind+readPlus
            read=ref[ind:end]
            if genos!='':
                read=mutate(read,ind,genos,refLen,1,readPlus,hd)
            read,quals=mkErrors(read,readLen,mx1,insD1,delD1,gQ,bQ,iQ,qual)
        elif dir==2:
            frag=ind-inter+1
            while frag<0:
                ind=random.randint(0,(refLen-1))
                frag=ind-inter+1
            end=ind-readPlus+1
            if end <0:
                end=0
            read=cRef[end:ind+1]
            if genos!='':
                read=mutate(read,end,genos,refLen,2,readPlus,hd)
            read=read[::-1]
            read,quals=mkErrors(read,readLen,mx1,insD1,delD1,gQ,bQ,iQ,qual)
    return read, ind, dir, quals

def readGen2(reference,cRef,pos,dir,readLen,genos,inter,mx2,insD2,delD2,gQ,bQ,iQ,qual,circular,hd):
    """Generates the 2nd read of a random pair of reads."""
    refLen=len(reference)
    readPlus=int(readLen*1.5)
    if circular:
        if dir==1:
            end=pos+inter
            start=end-readPlus
            if start<0 and end>=0:
                read=cRef[start:]+cRef[:end]
            elif end<=refLen:
                read=cRef[start:end]
            elif start<refLen and end >= refLen:
                read=cRef[start:]+cRef[:end-refLen]
                end-=refLen
            elif start>=refLen:
                read=cRef[start-refLen:end-refLen]
                end-=refLen
            if genos!='':
                if 0<=start<refLen:
                    read=mutate(read,start,genos,refLen,2,readPlus,hd)
                elif start<0:
                    read=mutate(read,start+refLen,genos,refLen,2,readPlus,hd)
                else:
                    read=mutate(read,start-refLen,genos,refLen,2,readPlus,hd)
            read=read[::-1]
            read, quals=mkErrors(read,readLen,mx2,insD2,delD2,gQ,bQ,iQ,qual)
        else:
            start=pos-inter+1
            end=start+readPlus
            if start>=0:
                read=reference[start:end] 
            elif start<0 and end >=0:
                read=reference[start:]+reference[:end]
                start+=refLen
            elif start < 0 and end <0:
                read=reference[start:end] 
                start+=refLen
            if genos!='':
                read=mutate(read,start,genos,refLen,1,readPlus,hd)
            read,quals=mkErrors(read,readLen,mx2,insD2,delD2,gQ,bQ,iQ,qual)
    else:
        if dir==1:
            end=pos+inter
            start=end-readPlus
            if start<0:
                start=0
            read=cRef[start:end]
            if genos!='':
                read=mutate(read,start,genos,refLen,2,readPlus,hd)
            read=read[::-1]
            read,quals=mkErrors(read,readLen,mx2,insD2,delD2,gQ,bQ,iQ,qual)
        else:
            start=pos-inter+1
            end=start+readPlus
            read=reference[start:end]
            if genos!='':
                read=mutate(read,start,genos,refLen,1,readPlus,hd)
            read,quals=mkErrors(read,readLen,mx2,insD2,delD2,gQ,bQ,iQ,qual)
    return read, quals

def mkInserts(mx,insD):
    """Returns a dictionary consisting of compiled functions to make inserts."""
    insertDict={}
    posKeys=insD.keys()
    posKeys.sort()
    for p in posKeys:
        indicies=p.split('.')
        tot=mx[int(indicies[0])][int(indicies[1])][int(indicies[2])][int(indicies[3])][int(indicies[4])][int(indicies[5])][5]
        insertKeys=insD[p].keys()
        insertKeys.sort() 
        insertList=[]
        iSum=0
        for i in insertKeys:
            insertList.append((i,insD[p][i])) 
            iSum+=0
        insertList.append(('',tot-iSum))
        insert=bisect_choiceTUP(insertList)
        insertDict[p]=insert
    return insertDict

def mkDels(mx,delD):
    """Returns a dictionary consisting of compiled functions to make deletiosn."""
    deletionDict={}
    posKeys=delD.keys()
    posKeys.sort()
    for p in posKeys:
        indicies=p.split('.')
        tot=mx[int(indicies[0])][int(indicies[1])][int(indicies[2])][int(indicies[3])][int(indicies[4])][int(indicies[5])][5]
        items=delD[p] 
        items.reverse()
        items.append(tot-sum(items))
        items.reverse()
        delete=bisect_choice(items)
        deletionDict[p]=delete
    return deletionDict

def bisect_choiceTUP(items):
    """Returns a function that makes a weighted random choice from a list of tuples."""
    added_weights = []
    last_sum = 0.0
    for item,weight in items:
        weight=float(weight)
        last_sum += weight
        added_weights.append(last_sum)
    def choice(rnd=random.random, bis=bisect.bisect):
        return items[bis(added_weights, rnd() * last_sum)][0]
    return choice

def genRef(ref):
    """Returns input as function"""
    def r():
        return ref
    return r

def ln(length):
    """Returns static length as a funtion."""
    def val():
        return length
    return val

def bisect_choice(items):
    """Returns a function that makes a weighted random choice from items."""
    added_weights = []
    last_sum = 0
    for weight in items:
        last_sum += weight
        added_weights.append(last_sum)
    def choice(rnd=random.random, bis=bisect.bisect):
        return bis(added_weights, rnd() * last_sum)
    return choice

def bisect_gens(items):
    added_weights = []
    last_sum = 0
    for weight,dict in items:
        last_sum +=weight
        added_weights.append(last_sum)
    def choice(rnd=random.random, bis=bisect.bisect):
        return items[bis(added_weights,rnd()*last_sum)][1]
    return choice

def mutate(read,ind,gens,refLen,dir,readLn,hd):
    """Adds predetermined mutations to reads."""
    d={'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c','N':'N','n':'n'}
    if gens=={}:
        return read    
    else:
        chroms=gens.keys()
        if hd not in chroms:
            return read
        else:
            posi=gens[hd].keys()
            if dir==1:
                for p in posi:
                    if p >ind and p<=(ind+readLn):
                        read1=read[:p-(ind+1)]+gens[hd][p]
                        read1=read1+read[p-ind:]
                        read=read1
                    elif p<=ind+readLn-refLen: 
                        read1=read[:refLen-ind+p-1]+gens[hd][p]
                        read1+=read[refLen-ind+p:]
                        read=read1
                return read
            elif dir==2:
                for p in posi:
                    if p >ind and p<=(ind+readLn):
                        read1=read[:p-(ind+1)]+d[gens[hd][p]]
                        read1=read1+read[p-ind:]
                        read=read1
                    elif p<=ind+readLn-refLen:
                        read1=read[:refLen-ind+p-1]+d[gens[hd][p]]
                        read1+=read[refLen-ind+p:]
                        read=read1
                return read

def parseModel(gzipFile,paired,readlen):
    """prepares error models for input to mkErrors."""
    rdlog.debug("Parsing error model file")
    try: 
        file=gzip.open(gzipFile,'rb')
    except IOError:
        rdlog.error("Cannot find input error model file "+gzipFile+". Please check pathname.") 
        sys.exit()
    if paired:
        modReadLen=cPickle.load(file)
        if readlen!='d' and readlen>modReadLen:
            rdlog.error("Inappropriate read length chosen for model. Maximum for this model: "+str(modReadLen))
            file.close()
            sys.exit() 
        mx1=cPickle.load(file)
        mx2=cPickle.load(file)
        insD1=cPickle.load(file)
        insD2=cPickle.load(file)
        delD1=cPickle.load(file)
        delD2=cPickle.load(file)
        intD=cPickle.load(file)
        gQualL=cPickle.load(file)
        bQualL=cPickle.load(file)
        iQualL=cPickle.load(file)
        mates=cPickle.load(file)
        rds=cPickle.load(file)
        rdLenD=cPickle.load(file)
        file.close()
        return mx1,mx2,insD1,insD2,delD1,delD2,intD,gQualL,bQualL,iQualL,mates,rds,rdLenD
    else:
        modReadLen=cPickle.load(file)
        if readlen!='d' and readlen>modReadLen:
            rdlog.error("Inappropriate read length chosen for model. Maximum for this model: "+str(modReadLen))
            file.close()
            sys.exit() 
        mx=cPickle.load(file)
        insD=cPickle.load(file)
        delD=cPickle.load(file)
        gQualL=cPickle.load(file)
        bQualL=cPickle.load(file)
        iQualL=cPickle.load(file)
        readCount=cPickle.load(file)
        rdLenD=cPickle.load(file)
        file.close()
        return mx,insD,delD,gQualL,bQualL,iQualL,readCount,rdLenD 
        
def mkErrors(read,readLen,mx,insD,delD,gQ,bQ,iQ,qual):
    """Adds random errors to read."""
    inds={'A':0,'T':1,'G':2,'C':3,'N':4,'a':0,'t':1,'g':2,'c':3,'n':4}
    pos=0
    quals=''
    index='0.4.4.4.4.'+str(inds[read[0]])
    if index in insD:
        insert=insD[index]()
        read='NNNN'+insert+read
        for i in insert:
            quals+=iQ[0]()
            pos+=1
    else:
        read='NNNN'+read
    pos+=1
    while pos<=readLen and pos<len(read)-4:                 
        prev=read[pos:pos+4]
        after=read[pos+4]
        d0=pos
        d1=inds[prev[3]]
        d2=inds[prev[2]]
        d3=inds[prev[1]]
        d4=inds[prev[0]]
        d5=inds[after]
        index=str(d0)+'.'+str(d1)+'.'+str(d2)+'.'+str(d3)+'.'+str(d4)+'.'+str(d5)
        tot=float(mx[d0][d1][d2][d3][d4][d5][5])

        # hack to force continue
        if tot == 0:
            E.warn("cannot process for %s read %s" % (mx, read))
            continue
        Mprobs=mx[d0][d1][d2][d3][d4][d5]/tot
        val=random.random()
        a=Mprobs[0]
        t=Mprobs[1]+a
        g=Mprobs[2]+t
        c=Mprobs[3]+g
        n=Mprobs[4]+c
        success=False
        if val>n or tot == 0:
            gPos=pos-1
            while gPos>=0:
                try:
                    quals+=gQ[gPos]()
                    success=True
                    break
                except:
                    gPos-=1
            if success==False:
                quals+=chr(30+qual) 
        elif val>c:
            read=read[:pos+3]+'N'+read[pos+4:]
            bPos=pos-1
            while bPos>=0:
                try:
                    quals+=bQ[bPos]()
                    success=True
                    break
                except:
                    bPos-1
                if success==False:
                    quals+=chr(2+qual)
        elif val>g:
            read=read[:pos+3]+'C'+read[pos+4:]
            bPos=pos-1
            while bPos>=0:
                try:
                    quals+=bQ[bPos]()
                    success=True
                    break
                except:
                    bPos-1
                if success==False:
                    quals+=chr(2+qual)
        elif val>t:
            read=read[:pos+3]+'G'+read[pos+4:]
            bPos=pos-1
            while bPos>=0:
                try:
                    quals+=bQ[bPos]()
                    success=True
                    break
                except:
                    bPos-1
                if success==False:
                    quals+=chr(2+qual)
        elif val>a:
            read=read[:pos+3]+'T'+read[pos+4:] 
            bPos=pos-1
            while bPos>=0:
                try:
                    quals+=bQ[bPos]()
                    success=True
                    break
                except:
                    bPos-1
                if success==False:
                    quals+=chr(2+qual)
        else:
            read=read[:pos+3]+'A'+read[pos+4:]
            bPos=pos-1
            while bPos>=0:
                try:
                    quals+=bQ[bPos]()
                    success=True
                    break
                except:
                    bPos-1
                if success==False:
                    quals+=chr(2+qual)
        if index in delD:
            delete=delD[index]()
            read=read[:pos+4]+read[pos+delete+4:]
        if index in insD:
            insert=insD[index]()
            read=read[:pos+4]+insert+read[pos+4:]
            for i in insert:
                iPos=pos-1
                while iPos>=0:
                    try:
                        quals+=iQ[iPos]()
                        success=True
                        break
                    except:
                        iPos-=1
                    if success==False:
                        quals+=chr(2+qual)
            pos+=len(insert)
        pos+=1
    quals+=quals[-1]
    read=read[4:readLen+4]
    quals=quals[:readLen]
    if len(quals)!=len(read):
        sys.exit()
    return read,quals      

def Comp(sequence):
    """ complements a sequence, preserving case."""
    d={'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c','N':'N','n':'n'}
    cSeq=''
    for s in sequence:
        cSeq+=d[s]
    return cSeq

def gParse(gens):
    try:
       gens=open(gens)
    except:
       rdlog.error('Cannot find file '+gens+'. Please check pathname')
       sys.exit('Cannot find file '+gens+'. Please check pathname')
    gL=[]
    for g in gens:
        gD={}
        g=g.rstrip()
        parts=g.split('\t')
        freq=float(parts[0])
        parts=parts[1:]
        for i in range(0,len(parts),3):
            chr=parts[i]
            pos=int(parts[i+1])
            var=parts[i+2]
            if chr in gD:
                gD[chr][pos]=var
            else:
                gD[chr]={pos:var}
        gL.append((freq,gD))
    gens.close()
    return gL

def usage():
    print '\n\n########################################################################'
    print '# GemSIM - Generic Error Model based SIMulator of N.G. sequencing data #'
    print '########################################################################\n'
    print '\nGemReads.py:\n'
    print 'Takes a reference genome, an empirical error model, and a haplotype file'
    print 'listing SNP locations and frequencies, and creates a simulated data set'
    print 'of random reads, as would be produced by a next-gen sequencing run.'
    print 'Output is in fastq format, suitable for input into popular alignment'
    print 'software.' 
    print '\nOptions:'
    print '      -h prints these instructions.'
    print '      -r reference genome, in fasta format.' 
    print '      -R Only for metagenome projects. Directory containing references.' 
    print '      -a Only for metagenome projects. Species-abundance file.'
    print '      -n number of reads to produce. For paired end reads, number of pairs.'
    print '      -g haplotype file, specifying location and frequency of snps.'
    print '      -G directory of haplotype files (metagenomics mode only).'
    print '      -l length of reads. Integer value, or -l d for empirical distribution.'
    print '      -m error model file *_single.gzip or *_paired.gzip.'
    print '      -c use this flag if you wish to draw reads from a circular genome.'
    print '      -q quality score offset. Usually 33 or 64 (see manual).' 
    print '      -o output file name prefix.'
    print '      -u Mean fragment length for paired end reads. -u d for empirical.'
    print '      -s standard deviation for fragment length. Use only with -u and -p.' 
    print '      -p use only to create paired end reads.\n\n'

def main(argv):
    refFile=''
    direct=''
    abund=''
    number=''
    length=''
    gens=''
    genDirect=''
    models=''
    circular=False
    qual=''
    out=''
    paired=False
    meta=False
    mean=''
    stdv=''
    try:
        opts, args = getopt.getopt(argv, "hr:R:a:n:g:G:l:m:ce:q:o:u:s:p")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '-r':
            refFile=arg
        elif opt == '-R':
            direct=arg
            meta=True
        elif opt =='-G':
            genDirect=arg
        elif opt == '-a':
            abund=arg
        elif opt == '-n':
            number = int(arg)
        elif opt == '-l':
            if arg=='d':
                length=arg
            else:
                length=int(arg)
        elif opt == '-u':
            if arg=='d':
                mean='emp'
            else:
                mean=int(arg) 
        elif opt == '-s':
            stdv=int(arg)
        elif opt == '-g':
            gens=arg
        elif opt =='-m':
            models=arg
        elif opt =='-c':
            circular=True
        elif opt =='-q':
            qual=int(arg)
        elif opt =='-o':
            out=arg
        elif opt =='-p':
            paired=True
    if number=='' or models=='' or length=='' or qual=='' or out=='':
        usage()
        sys.exit()
    if refFile=='' and direct=='':
        usage()
        sys.exit()
    if paired:
        if mean=='':
            print '\nPlease specify a fragment length\n'
            usage()
            sys.exit()
        rdlog.info('Generating paired end reads.')
        out1=open(out+'_fir.fastq','w')
    	out2=open(out+'_sec.fastq','w')
        rdlog.debug('Parsing model file.')
        mx1,mx2,insD1,insD2,delD1,delD2,intervals,gQualL,bQualL,iQualL,mates,rds,rdLenD=parseModel(models,paired,length)
        rdlog.debug('Model file parsed.')
    else:
        out1=open(out+'_single.fastq','w')
        rdlog.info('Generating single reads.')
        rdlog.debug('Parsing model file.')
        mx1,insD1,delD1,gQualL,bQualL,iQualL,readCount,rdLenD=parseModel(models,paired,length)
        rdlog.debug('Model file parsed.')
        #inserts
        insDict=mkInserts(mx1,insD1)
        #deletions
        delDict=mkDels(mx1,delD1)
    if paired:
        #choose insert length
        m0=float(mates[0])
        m1=float(mates[1])
        rd0=float(rds[0])
        rd1=float(rds[1])
        unAlign0=(m0*rd1-m1*m0)/(rd0*rd1-m1*m0)
        unAlign1=1.0-(unAlign0/(m0/rd0))
        keys=intervals.keys()
        keys.sort()
        if mean=='emp':
            inters=[]
            for k in keys:
                inters.append((k,intervals[k]))
            interval=bisect_choiceTUP(inters)
        #inserts1and2
        insDict1=mkInserts(mx1,insD1)
        insDict2=mkInserts(mx2,insD2)
        #deletions1and2
        delDict1=mkDels(mx1,delD1)
        delDict2=mkDels(mx2,delD2)
    #choose good quality bases
    gQList=[]             
    for i in (gQualL):
        gL=[]
        keys=i.keys()
        keys.sort()
        for k in keys:
            gL.append((chr(k+qual),i[k]))
        gQList.append(bisect_choiceTUP(gL))
    #choose bad quality bases
    bQList=[]
    for i in (bQualL):
        bL=[]
        keys=i.keys()
        keys.sort()
        for k in keys:
            bL.append((chr(k+qual),i[k]))
        bQList.append(bisect_choiceTUP(bL))
    #choose qualities for inserts
    iQList=[]
    for i in (iQualL):
        iL=[] 
        keys=i.keys()
        keys.sort()
        for k in keys:
            iL.append((chr(k+qual),i[k]))
        iQList.append(bisect_choiceTUP(iL))
    #choose read length
    if length=='d':
        rdlog.info('Using empirical read length distribution')
        lgth=[]
        keys=rdLenD.keys()
        keys.sort()
        for k in keys:
            lgth.append((k,rdLenD[k]))
        length=bisect_choiceTUP(lgth)
    else:
        length=ln(length)
    #choose reference
    genDict={}
    if meta:
        refDict,genDict,comDict,refList=getMet(direct,genDirect, abund)
    else:
        if gens!='':
            gens=gParse(gens)
            gens=bisect_gens(gens)
        else:
            gens=genRef('')
        refDict,comDict,refList=getRef(refFile)
        genDict[refFile]=gens
    reference=bisect_choiceTUP(refList)
    count=1
    #track read origin
    readOri={}
    for r,c in refList:
        readOri[r]=0    
    while count<=number:
        hd=reference()
        ref,refFile=refDict[hd]
        cRef=comDict[hd] 
        readOri[hd]+=1
        refLen=len(ref)
        if not paired:
            readLen=length()
            read1,pos,dir,quals1=readGen1(ref,cRef,refLen,readLen,genDict[refFile](),readLen,mx1,insDict,delDict,gQList,bQList,iQList,qual,circular,hd)
            head1='@'+'r'+str(count)+'_from_'+hd+'_#0/1\n'
        else:
            val=random.random()
            ln1=length()
            ln2=length()
            if mean=='emp':
                inter=interval()
            else:
                inter=int(random.normalvariate(mean,stdv))
            while (inter+ln1+ln2) > refLen*.9 or inter<ln1*1.5:
                if mean=='emp':
                    inter=interval()
                else:
                    inter=int(random.normalvariate(mean,stdv))
            if val > unAlign0+unAlign1:        
                read1,pos,dir,quals1=readGen1(ref,cRef,refLen,ln1,genDict[refFile](),inter,mx1,insDict1,delDict1,gQList,bQList,iQList,qual,circular,hd)
                read2,quals2=readGen2(ref,cRef,pos, dir, ln2, genDict[refFile](),inter,mx2,insDict2,delDict2,gQList,bQList,iQList,qual,circular,hd)
                p1=pos
                p2=pos+inter-ln2+1
            elif val > unAlign1:
                read1,pos,dir,quals1=readGen1(ref,cRef,refLen,ln1,genDict[refFile](),inter,mx1,insDict1,delDict1,gQList,bQList,iQList,qual,circular,hd)
                read2='N'*ln2
                quals2=chr(0+qual)*ln2
                p1=pos
                p2='*'
            else:
                read1,pos,dir,quals1=readGen1(ref,cRef,refLen,ln1,genDict[refFile](),inter,mx1,insDict1,delDict1,gQList,bQList,iQList,qual,circular,hd)
                read2,quals2=readGen2(ref,cRef,pos, dir, ln2,genDict[refFile](),inter,mx2,insDict2,delDict2,gQList,bQList,iQList,qual,circular,hd)
                read1='N'*ln1
                quals1=chr(0+qual)*ln1
                p1='*'
                p2=pos+inter-ln2+1
            head1='@'+'r'+str(count)+'_from_'+hd+'_ln'+str(inter)+'_#0/1\n'
            head2='@'+'r'+str(count)+'_from_'+hd+'_ln'+str(inter)+'_#0/2\n'
        out1.write(head1)
        out1.write(read1+'\n')
        out1.write('+\n')
        out1.write(quals1+'\n')
        if paired:
            out2.write(head2)
            out2.write(read2+'\n')
            out2.write('+\n')
            out2.write(quals2+'\n')
        count+=1
        if count%5000==0:
            if paired:
                rdlog.info('...simulated '+str(count)+' read pairs.')
            else:
                rdlog.info('...simulated '+str(count)+' single reads.')
    key=readOri.keys()
    key.sort()
    for k in key:
        rdlog.info('Simulated '+str(readOri[k])+' reads from reference '+k)

if __name__=="__main__":
    main(sys.argv[1:])
