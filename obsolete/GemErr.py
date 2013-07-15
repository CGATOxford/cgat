#!/usr/bin/python

# Copyright (C) 2011, Kerensa McElroy.
# kerensa@unsw.edu.au

# This file is part of the sequence simulator GemSIM. 
# It is used to calculate a platform- and run- specific
# error model for generating realistic sequencing reads.
# Alternatively, users may employ one of the precomputed
# error models distributed as part of the GemSIM package.

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
import cPickle
import gzip
import logging
import logging.handlers
import numpy as np

# Make a global logging object.
errlog=logging.getLogger("ErrLog")

# Set logging level, and write everything to a file
errlog.setLevel(logging.DEBUG)
LOG_FILENAME='./err.log'
h=logging.FileHandler(LOG_FILENAME,'w')
f=logging.Formatter("%(levelname)s %(asctime)s %(funcName)s %(lineno)d %(message)s")
h.setFormatter(f)
errlog.addHandler(h)


def rComp(sequence):
    """Reverse complements a sequence, preserving case."""
    d={'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c','N':'N','n':'n'}
    cSeq=''
    for s in sequence:
        cSeq+=d[s]
    cSeq=cSeq[::-1]
    return cSeq

def getRef(refFile):
    """Returns a genome reference."""
    refDict={}
    hdList=[]
    ref=''
    num=0
    try:
        f=open(refFile)
    except IOError:
        errlog.error('Cannot find reference file ' +refFile+'. Please check pathname.')
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
            refDict[head]=ref
            hdList.append(head)
            head=i[1:51].rstrip()
            i=f.readline()
            ref=''
    ref=ref.upper()
    for l in 'RYMKSWHBVD':
        ref=ref.replace(l,'N')
    refDict[head]=ref
    errlog.debug('Reference file successfully parsed.')
    return refDict


def parseFasta(file):
    """Returns sequence string from FASTA format."""
    f=open(file)
    ref=''
    for i in f:
        if i[0]!='>':
            ref+=i.rstrip()
    for l in 'RYLMKSWHBVD':
        ref=ref.replace(l,'N')
    return ref

def flip(refSlice,seq,qual,cigar):
    """Reverse complements a read."""
    flipSeq=''
    flipRS=''
    flipCig=[]
    comp={'A':'T','C':'G','G':'C','T':'A','N':'N','a':'t','t':'a','g':'c','c':'g','n':'n'}
    for i in seq:
        flipSeq+=comp[i]
    for i in refSlice:
        flipRS+=comp[i] 
    for i in range(0,len(cigar),2):
        flipCig.append(cigar[i+1])
        flipCig.append(cigar[i])
    flipCig.reverse()
    return flipRS[::-1],flipSeq[::-1],qual[::-1],flipCig

def parseMD(md):
    """Separates a cigar field into list of integers and character strings."""
    mdList=[]
    mdL=[]
    str=''
    if md[0].isdigit():
        before=True
    else:
        before=False
    for i in md:
        if i.isdigit():
            if before==True:
                str+=i
            else:
                mdList.append(str)
                str=i
                before=True
        else:
            if before==False:
                str+=i
            else:
                mdList.append(str)
                str=i
                before=False
    mdList.append(str)
    for i in mdList:
        if i.isdigit():
            mdL.append(int(i))
        else:
            mdL.append(i)
    return mdL


def updateM(ref,pos,seq,qual,cig,circ,mxNum,maxIndel,dir,readLen,excl):
    """Updates model with mutations, insertions, deletions in read."""
    swap={'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c','N':'N','n':'n'}
    inds={'A':0,'T':1,'G':2,'C':3,'N':4,'a':0, 't':1, 'g':2, 'c':3, 'n':4}
    RposL=pos-1                     #tracks leftmost pos of align against ref
    RposR=pos-1                     #tracks rightmost pos of align against ref
    Rpos=0                          #position within aligned reference slice
    Spos=0                          #sequence position within read
    bPos=pos-1                      #individual reference base position. 1 index
    refSlice=''
    cigH=cig
    if cig[1]=='H':
        cigH=cigH[2:]
    if cig[-1]=='H':
        cigH=cigH[:-2]
    if cigH[1]=='S':
        RposL-=cigH[0]
    for i in range(0,len(cig),2):   #slice alignment out of ref, excluding masked sections.
        if cig[i+1]=='M':
            RposR+=cig[i]
        elif cig[i+1]=='=':
            RposR+=cig[i]
        elif cig[i+1]=='X':
            RposR+=cig[i]
        elif cig[i+1]=='D':
            RposR+=cig[i]
        elif cig[i+1]=='N':
            refSlice+=ref[RposL:RposR] #cut before masked section.
            RposR+=cig[i]
            RposL=RposR 
    if cigH[-1]=='S':
        RposR+=cigH[-2]
    refSlice+=ref[RposL:RposR]
    refLen=len(ref)
    if dir=='f':
        if RposR<refLen:
            refSlice+=ref[RposR]     #+1 to allow consideration of base AFTER last read base.
        else:
            if circ:
                refSlice+=ref[0]         #+1 for reads ending at last reference base (circular).
            else:
                refSlice+='N'            #+1 for reads ending at last reference base (linear)
    elif dir=='r':
        if pos-2>0:
            refSlice=ref[pos-2]+refSlice
        else:
            if circ:
                refSlice=ref[-1]+refSlice
            else:
                refSlice='N'+refSlice
        refSlice,seq,qual,cig=flip(refSlice,seq,qual,cig)
        bPos=refLen-bPos-len(refSlice)  #so when we increment bpos it does in the right direction
    seq=seq[:readLen]           	#make sure fits in matrix
    seq=seq.upper()
    qual=qual[:readLen]
    d0=0
    d1=inds['N']
    d2=inds['N']
    d3=inds['N']
    d4=inds['N']
    d5=inds[refSlice[0]]
    d6=5                             #index for totals
    if cig[1]!='H':
        matrix[mxNum][d0][d1][d2][d3][d4][d5][d6]+=1
    for i in range(0,len(cig),2):
        if cig[i+1]=='H':
            seq=seq[:Spos]+'N'*cig[i]+seq[Spos:]
            Spos+=cig[i]            
        elif cig[i+1]=='M' or cig[i+1]=='S' or cig[i+1]=='X' or cig[i+1]=='=':
            matches=cig[i]
            count=0
            while count<matches:
                Spos+=1
                Rpos+=1
                bPos+=1
                count+=1
                refBase=refSlice[Rpos-1]
                mut=seq[Spos-1]
                after=refSlice[Rpos]
                qualIndex=ord(qual[Spos-1])-33
                if Spos>=4:
                    seq4=seq[Spos-4:Spos]
                else:
                    seq4='NNNN'+seq[:Spos]
                    seq4=seq4[-4:]
                d0=Spos
                d1=inds[refBase]
                d2=inds[seq4[2]]
                d3=inds[seq4[1]]
                d4=inds[seq4[0]]
                d5=inds[after]
                if mut!=refBase and refBase!='N':
                    snp=False
                    if dir=='f':
                        if str(bPos) in excl:
                            snp=True
                    else:
                        if (str(refLen-bPos)) in excl:
                            snp=True
                    if mut in 'ATCGatgcNn' and snp==False:
                        d6=inds[mut]
                        matrix[mxNum][d0][d1][d2][d3][d4][d5][d6]+=1
                        if qualIndex in bQualL[Spos-1]:
                            bQualL[Spos-1][qualIndex]+=1
                        else:
                            bQualL[Spos-1][qualIndex]=1
                    else:
                        if qualIndex in gQualL[Spos-1]:
                            gQualL[Spos-1][qualIndex]+=1
                        else:
                            gQualL[Spos-1][qualIndex]=1 
                else:
                    if qualIndex in gQualL[Spos-1]:
                        gQualL[Spos-1][qualIndex]+=1
                    else:
                        gQualL[Spos-1][qualIndex]=1
                matrix[mxNum][d0][d1][d2][d3][d4][d5][5]+=1
        elif cig[i+1]=='I':
            if cig[i]<=maxIndel:
                insert=seq[Spos:Spos+cig[i]]
                iQuals=qual[Spos:Spos+cig[i]]
                inDel=False
                if inDel==False: 
                    key=str(d0)+'.'+str(d1)+'.'+str(d2)+'.'+str(d3)+'.'+str(d4)+'.'+str(d5)
                    if key in insD[mxNum]:
                        if insert in insD[mxNum][key]:
                            insD[mxNum][key][insert]+=1
                        else:
                            insD[mxNum][key][insert]=1
                    else:
                        insD[mxNum][key]={insert:1}
                    for q in iQuals:
                        qualIndex=ord(q)-33
                        if qualIndex in iQualL[Spos]:
                            iQualL[Spos][qualIndex]+=1
                        else:
                            iQualL[Spos][qualIndex]=1
            Spos+=cig[i]
        elif cig[i+1]=='D':
            if cig[i]<=maxIndel:
                inDel=False
                if inDel==False: 
                    delete=cig[i]-1                                        #because of 0 index
                    key=str(d0)+'.'+str(d1)+'.'+str(d2)+'.'+str(d3)+'.'+str(d4)+'.'+str(d5)
                    if key in delD[mxNum]:
                        delD[mxNum][key][delete]+=1 
                    else:
			delD[mxNum][key]=[0]*maxIndel
                        delD[mxNum][key][delete]+=1
            Rpos+=cig[i]
            bPos+=cig[i]
        elif cig[i+1]=='N':
            bPos+=cig[i]

def kMers(reference,readLen,mxNum,maxIndel,minK):
    nucs= ['A','T','C','G']
    inds={'A':0,'T':1,'G':2,'C':3,'N':4,'a':0,'t':1,'g':2,'c':3,'n':4}
    kmers=[]
    reduce={}
    chrs=reference.values()
    ref=''
    for ch in chrs:
        ref+=ch+'.'
    for i in nucs:
        w1=i
        for j in nucs:
            w2=w1+j
            for k in nucs:
                w3=w2+k
                for l in nucs:
                    w4=w3+l
                    for m in nucs:
                        w5=w4+m
                        kmers.append(w5)
    for k in kmers:
        c=ref.count(k)
        cc=ref.count(rComp(k))
        if (c+cc)<minK:
            d1=inds[k[3]]
            d2=inds[k[2]]
            d3=inds[k[1]]
            d4=inds[k[0]]
            d5=inds[k[4]]
            new=k[1:]
            ikeys=insD[mxNum].keys()
            dkeys=delD[mxNum].keys()
            c=ref.count(new)
            cc=ref.count(rComp(new))
            if (c+cc)<minK:
                new=k[2:]
                c=ref.count(new)
                cc=ref.count(rComp(new))
                if (c+cc)<minK:
                    new=k[2:-1]
                    c=ref.count(new)
                    cc=ref.count(rComp(new))
                    if (c+cc)<minK:
                        for i in range(readLen+1):
                            tot=np.apply_over_axes(np.sum, matrix[mxNum][i], [1,2,3,4])[d1][0][0][0][0][5]
                            matrix[mxNum][i][d1][d2][d3][d4][d5][5]=tot
                            inserts={}
                            key=str(i)+'.'+str(d1)+'.'+str(d2)+'.'+str(d3)+'.'+str(d4)+'.'+str(d5)
                            for ik in ikeys:
                                if k[:-4]==''.join(ik.split('.')[:2]):
                                    ins=insD[mxNum][ik]
                                    ks=ins.keys()
                                    for s in ks:
                                        if s in inserts:
                                            inserts[s]+=ins[s]
                                        else:
                                            inserts[s]=ins[s]
                            if inserts!={}:
                                insD[mxNum][key]=inserts
                            dels=[0]*maxIndel
                            for dk in dkeys:
                                if k[:-4]==''.join(dk.split('.')[:2]):
                                    dele=delD[mxNum][dk]
                                    for e,d in enumerate(dele):
                                        dels[e]+=d
                            if dels!=[0]*maxIndel:
                                delD[mxNum][key]=dels
                            for n in range(5):
                                val=np.apply_over_axes(np.sum, matrix[mxNum][i], [1,2,3,4])[d1][0][0][0][0][n]
                                matrix[mxNum][i][d1][d2][d3][d4][d5][n]=val
                    else:
                        for i in range(readLen+1):
                            tot=np.apply_over_axes(np.sum, matrix[mxNum][i], [2,3,4])[d1][d2][0][0][0][5]
                            matrix[mxNum][i][d1][d2][d3][d4][d5][5]=tot
                            inserts={}
                            key=str(i)+'.'+str(d1)+'.'+str(d2)+'.'+str(d3)+'.'+str(d4)+'.'+str(d5)
                            for ik in ikeys:
                                if k[:-3]==''.join(ik.split('.')[:3]):
                                    ins=insD[mxNum][ik]
                                    ks=ins.keys()
                                    for s in ks:
                                        if s in inserts:
                                            inserts[s]+=ins[s]
                                        else:
                                            inserts[s]=ins[s]
                            if inserts!={}:
                                insD[mxNum][key]=inserts
                            dels=[0]*maxIndel
                            for dk in dkeys:
                                if k[:-3]==''.join(dk.split('.')[:3]):
                                    dele=delD[mxNum][dk]
                                    for e,d in enumerate(dele):
                                        dels[e]+=d
                            if dels!=[0]*maxIndel:
                                delD[mxNum][key]=dels
                            for n in range(5):
                                val=np.apply_over_axes(np.sum, matrix[mxNum][i], [2,3,4])[d1][d2][0][0][0][n]
                                matrix[mxNum][i][d1][d2][d3][d4][d5][n]=val
                else:
                    for i in range(readLen+1):
                        tot=np.apply_over_axes(np.sum, matrix[mxNum][i], [2,3])[d1][d2][0][0][d5][5]
                        matrix[mxNum][i][d1][d2][d3][d4][d5][5]=tot
                        inserts={}
                        key=str(i)+'.'+str(d1)+'.'+str(d2)+'.'+str(d3)+'.'+str(d4)+'.'+str(d5)
                        for ik in ikeys:
                           if k[:-3]==''.join(ik.split('.')[:3]) and k[-1]==ik[-1]:
                               ins=insD[mxNum][ik]
                               ks=ins.keys()
                               for s in ks:
                                   if s in inserts:
                                       inserts[s]+=ins[s]
                                   else:
                                       inserts[s]=ins[s]
                        if inserts!={}:
                            insD[mxNum][key]=inserts
                        dels=[0]*maxIndel
                        for dk in dkeys:
                            if k[:-3]==''.join(dk.split('.')[:3]) and k[-1]==dk[-1]: 
                                dele=delD[mxNum][dk]
                                for e,d in enumerate(dele):
                                    dels[e]+=d
                        if dels!=[0]*maxIndel:
                            delD[mxNum][key]=dels
                        for n in range(5):
                            val=np.apply_over_axes(np.sum, matrix[mxNum][i], [2,3])[d1][d2][0][0][d5][n]
                            matrix[mxNum][i][d1][d2][d3][d4][d5][n]=val 
            else:
                for i in range(readLen+1):
                    tot=np.apply_over_axes(np.sum, matrix[mxNum][i], [3])[d1][d2][d3][0][d5][5]
                    matrix[mxNum][i][d1][d2][d3][d4][d5][5]=tot
                    inserts={}
                    key=str(i)+'.'+str(d1)+'.'+str(d2)+'.'+str(d3)+'.'+str(d4)+'.'+str(d5)
                    for ik in ikeys:
                        if k[:-2]==''.join(ik.split('.')[:4]) and k[-1]==ik[-1]:
                            ins=insD[mxNum][ik]
                            ks=ins.keys()
                            for s in ks:
                                if s in inserts:
                                    inserts[s]+=ins[s]
                                else:
                                    inserts[s]=ins[s] 
                    if inserts!={}:
                        insD[mxNum][key]=inserts
                    dels=[0]*maxIndel
                    for dk in dkeys:
                        if k[:-2]==''.join(dk.split('.')[:4]) and k[-1]==dk[-1]:
                            dele=delD[mxNum][dk]
                            for e,d in enumerate(dele):
                                dels[e]+=d
                    if dels!=[0]*maxIndel:
                        delD[mxNum][key]=dels
                    for n in range(5):
                        val=np.apply_over_axes(np.sum, matrix[mxNum][i], [3])[d1][d2][d3][0][d5][n]
                        matrix[mxNum][i][d1][d2][d3][d4][d5][n]=val

def lowCov(readLen,mxNum):
    inds=[0,1,2,3,4]
    tot=np.apply_over_axes(np.sum,matrix[mxNum],[0,1,2,3,4,5])[0][0][0][0][0][0][5]
    muts=np.sum(matrix[mxNum])-tot
    avg=float(muts)/float(tot)
    for i in range(readLen+1):
        for j in inds:
            for k in inds:
                for l in inds:
                    for m in inds:
                        for n in inds:
                            if matrix[mxNum][i][j][k][l][m][n][5]<20:
                                for o in inds:
                                    matrix[mxNum][i][j][k][l][m][n][o]=int(avg*matrix[mxNum][i][j][k][l][m][n][5])

def mkMxSingle(readLen,ref,samFile,name,skip,circular,maxIndel,excl,minK):
    """Creates matrices of positional errors, insertions, deletions and bases in a sam file.""" 
    try:
        f=open(samFile)
    except:
        errlog.error('Cannot find samFile '+samFile+'. Please check pathname.')
        sys.exit('Cannot find samFile '+samFile+'. Please check pathname.')
    inds={'A':0,'T':1,'G':2,'C':3,'N':4,'a':0,'t':1,'g':2,'c':3,'n':4}
    global matrix, insD, delD, gQualL, bQualL, iQualL
    matrix=[np.zeros([readLen+1,5,5,5,5,5,6],dtype='int32'),None,None]
    insD=[{},None,None]
    delD=[{},None,None]
    gQualL=[]                              #tracks qualities for good bases
    for i in range(readLen):
        gQualL.append({})
    bQualL=[]                              #tracks qualities for bad (error) bases
    for i in range(readLen):
        bQualL.append({})
    iQualL=[]
    for i in range(readLen+1):
        iQualL.append({})                            #tracks average qualities for insertions
    readCount=0
    lineCount=0
    rdLenD={}
    line=f.readline()
    while line[0]=='@':
        line=f.readline()
    while line:
        lineCount+=1
        if skip==0 or lineCount%skip==0:       #take ith read
                parts=line.split('\t')
		flag=int(parts[1])
                if (flag & 0x04)==0:              #make sure read is aligned
                    #parse sam format
                    pos=int(parts[3])
                    seq=parts[9]
                    qual=parts[10]
                    cigar=parts[5]
                    cigList=parseMD(cigar)
                    chr=parts[2][:50]
                    #update read length dictionary
                    seqLen=len(seq)
                    if seqLen in rdLenD:
                        rdLenD[seqLen]+=1
                    else:
                        rdLenD[seqLen]=1
                    if flag & 0x10:          
                        #reverse complement
                        updateM(ref[chr],pos,seq,qual,cigList,circular,0,maxIndel,'r',readLen,excl)
                    else:
                        updateM(ref[chr],pos,seq,qual,cigList,circular,0,maxIndel,'f',readLen,excl)
                    if readCount%5000==0:
                        errlog.info('...parsed '+str(readCount)+' reads.')
                    readCount+=1
        line=f.readline()
    errlog.info('starting Kmers')
    if minK!=0:
        kMers(ref,readLen,0,maxIndel,minK)
    errlog.info('finished kmers')
    lowCov(readLen,0)
    errlog.debug('Finished parsing reads, writing matrices to files.')
    #write error models to files
    modelName=name+'_s.gzip'
    g=gzip.open(modelName,'wb')
    cPickle.dump(readLen,g)
    cPickle.dump(matrix[0], g)
    cPickle.dump(insD[0], g)
    cPickle.dump(delD[0], g)
    cPickle.dump(gQualL,g)
    cPickle.dump(bQualL,g)
    cPickle.dump(iQualL,g)
    cPickle.dump(readCount,g)
    cPickle.dump(rdLenD,g)
    g.close()
    errlog.info(str(lineCount)+' unpaired reads in total.')
    errlog.info('Parsed '+str(readCount)+'reads in total.')
    errlog.debug('Error models written to files.') 

def mkMxPaired(readLen,ref,samFile,name,skip,circular,maxIndel,excl,minK):
    """Creates matrices of positional errors, insertions, deletions and bases in a sam file."""
    try:
        f=open(samFile)
    except:
        errlog.error('Cannot find samFile '+samFile+'. Please check pathname.')
        sys.exit('Cannot find samFile '+samFile+'. Please check pathname.')
    inds={'A':0,'T':1,'G':2,'C':3,'N':4,'a':0,'t':1,'g':2,'c':3,'n':4}
    global matrix,insD,delD,gQualL,bQualL,iQualL
    matrix=[None,np.zeros([readLen+1,5,5,5,5,5,6],dtype='int32'),np.zeros([readLen+1,5,5,5,5,5,6],dtype='int32')]
    insD=[None,{},{}]
    delD=[None,{},{}]
    intD={}
    rds=[0,0]
    mates=[0,0]
    gQualL=[]                              #tracks qualities for good bases
    for i in range(readLen):
        gQualL.append({})
    bQualL=[]                              #tracks qualities for bad (error) bases
    for i in range(readLen):
         bQualL.append({})
    iQualL=[]                              #tracks average qualities for insertions
    for i in range(readLen+1):
         iQualL.append({})
    readCount=0
    rdLenD={}
    lineCount=0
    line=f.readline()
    lenD={}
    for i in ref.keys():
        lenD[i]=len(ref[i])
    while line[0]=='@':
        line=f.readline()
    while line:
        lineCount+=1 
        if skip==0 or lineCount%skip==0:                      #remove headers
                parts=line.split('\t')
                flag=int(parts[1])
                if (flag & 0x04)==0:              #make sure read is aligned
                    #parse sam format
                    pos=int(parts[3])
                    posMate=int(parts[7])
                    seq=parts[9]
                    qual=parts[10]
                    cigar=parts[5]
                    chr=parts[2][:50]
                    reflen=lenD[chr]
                    cigList=parseMD(cigar)
                    #update read length dictionary
                    if (readCount)%5000==0:
                        errlog.info('...parsed '+str(readCount)+' reads.')
                    seqLen=len(seq)
                    if seqLen in rdLenD:
                        rdLenD[seqLen]+=1
                    else:
                        rdLenD[seqLen]=1
                    insert=(int(parts[8]))
		    if insert < -reflen/2:
                        insert=reflen-pos+posMate+len(seq)
                    elif insert > reflen/2:
                        insert=reflen-posMate+pos+len(seq)
                    if (insert > 0):
                        if intD.has_key(insert):
                            intD[insert]+=1
                        else:
                            intD[insert]=1
                    if (flag & 0x40):      #reads unpaired or first in pair
                        if flag & 0x10:
                            updateM(ref[chr],pos,seq,qual,cigList,circular,1,maxIndel,'r',readLen,excl)
                        else:
                            updateM(ref[chr],pos,seq,qual,cigList,circular,1,maxIndel,'f',readLen,excl)
                        readCount+=1
                        rds[0]+=1
                        if (flag & 0x08): 
                            #track alignment of mates
                            mates[0]+=1
                    if (flag & 0x80):                     #matrices for 2nd read in pair
                        if flag & 0x10:
                            updateM(ref[chr],pos,seq,qual,cigList,circular,2,maxIndel,'r',readLen,excl)
                        else:
                            updateM(ref[chr], pos, seq,qual,cigList,circular,2,maxIndel,'f',readLen,excl)
                        readCount+=1
                        rds[1]+=1
                        if (flag & 0x08):
                            #track alignment of mates 
                            mates[1]+=1
        line=f.readline()
    mx=np.add.reduce(matrix[1],axis=0)
    if minK!=0:
        kMers(ref, readLen,1,maxIndel,minK)
        kMers(ref, readLen,2,maxIndel,minK)
    lowCov(readLen,1)
    lowCov(readLen,2)

    errlog.debug('Finished parsing reads, writing matrices to files.')
    #write error models to files
    modelName=name+'_p.gzip'
    g=gzip.open(modelName,'wb')
    cPickle.dump(readLen,g)
    cPickle.dump(matrix[1],g)
    cPickle.dump(matrix[2],g)
    cPickle.dump(insD[1],g)
    cPickle.dump(insD[2],g)
    cPickle.dump(delD[1],g)
    cPickle.dump(delD[2],g)
    cPickle.dump(intD,g)
    cPickle.dump(gQualL,g)
    cPickle.dump(bQualL,g)
    cPickle.dump(iQualL,g)
    cPickle.dump(mates,g)
    cPickle.dump(rds,g)
    cPickle.dump(rdLenD,g)
    g.close()
    errlog.info(str(lineCount)+' paired reads in total.')
    errlog.info('Parsed '+str(readCount)+' reads to create model.') 
    errlog.info(str(float(mates[0])/float(rds[0]))+'% first reads in pair with bad mates.')
    errlog.info(str(float(mates[1])/float(rds[1]))+'% second reads in pair with bad mates.')
                                                                                            
def usage():
    print '\n\n########################################################################'
    print '# GemSIM - Generic Error Model based SIMulator of N.G. sequencing data #'
    print '########################################################################\n'
    print '\nGemErr.py:\n'
    print 'Takes a sam file and catalogues all the mismatches, insertions, and deletions'
    print 'to create an error model for a particular sequencing run. Known true SNP'
    print 'positions may be excluded.' 
    print '\nOptions:'
    print '      -h prints these instructions.'
    print '      -r read length. Set to LONGEST read in dataset.'
    print '      -f reference genome in fasta format'
    print '      -s input file in sam format.'
    print '      -n desired output filename prefix.'
    print '      -c specifies reference genome is circular. Otherwise assumed linear.'
    print '      -i use only every ith read for model (optional, must be odd).' 
    print '      -m maximum indel size (optional, default=4).'
    print '      -p use only if your data contains paired end reads.'
    print '      -k minimum k-mer frequency in reference. (Default=0)'
    print '      -e comma separated list of reference positions to exclude e.g. 293, 342\n\n'

def main(argv):
    readLen=''
    samfile=''
    fasta=''
    circular=''
    name=''
    skip=0
    maxIndel=4
    excl=''
    paired=False
    circular=False
    minK=0
    try:
        opts, args = getopt.getopt(argv, "hr:f:s:n:ci:m:e:pk:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt =='-h':
            usage()
            sys.exit()
        elif opt =='-r':
            readLen=int(arg)
        elif opt =='-f':
            fasta=arg
        elif opt =='-s':
            samfile=arg
        elif opt =='-n':
            name=str(arg)
        elif opt =='-c':
            circular=True
        elif opt =='-i':
            skip=int(arg)
        elif opt =='-m':
            maxIndel=int(arg)
        elif opt =='-p':
            paired=True
        elif opt =='-e':
            excl=arg.split(',')
        elif opt =='-k':
            minK=int(arg)
    if readLen=='' or fasta=='' or samfile=='' or name=='':
        usage()
        sys.exit(2)
    reference=getRef(fasta)
    if skip!=0:
        if skip%2==0:
            usage()
            sys.exit(2)
    if circular:
        errlog.info('Treating reference genome as circular.')
    else:
        errlog.info('Treating reference genome as linear.')
    if paired:
        errlog.info('Treating reads as  paired.')
        mkMxPaired(readLen,reference,samfile,name,skip,circular,maxIndel,excl,minK)
    else:
        errlog.info('Treating reads as unpaired.')
        mkMxSingle(readLen,reference,samfile,name,skip,circular,maxIndel,excl,minK)
    
if __name__=="__main__":
    main(sys.argv[1:])
