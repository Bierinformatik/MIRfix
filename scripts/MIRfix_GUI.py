from __future__ import division
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
import os
import operator
import math
import sys
import re
from Bio.Align.Applications import ClustalwCommandline
from Tkinter import *
import Tkinter, tkFileDialog, Tkconstants
import Tkinter
import tkFileDialog
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import json
import gzip
import importlib
import multiprocessing, logging
import traceback as tb
from distutils.spawn import find_executable
sys.path.append('/scratch/fall/Ali_mirbase/lib/VRNA/249/lib64/python2.7/site-packages/')
import RNA
root=Tk()
#'/scratch/fall/Ali_mirbase/lib/VRNA/249/lib64/python2.7/site-packages/'

root.title("test")
root.geometry("500x400")

def getindex(sequence,specie,precID,precdesc,listnogenomes,listnotingenome,templong):#get the index of the original sequence in its genome
    try:
        specieitem=specie.split()
        listofgenomes=[]
        #lstgenomes = openfile(str(sys.argv[5]))
        lstgenomes = openfile(str(E1.get()))
        flaggenome=0#if there is a genomes for the specie, then it is 1
        flagseq=0#if the sequence found in its genome, then it is 1
        minusstrand=False
        for lines in lstgenomes:
            line=lines.strip()
            if specieitem[0].lower() in line.lower() and specieitem[1].lower() in line.lower():
                flaggenome=1
                listofgenomes.append(lines.strip())

        if flaggenome==0:
            listnogenomes.append(precID)
            returnlst=[]
            return returnlst,listnogenomes,listnotingenome,templong, minusstrand
        if len(listofgenomes)>0:
            returnlst=[]
            for gen in listofgenomes:
                filer=openfile(gen)
                fread = SeqIO.parse(filer,"fasta")
                for i in fread:
                    mixed=False
                    if "U" in str(i.seq).upper() and  "T" in str(i.seq).upper() :
                        mixed=True
                    precind = None
                    precind = str(i.seq).find(sequence)

                    if precind > 0:
                        printlog(["in genome",precID])
                        flagseq=1
                        gseq=str(i.seq)
                        cutlongbefore=250
                        cutlongafter=250
                        beforeseq=len(gseq[:precind])
                        afterseq=len(gseq[precind+len(sequence):])

                        if beforeseq<cutlongbefore:
                            cutlongbefore=beforeseq

                        if afterseq<cutlongafter:
                            cutlongafter=afterseq

                        longseq=gseq[precind-cutlongbefore:(precind+len(sequence)+cutlongafter)]
                        templong.append(precID.strip())
                        templong.append(str(longseq))
                        returnlst.append(precind)
                        returnlst.append(str(i.id))
                        returnlst.append(str(gen))
                        minusstrand=False #minus strand
                        return returnlst,listnogenomes,listnotingenome,templong,minusstrand #minus strand

                    else:       # We search for the reverse complement now
                        #precind = str(i.seq).find(str(Seq(sequence).reverse_complement()))
                        #if "U" in str(i.seq).upper() and  "T" in str(i.seq).upper() :
                        if not mixed:
                            precind =  str((i.seq).reverse_complement()).find(sequence) #minus strand
                            if precind > 0:
                                printlog(["in minus genome",precID])
                                flagseq=1
                                #gseq=str(i.seq)
                                gseq=str((i.seq).reverse_complement()) #minus strand
                                cutlongbefore=250
                                cutlongafter=250
                                beforeseq=len(gseq[:precind])
                                afterseq=len(gseq[precind+len(sequence):])

                                if beforeseq<cutlongbefore:
                                    cutlongbefore=beforeseq

                                if afterseq<cutlongafter:
                                    cutlongafter=afterseq

                                #longseq=str(Seq(gseq[precind-cutlongbefore:(precind+len(sequence)+cutlongafter)]).reverse_complement())  # we now search for the reverse complement and return this
                                longseq=str(gseq[precind-cutlongbefore:(precind+len(sequence)+cutlongafter)]) #minus strand
                                templong.append(precID.strip())
                                templong.append(str(longseq))
                                returnlst.append(precind)
                                returnlst.append(str(i.id))
                                returnlst.append(str(gen))
                                minusstrand=True #minus strand
                                return returnlst,listnogenomes,listnotingenome,templong,minusstrand #minus strand

        if flagseq==0 and flaggenome==1:
            listnotingenome.append(precID)
            returnlst=[]
            return returnlst,listnogenomes,listnotingenome,templong,minusstrand

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno
            root.destroy()
def getindex2mat(sequence,specie,precID,precdesc,listnogenomes,listnotingenome):#get the index of the original sequence in its genome
    try:
        specieitem=specie.split()
        listofgenomes=[]
        #lstgenomes = openfile(str(sys.argv[5]))
        lstgenomes = openfile(str(E1.get()))
        flaggenome=0#if there is a genomes for the specie, then it is 1
        flagseq=0#if the sequence found in its genome, then it is 1

        for lines in lstgenomes:
            line=lines.strip()
            if specieitem[0].lower() in line.lower() and specieitem[1].lower() in line.lower():
                flaggenome=1
                listofgenomes.append(lines.strip())

        if flaggenome==0:
            listnogenomes.append(precID)
            return "",listnogenomes,listnotingenome

        if len(listofgenomes)>0:
            printlog(["list of genomes>0:",sequence])
            for gen in listofgenomes:
                filer=openfile(gen)
                fread = SeqIO.parse(filer,"fasta")
                for i in fread:
                    mixed=False
                    if "U" in str(i.seq).upper() and  "T" in str(i.seq).upper() :
                        mixed=True
                    precind = None
                    precind = str(i.seq).find(sequence)

                    if precind > 0:
                        printlog(["in genome",precID])
                        flagseq=1
                        gseq=str(i.seq)
                        cutlongbefore=100
                        cutlongafter=100
                        beforeseq=len(gseq[:precind])
                        afterseq=len(gseq[precind+len(sequence):])

                        if beforeseq<cutlongbefore:
                            cutlongbefore=beforeseq

                        if afterseq<cutlongafter:
                            cutlongafter=afterseq

                        longseq=gseq[precind-cutlongbefore:(precind+len(sequence)+cutlongafter)]
                        return (str(longseq)),listnogenomes,listnotingenome
                    else:       # We search for the reverse complement now
                        #precind = str(i.seq).find(str(Seq(sequence).reverse_complement()))
                        if not mixed:
                            precind =  str((i.seq).reverse_complement()).find(sequence) #minus strand
                            if precind > 0:
                                printlog(["in minus genome",precID])
                                flagseq=1
                                #gseq=str(i.seq)
                                gseq=str((i.seq).reverse_complement()) #minus strand
                                cutlongbefore=100
                                cutlongafter=100
                                beforeseq=len(gseq[:precind])
                                afterseq=len(gseq[precind+len(sequence):])

                                if beforeseq<cutlongbefore:
                                    cutlongbefore=beforeseq

                                if afterseq<cutlongafter:
                                    cutlongafter=afterseq

                                #longseq=str(Seq(gseq[precind-cutlongbefore:(precind+len(sequence)+cutlongafter)]).reverse_complement())  # we now search for the reverse complement and return this
                                longseq=str(gseq[precind-cutlongbefore:(precind+len(sequence)+cutlongafter)]) #minus strand
                                return (str(longseq)),listnogenomes,listnotingenome

        if flagseq==0 and flaggenome==1:
            listnotingenome.append(precID)
            return "",listnogenomes,listnotingenome

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno
            
def flip(filename,filen,outdir,mappingfile,matfile,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,list2mat,listnogenomes,listnotingenome,templong,listgoodnew):#file name is the family name
    try:
        item=[]
        item1=[]
        precitem=[]
        returnlst=[]
        precfile = openfile(str(filen))
        printlog(["here prec", filen])
        for record in SeqIO.parse(precfile,'fasta'): #open the oriented file of a specific family
            precDes=record.description
            precseq=str(record.seq)
            precitem=precDes.split()
            specie=precitem[2].strip()+" "+precitem[3]
            precID=precitem[1].strip()
            flagprec=0
            flagmat=0

            if precID not in list2mat:
                mfi = openfile(mappingfile)
                for famline in mfi:#open the mapping file, to get the mature sequence ID, by the precursor ID in given family
                    famlinesplit=famline.split()
                    if precID in famline:# and not in 2mat list
                        printlog(["prec is here",precID])
                        flagprec=1
                        item=famline.split()
                        matID=item[4].strip()#mature sequence ID

                        mtf = openfile(matfile)
                        for mat in SeqIO.parse(mtf, "fasta"):#get the corresponding mature sequence and get position
                            if matID in str(mat.description):
                                matdesc=mat.description
                                flagmat=1
                                matseq=str(mat.seq).strip()
                                matseq=matseq.replace('T','U')#replace T by U because the RNAfold produces the sequences as U
                                spos=str(precseq).find(matseq)
#                                if not spos:  # again check minus strand
#                                    spos=str(precseq).find(reverse_complement(Seq(matseq)))
                                epos=spos+len(matseq)-1
                                break

                xcut=len(precseq[:spos])
                ycut=len(precseq[epos+1:])
                printlog(["cut", xcut, ycut])

                if xcut>ycut and ycut>10:#=> 3p cut from the end
                    precseq=precseq[:epos+11]
                elif xcut>ycut and ycut<=10:#=> 3p and no need to cut, already <=10
                    precseq=precseq
                elif xcut<ycut and xcut>10:#=> 5p cut at the top
                    precseq=precseq[spos-10:]
                elif xcut<ycut and xcut<=10:#=> 5p and no need to cut, already <=10
                    precseq=precseq

                spos=str(precseq).find(matseq)#spos after cut
                epos=spos+len(matseq)-1#epos after cut
                precseq=precseq.replace("U","T")#r
                returnlst,listnogenomes,listnotingenome,templong,minusstrand=getindex(precseq,specie,precID,precDes,listnogenomes,listnotingenome,templong)# returns 3 values received, the first is the index of the sequence, the ID where this sequence found in the genome and the genome filename

                if precID not in listnogenomes and precID not in listnotingenome:
                    listnewold=[]
                    precepos=returnlst[0]+len(precseq)# Get the end position of the prec by adding the length of the prec to the  index where the precursor found in the genome
                    x=0
                    m=0
                    sm=returnlst[0]+spos#start position of the mature sequence in the original precurson in the "GENOME"
                    em=sm+len(matseq)#end position of the mature sequence in the original precurson in the "GENOME"
                    x=sm-returnlst[0]#size/number of nucleotides before mature
                    y=precepos-em#size/number of nucleotides after mature
                    m=len(matseq)#mature sequence length
                    nx=sm-y#to replace X by Y
                    ny=em+x#to replace Y by X
                    newx=sm-nx #Calculate the new X, the new size/number of nucleotides before mature
                    newy=ny-em #Calculate the new Y, the new size/number of nucleotides after mature
                    newspos=newx #new start position of mature
                    newepos=newx+m-1 #new end position of mature
                    printlog(["mat= ",matdesc,matseq])

                    filer=openfile(returnlst[2]) #this loop to get the flipped sequence
                    fread = SeqIO.parse(filer,"fasta")

                    for i in fread:
                        if i.id==returnlst[1] and minusstrand is False:#
                            newseq=i.seq[nx:sm]+i.seq[sm:sm+m]+i.seq[em:ny]
                            newseq1 = newseq
                        elif i.id==returnlst[1] and minusstrand is True:#minus strand
                            revgenseq=(i.seq).reverse_complement() #reverse the genomic sequence
                            newseq=revgenseq[nx:sm]+revgenseq[sm:sm+m]+revgenseq[em:ny]
                            newseq1 = newseq


                    oldlstlstr=[]
                    oldlstlstl=[]
                    oldparts=0
                    finaloldcomp=[]
                    listnewold=dofold(listnewold,precDes,precseq,precDes,newseq)
                    printlog(["listnewold not in both: ",listnewold,filename,oldlstlstr,oldlstlstl,spos,epos,newspos,newepos,matdesc,matseq,outdir,oldparts,finaloldcomp,precDes,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew])
                    oldlstlstr,oldlstlstl,oldparts,finaloldcomp,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew=readfold(listnewold,filename,oldlstlstr,oldlstlstl,spos,epos,newspos,newepos,matdesc,matseq,outdir,oldparts,finaloldcomp,precDes,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew)#,newlstlstr,newlstlstl,oldlstlstr,oldlstlstl)#the orient sent to check the fold counts of the mature, the precseqold to print the full ID #The oldlstlstl oldlstlstr, are sent to save the list to the next step when checking the new sequence.

                elif precID in listnogenomes:
                    listnewold=[]
                    m=len(matseq)#mature sequence length
                    newspos=spos #new start position of mature
                    newepos=epos #new end position of mature
                    printlog(["mat= ",matdesc,matseq])
                    newseq=str(precseq)
                    newseq1=newseq
                    oldlstlstr=[]
                    oldlstlstl=[]
                    oldparts=0
                    finaloldcomp=[]
                    listnewold=dofold(listnewold,precDes,precseq,precDes,newseq)
                    printlog(["listnewold listnogenome: ",listnewold,filename,oldlstlstr,oldlstlstl,spos,epos,newspos,newepos,matdesc,matseq,outdir,oldparts,finaloldcomp,precDes,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew])
                    oldlstlstr,oldlstlstl,oldparts,finaloldcomp,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew=readfold(listnewold,filename,oldlstlstr,oldlstlstl,spos,epos,newspos,newepos,matdesc,matseq,outdir,oldparts,finaloldcomp,precDes,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew)

                elif precID in listnotingenome:
                    listnewold=[]
                    m=len(matseq)#mature sequence length
                    newspos=spos #new start position of mature
                    newepos=epos #new end position of mature
                    printlog(["mat= ",matdesc,matseq])
                    newseq=str(precseq)
                    newseq1=newseq
                    oldlstlstr=[]
                    oldlstlstl=[]
                    oldparts=0
                    finaloldcomp=[]
                    listnewold=dofold(listnewold,precDes,precseq,precDes,newseq)
                    printlog(["listnewold not in genome: ",listnewold,filename,oldlstlstr,oldlstlstl,spos,epos,newspos,newepos,matdesc,matseq,outdir,oldparts,finaloldcomp,precDes,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew])
                    oldlstlstr,oldlstlstl,oldparts,finaloldcomp,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew=readfold(listnewold,filename,oldlstlstr,oldlstlstl,spos,epos,newspos,newepos,matdesc,matseq,outdir,oldparts,finaloldcomp,precDes,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew)

        return listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listnogenomes,listnotingenome,templong,listgoodnew

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def dofold(listnewold,oldid,precseq,newid,newseq):
    try:
        md = RNA.md()
        md.dangles = 2 #int(sys.argv[10]) if sys.argv[10] else 3   JF: DO NOT CHANGE
        md.noLP=1
        fc = RNA.fold_compound(str(precseq), md)
        (oldstruct,oldscore) = fc.mfe()
        mdn = RNA.md()
        mdn.dangles = 2 #int(sys.argv[10]) if sys.argv[10] else 3   JF: DO NOT CHANGE
        mdn.noLP=1
        fcn = RNA.fold_compound(str(newseq), mdn)
        (newstruct,newscore)= fcn.mfe()
        printlog(["folding is here:",oldstruct,oldscore])
        printlog(["folding is here:",newstruct,newscore])
        listnewold.extend((oldid,str(precseq),oldstruct,oldscore,newid,str(newseq),newstruct,newscore))
        return listnewold

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def readfold(listnewold,filename,oldlstlstr,oldlstlstl,spos,epos,newspos,newepos,matdesc,matseq,outdir,oldparts,finaloldcomp,precDes,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew):
    try:
        printlog(["readfold",listnewold])
        for i in range(0,len(listnewold)):
            if i==0 or i==4:
                printlog(i)
                if i==0:
                    printlog("i old")
                    stat="old"
                    precid=listnewold[i]
                    oldprecseq=listnewold[i+1]
                    st=listnewold[i+2]
                    hairpin=listnewold[i+2]
                    oldscore=listnewold[i+3]
                    partslist=[]
                    parts=0
                    pos=0
                    realpos=0
                    left=[]
                    right=[]
                    lstlstr=[]#list of list for right ")"
                    lstlstl=[]#list of list for left "("
                    newlstlstr=[]#list of list for right ")" of the new sequence
                    newlstlstl=[]#list of list for left "(" of the new sequence
                    curr=0
                    nextl=0
                    cut=0
                    flag=0
                    countfold=0
                    newbroken=False
                    newloop=False
                    newparts=0
                    newncounts=0
                    warnlength=0

                else:
                    printlog("i new")
                    stat="new"
                    precid=listnewold[i]
                    newseq=listnewold[i+1]
                    newseq1=newseq
                    st=listnewold[i+2]
                    hairpin=listnewold[i+2]
                    newscore=listnewold[i+3]
                    partslist=[]
                    parts=0
                    pos=0
                    realpos=0
                    left=[]
                    right=[]
                    lstlstr=[]#list of list for right ")"
                    lstlstl=[]#list of list for left "("
                    newlstlstr=[]#list of list for right ")" of the new sequence
                    newlstlstl=[]#list of list for left "(" of the new sequence
                    curr=0
                    nextl=0
                    cut=0
                    flag=0
                    countfold=0
                    newbroken=False
                    newloop=False
                    newparts=0
                    newncounts=0
                    warnlength=0

                while pos<len(st): #to count the parts, loop in the string until we found another open "(" after a close ")", then means one more part.. in this case the pos assingned zero again, and the new st is the new part
                    parts=parts+1 #1 part: ..((.)).. #2parts: ..((.))..(((.)))..
                    tempright=[]  #positions of the ")"
                    templeft=[] ##positions of the "("
                    for i in st:
                        realpos=cut+pos #cut is where we find new "(" after closing ")",the real position is always the position in the original string (all string not when divided into new parts
                        if i=="(":
                            left.append(realpos)
                            curr=pos #current position, changes in each new string (new part, in case the main strign was more than one part)
                            templeft.append(realpos)
                        elif i==")":
                            right.append(realpos)
                            tempright.append(realpos)
                            nextl=pos

                        if nextl<curr and nextl>0: #here means that we have more one part, by finding "(" after ")"
                            endofcurrentpart=str(left[-1]-1)
                            startofnextpart=str(right[-1]+1)
                            left.pop()
                            templeft.pop()
                            realpos-=1
                            st=st[pos:]#make the next part as new string
                            partslist.append(int(endofcurrentpart))
                            partslist.append(int(startofnextpart))
                            cut=realpos # cut at the real position of thr whole original string
                            pos=0#position=zero for the new string (new part)
                            templst=right
                            templst=[]
                            nextl=0
                            break

                        pos=pos+1

                    if len(templeft)!=len(tempright):
                        warnlength=1

                    templst=right
                    lstlstr.append(tempright)
                    lstlstl.append(templeft)

                if (len(lstlstr) != len(lstlstl)) or ( (")" not in hairpin) and ('(' not in hairpin )) :#checking the parts
                    flag=1

                else:
                    for i in range(len(lstlstr)):
                        if len(lstlstr[i])!=len(lstlstl[i]):#check if the foldings not equal
                            flag=1

                if stat=="old":
                    x=spos
                    y=epos
                elif stat=="new":
                    x=newspos
                    y=newepos
                partslstlen=len(partslist)
                i=0
                #if the parts are overlapping, like: ..(((.(((((...)))))..((((..)))).)))...,,,in this case we take only the parts inside and ignore the outer foldings, in this section we correct the parts

                if warnlength==1:
                    printlog("here warn")
                    printlog(partslist)
                    #for i in range(0,partslstlen):
                    while i <=len(partslist)-1:
                        if i==0 :
                            warnleft=[]
                            warnright=[]
                            warnhairpin=hairpin[0:(int(partslist[0])+1)]
                            printlog(warnhairpin)
                            for k in range(0,len(warnhairpin)):
                                if hairpin[k]=="(":
                                    warnleft.append(k)
                                if hairpin[k]==")":
                                    warnright.append(k)
                            if len(warnleft)>len(warnright):
                                printlog("R")
                                newsindex=len(warnleft)-1-len(warnright)
                                partslist.insert(0,-1)
                                partslist.insert(1,warnleft[newsindex]+1)
                                if len(partslist)>4:
                                    i=4
                                elif len(partslist)<=4:
                                    i=3
                                printlog(partslist)
                                printlog(i)
                            if len(warnleft)<len(warnright):#which is impossible to happen
                                printlog("L")
                                diff=len(warnright)-len(warnleft)
                                neweindex=warnright[len(warnright)-diff]-1
                                printlog(neweindex)
                                partslist[0]=neweindex
                                if len(partslist)==2:
                                    i=1
                                elif len(partslist)>2:
                                    i=2
                            if len(warnleft)==len(warnright) and len(partslist)>2:
                                i=2
                            elif len(warnleft)==len(warnright) and len(partslist)==2:
                                i=1

                        if i==len(partslist)-1:
                            printlog("end")
                            printlog(partslist)
                            warnleft=[]
                            warnright=[]
                            warnhairpin=hairpin[int(partslist[i]):]
                            for k in range(partslist[i],len(hairpin)):

                                if hairpin[k]=="(":
                                    warnleft.append(k)
                                if hairpin[k]==")":
                                    warnright.append(k)

                            if len(warnleft)>len(warnright):
                                newsindex=len(warnleft)-1-len(warnright)
                                partslist[i]=warnleft[newsindex]+1
                                i=len(partslist)

                            if len(warnleft)<len(warnright):#which is impossible to happen
                                diff=len(warnright)-len(warnleft)
                                neweindex=warnright[len(warnright)-diff]-1
                                partslist.insert(i+1,-1)
                                partslist.insert(i+1,neweindex)
                                i=len(partslist)

                            if len(warnleft)==len(warnright):
                                i=len(partslist)

                        if i!=0 and i<len(partslist)-1 and i%2==0:
                            warnleft=[]
                            warnright=[]
                            printlog("mid")
                            printlog(partslist)
                            for k in range(int(partslist[i-1]),int(partslist[i])+1):
                                if hairpin[k]=="(":
                                    warnleft.append(k)
                                if hairpin[k]==")":
                                    warnright.append(k)

                            if len(warnleft)>len(warnright) and len(partslist)>4 and (int(len(partslist))-i)>2:
                                newsindex=len(warnleft)-1-len(warnright)
                                partslist[i-1]=warnleft[newsindex]+1
                                i=i+2
                            elif (len(warnleft)>len(warnright) and len(partslist)>4 and (int(len(partslist))-i)==2) or (len(warnleft)>len(warnright) and len(partslist)==4):
                                newsindex=len(warnleft)-1-len(warnright)
                                partslist[i-1]=warnleft[newsindex]+1
                                i=i+1

                            if len(warnleft)<len(warnright) and len(partslist)>4 and (int(len(partslist))-i)>2:#which is impossible to happen
                                diff=len(warnright)-len(warnleft)
                                neweindex=warnright[len(warnright)-diff]-1
                                partslist[i]=neweindex
                                i=i+2
                            elif (len(warnleft)<len(warnright) and len(partslist)>4 and (int(len(partslist))-i)==2) or (len(warnleft)<len(warnright) and len(partslist)==4):#which is impossible to happen
                                diff=len(warnright)-len(warnleft)
                                neweindex=warnright[len(warnright)-diff]-1
                                partslist[i]=neweindex
                                i=i+1

                            if len(warnleft)==len(warnright) and len(partslist)>4 and (int(len(partslist))-i)>2:
                                i=i+2
                            elif len(warnleft)==len(warnright) and len(partslist)>4 and (int(len(partslist))-i)==2:
                                i=i+1
                            elif len(warnleft)==len(warnright) and len(partslist)==4:
                                i=i+1
                        partslstlen=len(partslist)
                #at this section the the parts list is corrected ( if needed )
                if parts>1:
                    lstcompinside=[]
                    lstcompstart=[]
                    lstcompend=[]
                    finalcomp=[]
                    resbroken=None
                    resloop=None
                    respos=1000000000
                    rescounts=0

                    if stat=="old":
                        oldparts=parts
                    else:
                        newparts=parts
                    for i in range(0,len(partslist)):
                        if i==0 and partslist[0]!=-1:#if the partslist starting from position zero, i.e. the first part starts at position zero
                            loopstart=hairpin.rfind('(',0,int(partslist[0])+1)+1
                            loopend=hairpin.find(')',0,int(partslist[0])+1)-1
                            if (((x in range(0,loopstart)) and (y in range(0,loopstart))) or ((x  in range(loopend+1,int(partslist[0])+1)) and (y in range(loopend+1,int(partslist[0])+1)))) and stat=="old":
                                oldbroken=False
                                oldloop=False
                                oldncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                lstcompstart.append(oldbroken)
                                lstcompstart.append(oldloop)
                                lstcompstart.append(i)
                                lstcompstart.append(oldncounts)
                            elif ( ((x in range(0,loopstart)) or (x in range(loopend+1,int(partslist[0])+1)))  and  (y in range(loopstart,loopend+1))) and stat=="old":
                                oldloop=True
                                oldbroken=False
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                oldncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                lstcompstart.append(oldbroken)
                                lstcompstart.append(oldloop)
                                lstcompstart.append(i)
                                lstcompstart.append(oldncounts)
                            elif ( ((y in range(0,loopstart)) or (y in range(loopend+1,int(partslist[0])+1)))  and  (x in range(loopstart,loopend+1))) and stat=="old":
                                oldloop=True
                                oldbroken=False
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                oldncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                lstcompstart.append(oldbroken)
                                lstcompstart.append(oldloop)
                                lstcompstart.append(i)
                                lstcompstart.append(oldncounts)
                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))) and stat=="old":
                                oldloop=True
                                oldbroken=True
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                oldncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                            elif (((x in range (0,int(partslist[0])+1)) and (y not in range (0,int(partslist[0])+1))) or ((x not in range (0,int(partslist[0])+1)) and (y in range (0,int(partslist[0])+1))))and stat=="old" :#if x or y out of the  hairpin
                                oldbroken=True
                                oldncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                            elif ((x in range(0,loopstart)) and (y in range(loopend+1,int(partslist[0])+1))) or ((x in range(loopend+1,int(partslist[0])+1)) and ((y in range(0,loopstart)))) and stat=="old":
                                oldbroken=True
                                oldloop=True
                                oldncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))

                            elif (((x in range(0,loopstart)) and (y in range(0,loopstart))) or ((x  in range(loopend+1,int(partslist[0])+1)) and (y in range(loopend+1,int(partslist[0])+1)))) and stat=="new":
                                newbroken=False
                                newloop=False
                                newncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                lstcompstart.append(newbroken)
                                lstcompstart.append(newloop)
                                lstcompstart.append(i)
                                lstcompstart.append(newncounts)
                            elif ( ((x in range(0,loopstart)) or (x in range(loopend+1,int(partslist[0])+1)))  and  (y in range(loopstart,loopend+1))) and stat=="new":
                                newloop=True
                                newbroken=False
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                newncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                lstcompstart.append(newbroken)
                                lstcompstart.append(newloop)
                                lstcompstart.append(i)
                                lstcompstart.append(newncounts)
                            elif ( ((y in range(0,loopstart)) or (y in range(loopend+1,int(partslist[0])+1)))  and  (x in range(loopstart,loopend+1))) and stat=="new":
                                newloop=True
                                newbroken=False
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                newncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                lstcompstart.append(newbroken)
                                lstcompstart.append(newloop)
                                lstcompstart.append(i)
                                lstcompstart.append(newncounts)
                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))) and stat=="new":
                                newloop=True
                                newbroken=True
                                newncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                            elif (((x in range (0,int(partslist[0])+1)) and (y not in range (0,int(partslist[0])+1))) or ((x not in range (0,int(partslist[0])+1)) and (y in range (0,int(partslist[0])+1))))and stat=="new" :#if x or y out of the  hairpin
                                newbroken=True
                                newloop=True
                                newncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                            elif ((x in range(0,loopstart)) and (y in range(loopend+1,int(partslist[0])+1))) or ((x in range(loopend+1,int(partslist[0])+1)) and ((y in range(0,loopstart)))) and stat=="new":
                                newloop=True
                                newbroken=True
                                newncounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))
                                currcounts=int(hairpin[0:int(partslist[0])+1].count(")"))+int(hairpin[0:int(partslist[0])+1].count("("))


                        if i==len(partslist)-1 and partslist[i]!=-1:
                            loopstart=hairpin.rfind('(',int(partslist[i])+1)+1
                            loopend=hairpin.find(')',int(partslist[i])+1)-1

                            if (((x in range(int(partslist[i]),loopstart)) and (y in range(int(partslist[i]),loopstart) )) or ((x in range(loopend+1,len(hairpin))) and (y in range(loopend+1,len(hairpin))))) and stat=="old":
                                oldncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                oldbroken=False
                                oldloop=False
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                lstcompend.append(oldbroken)
                                lstcompend.append(oldloop)
                                lstcompend.append(i)
                                lstcompend.append(oldncounts)
                            elif ( ((x in range(int(partslist[i]),loopstart)) or (x in range(loopend+1,len(hairpin)))) and (y in range(loopstart,loopend+1)) ) and stat=="old":
                                oldbroken=False
                                oldncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                oldloop=True
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                lstcompend.append(oldbroken)
                                lstcompend.append(oldloop)
                                lstcompend.append(i)
                                lstcompend.append(oldncounts)
                            elif ( ((y in range(int(partslist[i]),loopstart)) or (y in range(loopend+1,len(hairpin)))) and (x in range(loopstart,loopend+1)) ) and stat=="old":
                                oldbroken=False
                                oldncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                oldloop=True
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                lstcompend.append(oldbroken)
                                lstcompend.append(oldloop)
                                lstcompend.append(i)
                                lstcompend.append(oldncounts)
                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))) and stat=="old":
                                oldbroken=True
                                oldncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                oldloop=True
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                            elif (((x not in range(int(partslist[i]),len(hairpin))) and (y in range(int(partslist[i]),len(hairpin)))) or ((y not in range(int(partslist[i]),len(hairpin))) and (x in range(int(partslist[i]),len(hairpin))))) and stat=="old":
                                oldncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                oldbroken=True
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                            elif (((x in range(int(partslist[i]),loopstart)) and (y in range(loopend+1,len(hairpin)))) or ((x in range(loopend+1,len(hairpin))) and (y in range(int(partslist[i]),loopstart)))) and stat=="old":
                                oldloop=True
                                oldbroken=True
                                oldncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                            elif (((x in range(int(partslist[i]),loopstart)) and (y in range(int(partslist[i]),loopstart) )) or ((x in range(loopend+1,len(hairpin))) and (y in range(loopend+1,len(hairpin))))) and stat=="new":
                                newbroken=False
                                newloop=False
                                newncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                currncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                lstcompend.append(newbroken)
                                lstcompend.append(newloop)
                                lstcompend.append(i)
                                lstcompend.append(newncounts)
                            elif ( ((x in range(int(partslist[i]),loopstart)) or (x in range(loopend+1,len(hairpin)))) and (y in range(loopstart,loopend+1)) ) and stat=="new":
                                newbroken=False
                                newloop=True
                                newncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                lstcompend.append(newbroken)
                                lstcompend.append(newloop)
                                lstcompend.append(i)
                                lstcompend.append(newncounts)
                            elif ( ((y in range(int(partslist[i]),loopstart)) or (y in range(loopend+1,len(hairpin)))) and (x in range(loopstart,loopend+1)) ) and stat=="new":
                                newbroken=False
                                newloop=True
                                newncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                lstcompend.append(newbroken)
                                lstcompend.append(newloop)
                                lstcompend.append(i)
                                lstcompend.append(newncounts)
                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))) and stat=="new":
                                newbroken=True
                                newloop=True
                                newncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                            elif (((x not in range(int(partslist[i]),len(hairpin))) and (y in range(int(partslist[i]),len(hairpin)))) or ((y not in range(int(partslist[i]),len(hairpin))) and (x in range(int(partslist[i]),len(hairpin))))) and stat=="new":
                                newbroken=True
                                newncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                            elif (((x in range(int(partslist[i]),loopstart)) and (y in range(loopend+1,len(hairpin)))) or ((x in range(loopend+1,len(hairpin))) and (y in range(int(partslist[i]),loopstart)))) and stat=="new":
                                newloop=True
                                newbroken=True
                                newncounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))
                                currcounts=int(hairpin[int(partslist[i]):].count("("))+int(hairpin[int(partslist[i]):].count(")"))

                        if i!=0 and i!=len(partslist)-1 and i%2==0:
                            loopstart=hairpin.rfind('(',int(partslist[i-1]),int(partslist[i])+1)+1
                            loopend=hairpin.find(')',int(partslist[i-1]),int(partslist[i])+1)-1
                            currflag=0

                            if (((x in range(int(partslist[i-1]),loopstart)) and (y in range(int(partslist[i-1]),loopstart))) or ((x in range(loopend+1,int(partslist[i])+1)) and (y in range(loopend+1,int(partslist[i])+1)))) and stat=="old":
                                oldncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                oldbroken=False
                                oldloop=False
                                currloop=oldloop
                                currbroken=oldbroken
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                            elif (((x in range(int(partslist[i-1]),loopstart)) or (x in range(loopend+1,int(partslist[i])+1))) and (y in range(loopstart,loopend+1))) and stat=="old":
                                oldbroken=False
                                oldncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                oldloop=True
                                currloop=oldloop
                                currbroken=oldbroken
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                            elif (((y in range(int(partslist[i-1]),loopstart)) or (y in range(loopend+1,int(partslist[i])+1))) and (x in range(loopstart,loopend+1))) and stat=="old":
                                oldbroken=False
                                oldncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                oldloop=True
                                currloop=oldloop
                                currbroken=oldbroken
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))) and stat=="old":
                                oldbroken=True
                                oldncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                oldloop=True
                                currloop=oldloop
                                currbroken=oldbroken
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                currflag=1
                            elif ((x not in range(int(partslist[i-1]),int(partslist[i])+1) and (y in range(int(partslist[i-1]),int(partslist[i])+1))) or ((y not in range(int(partslist[i-1]),int(partslist[i])+1)) and (x in range(int(partslist[i-1]),int(partslist[i])+1)))) and stat=="old":
                                oldbroken=True
                                oldloop=True
                                oldncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                currloop=oldloop
                                currbroken=oldbroken
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                            elif (((x in range(int(partslist[i-1]),loopstart)) and (y in range(loopend+1,int(partslist[i])+1))) or ((x in range(loopend+1,int(partslist[i])+1)) and (y in range(int(partslist[i-1]),loopstart)))) and stat=="old":
                                oldncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                oldloop=True
                                oldbroken=True
                                currloop=oldloop
                                currbroken=oldbroken
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                            elif (((x in range(int(partslist[i-1]),loopstart)) and (y in range(int(partslist[i-1]),loopstart))) or ((x in range(loopend+1,int(partslist[i])+1)) and (y in range(loopend+1,int(partslist[i])+1)))) and stat=="new":
                                newncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                newbroken=False
                                newloop=False
                                currbroken=newbroken
                                currloop=newloop
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                            elif (((x in range(int(partslist[i-1]),loopstart)) or (x in range(loopend+1,int(partslist[i])+1))) and (y in range(loopstart,loopend+1))) and stat=="new":
                                newbroken=False
                                newncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                newloop=True
                                currbroken=newbroken
                                currloop=newloop
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                            elif (((y in range(int(partslist[i-1]),loopstart)) or (y in range(loopend+1,int(partslist[i])+1))) and (x in range(loopstart,loopend+1))) and stat=="new":
                                newbroken=False
                                newncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                newloop=True
                                currbroken=newbroken
                                currloop=newloop
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))


                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))) and stat=="new":
                                newbroken=True
                                newncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                newloop=True
                                currbroken=newbroken
                                currloop=newloop
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                            elif ((x not in range(int(partslist[i-1]),int(partslist[i])+1) and (y in range(int(partslist[i-1]),int(partslist[i])+1))) or ((y not in range(int(partslist[i-1]),int(partslist[i])+1)) and (x in range(int(partslist[i-1]),int(partslist[i])+1)))) and stat=="new":
                                newncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                newbroken=True
                                newloop=True
                                currbroken=newbroken
                                currloop=newloop
                                currflag=1
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                            elif (((x in range(int(partslist[i-1]),loopstart)) and (y in range(loopend+1,int(partslist[i])+1))) or ((x in range(loopend+1,int(partslist[i])+1)) and (y in range(int(partslist[i-1]),loopstart)))) and stat=="new":
                                newncounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                currcounts=int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count("("))+int(hairpin[int(partslist[i-1]):int(partslist[i])+1].count(")"))
                                newloop=True
                                newbroken=True
                                currbroken=newbroken
                                currloop=newloop
                                currflag=1

                                if currflag==1 and i==2:
                                    prevbroken=currbroken
                                    prevloop=currbroken
                                    prevpos=i
                                    prevcounts=currcounts

                                    if currbroken==False:
                                        resbroken=currbroken
                                        resloop=currloop
                                        respos=i
                                        rescounts=currcounts
                                    else:
                                        resbroken=None
                                        resloop=None
                                        respos=1000000000
                                        rescounts=0

                                elif currflag==1 and i!=2 and currbroken==False and ((((currloop==False and resloop==False) or (currloop==True and resloop==True)) and currcounts>rescounts) or resbroken==None):
                                    resbroken=currbroken
                                    resloop=currloop
                                    respos=i
                                    rescounts=currcounts

                                if resbroken==False and respos not in lstcompinside:
                                    lstcompinside=[]
                                    lstcompinside.append(resbroken)
                                    lstcompinside.append(resloop)
                                    lstcompinside.append(respos)
                                    lstcompinside.append(rescounts)

                    #all if here and elif here, are not broken when len>0 then we compare the loop
                    if len(lstcompstart)>0 and len(lstcompinside)==0 and len(lstcompend)==0:#start only
                        finalcomp=lstcompstart
                        posi=finalcomp[2]
                        finalcomp.append(0)
                        finalcomp.append(int(partslist[0]))
                    elif len(lstcompstart)==0 and len(lstcompinside)>0 and len(lstcompend)==0:#middle only
                        finalcomp=lstcompinside
                        posi=finalcomp[2]#the index in partslist, that stores the real end position of a given sub-stem

                        finalcomp.append(int(partslist[posi-1]))
                        finalcomp.append(int(partslist[posi]))
                    elif len(lstcompstart)==0 and len(lstcompinside)==0 and len(lstcompend)>0:#end only
                        finalcomp=lstcompend
                        posi=finalcomp[2]
                        finalcomp.append(int(partslist[posi]))
                        finalcomp.append(len(hairpin)-1)
                    elif len(lstcompstart)>0 and len(lstcompinside)>0:#start + middle
                        if lstcompstart[1]==lstcompinside[1] and lstcompstart[3]>=lstcompinside[3]:#if both true or false loops we compare the number of foldings
                            finalcomp=lstcompstart
                            posi=finalcomp[2]
                            finalcomp.append(0)
                            finalcomp.append(int(partslist[0]))
                        elif lstcompstart[1]==lstcompinside[1] and lstcompstart[3]<lstcompinside[3]:
                            finalcomp=lstcompinside
                            posi=finalcomp[2]
                            finalcomp.append(int(partslist[posi-1]))
                            finalcomp.append(int(partslist[posi]))
                        elif lstcompstart[1]==True and lstcompinside==False:
                            finalcomp=lstcompinside
                            posi=finalcomp[2]
                            finalcomp.append(int(partslist[posi-1]))
                            finalcomp.append(int(partslist[posi]))
                        elif lstcompstart[1]==False and lstcompinside==True:
                            finalcomp=lstcompstart
                            posi=finalcomp[2]
                            finalcomp.append(0)
                            finalcomp.append(int(partslist[0]))

                    elif len(lstcompinside)>0 and len(lstcompend)>0:#middle + end
                        if lstcompinside[1]==lstcompend[1] and lstcompinside[3]>lstcompend[3]:
                            finalcomp=lstcompinside
                            posi=finalcomp[2]
                            finalcomp.append(int(partslist[posi-1]))
                            finalcomp.append(int(partslist[posi]))
                        elif lstcompinside[1]==lstcompend[1] and lstcompinside[3]<=lstcompend[3]:
                            finalcomp=lstcompend
                            posi=finalcomp[2]
                            finalcomp.append(int(partslist[posi]))
                            finalcomp.append(len(hairpin)-1)
                        elif lstcompinside[1]==False and lstcompend[1]==True:
                            finalcomp=lstcompinside
                            posi=finalcomp[2]
                            finalcomp.append(int(partslist[posi-1]))
                            finalcomp.append(int(partslist[posi]))
                        elif lstcompinside[1]==True and lstcompend[1]==False:
                            finalcomp=lstcompend
                            posi=finalcomp[2]
                            finalcomp.append(int(partslist[posi]))
                            finalcomp.append(len(hairpin)-1)

                    elif len(lstcompstart)>0 and len(lstcompend)>0:
                        if lstcompstart[1]==lstcompend[1] and lstcompstart[3]>=lstcompend[3]:
                            finalcomp=lstcompstart
                            posi=finalcomp[2]
                            finalcomp.append(0)
                            finalcomp.append(int(partslist[0]))
                        elif lstcompstart[1]==lstcompend[1] and lstcompstart[3]<lstcompend[3]:
                            finalcomp=lstcompend
                            posi=finalcomp[2]
                            finalcomp.append(int(partslist[posi]))
                            finalcomp.append(len(hairpin)-1)
                        elif lstcompstart[1]==True and lstcompend[1]==False:
                            finalcomp=lstcompend
                            posi=finalcomp[2]
                            finalcomp.append(int(partslist[posi]))
                            finalcomp.append(len(hairpin)-1)
                        elif lstcompstart[1]==False and lstcompend[1]==True:
                            finalcomp=lstcompstart
                            posi=finalcomp[2]
                            finalcomp.append(0)
                            finalcomp.append(int(partslist[0]))

                    if len(finalcomp)>0 and stat=="old":
                        currhairpin=hairpin[finalcomp[4]:finalcomp[5]+1]
                        oldlstlstr=[]
                        oldlstlstl=[]
                        for nuc in range(0,len(currhairpin)):
                            if currhairpin[nuc]==")":
                                oldlstlstr.append(nuc+finalcomp[4])
                            elif currhairpin[nuc]=="(":
                                oldlstlstl.append(nuc+finalcomp[4])

                    if len(finalcomp)>0 and stat=="new":
                        currhairpin=hairpin[finalcomp[4]:finalcomp[5]+1]
                        newlstlstr=[]
                        newlstlstl=[]
                        for nuc in range(0,len(currhairpin)):
                            if currhairpin[nuc]==")":
                                newlstlstr.append(nuc+finalcomp[4])
                            elif currhairpin[nuc]=="(":
                                newlstlstl.append(nuc+finalcomp[4])

                    if stat=="old":
                        finaloldcomp=[]
                        finaloldcomp=finalcomp

                    if stat=="new":
                        finalnewcomp=[]
                        finalnewcomp=finalcomp

                if parts==1:
                    loopstart=hairpin.rfind('(')+1
                    loopend=hairpin.find(')')-1
                    finalcomp1=[]

                    if (((x in range(loopstart,loopend+1)) and (y not in range(loopstart,loopend+1))) or ((x not in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1)))) and stat=="old":
                        oldbroken=False
                        inloopold=True
                        oldloop=True
                        oldncounts=hairpin.count('(')+hairpin.count(')')
                        finalcomp1.append(oldbroken)
                        finalcomp1.append(oldloop)
                        finalcomp1.append(0)
                        finalcomp1.append(oldncounts)

                    if ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))) and stat=="old":
                        oldbroken=True
                        inloopold=True
                        oldloop=True

                    elif (((x in range(0,loopstart)) and (y in range(0,loopstart))) or ((x in range(loopend+1,len(hairpin))) and (y in range(loopend+1,len(hairpin))))) and stat=="old":
                        oldloop=False
                        oldbroken=False
                        oldncounts=hairpin.count('(')+hairpin.count(')')
                        finalcomp1.append(oldbroken)
                        finalcomp1.append(oldloop)
                        finalcomp1.append(0)
                        finalcomp1.append(oldncounts)

                    elif ( ((x in range(0,loopstart)) and (y in range(loopend+1,len(hairpin)-1))) or ((y in range(0,loopstart)) and (x in range(loopend+1,len(hairpin)-1))))and stat=="old":
                        oldbroken=True
                        oldloop=True

                    elif (((x in range(loopstart,loopend+1)) and (y not in range(loopstart,loopend+1))) or ((x not in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1)))) and stat=="new":
                        newbroken=False
                        newloop=True
                        inloopnew=True
                        newncounts=hairpin.count('(')+hairpin.count(')')
                        finalcomp1.append(newbroken)
                        finalcomp1.append(newloop)
                        finalcomp1.append(0)
                        finalcomp1.append(newncounts)
                    elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))) and stat=="new":
                        newbroken=True
                        newloop=True
                        inloopnew=True
                    elif (((x in range(0,loopstart)) and (y in range(0,loopstart))) or ((x in range(loopend+1,len(hairpin))) and (y in range(loopend+1,len(hairpin)))))  and stat=="new":
                        newloop=False
                        newbroken=False
                        newncounts=hairpin.count('(')+hairpin.count(')')
                        finalcomp1.append(newbroken)
                        finalcomp1.append(newloop)
                        finalcomp1.append(0)
                        finalcomp1.append(newncounts)
                    elif ( ((x in range(0,loopstart)) and (y in range(loopend+1,len(hairpin)-1))) or ((y in range(0,loopstart)) and (x in range(loopend+1,len(hairpin)-1))))and stat=="new":
                        newbroken=True
                        oldloop=True

                    if stat=="old":
                        oldparts=parts
                        oldncounts=int(hairpin.count("("))+int(hairpin.count(")"))

                    if stat=="new":
                        newparts=parts
                        newncounts=int(hairpin.count("("))+int(hairpin.count(")"))

                    if len(finalcomp1)>0 and stat=="old":
                        oldlstlstr=[]
                        oldlstlstl=[]
                        for nuc in range(0,len(hairpin)):
                            if hairpin[nuc]==")":
                                oldlstlstr.append(nuc)
                            elif hairpin[nuc]=="(":
                                oldlstlstl.append(nuc)

                    if len(finalcomp1)>0 and stat=="new":
                        newlstlstr=[]
                        newlstlstr=[]
                        for nuc in range(0,len(hairpin)):
                            if hairpin[nuc]==")":
                                newlstlstr.append(nuc)
                            elif hairpin[nuc]=="(":
                                newlstlstl.append(nuc)

                    finalcomp1.append(0)
                    finalcomp1.append(len(hairpin)-1)

                    if stat=="old":
                        finaloldcomp=[]
                        finaloldcomp=finalcomp1

                    if stat=="new":
                        finalnewcomp=[]
                        finalnewcomp=finalcomp1

                done=False
                if stat=="new":
                    if len(finaloldcomp)>2:
                        oldbroken=False
                        oldloop=finaloldcomp[1]
                        oldhairpstart=finaloldcomp[4]
                        oldhairpend=finaloldcomp[5]
                        oldncounts=finaloldcomp[3]

                    else:
                        oldbroken=True
                        oldloop=True
                        oldhairpstart=-1
                        oldhairpend=-1
                        oldncounts=0

                    if len(finalnewcomp)>2:
                        newbroken=False
                        newloop=finalnewcomp[1]
                        newhairpstart=finalnewcomp[4]
                        newhairpend=finalnewcomp[5]
                        newncounts=finalnewcomp[3]
                    else:
                        newbroken=True
                        newloop=True
                        newhairpstart=-1
                        newhairpend=-1
                        newncounts=-1

                    printlog(["old: oldbroken= ",oldbroken,"/oldloop= ",oldloop,"start= ",oldhairpstart,"end= ",oldhairpend,"Counts= ",oldncounts, "old score= ",oldscore,"old parts= ",oldparts])
                    printlog(["new: newbroken= ",newbroken,"/newloop= ",newloop,"start= ",newhairpstart,"end= ",newhairpend,"Counts= ",newncounts, "new score= ",newscore,"new parts= ",newparts])
                    temcorrectsplit=precid.split()
                    precId=temcorrectsplit[0].strip()
                    precides=temcorrectsplit[1].strip()
                    if (oldbroken==False and newbroken==False):# I only want to save the old when are possibly flipped into new, so no need to add old when new is broken!
                        printlog("save old good")
                        mirstarcorr,mirstarsposcorr,mirstareposcorr,miroriencorr=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)

                        if mirstarsposcorr!=-1 and mirstareposcorr!=-1:
                            temcorrectsplit=precid.split()
                            precId=temcorrectsplit[0].strip()
                            precides=temcorrectsplit[1].strip()
                            listoldstatus.append(precId)
                            listoldstatus.append(precides)
                            listoldstatus.append(oldprecseq.replace('T','U'))
                            listoldstatus.append(matseq)
                            listoldstatus.append(mirstarcorr)
                            listoldstatus.append(str(miroriencorr))

                    done=True

                if done:
                    if oldbroken and newbroken:
                        printlog("remove old and no new")
                        listofboth.append(precid)
                        listofboth.append(oldprecseq)
                    elif (not oldbroken and not oldloop and oldparts==1) and ((oldscore<=newscore) or (newbroken)):
                        printlog("old is Better1")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofold.append(precid)
                            listofold.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (not newbroken and not newloop and newparts==1 and newscore<=-10.00) and ((oldscore>newscore) or (oldbroken)):
                        printlog("new is Better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)

                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(newseq[:mirstarspos])
                            ycut=len(newseq[mirstarepos+1:])
                            ##mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq=newseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq=newseq
                            elif mirorien=='3p' and xcut>10:
                                newseq=newseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq=newseq

                            listofnew.append(precid)
                            listofnew.append(newseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (not newbroken and not newloop and newparts==1 and newscore>-10.00) and ((oldscore>newscore) or (oldbroken)):
                        printlog("remove old and no new")
                        listofboth.append(precid)
                        listofboth.append(oldprecseq)
                    elif (oldscore<=newscore and not oldbroken and oldloop and not newbroken and not newloop and newparts==1 and newncounts> oldncounts and newscore<=-10.00):
                        printlog("new is Better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(newseq[:mirstarspos])
                            ycut=len(newseq[mirstarepos+1:])
                            ##mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq=newseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq=newseq
                            elif mirorien=='3p' and xcut>10:
                                newseq=newseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq=newseq

                            listofnew.append(precid)
                            listofnew.append(newseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (oldscore<=newscore and not oldbroken and oldloop and not newbroken and not newloop and newparts==1 and (newscore>-10.00 or newncounts<=oldncounts)) or (oldscore<newscore and not oldbroken and oldloop and(newloop or newbroken or newparts>1)):
                        printlog("old with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofoldloop.append(precid)
                            listofoldloop.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif oldscore<=newscore and not oldbroken and not oldloop and oldparts>1 and newparts>=1 and not newbroken and not newloop and newncounts>oldncounts:
                        printlog("new is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(newseq[:mirstarspos])
                            ycut=len(newseq[mirstarepos+1:])
                            ##mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq=newseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq=newseq
                            elif mirorien=='3p' and xcut>10:
                                newseq=newseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq=newseq

                            listofnew.append(precid)
                            listofnew.append(newseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (oldscore<=newscore and not oldbroken and not oldloop and oldparts>1) and ((newparts>=1 and not newbroken and not newloop and newncounts<=oldncounts) or (newbroken)):
                        printlog("old is better2")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])

                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofold.append(precid)
                            listofold.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)


                    elif ((newscore<=oldscore or oldscore<=newscore) and newscore<=-10.00 and not newbroken and newloop and oldbroken):
                        printlog("new is better with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(newseq[:mirstarspos])
                            ycut=len(newseq[mirstarepos+1:])
                            ##mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq=newseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq=newseq
                            elif mirorien=='3p' and xcut>10:
                                newseq=newseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq=newseq

                            listofnewloop.append(precid)
                            listofnewloop.append(newseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif ((newscore<=oldscore or oldscore<=newscore) and newscore>-10.00 and not newbroken and newloop and oldbroken):
                        printlog("remove old and no news")
                        listofboth.append(precid)
                        listofboth.append(oldprecseq)

                    elif (newscore<oldscore and not newbroken and newloop and not oldbroken and not oldloop):
                        printlog("old is better3")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofold.append(precid)
                            listofold.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (newscore<oldscore and not newbroken and newloop and not oldbroken and oldloop):
                        printlog("old with warnings")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofoldloop.append(precid)
                            listofoldloop.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (oldscore<=newscore and oldparts>1 and not oldbroken and not oldloop) and (newbroken or newloop):
                        printlog("old is better4")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofold.append(precid)
                            listofold.append(oldprecseq)
                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (newscore<=oldscore and  newparts>1 and not newbroken and not newloop and newscore<=-10.00) and ((not oldbroken and not oldloop and oldparts>1 and newncounts>oldncounts) or oldbroken):
                        printlog("new is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(newseq[:mirstarspos])
                            ycut=len(newseq[mirstarepos+1:])
                            ##mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq=newseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq=newseq
                            elif mirorien=='3p' and xcut>10:
                                newseq=newseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq=newseq

                            listofnew.append(precid)
                            listofnew.append(newseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (newscore<=oldscore and  newparts>1 and not newbroken and not newloop and newscore<=-10.00) and ((not oldbroken and not oldloop and oldparts>1 and newncounts<=oldncounts)):
                        printlog("old is better5")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""

                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])

                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofold.append(precid)
                            listofold.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (newscore<oldscore and not newbroken and not newloop and newparts>1 and newscore>-10.00) and (oldparts>1 and not oldbroken and not oldloop):
                        printlog("old is better6")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofold.append(precid)
                            listofold.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (newscore<oldscore and not newbroken and not newloop and not oldbroken and not oldloop and oldparts==1 and newparts>1):
                        printlog("old is Better7")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""

                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofold.append(precid)
                            listofold.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (newscore<oldscore and not newbroken and not newloop and newparts>1 and newscore<=-10.00) and ( (not oldbroken and oldloop)):
                        printlog("new is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(newseq[:mirstarspos])
                            ycut=len(newseq[mirstarepos+1:])
                            ##mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq=newseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq=newseq
                            elif mirorien=='3p' and xcut>10:
                                newseq=newseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq=newseq

                            listofnew.append(precid)
                            listofnew.append(newseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (newscore<oldscore and not newbroken and not newloop and newparts>1) and (newscore>-10.00 and not oldbroken and oldloop):
                        printlog("old with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofoldloop.append(precid)
                            listofoldloop.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (newscore<oldscore or oldscore<newscore) and (not newbroken and not newloop and newparts>1 and oldbroken and newscore<=-10.00):
                        printlog("new is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(newseq[:mirstarspos])
                            ycut=len(newseq[mirstarepos+1:])
                            ##mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq=newseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq=newseq
                            elif mirorien=='3p' and xcut>10:
                                newseq=newseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq=newseq

                            listofnew.append(precid)
                            listofnew.append(newseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (newscore<oldscore or oldscore<newscore) and (not newbroken and not newloop and newparts>1 and oldbroken and newscore>-10.00):
                        printlog("remove old and no new")
                        listofboth.append(precid)
                        listofboth.append(oldprecseq)

                    elif newscore<oldscore and (newbroken or newloop) and newparts>=1 and not oldbroken and not oldloop and oldparts>=1:
                        printlog("old is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofold.append(precid)
                            listofold.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif newscore<oldscore and (newbroken or newloop) and newparts>=1 and not oldbroken and oldloop and oldparts>=1:
                        printlog("old with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofoldloop.append(precid)
                            listofoldloop.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif newscore<=-10.00 and not newbroken and not newloop and oldbroken:
                        printlog("new is better here")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(newseq[:mirstarspos])
                            ycut=len(newseq[mirstarepos+1:])
                            ##mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq=newseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq=newseq
                            elif mirorien=='3p' and xcut>10:
                                newseq=newseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq=newseq

                            listofnew.append(precid)
                            listofnew.append(newseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif newscore<=-10.00 and not newbroken and newloop and oldbroken:
                        printlog("new warning better here")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(newseq[:mirstarspos])
                            ycut=len(newseq[mirstarepos+1:])
                            #mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq=newseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq=newseq
                            elif mirorien=='3p' and xcut>10:
                                newseq=newseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq=newseq

                            listofnewloop.append(precid)
                            listofnewloop.append(newseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    elif (oldscore==newscore and not oldbroken and oldloop and not newbroken and newloop and newncounts==oldncounts) and ( (oldparts>1 and newparts>1) or (oldparts==1 and newparts==1)):
                        printlog("old with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog([mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listofmirstar.append(matstardesc)
                            listofmirstar.append(str(mirstar))
                            listofmirstar.append(str(mirstarspos)+".."+str(mirstarepos))
                            precsplit=precid.split()
                            listofmirstar.append(str(precsplit[1]))
                            xcut=len(oldprecseq[:mirstarspos])
                            ycut=len(oldprecseq[mirstarepos+1:])
                            if mirorien=='5p' and ycut>10:
                                oldprecseq=oldprecseq[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                oldprecseq=oldprecseq
                            elif mirorien=='3p' and xcut>10:
                                oldprecseq=oldprecseq[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                oldprecseq=oldprecseq

                            listofoldloop.append(precid)
                            listofoldloop.append(oldprecseq)

                        else:
                            listofboth.append(precid)
                            listofboth.append(oldprecseq)

                    if (newbroken==False and newncounts>0 and newscore<=-10.00):
                        splitgoodprec=precid.split()
                        listgoodnew.append(str(splitgoodprec[0]).strip())
                        listgoodnew.append(str(splitgoodprec[1]).strip())
                        printlog(["here look test:",newspos,newepos,matseq,newlstlstl,newlstlstr,newseq1,newhairpstart,newhairpend])
#                        for i in len(range(newlstlstl)):
#                            newlstlstl[i] -= newhairpstart
#
#                        for i in len(range(newlstlstr)):
#                            newlstlstr[i] -= newhairpstart
#
#                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,0,newhairpend-newhairpstart)
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq1,newhairpstart,newhairpend)
                        printlog(["test 1515",mirstar,mirstarspos,mirstarepos,mirorien])
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            printlog(["also here",mirstar,mirstarspos,mirstarepos])
                            matdescsplit=matdesc.split()
                            matdescsplit[1]=matdescsplit[1]+"-star"
                            matstardesc=""
                            for i in matdescsplit:
                                matstardesc=str(matstardesc)+" "+str(i)

                            listgoodnew.append(matstardesc)
                            listgoodnew.append(str(matseq))
                            listgoodnew.append(str(mirstar))
                            listgoodnew.append(str(mirstarspos)+".."+str(mirstarepos))
                            listgoodnew.append(str(mirorien))

                            xcut=len(newseq1[:mirstarspos])
                            ycut=len(newseq1[mirstarepos+1:])
                            ##mirorien
                            if mirorien=='5p' and ycut>10:
                                newseq1=newseq1[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq1=newseq1
                            elif mirorien=='3p' and xcut>10:
                                newseq1=newseq1[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq1=newseq1
                            listgoodnew.append(newseq1)#[precID[0](x-mir-y),precid(MI),matstardesc,mirseq,mirstar,start..end,newseq]
                            if str(splitgoodprec[0]).strip()=='ola-mir-30c' or str(splitgoodprec[1]).strip()=='MI0019480':
                                printlog("for test final here")
                                printlog([str(splitgoodprec[0]).strip(),str(splitgoodprec[1]).strip()])
                                printlog([mirstar,mirstarspos,mirstarepos,mirorien])
                                printlog(listgoodnew)

                        else:
                            del listgoodnew[-2:]

        return oldlstlstr,oldlstlstl,oldparts,finaloldcomp,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listgoodnew

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def comp_5p(ll, lr, sm, ss, precursor, run):
    try:
        printlog(['comp5p', ll, lr, sm, ss, precursor, run])
        if run > 3:
            return ss-run+1
        if sm in ll:
            if any(x in precursor[sm].upper() for x in ['G', 'C']):
                return ss-run+1
            if ss in lr :
                if any(x in precursor[sm].upper() for x in ['A', 'U']) and any(x in precursor[ss].upper() for x in ['G', 'C']):
                    return ss
                else:
                    return comp_5p(ll, lr, sm, ss+1, precursor, run+1)
            else:
                return comp_5p(ll, lr, sm, ss+1, precursor, run+1)
        else:
            if ss in lr:
                return ss
            else:
                return comp_5p(ll, lr, sm, ss+1, precursor, run+1)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno


def getmirstar(spos,epos,mature,lstl,lstr,precursor,hairpstart,hairpend):
    try:
        printlog("get mirstar here 18")
        mirflag=False
#        cutprec=precursor[hairpstart:hairpend]
        if epos>=lstr[0]:
            orien="3p"
        elif epos<lstr[0]:
            orien="5p"

        if ((orien == '5p' and (spos > lstl[-1] or epos>=lstr[0])) or (orien == '3p' and (spos < lstl[-1] or epos < lstr[0])) or epos<=spos):
             printlog("-1 wrong")
             mirflag=True
             return "",-1,-1,'p'

        if orien=="5p":
            tempr=lstr
            templ=lstl
            rev=tempr[::-1]
            printlog("get mirstar here 19")
            if spos in templ:
                sind=templ.index(spos)
                sposstar=rev[sind]
                mirstarepos=sposstar+2
            if epos in templ:
                eind=templ.index(epos)
                enposstar=rev[eind]
                mirstarspos=enposstar+2
            if spos not in templ:
                for i in templ:
                    if i>int(spos):
                        sind=templ.index(i)#first ( after spos
                        diff=i-spos
                        mirstarepos=rev[sind]+diff+2
                        break
            if epos not in templ:
                if max(lstl)<epos:
                    sind=templ.index(max(lstl))
                    sind1=max(lstl)
                    diff=epos-sind1
                    mirstarspos=rev[sind]-diff+2

                else:
                    for i in templ:
                        if i>int(epos):
                            sind=templ.index(i)#after the end of the mat
                            sind1=templ[(sind-1)]#before the end of the mat
                            diff=epos-sind1
                            mirstarspos=rev[sind-1]-diff+2
                            break
            if mirstarspos<=epos:#To avoid overlapping between mir and mir*, this would happen in case the mir is in the loop
                mirstarepos=mirstarepos+(epos-mirstarspos+1)
                mirstarspos=mirstarspos+(epos-mirstarspos+1)

            if mirstarepos>len(precursor)-1:#To avoid mir* end position going outside the precursor
                mirstarepos=len(precursor)-1
            mirstarspos = comp_5p(lstl, lstr, spos, mirstarspos, precursor, 1)  # Comparing 5'ends of mir and mir*
            mirstar=precursor[mirstarspos:mirstarepos+1]
            mirstar= mirstar.replace("T","U")#here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
            printlog(["get mirstar here 20",mirstar,mirstarspos,orien])
            mirflag=True
            return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))
        #break

        elif orien=="3p":
            tempr=lstr
            templ=lstl
            rev=templ[::-1]

            if spos in tempr:
                sind=tempr.index(spos)
                sposstar=rev[sind]
                mirstarepos=sposstar+2
            if epos in tempr:
                eind=tempr.index(epos)
                enposstar=rev[eind]
                mirstarspos=enposstar+2
            if spos not in tempr:
                for i in tempr:
                    if i>int(spos):
                        sind=tempr.index(i)#first ( after spos
                        diff=i-spos
                        mirstarepos=rev[sind]+diff+2
                        break
            if epos not in tempr:
                if max(lstr)<epos:
                    sind=tempr.index(max(lstr))
                    sind1=max(lstr)
                    diff=epos-sind1
                    mirstarspos=rev[sind]-diff+2
                else:
                    for i in tempr:
                        if i>int(epos):
                            sind=tempr.index(i)#after the end of the mat
                            sind1=tempr[sind-1]#before the end of the mat
                            diff=epos-sind1
                            mirstarspos=rev[sind-1]-diff+2
                            break
            if mirstarepos>=spos:#To avoid overlapping between mir and mir*, this would happen in case the mir is in the loop
                mirstarspos=mirstarspos-(mirstarepos-spos+1)
                mirstarepos=mirstarepos-(mirstarepos-spos+1)

            if mirstarspos<=0:
                mirstarspos=0

            mirstarspos = comp_5p(lstr, lstl, spos, mirstarspos, precursor, 1)  # Comparing 5'ends of mir and mir*
            mirstar=precursor[mirstarspos:mirstarepos+1]
            mirstar=mirstar.replace("T","U")#here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
            printlog(["get mirstar here 24",mirstar,mirstarspos,orien])
            mirflag=True
            return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))

        if not mirflag:
            printlog("no predicted mir*")
            return "",-1,-1,'p'

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def getmirstarbak(spos,epos,mature,lstl,lstr,precursor,hairpstart,hairpend):
    try:
        printlog("get mirstar here 18")
        mirflag=False
        cutprec=precursor[hairpstart:hairpend]
        if epos>=lstr[0]:
            orien="3p"
        elif epos<lstr[0]:
            orien="5p"

        if (spos in range(hairpstart,lstr[0]) and (epos>=lstr[0] or epos<spos)) or (spos in range(lstr[0],hairpend+1) and (epos<spos or epos>hairpend)):
             printlog("-1 wrong")
             mirflag=True
             return "",-1,-1,'p'

        if orien=="5p":
            tempr=lstr
            templ=lstl
            rev=tempr[::-1]
            printlog("get mirstar here 19")
            for i in templ:
                if i==int(spos):
                    sind=templ.index(i)#list index(in the list of foldings), of the nucleotide representing the first fold, after the first nucleotide of the mature
                    sx=templ[sind]#get the item at sind, which is the position in the precursor, of the first folding nucleotide after the first nucleotide of mature
                    sposstar=rev[sind]#start position of the mature*
                    mirstarspos=sposstar-len(mature)+1
                    mirstarepos=sposstar
                    mirstar=precursor[mirstarspos+2:mirstarepos+3]
                    if mirstarspos<=epos:#To avoid overlapping between mir and mir*, this would happen in case the mir is in the loop
                        mirstarepos=mirstarepos+(epos-mirstarspos+1)
                        mirstarspos=mirstarspos+(epos-mirstarspos+1)

                    if mirstarepos>len(precursor)-1:#To avoid mir* end position going outside the precursor
                        mirstarepos=len(precursor)-1

                    mirstarspos = comp_5p(lstl, lstr, spos, mirstarspos+2, precursor, 1)  # Comparing 5'ends of mir and mir*
                    mirstar=precursor[mirstarspos:mirstarepos+1]
                    mirstar= mirstar.replace("T","U")#here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
                    printlog(["get mirstar here 20",mirstar,mirstarspos,orien])
                    mirflag=True
                    return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))
                    break

                if i>int(spos):
                    printlog("get mirstar here 21")
                    sind=templ.index(i)#list index(in the list of foldings), of the nucleotide representing the first fold, after the first nucleotide of the mature
                    sx=templ[sind]#get the item at sind, which is the position in the precursor, of the first folding nucleotide after the first nucleotide of mature
                    snnuc=sx-spos#get the number of nucloetides between the first nucleotide of mature and the first folding
                    #TODO: check if the number less than the length of the mature
                    sx1=rev[sind]#get the item folding with sx in the other side, means in the other list (when it is already reversed)
                    sposstar=sx1+snnuc#start position of the mature*, is the first folding before the mature* + number of nucleotides found between first folding after mature
                    mirstarspos=sposstar-len(mature)+1
                    mirstarepos=sposstar

                    if mirstarspos<=epos:#To avoid overlapping between mir and mir*, this would happen in case the mir is in the loop
                        mirstarepos=mirstarepos+(epos-mirstarspos+1)
                        mirstarspos=mirstarspos+(epos-mirstarspos+1)

                    if mirstarepos>len(precursor)-1:#To avoid mir* end position going outside the precursor
                        mirstarepos=len(precursor)-1

                    mirstar=precursor[mirstarspos:mirstarepos+1]
                    mirstar=mirstar.replace("T","U")#here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
                    printlog(["get mirstar here 22",mirstar,mirstarspos,orien])
                    mirflag=True
                    return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))
                    break

        elif orien=="3p":
            tempr=lstr
            templ=lstl
            rev=templ[::-1]
            printlog("get mirstar here 23")
            printlog(rev)
            for i in tempr:
                printlog("i "+str(i))
                if i==int(spos):
                    printlog("i=spos "+str(spos))
                    sind=tempr.index(i)#list index(in the list of foldings), of the nucleotide representing the first fold, after the first nucleotide of the mature
                    sx=tempr[sind]#get the item at sind, which is the position in the precursor, of the first folding nucleotide after the first nucleotide of mature,WHICH IS i HERE
                    sposstar=rev[sind]#start position of the mature*
                    mirstarspos=sposstar-len(mature)+1
                    mirstarepos=sposstar
                    if mirstarepos>=spos:#To avoid overlapping between mir and mir*, this would happen in case the mir is in the loop
                        mirstarspos=mirstarspos-(mirstarepos-spos+1)
                        mirstarepos=mirstarepos-(mirstarepos-spos+1)

                    if mirstarspos<=0:
                        mirstarspos=0

                    mirstarspos = comp_5p(lstl, lstr, spos, mirstarspos+2, precursor, 1)  # Comparing 5'ends of mir and mir*
                    mirstar=precursor[mirstarspos:mirstarepos+1]
                    mirstar=mirstar.replace("T","U")#here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
                    printlog(["get mirstar here 24",mirstar,mirstarspos,orien])
                    mirflag=True
                    return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))
                    break

                if i>int(spos):
                    printlog("i>spos")
                    sind=tempr.index(i)#list index(in the list of foldings), of the nucleotide representing the first fold, after the first nucleotide of the mature
                    sx=tempr[sind]#get the item at sind, which is the position in the precursor, of the first folding nucleotide after the first nucleotide of mature. WHICH IS i HERE
                    snnuc=sx-spos#get the number of nucloetides between the first nucleotide of mature and the first folding
                    #TODO: check if the number less than the length of the mature
                    sx1=rev[sind]#get the item folding with sx in the other side, means in the other list (when it is already reversed)
                    sposstar=sx1+snnuc#start position of the mature*, is the first folding before the mature* + number of nucleotides found between first folding after mature

                    mirstarspos=sposstar-len(mature)+1
                    mirstarepos=sposstar
                    if mirstarepos>=spos:#To avoid overlapping between mir and mir*, this would happen in case the mir is in the loop
                        mirstarspos=mirstarspos-(mirstarepos-spos+1)
                        mirstarepos=mirstarepos-(mirstarepos-spos+1)

                    if mirstarspos<=0:#can be edited here
                        mirstarspos=0

                    printlog("mir*")
                    printlog([mirstarspos,mirstarepos])
                    mirstar=precursor[mirstarspos:mirstarepos+1]#
                    mirstar=mirstar.replace("T","U")#here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
                    printlog(["get mirstar here 25",mirstar,mirstarspos,orien])
                    mirflag=True
                    return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))
                    break

        if not mirflag:
            printlog("no predicted mir*")
            return "",-1,-1,'p'

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def foldnomat(inputfasta,outputfasta):#fold the new and the old sequences, using temp file every time I get the new sequence from the original
    try:
        f=os.popen("RNAfold -d3 --noPS --noLP <"+inputfasta)
        fi=f.read()
        wr=open(outputfasta,"w")
        f.close()
        wr.write(fi)
        wr.close
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def alignTostock(align):
    try:
        reads=openfile(align)
        stkfile=align+".stk"
        writes=open(stkfile,'a')
        numberlines=0#count lines
        items=[]
        listofids=[]
        maxspaces=0
        minlen=100000
        maxlen=0
        for line in reads:
            item=line.split()
            if len(item) > 0 and item[0]!="CLUSTAL" and ("*" not in line):
                numberlines=numberlines+1
                listofids.append(item[0])
                items.append(line.strip())
                numberofspaces=items[0].count(" ")
                if maxspaces<numberofspaces:
                    maxspaces=numberofspaces
        for i in listofids:
            s=len(i)
            if minlen>s:
                minlen=s
            if maxlen<s:
                maxlen=s
        if maxlen>len('#=GC SS_cons'):
            maxlen=maxlen
        if maxlen<=len('#=GC SS_cons'):
            maxlen=len('#=GC SS_cons')

        numberofblocks=listofids.count(listofids[0])
        numberofseqs=int(numberlines/numberofblocks)
        ranges=int(len(items)/2)
        writes.write('# STOCKHOLM 1.0'+"\n"+"\n")
        offset=0
        for i in range (0,numberofseqs):
            offset=0

            seqitem=items[i].split()
            seqid=seqitem[0]
            seq=seqid.strip()+" "*(maxlen-len(seqid)+1)+seqitem[1]

            for k in range(0,numberofblocks-1):
                offset=offset+numberofseqs
                item2=items[i+offset].split()
                seq=seq+item2[1]

            writes.write(seq+"\n")
        tempitem=seq.split()
        writes.write('#=GC SS_cons'+" "*(maxlen-len('#=GC SS_cons')+1)+"."*len(tempitem[1])+"\n"+"//")
        writest=str(writes)
        writeitem=writest.split("'")
        writes.close()
        return str(writeitem[1])

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def predict(align,matId,newmatID,matfile,filename,precdescrip,mapfile,directory,listremovedbroken,listremovedscore):
    try:
        if directory[-1]!="/":# to work in both cases, user puts / at the end of the directory or not
            directory=str(directory)+"/"
        stockfile=alignTostock(align)
        alignment = AlignIO.read(stockfile, "stockholm")
        listrecords=[]
        for record in alignment:
            if matId in record.id:
               originalmature =str(record.seq)
               alnmat=str(record.seq)
            else:
                listrecords.append(record.id)
                listrecords.append(str(record.seq))

        for k in range(0,len(listrecords)):
            if k%2!=0:
                originalseq=str(listrecords[k])#oprecursors with gaps
                originalID=str(listrecords[k-1])
                originalseq=re.sub('-', '',originalseq)#remove the gaps and get the original prec
                originalmature=re.sub('-', '',originalmature)#remove the gaps and get the original mature
                alnseq=str(listrecords[k])
                lastindex=len(alnseq)-len(originalmature)#get the last index where I should stop searching for mature in the prec, to avoid going out of range
                nuc=['A','a','C','c','G','g','T','t','u','U']
                maxstar=0#max exact matching
                starlist=[]
                maxnuc=0#max number of nucleotides in each window (the window is the size of the real mature)
                nuclist=[]
                maxindex=1000000
                maxindexlist=[]
                dup=[]
                # the window is taking substring from the aligned mature line with szie of  real mature, and move nuc by nuc on the aligned prec to see how many nucleotides matching
                for i in range (0,lastindex+1):#moving alon the precurson, and stop at the last index, where after that will be out of range
                    countnuc=0
                    countstar=0
                    itemmature=str(alnmat[i:i+len(originalmature)])#the window with size of the mature sequence, moving along the

                    for char in itemmature:#move in the window to check how many nucleotides we have and not gaps
                        if char in nuc:
                            countnuc=countnuc+1
                    nuclist.append(countnuc)#save number of nucleotides and not gaps in each window

                    itemseq=alnseq[i:i+len(originalmature)]
                    for chars in range(0,len(originalmature)):#count the stars in each window (how many exact matching between mature aligned and the prec aligned)
                        if itemseq[chars].upper()==itemmature[chars].upper():
                            countstar=countstar+1
                    starlist.append(countstar)#save the number of stars at each window

                    if countstar>maxstar:#get the index where the window with max exact matching starts,and store the max number of stars

                        maxstar=countstar
                        oldmaxindex=maxindex
                        maxindex=i

                    if countnuc>maxnuc:
                        maxnuc=countnuc

                    maxnuctemp=0
                    for k in range (0,len(starlist)):#move in number of windows, where the number of stars in a given window, at specific position is stored. the index of the list is the position number in the prec
                    #in case we have more than one window with maximum stars, we chooe the one with more nucleotides
                        if maxstar==starlist[k]:#check all the positions with maximum number of stars
                            if nuclist[k]>maxnuctemp:#get the max number of nuc, in the max star position
                                startindex=k
                                maxnuctemp=nuclist[k]

                updatematfile=open(matfile,'a')
                updatemapfile=open(mapfile,'a')
                predictedspos=startindex
                predictedepos=startindex+len(originalmature)-1
                predictedtofoldfile=directory+'temptofold.fa'
                tempredictedtofold=open(predictedtofoldfile,'w')
                tempredictedtofold.write(">"+originalID+"\n"+originalseq+"\n")
                tempredictedtofold.close()
                tempredictfold=str(directory+'tempfold.fa')

                foldnomat(predictedtofoldfile,tempredictfold)
                finalpredspos,finalpredepos=readfoldpredict(tempredictfold,predictedspos,predictedepos)

                if finalpredspos!=-1:
                    predictedmat=originalseq[finalpredspos:finalpredepos]
                    updatematfile.write(">"+newmatID+"\n"+predictedmat+"\n")
                    newmatIDsplit=newmatID.split()
                    famsplit=precdescrip.split()
                    tempfamname=str(famsplit[0])
                    famname=tempfamname[tempfamname.find('-')+1:]
                    pos=str(finalpredspos)+".."+str(finalpredepos-1)#to remove the plus added before (in the readfold function) we did that there to make it included when substringing
                    updatemapfile.write(filename+" "+famname+" "+ famsplit[1]+" "+famsplit[0]+" "+newmatIDsplit[1]+" "+pos+" "+newmatIDsplit[0]+"\n")
                elif finalpredspos==-1:
                    famsplit=precdescrip.split()
                    printlog(["broken",famsplit[1]])
                    if famsplit[1] not in listremovedbroken:
                        listremovedbroken.append(famsplit[1].strip())
        return listremovedbroken,listremovedscore

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def checknomat(precfile,mapfile,matfile,directory,precfilename,listremovedbroken,listremovedscore,nomats,listnomat):
    try:
        flagnomatexists=False
        if directory[-1]!="/":# to work in both cases, user puts / at the end of the directory or not
            directory=str(directory)+"/"
        else:
            directory=directory
        matIDs=[]
        precnomat=[]
        countnomat=0#count number of precs without mature

        pf = openfile(precfile)
        for record in SeqIO.parse(pf, 'fasta'):
            precAllID=record.description
            precID=record.id
            precseq=record.seq
            flag=0
            mf = openfile(mapfile)
            for line in mf:
                item=line.split()
                if precID in line:#and not 2 mat
                    flag=1
                    numofmat=((len(item)+1)//3)-1
                    for i in range (0,numofmat):
                        matIDs.append(item[4+i])
            if flag==0:
                subprecfile=open(directory+'nomat-'+precfilename+'.fa','a')
                countnomat=countnomat+1
                precnomat.append(precID)

                subprecfile.write(">"+precAllID+"\n"+str(precseq.strip())+"\n")
                subprecfile.close()

        if countnomat>0 and len(matIDs)>0:
            nomats=countnomat
            flagnomatexists=True
            printlog("here is True")
            for i in matIDs:
                mtf = openfile(matfile)
                for record in SeqIO.parse(mtf, 'fasta'):
                    if i.strip() in record.description:
                        tempmaturefile=open(directory+'tempmat.fa','a')
                        tempmaturefile.write(">"+record.description+"\n"+str(record.seq.strip())+"\n")

            tempmaturefile.close()
            nomatfilefiller = openfile(directory+'nomat-'+precfilename+'.fa')
            nomatfile=SeqIO.parse(nomatfilefiller,'fasta')
            for prec in nomatfile:
                listofmat=[]
                listofscore=[]
                splitprecdisc=prec.description.split()#because the clustal reads only id, so it will only write in thr result the id without MI00....
                prectempid=splitprecdisc[0]+"-"+splitprecdisc[1]
                listnomat.append((splitprecdisc[1]).strip())
                for mat in SeqIO.parse(openfile(directory+'tempmat.fa'),'fasta'):

                    temptoalign=open(directory+'temptoalign.fa','w')

                    splitmatdisc=mat.description.split()
                    mattempid=splitmatdisc[0]+"-"+splitmatdisc[1]#because the clustal reads only id, so it will only write in thr result the id without MIMAT...
                    temptoalign.write(">"+prectempid+"\n"+str(prec.seq)+"\n"+">"+mattempid+"\n"+str(mat.seq)+"\n")
                    temptoalign.close()
                    infile=directory+'temptoalign.fa'
                    outfile=directory+'temptoalign.aln'
                    cline = ClustalwCommandline("clustalw", infile=infile, outfile=outfile)

                    stdout,stder=cline()
                    fi=open(directory+'temptoalign.txt','w')
                    fi.write(stdout)
                    fi.close()

                    for line in openfile(directory+'temptoalign.txt'):
                        if 'Sequence 2:' in line:
                            lineitem=line.split()
                            lineitem1=lineitem[2].split('-')
                            tempmatID=lineitem1[-1]
                            listofmat.append(tempmatID)

                        if 'Sequences (1:2)' in line:
                            lineitem=line.split()
                            score=float(lineitem[4])
                            listofscore.append(score)
                maxscore=max(listofscore)
                indexofbest=listofscore.index(maxscore)
                bestMatID=listofmat[indexofbest]
                prectempsplit=prec.description.split(" ")
                printlog(["best for ",prectempsplit[1]," is ",bestMatID, "of ", listofscore, listofmat])

                if maxscore>=21:
                    for mattopredict in SeqIO.parse(openfile(directory+'tempmat.fa'),'fasta'):
                        if bestMatID in mattopredict.description:
                            matpredictid=mattopredict.id#the first part of the ID which be read and mentioned in the alignment file
                            matpredictID=str(mattopredict.description)#the full ID (with description to write ti full in the fasta)
                            matpredictseq=str(mattopredict.seq)
                            break
                    temptopredictfa=open(directory+'temptopredict.fa','w')
                    tempinfilepredict=directory+'temptopredict.fa'
                    tempoutfilepredict=directory+'temptopredict.aln'
                    temptopredictfa.write(">"+matpredictID+"\n"+matpredictseq+"\n"+">"+prec.description+"\n"+str(prec.seq))
                    temptopredictfa.close()
                    predictcline = ClustalwCommandline("clustalw", infile=tempinfilepredict, outfile=tempoutfilepredict)
                    predictcline()

                    matId=matpredictid
                    prectempsplit=prec.description.split()
                    mattempsplit=matpredictID.split()
                    newmatID=prectempsplit[0]+"-mat "+mattempsplit[1]+"/"+prectempsplit[1]+" "+prectempsplit[2]+" "+prectempsplit[3]
                    filename=getfilename(precfile)

                    listremovedbroken,listremovedscore=predict(tempoutfilepredict,matId,newmatID,matfile,filename,prec.description,mapfile,directory,listremovedbroken,listremovedscore)
                    stktemptopredict=tempoutfilepredict+".stk"
                    f=os.popen("rm "+str(stktemptopredict))

                elif maxscore<21:
                    prectempsplit=prec.description.split()
                    printlog(["if not score",prectempsplit[1]])
                    if prectempsplit[1] not in listremovedscore:
                        listremovedscore.append(prectempsplit[1].strip())
            nomatfile.close()
        elif countnomat>0 and len(matIDs)==0:
            flagnomatexists=True
            nomats=-1
            return flagnomatexists,nomats,listremovedbroken,listremovedscore,listnomat
        elif countnomat<=0:
            printlog("do the normal procedure")
            #flip here
            nomats=0
            flagnomatexists=False
            printlog("all precursors has mature, at least one")
        return flagnomatexists,nomats,listremovedbroken,listremovedscore,listnomat
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def getfilename(dirfile):#to get the name of the file without directory or extension, in principal was used to get famname
    try:
        if '/' in dirfile and "." in dirfile:
            dirf=dirfile[dirfile.rfind('/')+1:]
            filename=dirf[:dirf.rfind('.')]
        elif '/' in dirfile and "." not in dirfile:
            filename=dirfile[dirfile.rfind('/')+1:]
        elif '/' not in dirfile and "." in dirfile:
            filename=dirfile[:dirfile.rfind('.')]
        elif '/' not in dirfile and "." not in dirfile:
            filename=dirfile
        return str(filename)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def readfoldpredict(foldfile,predictedspos,predictedepos):#,predictedspos,predictedepos):
    try:
        fofi = openfile(foldfile)
        for line in SeqIO.parse(fofi,"fasta"):
            string=str(line.seq)
            if ((")") or ("(") or (".")) in str(line.seq):#to get where the folding starts, ex1: ACGTtgatagt..((..))(score) ex2: acgtatgat(((.)))(socre)
                rindex=string.find(")")
                lindex=string.find("(")
                pindex=string.find(".")
                splitindex=min(rindex,lindex,pindex)
                stall=string[splitindex:]
                st=stall[0:stall.rfind('(')]
                hairpin=stall[0:stall.rfind('(')]
                precsequence=string[0:splitindex]
                partslist=[]
                parts=0
                pos=0
                realpos=0
                left=[]
                right=[]
                lstlstr=[]#list of list for right ")"
                lstlstl=[]#list of list for left "("
                newlstlstr=[]#list of list for right ")" of the new sequence
                newlstlstl=[]#list of list for left "(" of the new sequence
                curr=0
                next=0
                cut=0
                flag=0
                countfold=0
                newbroken=False
                newloop=False
                newparts=0
                newncounts=0
                while pos<len(st): #to count the parts, loop in the string until we found another open "(" after a close ")", then means one more part.. in this case the pos assingned zero again, and the new st is the new part
                    parts=parts+1 #1 part: ..((.)).. #2parts: ..((.))..(((.)))..
                    tempright=[]  #positions of the ")"
                    templeft=[] ##positions of the "("
                    for i in st:
                        realpos=cut+pos #cut is where we find new "(" after closing ")",the real position is always the position in the original string (all string not when divided into new parts
                        if i=="(":
                            left.append(realpos)
                            curr=pos #current position, changes in each new string (new part, in case the main strign was more than one part)
                            templeft.append(realpos)
                        elif i==")":
                            right.append(realpos)
                            tempright.append(realpos)
                            next=pos
                        if next<curr and next>0: #here means that we have more one part, by finding "(" after ")"
                            endofcurrentpart=str(left[-1]-1)
                            startofnextpart=str(right[-1]+1)
                            left.pop()
                            templeft.pop()
                            realpos-1
                            st=st[pos:]#make the next part as new string
                            partslist.append(endofcurrentpart)
                            partslist.append(startofnextpart)
                            cut=realpos # cut at the real position of thr whole original string
                            pos=0#position=zero for the new string (new part)
                            templst=right
                            templst=[]
                            next=0
                            break

                        pos=pos+1

                    templst=right
                    lstlstr.append(tempright)
                    lstlstl.append(templeft)

                if (len(lstlstr) != len(lstlstl)) or ( (")" not in hairpin) and ('(' not in hairpin)) :#checking the parts
                    flag=1
                else:
                    for i in range(len(lstlstr)):
                        if len(lstlstr[i])!=len(lstlstl[i]):#check if the foldings not equal
                            flag=1

                x=predictedspos
                y=predictedepos
                finalpredictedspos=-1
                finalpredictedepos=-1
                if parts>1:
                    flagvalid=0 #flag if x and y are valid at least once in the parts in between the first hairpin and last hairpin
                    for i in range(0,len(partslist)):
                        if i==0:
                            loopstart=hairpin.rfind('(',0,int(partslist[0])+1)+1
                            loopend=hairpin.find(')',0,int(partslist[0])+1)-1
                            if (((x in range(0,loopstart)) and (y in range(0,loopstart))) or ((x  in range(loopend+1,int(partslist[0])+1)) and (y in range(loopend+1,int(partslist[0])+1)))):
                                printlog("a")
                                flagvalid=1
                                finalpredictedspos=predictedspos
                                finalpredictedepos=predictedepos

                            elif (x in range(0,loopstart))  and  (y in range(loopstart,loopend+1)):
                                flagvalid=1
                                shiftby=y-(loopstart-1)
                                if x<=shiftby:
                                    printlog("b")
                                    x=0
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    printlog("c")
                                    x=x-shiftby
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif (x in range(loopstart,loopend+1)) and (y in range(loopend+1,int(partslist[0])+1)):
                                flagvalid=1
                                shiftby=(loopend+1)-x
                                if (y+shiftby)>=(len(hairpin)-1):
                                    printlog("d")
                                    y=len(hairpin)-1
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                    break
                                else:
                                    printlog("e")
                                    y=y+shiftby
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))):
                                printlog("f")
                                finalpredictedspos=-1
                                finalpredictedepos=-1
                            elif (((x in range (0,int(partslist[0])+1)) and (y not in range (0,int(partslist[0])+1))) or ((x not in range (0,int(partslist[0])+1)) and (y in range (0,int(partslist[0])+1)))):#if x or y out of the  hairpin
                                printlog("g")
                                finalpredictedspos=-1
                                finalpredictedepos=-1
                            elif ((x in range(0,loopstart)) and (y in range(loopend+1,int(partslist[0])+1))) or ((x in range(loopend+1,int(partslist[0])+1)) and ((y in range(0,loopstart)))):
                                printlog("h")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                        if i==len(partslist)-1 and flagvalid!=1:
                            loopstart=hairpin.rfind('(',int(partslist[i])+1)+1
                            loopend=hairpin.find(')',int(partslist[i])+1)-1

                            if (((x in range(int(partslist[i]),loopstart)) and (y in range(int(partslist[i]),loopstart) )) or ((x in range(loopend+1,len(hairpin))) and (y in range(loopend+1,len(hairpin))))):
                                printlog("1")
                                flagvalid=1
                                finalpredictedspos=predictedspos
                                finalpredictedepos=predictedepos
                            elif (x in range(int(partslist[i]),loopstart)) and (y in range(loopstart,loopend+1)):
                                flagvalid=1
                                shiftby=y-(loopstart-1)
                                if x<=shiftby:
                                    printlog("2")
                                    x=0
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    printlog("3")
                                    x=x-shiftby
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif  (x in range(loopstart,loopend+1)) and (y in range(loopend+1,len(hairpin))):
                                shiftby=(loopend+1)-x
                                if (y+shiftby)>=(len(hairpin)-1):
                                    printlog("4")
                                    flagvalid=1
                                    y=len(hairpin)-1
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    printlog("5")
                                    flagvalid=1
                                    y=y+shiftby
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))):
                                printlog("6")
                                finalpredictedspos=-1
                                finalpredictedepos=-1
                            elif (((x not in range(int(partslist[i]),len(hairpin))) and (y in range(int(partslist[i]),len(hairpin)))) or ((y not in range(int(partslist[i]),len(hairpin))) and (x in range(int(partslist[i]),len(hairpin))))):
                                printlog([x,y])
                                printlog("7")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                            elif (((x in range(int(partslist[i]),loopstart)) and (y in range(loopend+1,len(hairpin)))) or ((x in range(loopend+1,len(hairpin))) and (y in range(int(partslist[i]),loopstart)))):
                                printlog("8")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                        if i!=0 and i!=len(partslist)-1 and i%2==0 and flagvalid!=1:
                            loopstart=hairpin.rfind('(',int(partslist[i-1]),int(partslist[i])+1)+1
                            loopend=hairpin.find(')',int(partslist[i-1]),int(partslist[i])+1)-1
                            if (((x in range(int(partslist[i-1]),loopstart)) and (y in range(int(partslist[i-1]),loopstart))) or ((x in range(loopend+1,int(partslist[i])+1)) and (y in range(loopend+1,int(partslist[i])+1)))):
                                printlog("9")
                                flagvalid=1
                                finalpredictedspos=predictedspos
                                finalpredictedepos=predictedepos

                            elif (x in range(int(partslist[i-1]),loopstart)) and (y in range(loopstart,loopend+1)):
                                flagvalid=1
                                shiftby=y-(loopstart-1)
                                if x<=shiftby:
                                    printlog("10")
                                    x=0
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    printlog("11")
                                    x=x-shiftby
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif  (y in range(loopend+1,int(partslist[i])+1)) and (x in range(loopstart,loopend+1)):
                                flagvalid=1
                                shiftby=(loopend+1)-x
                                if (y+shiftby)>=(len(hairpin)-1):
                                    printlog("12")
                                    y=len(hairpin)-1
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    printlog("13")
                                    y=y+shiftby
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))):
                                printlog("14")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                            elif ((x not in range(int(partslist[i-1]),int(partslist[i])+1) and (y in range(int(partslist[i-1]),int(partslist[i])+1))) or ((y not in range(int(partslist[i-1]),int(partslist[i])+1)) and (x in range(int(partslist[i-1]),int(partslist[i])+1)))):
                                printlog("15")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                            elif (((x in range(int(partslist[i-1]),loopstart)) and (y in range(loopend+1,int(partslist[i])+1))) or ((x in range(loopend+1,int(partslist[i])+1)) and (y in range(int(partslist[i-1]),loopstart)))):
                                printlog("16")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                if parts==1:
                    loopstart=hairpin.rfind('(')+1
                    loopend=hairpin.find(')')-1
                    if (((x in range(0,loopstart)) and (y in range(0,loopstart))) or ((x in range(loopend+1,len(hairpin))) and (y in range(loopend+1,len(hairpin))))):
                        finalpredictedspos=x
                        finalpredictedepos=y

                    elif (x in range(0,loopstart)) and (y in range(loopstart,loopend+1)):
                        shiftby=y-(loopstart-1)
                        if x<=shiftby:
                            x=0
                            y=y-shiftby
                            finalpredictedspos=x
                            finalpredictedepos=y
                        else:
                            x=x-shiftby
                            y=y-shiftby
                            finalpredictedspos=x
                            finalpredictedepos=y

                    elif (x in range(loopstart,loopend+1))  and (y in range(loopend+1,len(hairpin))):
                        shiftby=(loopend+1)-x
                        if (y+shiftby)>=(len(hairpin)-1):
                            y=len(hairpin)-1
                            x=x+shiftby
                            finalpredictedspos=x
                            finalpredictedepos=y
                        else:
                            y=y+shiftby
                            x=x+shiftby
                            finalpredictedspos=x
                            finalpredictedepos=y

                    elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))):
                        finalpredictedspos=-1
                        finalpredictedepos=-1

                    elif ( ((x in range(0,loopstart)) and (y in range(loopend+1,len(hairpin)-1))) or ((y in range(0,loopstart)) and (x in range(loopend+1,len(hairpin)-1)))):
                        finalpredictedspos=-1
                        finalpredictedepos=-1
        printlog([finalpredictedspos,(finalpredictedepos+1)])
        flagfinalparts=0
        if (finalpredictedspos!=-1 and finalpredictedepos!=-1 and partslist>1):
            for i in range(0,len(partslist)):
                if i==0:
                    if (finalpredictedspos in range(0,int(partslist[0])+1)) and (finalpredictedepos  in range(0,int(partslist[0])+1)):
                        flagfinalparts=1
                        break

                if i!=0 and i!=len(partslist)-1 and i%2==0:
                    if  (finalpredictedspos in range(int(partslist[i-1]),int(partslist[i])+1)) and (finalpredictedepos in range(int(partslist[i-1]),int(partslist[i])+1)):
                        printlog("here 2")
                        flagfinalparts=1
                        break

                if i==len(partslist)-1:
                    if  (finalpredictedspos in range(int(partslist[i]),len(hairpin))) and (finalpredictedepos in range(int(partslist[i]),len(hairpin))):
                        flagfinalparts=1
                        break

            if flagfinalparts==0 and parts>1:
                finalpredictedspos=-1
                finalpredictedepos=-1

        return finalpredictedspos,(finalpredictedepos+1)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def doalifold(alnfile,outdir):
    try:
        f=os.popen("RNAalifold --noPS "+alnfile)
        alifoldtemp=outdir+'alifoldtemp.txt'
        foldrestemp=open(alifoldtemp,'w')
        fi=f.read()
        foldrestemp.write(fi)
        foldrestemp.close()
        f.close()
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def getstructure(alifoldtemp):
    try:
        for line in (openfile(alifoldtemp)):
            item=line.split()
            hairpin = ''
            if len(item)>0 and 'A' not in item[0] and 'C' not in item[0] and 'G' not in item[0] and 'T' not in item[0] and 'a' not in item[0] and 'c' not in item[0] and 'g' not in item[0] and 't' not in item[0] and 'N' not in item[0] and 'n' not in item[0] and ('.' in item[0] or '(' in item[0] or ')' in item[0]):
                hairpin=item[0]
                return str(hairpin)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def cal_ent(s):
    try:
        listofchars=[0,0,0,0,0]

        seqi=s
        lengseqi=len(seqi)
        listofsites = []
        for i in seqi:#loop on the sequences on each site
            if i == 'A'or i =='a':
                listofchars[0]=(listofchars[0]+1)#/numofseqs
            if i == 'C'or i =='c':
                listofchars[1]=(listofchars[1]+1)#/numofseqs
            if i == 'G'or i =='g':
                listofchars[2]=(listofchars[2]+1)#/numofseqs
            if i == 'T'or i =='t' or i =='u' or i =='U'  :
                listofchars[3]=(listofchars[3]+1)#/numofseqs
            if i == '-'or i =='.':
                listofchars[4]=(listofchars[4]+1)#/numofseqs

        listofsites.append(listofchars)#list of sites adding the ratio for each of them
        listofchars = [0,0,0,0,0]
        for k in range (0,5):#find ratio of each nucleotide at a specific site
            listofchars[k]=listofchars[k]/lengseqi

        listofsites.append(listofchars)#list of sites adding the ratio for each of them
        PA=0
        PC=0
        PG=0
        PT=0
        PGA=0
        PA=(listofsites[0][0])/lengseqi
        PC=(listofsites[0][1])/lengseqi
        PG=(listofsites[0][2])/lengseqi
        PT=(listofsites[0][3])/lengseqi
        PGA=(listofsites[0][4])/lengseqi
        EA=0
        EC=0
        EG=0
        ET=0
        EGA=0

        if PA > 0:
            EA=(PA*(math.log(PA,2)))
        else:
            EA=0

        if PC > 0:
            EC=(PC*(math.log(PC,2)))
        else:
            EC=0

        if PG > 0:
            EG=(PG*(math.log(PG,2)))
        else:
            EG=0

        if PT > 0:
            ET=(PT*(math.log(PT,2)))
        else:
            ET=0

        if PGA > 0:
            EGA= (PGA*(math.log(PGA,2)))
        else:
            EGA=0

        #Entropy=-1* ((PA*(math.log(PA,2)))+ (PC*(math.log(PC,2))) + (PG*(math.log(PG,2))) + (PT*(math.log(PT,2))) + (PGA*(math.log(PGA,2))))
        Entropy=-1* (EA+EC+EG+ET+EGA)

        return Entropy
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def CalShanon(stkfile):
    try:
        aligEnt=0
        nameOffile=stkfile
        listofchars = [] #use as temp to count nucleotides at each site
        listofsites = [] #use to store the ratio of each nucleotide at each site, i filled it as list of list
        lengseq=0  #lenght of sequence, NOTE: all should be same
        numofseqs=0
        #loop to initialize nucleotides' list to zero
        for i in range (0,5):
            listofchars.append(0)

            #Fill the list of sites with the ratios of nucleotides
        alignment = AlignIO.read(nameOffile, "stockholm")

        for record1 in  alignment:
            numofseqs+=1
            recid1= record1.id
            seqi=str(record1.seq)
            lengseq=len(seqi)

        SeqEachSite=[]
        EntEachSite=[]

        for i in  range (0,lengseq):
            EntEachSite.append("")
            SeqEachSite.append("")

        for record in alignment:
            seqtemp=str(record.seq)
            for i in range(0,lengseq):
                SeqEachSite[i]=SeqEachSite[i]+seqtemp[i]

        for i in range(0,lengseq):
            EntEachSite=cal_ent(SeqEachSite[i])
            aligEnt=aligEnt+EntEachSite

        return aligEnt

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def correct(corid,flanking,countcorrected,countcorrectedTonew,listofnew,listofnewloop,listoldstatus,templong,listmisalignedcorr,listcorrected,listcorrectedori,listgoodnew):
    try:
        printlog("def correct")
        listidsnewnewloop=[]
        for i in range(0,len(listofnew),2):
            idsplit=listofnew[i].split()
            Id=idsplit[0].strip()#to search with id come from the final stk
            IdDes=idsplit[1].strip()#to search in longmat
            listidsnewnewloop.append(Id)
            listidsnewnewloop.append(IdDes)
        for i in range(0,len(listofnewloop),2):
            idsplit=listofnewloop[i].split()
            Id=idsplit[0].strip()#to search with id come from the final stk
            IdDes=idsplit[1].strip()#to search in longmat
            listidsnewnewloop.append(Id)
            listidsnewnewloop.append(IdDes)
        printlog(listidsnewnewloop)

        if corid in listoldstatus:
            corrind=int(listoldstatus.index(corid))#the x-mir-y
            corriddes=str(listoldstatus[corrind+1])#the MI00001 id
            if corid in listidsnewnewloop and corriddes in templong:
                indexlongmat=int(templong.index(corriddes.strip()))
                longseq=str(templong[indexlongmat+1]).replace('T','U')
                printlog(["long seq",longseq,templong[indexlongmat]])
                corrind=int(listoldstatus.index(corid))
                correctseq=str(listoldstatus[corrind+2])
                correctmir=str(listoldstatus[corrind+3])
                correctmirstar=str(listoldstatus[corrind+4])
                orien=str(listoldstatus[corrind+5])
                coor1=int(longseq.find(correctmir))
                coor2=int(longseq.find(correctmirstar))

                if coor2<coor1:
                    tempcor=correctmir[:]
                    correctmir=correctmirstar[:]
                    correctmirstar=tempcor[:]

                startmatlong=int(longseq.find(correctmir))#position of mir in the long seq
                startfinalseq=startmatlong-flanking#the start position in the long seq, based on user flanking
                templongseq=longseq[startfinalseq:]#cut the long, with userflanking number of nucleotides
                printlog(["1st temp",startmatlong,startfinalseq,templongseq])
                startmatstarlong=int(templongseq.find(correctmirstar))
                endmatstarlong=int(startmatstarlong+len(correctmirstar)-1)
                endfinalseq=endmatstarlong+flanking
                printlog(["2nd temp",startmatstarlong,endmatstarlong,endfinalseq])
                correctfinalseq=str(templongseq[:endfinalseq+1])
                listmisalignedcorr.append(corriddes)
                listmisalignedcorr.append(correctfinalseq)
                listmisalignedcorr.append(int(correctfinalseq.find(correctmir)))
                listmisalignedcorr.append(int(correctfinalseq.find(correctmir))+len(correctmir))
                listmisalignedcorr.append(int(correctfinalseq.find(correctmirstar)))
                listmisalignedcorr.append(int(correctfinalseq.find(correctmirstar))+len(correctmirstar))
                listmisalignedcorr.append(str(orien))
                countcorrected=countcorrected+1
                listcorrected.append(corriddes.strip())
                printlog(["corrected seq",correctfinalseq])

        if corid.strip() in listgoodnew and corid.strip() not in listidsnewnewloop and len(listgoodnew)>2:
            corrind=int(listgoodnew.index(corid.strip()))#x-mir-y
            corriddes=listgoodnew[corrind+1]#MI00
            if corriddes.strip() in templong:
                indexlongmat=int(templong.index(corriddes.strip()))
                longseq=str(templong[indexlongmat+1]).replace('T','U')
                printlog(["long seq",longseq,templong[indexlongmat]])
                correctseq=str(listgoodnew[corrind+6])
                correctmir=str(listgoodnew[corrind+3])
                correctmirstar=str(listgoodnew[corrind+4])
                orien=str(listgoodnew[corrind+6])
                coor1=int(longseq.find(correctmir))
                coor2=int(longseq.find(correctmirstar))

                if coor2<coor1:
                    tempcor=correctmir[:]
                    correctmir=correctmirstar[:]
                    correctmirstar=tempcor[:]

                startmatlong=int(longseq.find(correctmir))#position of mir in the long seq
                startfinalseq=startmatlong-flanking#the start position in the long seq, based on user flanking
                templongseq=longseq[startfinalseq:]#cut the long, with userflanking number of nucleotides
                printlog(["1st temp",startmatlong,startfinalseq,templongseq])
                startmatstarlong=int(templongseq.find(correctmirstar))
                endmatstarlong=int(startmatstarlong+len(correctmirstar)-1)
                endfinalseq=endmatstarlong+flanking
                printlog(["2nd temp",startmatstarlong,endmatstarlong,endfinalseq])
                correctfinalseq=str(templongseq[:endfinalseq+1])
                listmisalignedcorr.append(corriddes)
                listmisalignedcorr.append(correctfinalseq)
                listmisalignedcorr.append(int(correctfinalseq.find(correctmir)))
                listmisalignedcorr.append(int(correctfinalseq.find(correctmir))+len(correctmir))
                listmisalignedcorr.append(int(correctfinalseq.find(correctmirstar)))
                listmisalignedcorr.append(int(correctfinalseq.find(correctmirstar))+len(correctmirstar))
                listmisalignedcorr.append(str(orien))
                countcorrectedTonew=int(countcorrectedTonew)+1
                listcorrectedori.append(corriddes.strip())
                printlog(["corrected seq",correctfinalseq])
                if corid.strip()=='la-mir-30c' or corriddes.strip()=='MI0019480':
                    printlog("IN correct test here")
                    printlog(["long seq",longseq,templong[indexlongmat]])
                    printlog(["1st temp",startmatlong,startfinalseq,templongseq])
                    printlog(["2nd temp",startmatstarlong,endmatstarlong,endfinalseq])
                    printlog(["list misaligned",listmisalignedcorr])

        return int(countcorrected),int(countcorrectedTonew),listmisalignedcorr,listcorrected,listcorrectedori

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def outdirectory():
    try:
        E3.delete(0, END)
        outdir = tkFileDialog.askdirectory(parent=root,  title='Select output directory')
        dirbut["text"] = str(outdir) if outdir else 'Please enter your output directory'
        dirbut.config(fg='green')
        if dirbut["text"] =="Please enter your output directory":
            dirbut.config(fg='red')
        else:
            dirbut.config(fg='green')
        E3.insert(0,str(outdir)+'/')

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def filesdirectory():
    try:
        E6.delete(0, END)
        fdir = tkFileDialog.askdirectory(parent=root,  title='Select the files directory')
        fdirbut["text"] = str(fdir) if fdir else 'Please select the directory of the fasta files, of your list'
        if fdirbut["text"] =="Please select the directory of the fasta files, of your list":
            fdirbut.config(fg='red')
        else:
            fdirbut.config(fg='green')

        E6.insert(0,str(fdir)+'/')

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def browse_matrices():
    try:
        E7.delete(0, END)
        mrdir = tkFileDialog.askdirectory(parent=root,  title='Select the matrices directory')
        matricesbutton["text"] = str(mrdir) if mrdir else 'Please select the directory of the matrices for the dialign tool'
        if matricesbutton["text"] =="Please select the directory of the matrices for the dialign tool":
            matricesbutton.config(fg='red')
        else:
            matricesbutton.config(fg='green')

        E7.insert(0,str(mrdir))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def browse_mappingfile():
    try:
        E4.delete(0, END)
        E4.insert(0,"")#to avoid reading wrong directory if the user uploaded a file twice!
        mapfile = tkFileDialog.askopenfilename(title="Selecting the mapping file", filetypes = (("Template files", "*.type"), ("Txt files", "*.txt"), ("All files", "*")))
        mapbutton["text"] = str(mapfile) if mapfile else 'Please browse the mapping file .txt'
        if mapbutton["text"] =="Please browse the mapping file .txt":
            mapbutton.config(fg='red')
        else:
            mapbutton.config(fg='green')
        E4.insert(0,str(mapfile))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def browse_maturefile():
    try:
        E5.delete(0, END)
        E5.insert(0,"")#to avoid reading wrong directory if the user uploaded a file twice!
        matfile = tkFileDialog.askopenfilename(title="Selecting the mature fasta file",filetypes = (("Template files", "*.type"), ("Fasta files", "*.fa"), ("All files", "*")))
        maturebutton["text"] = str(matfile) if matfile else 'Please browse the mature fasta file'
        if maturebutton["text"] =="Please browse the mature fasta file":
            maturebutton.config(fg='red')
        else:
            maturebutton.config(fg='green')

        E5.insert(0,str(matfile))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def browse_file_list():
    try:
        E2.delete(0, END)#to avoid reading wrong directory if the user uploaded a file twice!
        fname = tkFileDialog.askopenfilename(title="Selecting list of fasta files",filetypes = (("Template files", "*.type"),("Fasta files", "*.fa") ,("Txt files", "*.txt"),("All files", "*")))
        broButton["text"] = str(fname) if fname else 'Please browse one txt file as a list of fasta files'
        if broButton["text"] =="Please browse one txt file as a list of fasta files":
            broButton.config(fg='red')
        else:
            broButton.config(fg='green')
        E2.insert(0,str(fname))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def browse_genomes():
    try:
        E1.delete(0, END)#to avoid reading wrong directory if the user uploaded a file twice!
        gname = tkFileDialog.askopenfilename(title="Selecting list of related genomes",filetypes = (("Template files", "*.type"), ("Txt files", "*.txt"), ("All files", "*")))
        genbutton["text"] = str(gname) if gname else 'Please browse the list of related genomes'
        if genbutton["text"] =="Please browse the list of related genomes":
            genbutton.config(fg='red')
        else:
            genbutton.config(fg='green')

        E1.insert(0,gname)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def sublist(filename):
    try:
        #filesdir=str(sys.argv[3])
        filesdir=str(E6.get())
        command="list"
        #mapfile=str(sys.argv[6])
        mapfile=str(E4.get())
        #matfile=str(sys.argv[7])
        matfile=str(E5.get())
        #matrdir=None            # Not needed of dialign run from bioconda
        #if len(sys.argv) > 9:
        #matrdir=str(sys.argv[9])
        matrdir=str(E7.get())
        tProcessed=0#all start with 't', are for the total of all families
        tRemoved=0
        tTotalnumberofSequences=0
        tpredicted=0
        tnoprediction=0
        tflippednotcorrected=0
        tflippedcorrected=0
        twith2mats=0
        twith0mats=0
        twith1mats=0
        toldshanon=0
        tnewshanon=0
        numberoffamilies=0
        tcountcorrected=0
        #moved from global:
        listofnew=[]
        listofnewloop=[]
        listoldstatus=[]
        listofoldloop=[]#add seq id and old seq
        listofold=[]#add seq id and old seq
        listofboth=[]#add the id, that shouldn't be added to the new file
        listremovedbroken=[]
        listremovedscore=[]
        listofmirstar=[]
        listnomat=[]
        list2mat=[]
        list2matcoor=[]
        list1matcoor=[]
        listmatcoor=[]
        listnogenomes=[]#list of IDs that no existing genome to search, for their species
        listnotingenome=[]#list of IDs that are not found in the genome searched, for their species
        nomats=0
        templong=[]
        userflanking=0
        listmisalignedcorr=[]
        listgoodnew=[]
        listcorrected=[]
        listcorrectedori=[]

        printlog(filename)
        del listofnew[:]#add seq id and new sequence
        del listofnewloop[:]#add seq id and new sequence
        del listofoldloop[:]#add seq id and old seq
        del listofold[:]#add seq id and old seq
        del listofboth[:]#add the id, that shouldn't be added to the new file
        del listremovedbroken[:]
        del listremovedscore[:]
        del listofmirstar[:]
        del listnomat[:]
        del list2mat[:]
        del listnogenomes[:]#list of IDs that no existing genome to search, for their species
        del listnotingenome[:]
        del list2matcoor[:]
        del list1matcoor[:]
        del templong[:]
        del listoldstatus[:]
        del listmisalignedcorr[:]
        del listgoodnew[:]
        del listcorrected[:]
        del listcorrectedori[:]
        numberoffamilies=numberoffamilies+1
        countcorrected=0
        countcorrectedTonew=0
        Processed=0
        Removed=0
        TotalnumberofSequences=0
        predicted=0
        noprediction=0
        flippednotcorrected=0
        flippedcorrected=0
        with2mats=0
        with0mats=0
        with1mats=0
        oldshanon=0
        newshanon=0
        Totalproc=[]
        Totalrem=[]
        Nomatspred=[]
        Nomatsnotpred=[]
        flippednotcorr=[]
        flippedcorr=[]
        seq0mat=[]
        seq1mat=[]
        seq2mat=[]
        newshan=[]
        oldshan=[]
        suma=[]
        sumall=[]
        nomats=0#in case no mat exists and no mats at all
        userflanking=0
        listnomatremoved=[]
        flagnomatexists=False

        if ".fa" in filename:
            filen=filesdir+str(filename).strip()
            #outdir=str(sys.argv[2])+filename.strip()+".out/"
            outdir=str(E3.get())+filename.strip()+".out/"
            f=os.popen("mkdir "+outdir)
            f.close()
            familyfileres=open(outdir+filename.strip()+"-res.fa","a")
        else:
            filename=str(filename).strip()
            #outdir=str(sys.argv[2])+filename.strip()+".out/"
            outdir=str(E3.get())+filename.strip()+".out/"
            f=os.popen("mkdir "+outdir)
            f.close()
            familyfileres=open(outdir+filename.strip()+"-res.fa","a")
            filen=filesdir+str(filename.strip())+".fa"

        OldShanon=0
        NewShanon=0
        familyfileresfinal=open(outdir+filename.strip()+"-Final.fasta","a")
        summaryfile=open(outdir+filename.strip()+"-summ.txt","a")
        anchorcoorfile=open(outdir+filename.strip()+"-Final.anc","a")
        tempcountsucnomat=0
        infile=""
        outfile=""
        infile=filen[:]
        outfile=outdir+filename.strip()+"-tempshan.aln"
        clustaline = ClustalwCommandline("clustalw", infile=infile, outfile=outfile)
        stdoutshan,stdershan=clustaline()
        alignTostock(outfile)
        OldShanon=CalShanon(outfile+'.stk')
        printlog(["OldShanon",OldShanon])
        fshan=os.popen("rm "+outfile)
        fshan.close()
        fshan=os.popen("rm "+filesdir+filename.strip()+".dnd")
        fshan.close()
        fshan=os.popen("rm "+outfile+".stk")
        fshan.close()
        infile=""
        outfile=""
        #userflanking=int(sys.argv[8])
        if flanking.get()=="Please Choose The number of Flanking Nucleotides (upstream/downstream the 5'/3' matures)":
                
                userflanking=int(10)
        else:
                userflanking=int(flanking.get())
        fl = openfile(filen)
        for rec in SeqIO.parse(fl,'fasta'):
            pidsplit=rec.description.split()
            pid=str(pidsplit[1])
            mat2seq=str(rec.seq)
            coorflag=0
            smat=""
            emat=""
            mf = openfile(mapfile)
            for line in mf:
                linesplit=line.split()
                if len(linesplit)>8 and pid in line:
                    list2mat.append(linesplit[2].strip())
                    numofmat=((len(linesplit)+1)//3)-1
                    smat=linesplit[4]#first mat ID
                    emat=linesplit[4+numofmat-1].strip()#last mat ID
                    mtf = openfile(matfile)
                    first = None
                    second = None
                    for k in SeqIO.parse(mtf,'fasta'):
                        if smat in k.description:
                            first = str(k.seq)
                        if emat in k.description:
                            second = str(k.seq)
                        if first and second:
                            list2mat.append(first)#add the first mat seq
                            list2mat.append(second)#add the second mat seq
                            break
                    if first and not second:
                        list2mat.append(first)#add the first mat seq
                        list2mat.append('')#empty string for second
                    if second and not first:
                        list2mat.append('')#empty string for first
                        list2mat.append(second)#add the second mat seq

        printlog(["list 2 mat is here",list2mat])

        if os.path.isfile(outdir+'nomat-'+filename.strip()+'.fa'):#before calling checknomat in Submit
            printlog(["The file "+filename.strip()+" already processed, will be done again"])
            f=os.popen("rm "+outdir+'nomat-'+filename+'.fa')

        flagnomatexists,nomats,listremovedbroken,listremovedscore,listnomat=checknomat(filen,mapfile,matfile,outdir,filename,listremovedbroken,listremovedscore,nomats,listnomat)
        printlog(["test A"])

        if flagnomatexists and nomats!=-1:
            printlog("flagmathere")
            if os.path.isfile(outdir+filename.strip()+"-new.fa"):
                f=os.popen("rm "+outdir+filename.strip()+"-new.fa")

            listnomatremoved=listremovedbroken+listremovedscore
            printlog(["all removed",listnomatremoved])

            if len(listnomatremoved)>0:
                printlog("all list")
                fl = openfile(filen)
                for record in SeqIO.parse(fl, 'fasta'):
                    tempdes=record.description
                    tempdeslst=record.description.split()
                    if tempdeslst[1].strip() not in listremovedscore and tempdeslst[1].strip() not in listremovedbroken:# and tempdeslst[1].strip() not in lst2mat:
                        newprecfile=open(outdir+filename+"-new.fa",'a')
                        newprecfile.write(">"+str(tempdes)+"\n"+str(record.seq)+"\n")
                        newprecfile.close()

            f=os.popen("rm "+outdir+"tempfold.fa "+outdir+"tempmat.fa "+outdir+"temptoalign.aln "+outdir+"temptoalign.dnd "+outdir+"temptoalign.fa "+outdir+"temptoalign.txt "+outdir+"temptofold.fa "+outdir+"temptopredict.aln "+outdir+"temptopredict.dnd "+outdir+"temptopredict.fa ")

            if len(listnomatremoved)>0:
                listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listnogenomes,listnotingenome,templong,listgoodnew=flip(filename.strip(),outdir+filename+"-new.fa",outdir,mapfile,matfile,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,list2mat,listnogenomes,listnotingenome,templong,listgoodnew)
            else:
                listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listnogenomes,listnotingenome,templong,listgoodnew=flip(filename.strip(),filen,outdir,mapfile,matfile,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,list2mat,listnogenomes,listnotingenome,templong,listgoodnew)

        elif flagnomatexists and nomats==-1:
            summaryfile.write("no matures for the sequences and no related matures in the mapping file\n")

        elif not flagnomatexists:
            printlog(["this flip",flagnomatexists])
            listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,listnogenomes,listnotingenome,templong,listgoodnew=flip(filename.strip(),filen,outdir,mapfile,matfile,listofnew,listofnewloop,listoldstatus,listofoldloop,listofold,listofboth,listofmirstar,listnomat,list2mat,listnogenomes,listnotingenome,templong,listgoodnew)#filename: filename/family, filen: the file itself(with the directory)

        printlog(["test B"])
        if os.path.isfile(outdir+filename+"-new.fa"):
            f=os.popen("rm "+outdir+filename+"-new.fa")

        if len(listofnew)>0:
            for i in range(0,len(listofnew)):
                if i%2==0:
                    printlog("new")
                    printlog([">"+listofnew[i],listofnew[i+1]])
                    familyfileres.write(">"+str(listofnew[i])+"\n"+str(listofnew[i+1].replace('T','U'))+"\n")

        if len(listofnewloop)>0:
            for i in range(0,len(listofnewloop)):
                if i%2==0:
                    printlog("new loop")
                    printlog([">"+listofnewloop[i],listofnewloop[i+1]])
                    familyfileres.write(">"+str(listofnewloop[i])+"\n"+str(listofnewloop[i+1].replace('T','U'))+"\n")

        if len(listofold)>0:
            for i in range(0,len(listofold)):
                if i%2==0:
                    printlog("old")
                    printlog([">"+listofold[i],listofold[i+1]])
                    familyfileres.write(">"+str(listofold[i])+"\n"+str(listofold[i+1].replace('T','U'))+"\n")

        if len(listofoldloop)>0:
            for i in range(0,len(listofoldloop)):
                if i%2==0:
                    printlog("old loop")
                    printlog([">"+listofoldloop[i],listofoldloop[i+1]])
                    familyfileres.write(">"+str(listofoldloop[i])+"\n"+str(listofoldloop[i+1].replace('T','U'))+"\n")

        listnomatbroken=listremovedbroken
        listnomatscore=listremovedscore

        if len(list2mat)>0:
            printlog("list2mat: "+str(list2mat))
            for i in range(0,len(list2mat),3):
                precID=list2mat[i]
                fmatseq=list2mat[i+1].strip()#first mat id
                ematseq=list2mat[i+2].strip()#last mat id
                fl = openfile(filen)
                for record in SeqIO.parse(fl,'fasta'):
                    pdescsplit=record.description.split()
                    pid=str(pdescsplit[1].strip())
                    pseq=str(record.seq)
                    xcut=userflanking
                    ycut=userflanking
                    precDes=record.description
                    precitem=precDes.split()
                    specie=precitem[2].strip()+" "+precitem[3]
                    #ALI CHECK THIS
                    mspos=None
                    mepos=None
                    mspos=pseq.find(fmatseq)
                    mepos=pseq.rfind(ematseq)  # rfind returns the last index of the match

#                    if not mspos and mepos:
#                        mspos=pseq.rfind(str(Seq(fmatseq).reverse_complement()))
#                        mepos=pseq.find(str(Seq(ematseq).reverse_complement()))
                    # Now we have found start and end position on plus and minus strands

                    xcutseq=len(pseq[:mspos])#used in case the seq not found in the genome
                    ycutseq=len(pseq[mepos+1:])

                    if precID==pid:
                        long2matseq=""
                        long2matseq,listnogenomes,listnotingenome=getindex2mat(pseq.replace("U","T"),specie,precID,precDes,listnogenomes,listnotingenome)  # from here on we have the reverse complement if seq on minus strand
                        printlog("here get")
                        printlog(long2matseq)

                        if long2matseq!="":
                            if xcut>=0 and xcut<=50:
                                fspos=(mspos+100)-xcut
                            elif xcut>50 or xcut<0:
                                fspos=(mspos+100)-10

                            if ycut>=0 and ycut<=50:
                                fepos=(mepos+100)+ycut
                            elif ycut>50 or ycut<0:
                                fepos=(mepos+100)+10

                            #cutpseq=long2matseq[fspos:fepos+1]
                            cutpseq=long2matseq[fspos:fepos+len(ematseq)]
                            familyfileres.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")
                            familyfileresfinal.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")

                        elif long2matseq=="":
                            if xcutseq>=xcut and xcut>=0 and xcut<=50:
                                fspos=mspos-xcut
                            elif xcutseq<xcut and xcut>=0 and xcut<=50:#if requested flankings more than available for this seq, because not found in genome,then take all available
                                fspos=0
                            elif xcutseq>10 and (xcut>50 or xcut<0):
                                fspos=mspos-10
                            elif xcutseq<=10 and (xcut50 or xcut<0):
                                fspos=0

                            if ycutseq>=ycut and ycut>=0 and ycut<=50:
                                fepos=mepos+ycut+1
                            elif ycutseq<ycut and ycut>=0 and ycut<=50:
                                fepos=len(pseq)
                            elif ycutseq>10 and (ycut>50 or ycut<0):
                                fepos=mepos+11
                            elif ycutseq<=10 and (ycut50 or ycut<0):
                                fepos=len(pseq)

                            #cutpseq=pseq[fspos:fepos]
                            cutpseq=pseq[fspos:fepos+len(ematseq)]
                            familyfileres.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")
                            familyfileresfinal.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")

        familyfileres.close()

        if len(listofmirstar)>0:
            mirstarfile=open(outdir+filename.strip()+"-mirstar.fa","a")
            mirstarmapfile=open(outdir+filename.strip()+"-mirstar-map.txt","a")

            for star in range(0,len(listofmirstar),4):#devided into 4s..#1:star desc 2:mirstar seq 3:positions 4:related prec ID
                mirstarfile.write(">"+str(listofmirstar[star].strip())+" "+str(listofmirstar[star+3])+"\n"+str(listofmirstar[star+1])+"\n")
                tempsplit=listofmirstar[star].split()#to get only the ID of the mirstar
                starsplit=listofmirstar[star].split()
                mirstarmapfile.write(str(listofmirstar[star+3])+" "+str(starsplit[1])+" "+str(listofmirstar[star+2])+"\n")
            mirstarmapfile.close()
            mirstarfile.close()
            del listofmirstar[:]

        if ".fa" in filename:
            filename=filename[:filename.find(".fa")]
            resultfastafile=outdir+filename.strip()+"-res.fa"
        else:
            resultfastafile=outdir+filename.strip()+"-res.fa"

        rff = openfile(resultfastafile)
        for frec in SeqIO.parse(rff, 'fasta'):
            resprecdesc=str(frec.description)
            resfilesplit=(frec.description).split()
            resprecid=str(resfilesplit[1].strip())
            mat2seq=str(frec.seq).replace('T','U')
            mat1seq=str(frec.seq).replace('T','U')
            printlog(["res",resprecid,mat2seq,mat1seq])

            mf = openfile(mapfile)
            for line in mf:
                coorflag =0   #
                linesplit=line.split()
                if len(linesplit)>8 and resprecid.strip() in linesplit:
                    numofmat=((len(linesplit)+1)//3)-1
                    list2matcoor.append(linesplit[2].strip())
                    firstmat=linesplit[4]
                    lastmat=linesplit[4+numofmat-1]

                    mtf = openfile(matfile)
                    first = None
                    second = None
                    for record in SeqIO.parse(mtf, 'fasta'):
                        curmatseq=str(record.seq)
                        if firstmat.strip() in record.description:
                            startmat=int(mat2seq.find(curmatseq))
                            endmat=startmat+len(curmatseq)-1
                            first = 'Found'
                            printlog(["heres new",mat2seq,startmat,endmat,curmatseq])

                        if lastmat.strip() in record.description:
                            startmatstar=int(mat2seq.find(curmatseq))
                            endmatstar=startmatstar+len(curmatseq)-1
                            second = 'Found'
                            printlog(["heres new",startmatstar,endmatstar,curmatseq])

                        if first and second:
                            coorflag=1
                            break

                    if coorflag==1:
                        list2matcoor.append(firstmat.strip())
                        list2matcoor.append(lastmat.strip())
                        list2matcoor.append(startmat)
                        list2matcoor.append(endmat)
                        list2matcoor.append(startmatstar)
                        list2matcoor.append(endmatstar)

                elif len(linesplit)<8 and resprecid in linesplit and resprecid.strip() not in templong:
                    star=False
                    nstar=False
                    curmatID=linesplit[4]
                    startmattemp=0
                    endmattemp=0
                    startmat=0
                    endmat=0
                    startmatstar=0
                    endmatstar=0
                    list1matcoor.append(resprecid)
                    finalseq="A"
                    tempseq="A"
                    numberendflank=0
                    startfinalseq=0
                    endfinalseq=0
                    mtf = openfile(matfile)
                    for reco in SeqIO.parse(mtf, 'fasta'):
                        splitreco=reco.description.split()

                        if curmatID == splitreco[1]:
                            curmatseq=str(reco.seq)
                            nstar=True
                            break

                    mtfs = openfile(outdir+filename.strip()+"-mirstar.fa")
                    for starrec in SeqIO.parse(mtfs,'fasta'):
                        curmatsplit=(starrec.description).split()
                        curmatsplit1=(curmatsplit[1]).split('-')
                        starrecID=curmatsplit1[0]

                        if (curmatID).strip()==(starrecID).strip() and (resprecid.strip() in starrec.description):
                            curmatstar=str(starrec.seq)
                            star=True
                            break

                    coortemp1=int(mat1seq.index(curmatseq))
                    coortemp2=int(mat1seq.index(curmatstar))

                    if coortemp2<coortemp1:
                        tempseqex=curmatseq[:]
                        curmatseq=curmatstar[:]
                        curmatstar=tempseqex[:]

                    printlog(["no long",coortemp1,coortemp2,curmatseq,curmatstar])

                    startmattemp=int(mat1seq.index(curmatseq))

                    if startmattemp<=userflanking:
                        startmat=startmattemp
                        startfinalseq=0
                    elif startmattemp>userflanking:
                        startmat=userflanking
                        startfinalseq=startmattemp-(startmat-userflanking)

                    endmat=startmat+len(curmatseq)-1
                    tempseq=mat1seq[startfinalseq:]

                    printlog(["curmat not",mat1seq,curmatseq])#,longseq)
                    printlog(['coor not',startmat,endmat])

                    if star:
                        startmatstar=int(tempseq.find(curmatstar))
                        endmatstar=startmatstar+int(len(curmatstar)-1)
                        numberendflank=int(len(tempseq)-endmatstar-1)

                        if numberendflank<=userflanking:
                            finalseq=tempseq[0:]
                        elif numberendflank>userflanking:
                            endfinalseq=endmatstar+userflanking
                            finalseq=tempseq[0:endfinalseq+1]

                        printlog(['coor star not',startmatstar,endmatstar,finalseq])#,startmatstarlong,endmatstarlong

                    if star and nstar and  startmatstar!=-1 and endmatstar!=-1 and startmatstar>endmat:
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break
                    elif star and nstar and  startmatstar!=-1 and endmatstar!=-1 and startmatstar<endmat:#to put them in order
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif (not star and nstar) or startmatstar==-1 or endmatstar==-1:
                        startmat=int(mat1seq.find(curmatseq))
                        endmat=startmat+len(curmatseq)-1
                        finalseq=mat1seq
                        startmstar=0
                        endmatstar=len(mat1seq)-1
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif (not star and not nstar) or startmatstar==-1 or endmatstar==-1:#impossible but just in case
                        startmat=0
                        endmat=len(mat1seq)-1
                        startmstar=0
                        endmatstar=len(mat1seq)-1
                        finalseq=mat1seq
                        list1matcoor.append("NULL")#mat original/here no mir
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                elif len(linesplit)<8 and resprecid in linesplit and resprecid.strip() in templong:
                    star=False
                    nstar=False
                    curmatID=linesplit[4]
                    startmat=0
                    endmat=0
                    startmatstar=0
                    endmatstar=0
                    longseq=""
                    startmatlong=0
                    startfinalseq=0
                    startmat=0
                    templongseq=""
                    startmatstarlong=0
                    endmatstarlong=0
                    endfinalseq=0
                    finalseq="A"
                    startmatstar=0
                    endmatstar=0
                    list1matcoor.append(resprecid)
                    indexlongmat=int(templong.index(resprecid.strip()))
                    longseq=str(templong[indexlongmat+1]).replace('T','U')

                    mtf = openfile(matfile)
                    for reco in SeqIO.parse(mtf, 'fasta'):
                        splitreco=reco.description.split()
                        if curmatID == splitreco[1]:
                            curmatseq=str(reco.seq)
                            nstar=True
                            break

                    for starrec in SeqIO.parse(openfile(outdir+filename.strip()+"-mirstar.fa"),'fasta'):
                        curmatsplit=(starrec.description).split()
                        curmatsplit1=(curmatsplit[1]).split('-')
                        starrecID=curmatsplit1[0]
                        if (curmatID).strip()==(starrecID).strip() and (resprecid.strip() in starrec.description):
                            curmatstar=str(starrec.seq)
                            star=True
                            break

                    printlog(["coor1temp",longseq,curmatseq])
                    coortemp1=int(longseq.index(curmatseq))
                    coortemp2=int(longseq.index(curmatstar))
                    if coortemp2<coortemp1:
                        tempseqex=curmatseq[:]
                        curmatseq=curmatstar[:]
                        curmatstar=tempseqex[:]

                    startmatlong=int(longseq.find(curmatseq))
                    startfinalseq=startmatlong-userflanking
                    startmat=userflanking
                    endmat=startmat+len(curmatseq)-1
                    templongseq=longseq[startfinalseq:]
                    printlog(["curmat long",mat1seq,curmatseq,longseq])
                    printlog(['coor long',startmat,endmat])

                    startmatstarlong=int(templongseq.find(curmatstar))
                    endmatstarlong=int(startmatstarlong+len(curmatstar)-1)
                    endfinalseq=endmatstarlong+userflanking
                    finalseq=str(templongseq[:endfinalseq+1])

                    printlog(["final seq",finalseq,curmatstar,endfinalseq+1])
                    startmatstar=int(finalseq.find(curmatstar))
                    endmatstar=int(startmatstar+len(curmatstar)-1)
                    printlog(['coor',startmatstar,endmatstar,startmatstarlong,endmatstarlong])

                    if star and nstar and  startmatstar!=-1 and endmatstar!=-1 and startmatstar>endmat:
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif star and nstar and  startmatstar!=-1 and endmatstar!=-1 and startmatstar<endmat:#to put them in order
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif (not star and nstar) or startmatstar==-1 or endmatstar==-1:
                        startmat=int(mat1seq.find(curmatseq))
                        endmat=startmat+len(curmatseq)-1
                        finalseq=mat1seq
                        startmstar=0
                        endmatstar=len(mat1seq)-1
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif (not star and not nstar) or startmatstar==-1 or endmatstar==-1:
                        startmat=0
                        endmat=len(mat1seq)-1
                        startmstar=0
                        endmatstar=len(mat1seq)-1
                        finalseq=mat1seq
                        list1matcoor.append("NULL")#mat original/here no mir
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

        familyfileresfinal.close()
        printlog([len(list2matcoor), len(list1matcoor), list2matcoor, list1matcoor])
        listmatcoor=list2matcoor+list1matcoor

        mi=0
        mk=0
        printlog(len(listmatcoor))
        printlog(len(listmatcoor)/7)

        while mi<=int(len(listmatcoor))-7:
            mk=mk+1
            r=0

            for n in range(mk+1,int(len(listmatcoor)/7)+1):
                anchorcoorfile.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+3])+" "+str(listmatcoor[mi+10+r])+" "+str(22)+" "+str(1)+"\n")
                anchorcoorfile.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+5])+" "+str(listmatcoor[mi+12+r])+" "+str(22)+" "+str(1)+"\n")
                r=r+7

            mi=mi+7

        printlog(["here listmatcoor 1 ",listmatcoor])
        maxidesc=0

        finalstk=open(outdir+filename.strip()+'.stk','a')
        finalstk.write('# STOCKHOLM 1.0\n')
        if matrdir:
            fs=os.environ["DIALIGN2_DIR"]=matrdir
            printlog(matrdir)
        f1=os.popen("dialign2-2 -n -fa  "+outdir+filename.strip()+'-Final.fasta')
        printlog(f1)
        f1.close()

        if os.path.isfile(outdir+filename.strip()+'-Final.fa'):
            doalifold(outdir+filename.strip()+"-Final.fa",outdir)
            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-Final.fa'),'fasta'):
                if len(rec.description.strip())>maxidesc:
                    maxidesc=len(rec.description.strip())

                if maxidesc<len('#=GC SS_cons'):
                    maxidesc=len('#=GC SS_cons')

            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-Final.fa'),'fasta'):
                finalstk.write(str(rec.description.strip())+" "*(maxidesc-len(rec.description.strip())+2)+str(rec.seq.strip())+"\n")

            struct=getstructure(outdir+'alifoldtemp.txt')
            f2=os.popen("rm "+outdir+'alifoldtemp.txt')
            f4=os.popen("rm "+outdir+filename.strip()+"-Final.ali")
            f2.close()
            f4.close()
            finalstk.write('#=GC SS_cons'+" "*(maxidesc-len('#=GC SS_cons')+2)+str(struct))
            finalstk.close()

            stki = openfile(outdir+filename.strip()+'.stk')
            alignment=AlignIO.read(stki,"stockholm")
            countcorrected=0
            countcorrectedTonew=0
            sthalfgaps=0
            ndhalfgaps=0
            sthalfsum=0
            ndhalfsum=0
            stnucnum=0
            ndnucnum=0
            numofseqs=0
            totalstnucnum=0
            totalndnucnum=0

            for record in alignment:
                lenseq=len(str(record.seq))
                seq=str(record.seq)
                sthalf=seq[0:int((lenseq/2))]
                ndhalf=seq[int((lenseq/2)):]
                sthalfgaps=sthalf.count('-')
                stnucnum=len(sthalf)-sthalfgaps
                totalstnucnum=totalstnucnum+stnucnum
                ndhalfgaps=ndhalf.count('-')
                ndnucnum=len(ndhalf)-ndhalfgaps
                totalndnucnum=totalndnucnum+ndnucnum
                sthalfsum=sthalfgaps+sthalfsum
                ndhalfsum=ndhalfgaps+ndhalfsum
                numofseqs=numofseqs+1
            printlog([sthalfsum,ndhalfsum,totalstnucnum,totalndnucnum])

            stnucavg=totalstnucnum/numofseqs
            ndnucavg=totalndnucnum/numofseqs
            alignment=AlignIO.read(outdir+filename.strip()+'.stk',"stockholm")

            for record in alignment:
                lenseq=len(str(record.seq))
                seq=str(record.seq)
                corid=str(record.id)
                sthalf=seq[0:int((lenseq/2))]
                ndhalf=seq[int((lenseq/2)):]
                sthalfgaps=sthalf.count('-')
                stnucnum=len(sthalf)-sthalfgaps
                ndhalfgaps=ndhalf.count('-')
                ndnucnum=len(ndhalf)-ndhalfgaps
                nucnum=stnucnum+ndnucnum

                if (stnucnum>(stnucavg*1.5) and (ndnucnum<(ndnucavg/2) or ndnucnum==0) ) or (stnucnum>0.7*nucnum and (ndnucnum<(ndnucavg/2) or ndnucnum==0)):#-(stnucavg*0.1)):
                    printlog([record.id,'st'])
                    printlog([stnucnum,ndnucnum,nucnum])
                    printlog(["list new good",corid,listgoodnew])
                    countcorrected,countcorrectedTonew,listmisalignedcorr,listcorrected,listcorrectedori=correct(corid.strip(),userflanking,countcorrected,countcorrectedTonew,listofnew,listofnewloop,listoldstatus,templong,listmisalignedcorr,listcorrected,listcorrectedori,listgoodnew)

                if (ndnucnum>(ndnucavg*1.5) and (stnucnum<(stnucavg/2) or stnucnum==0)) or (ndnucnum>0.7*nucnum and (stnucnum<(stnucavg/2) or stnucnum==0)):#-(ndnucavg*0.1)):
                    printlog([record.id,'nd'])
                    printlog([stnucnum,ndnucnum,nucnum])
                    printlog(["list new good",corid,listgoodnew])
                    countcorrected,countcorrectedTonew,listmisalignedcorr,listcorrected,listcorrectedori=correct(corid.strip(),userflanking,countcorrected,countcorrectedTonew,listofnew,listofnewloop,listoldstatus,templong,listmisalignedcorr,listcorrected,listcorrectedori,listgoodnew)

        if len(listmisalignedcorr)>0:
            familyfilerescorrected=open(outdir+filename+"-corrected.fasta","a")
            mirstarcorrectedfile=open(outdir+filename+"-mirstar-corrected.fa","a")

            for corrrecord in SeqIO.parse(openfile(outdir+filename.strip()+'-Final.fasta'),'fasta'):
                splitdes=(corrrecord.description).split()
                splitID=splitdes[1].strip()

                if splitID.strip() in listmisalignedcorr:
                    coind=listmisalignedcorr.index(splitID.strip())
                    corcoorind=listmatcoor.index(splitID.strip())
                    listmatcoor[corcoorind+3]=listmisalignedcorr[coind+2]#set again the coordinates of the corrected sequence, in the liost related to build the anc file
                    listmatcoor[corcoorind+4]=listmisalignedcorr[coind+3]
                    listmatcoor[corcoorind+5]=listmisalignedcorr[coind+4]
                    listmatcoor[corcoorind+6]=listmisalignedcorr[coind+5]
                    familyfilerescorrected.write(">"+str(corrrecord.description)+"\n"+str(listmisalignedcorr[coind+1])+"\n")

                else:
                    familyfilerescorrected.write(">"+str(corrrecord.description)+"\n"+str(corrrecord.seq)+"\n")

            familyfilerescorrected.close()
            printlog(["here listmatcoor 2 ",listmatcoor])

            for mstarrec in SeqIO.parse(openfile(outdir+filename.strip()+'-mirstar.fa'),'fasta'):
                splitdes=(mstarrec.description).split()
                descstar=str(mstarrec.description)
                oriseq=str(mstarrec.seq)
                splitID=str(splitdes[-1]).strip()

                if splitID.strip() in listmisalignedcorr:
                    newcorrind=listmisalignedcorr.index(splitID.strip())
                    orien=str(listmisalignedcorr[newcorrind+6])

                    if orien=='3p':
                        newsstar=listmisalignedcorr[newcorrind+2]
                        newestar=listmisalignedcorr[newcorrind+3]
                    else:
                        newsstar=listmisalignedcorr[newcorrind+4]
                        newestar=listmisalignedcorr[newcorrind+5]

                    correctseq=str(listmisalignedcorr[newcorrind+1])
                    newmatstar=correctseq[newsstar:newestar+1]
                    mirstarcorrectedfile.write(">"+descstar.strip()+"\n"+newmatstar+"\n")
                else:
                    mirstarcorrectedfile.write(">"+descstar.strip()+"\n"+oriseq+"\n")

            mirstarcorrectedfile.close()

            anchorcoorfilecorrected=open(outdir+filename.strip()+"-corrected.anc","a")
            printlog(["here listmatcoor 3 ",listmatcoor])
            mi=0
            mk=0
            while mi<=int(len(listmatcoor))-7:
                mk=mk+1
                r=0
                for n in range(mk+1,int(len(listmatcoor)/7)+1):
                    anchorcoorfilecorrected.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+3])+" "+str(listmatcoor[mi+10+r])+" "+str(22)+" "+str(1)+"\n")
                    anchorcoorfilecorrected.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+5])+" "+str(listmatcoor[mi+12+r])+" "+str(22)+" "+str(1)+"\n")
                    r=r+7

                mi=mi+7

            maxidesc=0
            anchorcoorfilecorrected.close()
            finalstkcorrected=open(outdir+filename.strip()+'corrected.stk','a')
            finalstkcorrected.write('# STOCKHOLM 1.0\n')
            if matrdir:
                fe=os.environ["DIALIGN2_DIR"]=matrdir
            f11=os.popen("dialign2-2 -n -fa  "+outdir+filename.strip()+'-corrected.fasta')
            f11.close()

            doalifold(outdir+filename.strip()+"-corrected.fa",outdir)

            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-corrected.fa'),'fasta'):
                if len(rec.description.strip())>maxidesc:
                    maxidesc=len(rec.description.strip())

                if maxidesc<len('#=GC SS_cons'):
                    maxidesc=len('#=GC SS_cons')

            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-corrected.fa'),'fasta'):
                finalstkcorrected.write(str(rec.description.strip())+" "*(maxidesc-len(rec.description.strip())+2)+str(rec.seq.strip())+"\n")

            struct=getstructure(outdir+'alifoldtemp.txt')
            f22=os.popen("rm "+outdir+'alifoldtemp.txt')
            f33=os.popen("rm "+outdir+filename.strip()+"-corrected.fa")
            f44=os.popen("rm "+outdir+filename.strip()+"-corrected.ali")
            f22.close()
            f33.close()
            f44.close()
            finalstkcorrected.write('#=GC SS_cons'+" "*(maxidesc-len('#=GC SS_cons')+2)+str(struct))
            finalstkcorrected.close()
            NewShanon=CalShanon(outdir+filename.strip()+'corrected.stk')

        else:
            if os.path.isfile(outdir+filename.strip()+'-Final.fa'):
                f3=os.popen("rm "+outdir+filename.strip()+"-Final.fa")
                f3.close()
                NewShanon=CalShanon(outdir+filename.strip()+'.stk')
                printlog(["stk file studied is: "+outdir+filename.strip()+'.stk'])
                printlog(["final.fa exists",NewShanon])

            else:
                NewShanon=OldShanon
                printlog(["final.fa DO NOT exists",NewShanon])

        printlog(["new shanon",NewShanon])
        finalcoor=open(outdir+filename+"-FinalCoor.txt","a")
        printlog(["here listmatcoor 4 ",len(listmatcoor),listmatcoor])

        for l in range(0,len(listmatcoor),7):
            finalcoor.write(str(listmatcoor[l]).strip()+" "+str(listmatcoor[l+1]).strip()+" "+str(listmatcoor[l+2]).strip()+" "+str(listmatcoor[l+3]).strip()+" "+str(listmatcoor[l+4]).strip()+" "+str(listmatcoor[l+5]).strip()+" "+str(listmatcoor[l+6]).strip()+"\n")

        finalcoor.close()

        summaryfile.write("---------------------------Original precursors with bad positioned matures ----------------------------\n")
        if len(listofoldloop)>0:
            for k in range(0,len(listofoldloop)):
                if k%2==0:
                    tempsplit=listofoldloop[k].split()
                    summaryfile.write(str(tempsplit[1].strip())+"\n")
        else:
            summaryfile.write("---> NO Original precursors with bad positioned matures\n")
        summaryfile.write("\n")
        summaryfile.write("---------------------------flipped/changed precursors ----------------------------\n")

        if len(listofnew)>0:
            for k in range(0,len(listofnew)):
                if k%2==0:
                    tempsplit=listofnew[k].split()
                    summaryfile.write(str(tempsplit[1].strip())+"\n")
        else:
            summaryfile.write("--->NO flipped/changed precursors\n")
        summaryfile.write("\n")
        summaryfile.write("---------------------------flipped/changed precursors; inloop ----------------------------\n")

        if len(listofnewloop)>0:
            for k in range(0,len(listofnewloop)):
                if k%2==0:
                    tempsplit=listofnewloop[k].split()
                    summaryfile.write(str(tempsplit[1].strip())+"\n")
        else:
            summaryfile.write("---> NO flipped/changed precursors; inloop\n")
        summaryfile.write("\n")

        summaryfile.write("---------------------------Flipped/changed precursors don't fit with the final alignment (changed back to original)----------------------------\n")
        if len(listcorrected)>0:
            for k in (listcorrected):
                summaryfile.write(k.strip()+"\n")
        else:
            summaryfile.write("--->None of the changed/flipped precursors was changed back to original\n")

        summaryfile.write("\n")

        if len(listcorrectedori)>0:
            summaryfile.write("---------------------------Original precursors, Change at the Alignment level (to fit the alignment----------------------------\n")
            for k in (listcorrectedori):
                summaryfile.write(k.strip()+"\n")

        summaryfile.write("\n")
        summaryfile.write("---------------------------Precursors without annotated mature----------------------------:\n")

        if len(listnomat)>0:
            summaryfile.write("#Precursors without annotated mature; Successfully predicted mature---------------------------\n")
            tempcountsucnomat=0
            for k in range(0,len(listnomat)):
                if listnomat[k] not in listnomatbroken and listnomat[k] not in listnomatscore:
                    tempcountsucnomat=tempcountsucnomat+1
                    summaryfile.write(str(listnomat[k].strip())+"\n")

            if tempcountsucnomat==0:
                summaryfile.write("---> NO Precursors without annotated mature; Successfully predicted mature\n")

            summaryfile.write("\n")
            summaryfile.write("#Precursors without annotated mature; bad positioned predicted mature---------------------------\n")

            if len(listnomatbroken)>0:
                for k in listnomatbroken:
                    summaryfile.write(k+"\n")

            else:
                summaryfile.write("--->NO Precursors without annotated mature; bad positioned predicted mature:\n")

            summaryfile.write("\n")
            summaryfile.write("#Precursors without annotated mature; no similar mature---------------------------\n")

            if len(listnomatscore)>0:
                for k in listnomatscore:
                    summaryfile.write(k+"\n")

            else:
                summaryfile.write("--->NO Precursors without annotated mature; no similar mature\n")

        else:
            summaryfile.write("--->All precurssors have annotated matures\n")

        summaryfile.write("\n")
        summaryfile.write("---------------------------Original precursors totally removed ----------------------------\n")

        if len(listofboth)>0:
            for i in range(0,len(listofboth)):
                if i%2==0:
                    tempsplit=listofboth[i].split()
                    summaryfile.write(str(tempsplit[1].strip())+"\n")

        else:
            summaryfile.write("---> NO Original precursors totally removed ----------------------------\n")

        summaryfile.write("\n")
        summaryfile.write("---------------------------Precursors without given genome file(s)----------------------------\n")

        if len(listnogenomes)>0:
            for i in listnogenomes:
                summaryfile.write(i.strip()+"\n")

        else:
            summaryfile.write("---> All the precursors have given genomes to search----------------------------\n")

        summaryfile.write("\n")
        summaryfile.write("---------------------------Precursors NOT found in their given genomes file(s)----------------------------\n")

        if len(listnotingenome)>0:
            for i in listnotingenome:
                summaryfile.write(i.strip()+"\n")

        else:
            summaryfile.write("--->All the precursors were found in their given genomes file(s)----------------------------\n")


        #Total
        Processed=(len(listofold)/2)+(len(list2mat)/3)+(len(listofoldloop)/2)+(len(listofnew)/2)+(len(listofnewloop)/2)#+len(listnomat)
        Removed=len(listofboth)/2
        TotalnumberofSequences=Processed+Removed
        #no mats/predicted-no prediction
        predicted=tempcountsucnomat
        noprediction=len(listnomat)-tempcountsucnomat

        #Corrected/notcorrected
        flippednotcorrected=int((len(listofnew)/2))+int((len(listofnewloop)/2))-int(countcorrected)
        flippedcorrected=int(countcorrected)+int(countcorrectedTonew)


        #sequences distribution
        with2mats=len(list2mat)/3
        with0mats=len(listnomat)
        with1mats=TotalnumberofSequences-with2mats-with0mats

        #Shanon
        oldshanon=OldShanon/10
        newshanon=NewShanon/10

        pdf=plt.figure(figsize=(6,6))
        N = 5

        tProcessed+=Processed#all start with 't', are for the total of all families
        tRemoved+=Removed
        tTotalnumberofSequences+=TotalnumberofSequences
        tpredicted+=predicted
        tnoprediction+=noprediction
        tflippednotcorrected+=flippednotcorrected
        tflippedcorrected+=flippedcorrected
        twith2mats+=with2mats
        twith0mats+=with0mats
        twith1mats+=with1mats
        toldshanon+=oldshanon
        tnewshanon+=newshanon
        tcountcorrected=int(countcorrected)+int(countcorrectedTonew)+tcountcorrected
        #Total
        Totalproc = [Processed, 0, 0, 0, 0]
        Totalrem = [Removed, 0, 0, 0, 0]
        #nomats
        Nomatspred = [0,0,0,predicted,0]
        Nomatsnotpred = [0,0,0,noprediction,0]
        #flipped
        flippednotcorr=[0,0,flippednotcorrected,0,0]
        flippedcorr=[0,0,flippedcorrected,0,0]
        #sequences
        seq0mat=[0,with0mats,0,0,0]
        seq1mat=[0,with1mats,0,0,0]
        seq2mat=[0,with2mats,0,0,0]

        #Shanon
        newshan=[0,0,0,0,newshanon]
        oldshan=[0,0,0,0,oldshanon]
        suma=[]
        sumall=[]

        for i in range(0,len(Totalproc)):
            suma.append(seq0mat[i]+seq1mat[i])

        sumall.append(Totalproc[0])
        sumall.append(Totalrem[0])
        sumall.append(Nomatspred[3])
        sumall.append(Nomatsnotpred[3])
        sumall.append(flippednotcorr[2])
        sumall.append(flippedcorr[2])
        sumall.append(seq0mat[1])
        sumall.append(seq1mat[1])
        sumall.append(seq2mat[1])
        sumall.append(newshan[4])
        sumall.append(oldshan[4])

        ind = np.arange(N)    # the x locations for the groups
        width = 0.35       # the width of the bars: can also be len(x) sequence

        pdf, ax = plt.subplots()

        if Totalrem[0]==0:
            Totalrem[0]=0.15
            p1 = ax.bar(ind+width, Totalrem, width, color='#FFFFFF',edgecolor='black')#,tick_label='Total')
        else:
            p1 = ax.bar(ind+width, Totalrem, width, color='#FF0000')#,tick_label='Total')
        if Totalproc[0]==0:
            Totalproc[0]=0.15
            p2 = ax.bar(ind+width*2, Totalproc, width, color='#FFFFFF',edgecolor='black')
        else:
            p2 = ax.bar(ind+width*2, Totalproc, width,bottom=0,color='#BB0000')#,'grey','blue','yellow','lime'])

        #COL2/seq-distribution
        if seq0mat[1]==0:
            seq0mat[1]=0.15
            p3 = ax.bar(ind+width, seq0mat, width,color='#FFFFFF',edgecolor='black')
        else:
            p3 = ax.bar(ind+width, seq0mat, width,color='#660000')

        if seq1mat[1]==0:
            seq1mat[1]=0.15
            p4 = ax.bar(ind+width*2, seq1mat, width,bottom=0,color='#FFFFFF',edgecolor='black')
        else:
            p4 = ax.bar(ind+width*2, seq1mat, width,bottom=0,color='#669900')#,tick_label='Sequences\n Distribution')

        if seq2mat[1]==0:
            seq2mat[1]=0.15
            p5 = ax.bar(ind+width*3, seq2mat, width,bottom=0,color='#FFFFFF',edgecolor='black')
        else:
            p5 = ax.bar(ind+width*3, seq2mat, width,bottom=0,color='#996600')

        #col3/flipped
        if flippednotcorr[2]==0:
            flippednotcorr[2]=0.15
            p6 = ax.bar(ind+width*2, flippednotcorr, width,color='#FFFFFF',edgecolor='black')#,tick_label='Changed\n sequences')
        else:
            p6 = ax.bar(ind+width*2, flippednotcorr, width,color='#0000FF')#,tick_label='Changed\n sequences')

        if flippedcorr[2]==0:
            flippedcorr[2]=0.15
            p7 = ax.bar(ind+width*3, flippedcorr, width,bottom=0,color='#FFFFFF',edgecolor='black')
        else:
            p7 = ax.bar(ind+width*3, flippedcorr, width,bottom=0,color='#00CCFF')

        #col4/nomats
        if Nomatspred[3]==0:
            Nomatspred[3]=0.15
            p8 = ax.bar(ind+width*2, Nomatspred, width,color='#FFFFFF',edgecolor='black')#,tick_label='Precursors\n without \n mature(s)')
        else:
            p8 = ax.bar(ind+width*2, Nomatspred, width,color=['#7D7D7D'])#,tick_label='Precursors\n without \n mature(s)')

        if Nomatsnotpred[3]==0:
            Nomatsnotpred[3]=0.15
            p9 = ax.bar(ind+width*3, Nomatsnotpred, width,bottom=0,color='#FFFFFF',edgecolor='black')
        else:
            p9 = ax.bar(ind+width*3, Nomatsnotpred, width,bottom=0,color='#0D0D0D')#,tick_label='Total')

        #col5/shanon
        if newshan[4]==0:
            newshan[4]=0.15
            p10 = ax.bar(ind+width*2, newshan, width,color='#FFFFFF',edgecolor='black')
        else:
            p10 = ax.bar(ind+width*2, newshan, width,color='#003300')#,tick_label='Alignment\n Entropy')

        if oldshan[4]==0:
            oldshan[4]=0.15
            p11 = ax.bar(ind+width*3, oldshan, width,bottom=0,color='#FFFFFF',edgecolor='black')
        else:
            p11 = ax.bar(ind+width*3, oldshan, width,bottom=0,color='#009900')

        ax.set_title('Statistics for file (family)-'+str(filename).strip())
        ax.set_xticks(ind + width+(width))
        ax.set_ylabel('Number / Entropy')
        ax.set_title('Statistics for file (family)-'+str(filename).strip())
        ax.set_xticklabels(('Total', 'Sequences\n Distribution', 'Changed\n sequences', 'Precursors\n without \n mature(s)', 'Alignment\n Entropy'))
        if max(sumall)>25:
                offset=int(max(sumall)/25)
                printlog(['maxsum',max(sumall),sumall])
        else:
                offset=1

        ax.set_yticks(np.arange(0, max(sumall)+5, offset))
        ax.legend((p1[0], p2[0],p3[1],p4[1],p5[2],p6[2],p7[3],p8[3],p9[3],p10[4],p11[4]), ('Removed', 'Processed *','Without mature','With 1 mature','With 2 matures','Changed Not Corrected','Corrected misaligned','Predicted','No mature Predicted','New Entropy','Old Entropy'))#,title="White bars refer to zero/doesn't exist")

        pdf.savefig(outdir+filename.strip()+'statistics.pdf')
        plt.close('all')
        summaryfile.write("---------------------------Results In Numbers----------------------------\n")
        summaryfile.write("*Number of remained precursors= "+str(int((len(listofold)/2)+(len(list2mat)/3)))+"\n")
        summaryfile.write("*Number of remained precursors with bad positioned matures= "+str(int((len(listofoldloop)/2)))+"\n")
        summaryfile.write("*Number of flipped(changed) precursors= "+str(int((len(listofnew)/2)))+"\n")
        summaryfile.write("*Number of flipped(changed) precursors with bad positioned matures= "+str(int((len(listofnewloop)/2)))+"\n")
        summaryfile.write("*Number of removed precursors= "+str(int(len(listofboth)/2))+"\n")
        summaryfile.write("*Number of precursors without a given matures= "+str(int(len(listnomat)))+"\n")
        summaryfile.write("*Number of precursors with successfully predicted matures= "+str(int(len(listnomat)-len(listremovedbroken)-len(listremovedscore)))+"\n")
        summaryfile.write("*Number of precursors without a given genome file= "+str(int(len(listnogenomes)))+"\n")
        summaryfile.write("*Number of precursors not found in their given genomes= "+str(int(len(listnotingenome)))+"\n")
        summaryfile.write("---------------------------Numbers Used For The Graph----------------------------\n")
        summaryfile.write("Total number of sequences="+str(int(TotalnumberofSequences))+"\n")
        summaryfile.write("Processed="+str(int(Processed))+"\n")
        summaryfile.write("Removed="+str(int(Removed))+"\n")
        summaryfile.write("Predicted="+str(int(predicted))+"\n")
        summaryfile.write("without given mature/no predicted mature="+str(int(noprediction))+"\n")
        summaryfile.write("Changed(flipped)="+str(int(flippednotcorrected))+"\n")
        summaryfile.write("Misaligned (shifted) precursors corrected at the end="+str(int(flippedcorrected))+"\n")
        summaryfile.write("Number of precursors without annotated mature(s)="+str(int(with0mats))+"\n")
        summaryfile.write("Number of precursors with one annotated mature="+str(int(with1mats)-int(with0mats))+"\n")
        summaryfile.write("Number of precursors with two annotated matures="+str(int(with2mats))+"\n")
        summaryfile.write("Old Entropy="+str(oldshanon)+"\n")
        summaryfile.write("New Entropy="+str(newshanon)+"\n")
        summaryfile.close()

        listofoldloopjson=[]
        for k in range(0,len(listofoldloop)):
            if k%2==0:
                tempsplit=listofoldloop[k].split()
                listofoldloopjson.append(str(tempsplit[1].strip()))

        listofnewjson=[]
        for k in range(0,len(listofnew)):
            if k%2==0:
                tempsplit=listofnew[k].split()
                listofnewjson.append(str(tempsplit[1].strip()))

        listofnewloopjson=[]
        for k in range(0,len(listofnewloop)):
            if k%2==0:
                tempsplit=listofnewloop[k].split()
                listofnewloopjson.append(str(tempsplit[1].strip()))

        listsuccpred=[]
        for k in range(0,len(listnomat)):
            if listnomat[k] not in listnomatbroken and listnomat[k] not in listnomatscore:
                listsuccpred.append(str(listnomat[k].strip()))

        listremovedjson=[]
        for i in range(0,len(listofboth)):
            if i%2==0:
                tempsplit=listofboth[i].split()
                listremovedjson.append(str(tempsplit[1].strip()))

        data = {
            "Processed":int(Processed),
            "Precursors totally removed":{
                "Number":int(len(listremovedjson)),
                "IDs":listremovedjson,
            },
            "Remained precursors":int((len(listofold)/2)+(len(list2mat)/3)),
            "Remained precursors with bad positioned matures":{
                "Number":int(len(listofoldloopjson)),#int((len(listofoldloop)/2)),
                "IDs":listofoldloopjson,
            },
            "Flipped(changed) precursors":{
                "Number":int(len(listofnewjson)),
                "IDs":listofnewjson,
            },
            "Flipped(changed) precursors with bad positioned matures":{
                "Number":int(len(listofnewloopjson)),
                "IDs":listofnewloopjson,
            },
            "Flipped/changed precursors changed back to original":{
                "Number":int(len(listcorrected)),
                "IDs":listcorrected,
            },
            "Flipped/changed precursors changed back to original":{
                "Number":int(len(listcorrected)),
                "IDs":listcorrected,
            },
            "Original precursors, Change at the Alignment level":{
                "Number":int(len(listcorrectedori)),
                "IDs":listcorrectedori,
            },
            "Precursors without a given matures":{
                "Number":int(len(listnomat)-len(listremovedbroken)-len(listremovedscore)),
                "IDs":listnomat,
            },
            "Precursors with successfully predicted matures":{
                "Number":int(len(listsuccpred)),
                "IDs":listsuccpred,
            },
            "Precursors with bad positioned predicted matures":{
                "Number":int(len(listnomatbroken)),
                "IDs":listnomatbroken,
            },
            "Precursors without annotated and without predicted mature":{
                "Number":int(len(listnomatscore)),
                "IDs":listnomatscore,
            },
            "Precursors without a given genome":{
                "Number":int(len(listnogenomes)),
                "IDs":listnogenomes,
            },
            "Precursors not found in their genomes":{
                "Number":int(len(listnotingenome)),
                "IDs":listnotingenome,
            },
            "Old Entropy":oldshanon,
            "New Entropy":newshanon,
        }
        jsonData = json.dumps(data)

        with open(outdir+filename.strip()+'.json', 'w') as f:
            json.dump(jsonData, f)

        Processed=0
        Removed=0
        TotalnumberofSequences=0
        predicted=0
        noprediction=0
        flippednotcorrected=0
        flippedcorrected=0
        with2mats=0
        with0mats=0
        with1mats=0
        oldshanon=0
        newshanon=0
        countcorrected=0
        countcorrectedTonew=0
        del Totalproc[:]
        del Totalrem[:]
        del Nomatspred[:]
        del Nomatsnotpred[:]
        del flippednotcorr[:]
        del flippedcorr[:]
        del seq0mat[:]
        del seq1mat[:]
        del seq2mat[:]
        del newshan[:]
        del oldshan[:]
        del suma[:]
        del sumall[:]
        del listofnew[:]#add seq id and new sequence
        del listofnewloop[:]#add seq id and new sequence
        del listofoldloop[:]#add seq id and old seq
        del listofold[:]#add seq id and old seq
        del listofboth[:]#add the id, that shouldn't be added to the new file
        del listremovedbroken[:]
        del listremovedscore[:]
        del listnomat[:]
        del list2mat[:]
        del listmatcoor[:]
        del listnogenomes[:]#list of IDs that no existing genome to search, for their species
        del listnotingenome[:]

        familyfileres.close()
        printlog("done")
        printlog([listofold,listofnew,listofnewloop,listofoldloop,listremovedbroken,listremovedscore,listofmirstar])
        printlog(list2mat)
        if os.path.isfile(outdir+filename.strip()+"-res.fa"):
                fr1=os.popen("rm "+outdir+filename.strip()+"-res.fa")
                
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno
            print >>h, "*** Family:", filename

def eprint(msg):
    with open('error', 'a') as e:
        print >>e, 'ERROR: '+str(msg)

def printlog(msg):
    with open('log', 'a') as l:
        print >>l, str(msg)

def openfile(f):
    try:
        return open(f,'r') if not '.gz' in f[-4:] else gzip.open(f,'rb')
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno

def runprogram():
    try:
        #nthreads=int(sys.argv[1])
        if threads.get()=="Please Choose The number of Flanking Nucleotides (upstream/downstream the 5'/3' matures)":
                
            nthreads=int(1)
        else:
            nthreads=int(threads.get())
        #filelist=str(sys.argv[4])
        filelist=str(E2.get())
        lstfams = open(filelist,'r') if not '.gz' in filelist else gzip.open(filelist,'r')
        #outd=str(sys.argv[2])
        outd=str(E3.get())
        lfams=[]

        
        
        for line in lstfams:
            lfams.append(line.strip())

        printlog(lfams)
        pool = multiprocessing.Pool(processes=nthreads, maxtasksperchild=1)

        global RNA
        RNA = importlib.import_module('RNA')
        globals().update(
            {n: getattr(RNA, n) for n in RNA.__all__}
            if hasattr(RNA, '__all__')
            else {k: v for (k, v) in RNA.__dict__.items() if not k.startswith('_')
            })
        md = RNA.md()
        md = None

        find_executable('clustalw') or sys.exit('Please install clustalw to run this')
        find_executable('dialign2-2') or sys.exit('Please install dialign2-2 to run this')

        for fam in lfams:
            pool.apply_async(sublist, args=(fam,))

        pool.close()
        pool.join()

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        with open('error','a') as h:
            tb.print_exception(exc_type, exc_value, exc_tb,
                               limit=20, file=h)
            print >>h, "*** tb_lineno:", exc_tb.tb_lineno
            root.destroy()
    root.destroy()
##############################MAIN##############################

if __name__ == '__main__':
    
    button_opt = {'fill': Tkconstants.BOTH, 'padx': 5, 'pady': 5}
    dirbut= Button(root, text = 'Select output directory', fg = 'blue', command= outdirectory)
    dirbut.pack(**button_opt) ## must pack separately to get the value to dirbut   

    fdirbut= Button(root, text = 'Please select the directory of the fasta file(s) in the list', fg = 'blue', command= filesdirectory)
    fdirbut.pack(**button_opt) ## must pack separately to get the value to dirbut 

    broButton = Button(master = root, text = 'Browse the fasta file(s) list',fg='blue' ,  command=browse_file_list)
    broButton.pack(**button_opt)

    genbutton = Button(master = root, text = 'Browse Genomes',fg='blue' ,  command=browse_genomes)
    genbutton.pack(**button_opt)

    mapbutton = Button(master = root, text = 'Browse Mapping File',fg='blue' ,  command=browse_mappingfile)
    mapbutton.pack(**button_opt)

    maturebutton = Button(master = root, text = 'Browse Mature File',fg='blue' ,  command=browse_maturefile)
    maturebutton.pack(**button_opt)

    matricesbutton = Button(master = root, text = 'Browse the Matrices directory',fg='blue' ,  command=browse_matrices)
    matricesbutton.pack(**button_opt)

    #submit_file = Button(root, text ="Submit file", command = subfile)
    #submit_file.pack(side=LEFT, padx = 6, pady=6)

    submit_list = Button(root, text ="Submit list", command = runprogram)
    submit_list.pack(side=LEFT, padx = 6, pady=6)

    OPTIONS = []
    for i in range(0,51):
        OPTIONS.append(i)
    flanking = StringVar(root)
    #flanking.set(10) # default value
    flanking.set("Please Choose The number of Flanking Nucleotides (upstream/downstream the 5'/3' matures)")
    droplist = OptionMenu(root,flanking, *OPTIONS)
    droplist.config(fg="BLUE")
    #w.pack(side=TOP, padx = 0, pady=0)
    droplist.place(x = 0, y=280)
    OPTION = []
    for i in range(0,51):
        OPTION.append(i)
    threads = StringVar(root)
    #flanking.set(10) # default value
    threads.set("Please Choose The number of files to run in parallel")
    droplist1 = OptionMenu(root,threads, *OPTION)
    droplist1.config(fg="BLUE")
    droplist1.place(x = 750, y=280)
    E1 = Entry(root, bd =5)#genomes 
    E2 = Entry(root, bd =5)#file or list
    E3 = Entry(root, bd =5)#output directory
    E4 = Entry(root, bd =5)#mapping file
    E5 = Entry(root, bd =5)#matrure file
    E6 = Entry(root, bd =5)#files directory
    E7 = Entry(root, bd =5)#matrices directory
    root.mainloop()
