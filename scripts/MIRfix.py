#!/usr/bin/env python3

# import modules
import os
import argparse
import sys
import re
import shlex
import multiprocessing
import json
import traceback as tb
from distutils.spawn import find_executable
# import Bio modules
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
# Indexing
from pyfaidx import Fasta
# Logging
import logging
from lib.logger import makelogdir, makelogfile, listener_process, listener_configurer, worker_configurer
# load own modules
from lib.Collection import *

# Load in future devs
# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
# matplotlib.use('Agg')

log = logging.getLogger(__name__)  # use module name
scriptname = os.path.basename(__file__).replace('.py', '')


def getindex(sequence, specie, precID, precdesc, listnogenomes, listnotingenome, templong, args):#get the index of the original sequence in its genome
    logid = scriptname+'.getindex: '
    try:
        specieitem=specie.split()
        listofgenomes=[]
        lstgenomes = openfile(args.genomes)
        flaggenome=0  # if there is a genomes for the specie, then it is 1
        flagseq=0  # if the sequence found in its genome, then it is 1
        minusstrand=False
        for lines in lstgenomes:
            line=lines.strip()
            if specieitem[0].lower() in line.lower() and specieitem[1].lower() in line.lower():
                flaggenome=1
                listofgenomes.append(lines.strip())

        if flaggenome==0:
            listnogenomes.append(precID)
            returnlst=[]
            return returnlst, listnogenomes, listnotingenome, templong, minusstrand
        if len(listofgenomes)>0:
            returnlst=[]
            for gen in listofgenomes:  # List of available genomes
                extension = 250  # Extension
                if specieitem[0] in gen:
                    log.debug(["Location hairpin in genome", specieitem, specieitem[0], gen])
                    (longseq, precind, chr, minusstrand) = find_precursor_genome(precID, sequence, gen, extension, args.outdir)
                    if precind > 0:
                        log.debug(["in genome", precID])
                        flagseq = 1
                        templong.append(precID.strip())
                        templong.append(str(longseq))
                        returnlst.append(precind)
                        returnlst.append(str(chr))
                        returnlst.append(str(gen))
                        return returnlst, listnogenomes, listnotingenome, templong, minusstrand

        if flagseq==0 and flaggenome==1:
            # Sequence not found on genomes
            log.debug(["Not found in available genomes", precID])
            listnotingenome.append(precID)
            returnlst=[]
            return returnlst, listnogenomes, listnotingenome, templong, minusstrand

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def getindex2mat(sequence, specie, precID, precdesc, listnogenomes, listnotingenome, args):  # get the index of the original sequence in its genome
    logid = scriptname+'.getindex2mat: '
    try:
        specieitem=specie.split()
        listofgenomes=[]
        lstgenomes = openfile(args.genomes)
        flaggenome=0  # if there is a genomes for the specie, then it is 1
        flagseq=0  # if the sequence found in its genome, then it is 1

        for lines in lstgenomes:
            line=lines.strip()
            if specieitem[0].lower() in line.lower() and specieitem[1].lower() in line.lower():
                flaggenome=1
                listofgenomes.append(lines.strip())

        if flaggenome==0:
            listnogenomes.append(precID)
            return "", listnogenomes, listnotingenome

        if len(listofgenomes)>0:
            log.debug(["list of genomes>0:", sequence])
            for gen in listofgenomes:
                extension = 250
                if specieitem[0] in gen:
                    log.debug(["Location hairpin in genome", specieitem, specieitem[0], gen])
                    (longseq, precind, chr, minusstrand) = find_precursor_genome(precID, sequence, gen, extension, args.outdir)
                    if precind > 0:
                        log.debug(["in genome", precID])
                        flagseq=1
                        return (str(longseq)), listnogenomes, listnotingenome

        if flagseq==0 and flaggenome==1:
            log.debug(["Not found in available genomes", precID])
            listnotingenome.append(precID)
            return "", listnogenomes, listnotingenome

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def extend_sequence(precursor, indexed_genomes, extension, sense, workfolder, specie):
    temp_fasta = workfolder+"temporal_seq.fa"
    temporal = open(temp_fasta, "w")
    print(">target_sequence\n"+precursor, file=temporal)
    temporal.close()
    short_specie = specie.split()[0]
    for genome in indexed_genomes:
        genome = genome.rstrip()
        # Check only genomes that are related with precursor specie
        if short_specie in genome:
            log.debug(["validation genomes", short_specie, specie, genome])
            out_blast = workfolder+"temporal_blast.txt"
            cline = NcbiblastnCommandline(query=temp_fasta, db=genome,
                                          perc_identity=100.0,
                                          qcov_hsp_perc=100.0,
                                          word_size=40, max_target_seqs=5,
                                          out=out_blast, outfmt=6,
                                          num_threads=4)
            cline()
            filesize = os.path.getsize(out_blast)
            if filesize > 0:
                with open(out_blast) as filename:
                    #  target_sequence    chr12    100.000    57    0    0    1    57    8956681    8956737    3.96e-23    106
                    first_line = filename.readline().rstrip()
                    fields = first_line.split()
                    chr = fields[1]
                    start = fields[8]
                    end = fields[9]
                    if start <= end:
                        strand = 1
                    else:
                        start = fields[9]
                        end = fields[8]
                        strand = -1
                # Add additional nt to complete adjacent nt
                if sense == "5p":
                    start = int(start) - int(extension)
                    end = int(end) + int(extension)
                else:
                    start = int(start) - int(extension)
                    end = int(end) + int(extension)
                if start < 1:
                    start = 1
                # index sequence
                get_genome = Fasta(genome)
                # Get forward precursor
                if strand == 1:
                    new_sequence = get_genome[chr][start-1:end].seq
                else:
                    # Here reverse one
                    new_sequence = get_genome[chr][start-1:end].reverse.complement.seq
                os.remove(out_blast)
                os.remove(temp_fasta)
                new_sequence = new_sequence.replace("T", "U")
                return new_sequence


def find_precursor_genome(id, precursor, genome, extension, workfolder):
    temp_fasta = workfolder+"temporal_search_seq.fa"
    temporal = open(temp_fasta, "w")
    id = id.strip()
    print(">"+str(id)+"\n"+precursor, file=temporal)
    temporal.close()
    genome = genome.rstrip()
    out_blast = workfolder+"temporal_blast_genome.txt"
    cline = NcbiblastnCommandline(query=temp_fasta, db=genome, strand="both",
                                  perc_identity=100.0, qcov_hsp_perc=100.0,
                                  word_size=50, max_target_seqs=5,
                                  out=out_blast, outfmt=6,
                                  num_threads=4)
    cline()
    filesize = os.path.getsize(out_blast)
    minusstrand = False
    if filesize == 0:
        return (0, 0, 0, minusstrand)
    else:
        num_lines = sum(1 for line in open(out_blast))
        if num_lines > 1:
            log.warning("Sequence "+id+" has paralogs on "+genome+" genome. I will process first result")
        with open(out_blast) as filename:
            #  target_sequence    chr12    100.000    57    0    0    1    57    8956681    8956737    3.96e-23    106
            first_line = filename.readline().rstrip()
            fields = first_line.split()
            chr = fields[1]
            start = int(fields[8])
            start_genome = int(fields[8])
            end = int(fields[9])
            if start < end:
                strand = 1
            else:
                start = fields[9]
                start_genome = int(fields[9])
                end = fields[8]
                strand = -1
        # Add additional nt to complete adjacent nt
        start = int(start) - int(extension)
        end = int(end) + int(extension)
        if start < 1:
            start = 1
        # index sequence
        get_genome = Fasta(genome)
        # Get forward precursor
        if strand == 1:
            new_sequence = get_genome[chr][start-1:end].seq
            minusstrand = False
        else:
            # Here reverse one
            new_sequence = get_genome[chr][start-1:end].reverse.complement.seq
            minusstrand = True
        os.remove(out_blast)
        os.remove(temp_fasta)
        return (new_sequence, start_genome, chr, minusstrand)


def flip(filename, filen, outdir, mappingfile, matfile, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, list2mat, listnogenomes, listnotingenome, templong, listgoodnew, indexed_genomes, args):  # file name is the family name
    logid = scriptname+'.flip: '
    try:
        item=[]
        item1=[]
        precitem=[]
        returnlst=[]
        precfile = openfile(str(filen))
        log.debug(["here prec", filen])
        for record in SeqIO.parse(precfile, 'fasta'):  # open the oriented file of a specific family
            precDes=record.description
            precseq=str(record.seq)
            precitem=precDes.split()
            specie=precitem[2].strip()+" "+precitem[3]
            precID=precitem[1].strip()
            flagprec=0
            flagmat=0

            if precID not in list2mat:
                mfi = openfile(mappingfile)
                for famline in mfi:  # open the mapping file, to get the mature sequence ID, by the precursor ID in given family
                    famlinesplit=famline.split()
                    if precID in famlinesplit:  # and not in 2mat list
                        log.debug(["prec is here", precID])
                        flagprec=1
                        matID=famlinesplit[4].strip()  # mature sequence ID
                        log.debug(["mat is here", matID])
                        mtf = openfile(matfile)
                        for mat in SeqIO.parse(mtf, "fasta"):  # get the corresponding mature sequence and get position
                            if matID in str(mat.description):
                                log.debug(["mat is here1", matID])
                                matdesc=mat.description
                                flagmat=1
                                log.debug(["mat desc is here", matdesc])
                                matseq=str(mat.seq).strip()
                                matseq=matseq.replace('T', 'U')  # replace T by U because the RNAfold produces the sequences as U
                                spos=str(precseq).find(matseq)
                                epos=spos+len(matseq)-1
                                break

                xcut=len(precseq[:spos])
                ycut=len(precseq[epos+1:])
                log.debug(["cut", xcut, ycut])

                # Here be sure that the surrounding is 10nt always.
                precseqOLD = precseq  # Save a reference in case it is not located on the genome
                if xcut>ycut and ycut>10:  # => 3p cut from the end
                    precseq=precseq[:epos+11]
                elif xcut>ycut and ycut==10:  # => 3p and no need to cut, already <=10
                    precseq=precseq
                elif xcut>ycut and ycut<10:  # => 3p and no need to cut, already <=10
                    # Extent to 10 the sequence
                    diff = 10 - ycut
                    precseq=precseq.replace("U", "T")  # r
                    precseq = extend_sequence(precseq, indexed_genomes, diff, "3p", outdir, specie)
                    if precseq is None:
                        precseq = precseqOLD
                elif xcut<ycut and xcut>10:  # => 5p cut at the top
                    precseq=precseq[spos-10:]
                elif xcut<ycut and xcut==10:  # => 5p and no need to cut, already <=10

                    precseq=precseq
                elif xcut<ycut and xcut<10:  # => 5p and no need to cut, already <=10
                    # Extend to 10 the sequence
                    diff = 10 - xcut
                    precseq=precseq.replace("U", "T")  # r
                    precseq = extend_sequence(precseq, indexed_genomes, diff, "5p", outdir, specie)
                    if precseq is None or precseq == "":
                        precseq = precseqOLD

                precseq=precseq.replace("U","T") # r
                matseq=matseq.replace("U","T") # r
                spos=str(precseq).find(matseq)  # spos after cut
                epos=spos+len(matseq)-1  # epos after cut
                matseq=matseq.replace('T','U')
                returnlst, listnogenomes, listnotingenome, templong, minusstrand=getindex(precseq, specie, precID, precDes, listnogenomes, listnotingenome, templong, args)  # returns 3 values received, the first is the index of the sequence, the ID where this sequence found in the genome and the genome filename

                if precID not in listnogenomes and precID not in listnotingenome:
                    listnewold=[]
                    precepos=returnlst[0]+len(precseq)  # Get the end position of the prec by adding the length of the prec to the  index where the precursor found in the genome
                    x=0
                    m=0
                    sm=returnlst[0]+spos  # start position of the mature sequence in the original precurson in the "GENOME"
                    em=sm+len(matseq)  # end position of the mature sequence in the original precurson in the "GENOME"
                    x=sm-returnlst[0]  # size/number of nucleotides before mature
                    y=precepos-em  # size/number of nucleotides after mature
                    m=len(matseq)  # mature sequence length
                    nx=sm-y  # to replace X by Y
                    ny=em+x  # to replace Y by X
                    newx=sm-nx  # Calculate the new X, the new size/number of nucleotides before mature
                    newy=ny-em  # Calculate the new Y, the new size/number of nucleotides after mature
                    newspos=newx  # new start position of mature
                    newepos=newx+m-1  # new end position of mature
                    log.debug(["mat= ", matdesc, matseq])
                    if nx < 1:
                        nx = 1
                    genome_ref = Fasta(returnlst[2])
                    chr_ref = returnlst[1]
                    if minusstrand is False:
                        newseq = genome_ref[chr_ref][nx-1:ny].seq
                        newseq1 = newseq
                    elif minusstrand is True:
                        newseq = genome_ref[chr_ref][nx-1:ny].reverse.complement.seq
                        newseq1 = newseq
                    oldlstlstr=[]
                    oldlstlstl=[]
                    oldparts=0
                    finaloldcomp=[]
                    listnewold=dofold(listnewold, precDes, precseq, precDes, newseq)
                    log.debug(["listnewold not in both: ", listnewold, filename, oldlstlstr, oldlstlstl, spos, epos, newspos, newepos, matdesc, matseq, outdir, oldparts, finaloldcomp, precDes, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew])
                    oldlstlstr, oldlstlstl, oldparts, finaloldcomp, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew=readfold(listnewold, filename, oldlstlstr, oldlstlstl, spos, epos, newspos, newepos, matdesc, matseq, outdir, oldparts, finaloldcomp, precDes, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew)  # newlstlstr, newlstlstl, oldlstlstr, oldlstlstl)#the orient sent to check the fold counts of the mature, the precseqold to print the full ID #The oldlstlstl oldlstlstr, are sent to save the list to the next step when checking the new sequence.

                elif precID in listnogenomes:
                    listnewold=[]
                    m=len(matseq)  # mature sequence length
                    newspos=spos   # new start position of mature
                    newepos=epos   # new end position of mature
                    log.debug(["mat= ", matdesc, matseq])
                    newseq=str(precseq)
                    newseq1=newseq
                    oldlstlstr=[]
                    oldlstlstl=[]
                    oldparts=0
                    finaloldcomp=[]
                    listnewold=dofold(listnewold, precDes, precseq, precDes, newseq)
                    log.debug(["listnewold listnogenome: ", listnewold, filename, oldlstlstr, oldlstlstl, spos, epos, newspos, newepos, matdesc, matseq, outdir, oldparts, finaloldcomp, precDes, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew])
                    oldlstlstr, oldlstlstl, oldparts, finaloldcomp, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew=readfold(listnewold, filename, oldlstlstr, oldlstlstl, spos, epos, newspos, newepos, matdesc, matseq, outdir, oldparts, finaloldcomp, precDes, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew)

                elif precID in listnotingenome:
                    listnewold=[]
                    m=len(matseq) # mature sequence length
                    newspos=spos  # new start position of mature
                    newepos=epos  # new end position of mature
                    log.debug(["mat= ", matdesc, matseq])
                    newseq=str(precseq)
                    oldlstlstr=[]
                    oldlstlstl=[]
                    oldparts=0
                    finaloldcomp=[]
                    listnewold=dofold(listnewold, precDes, precseq, precDes, newseq)
                    log.debug(["listnewold not in genome: ", listnewold, filename, oldlstlstr, oldlstlstl, spos, epos, newspos, newepos, matdesc, matseq, outdir, oldparts, finaloldcomp, precDes, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew])
                    oldlstlstr, oldlstlstl, oldparts, finaloldcomp, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew=readfold(listnewold, filename, oldlstlstr, oldlstlstl, spos, epos, newspos, newepos, matdesc, matseq, outdir, oldparts, finaloldcomp, precDes, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew)

        return listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listnogenomes, listnotingenome, templong, listgoodnew

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def readfold(listnewold, filename, oldlstlstr, oldlstlstl, spos, epos, newspos, newepos, matdesc, matseq, outdir, oldparts, finaloldcomp, precDes, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew):
    logid = scriptname+'.readfold: '
    try:
        log.debug(["readfold", listnewold])
        for i in range(0, len(listnewold)):
            if i==0 or i==4:
                log.debug(i)
                if i==0:
                    log.debug("i old")
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
                    lstlstr=[]  # list of list for right ")"
                    lstlstl=[]  # list of list for left "("
                    newlstlstr=[]  # list of list for right ")" of the new sequence
                    newlstlstl=[]  # list of list for left "(" of the new sequence
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
                    log.debug("i new")
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
                    lstlstr=[]  # list of list for right ")"
                    lstlstl=[]  # list of list for left "("
                    newlstlstr=[]  # list of list for right ")" of the new sequence
                    newlstlstl=[]  # list of list for left "(" of the new sequence
                    curr=0
                    nextl=0
                    cut=0
                    flag=0
                    newbroken=False
                    newloop=False
                    newparts=0
                    newncounts=0
                    warnlength=0

                while pos<len(st):  # to count the parts, loop in the string until we found another open "(" after a close ")", then means one more part.. in this case the pos assingned zero again, and the new st is the new part
                    parts=parts+1  # 1 part: ..((.)).. #2parts: ..((.))..(((.)))..
                    tempright=[]  # positions of the ")"
                    templeft=[]  # positions of the "("
                    for i in st:
                        realpos=cut+pos  # cut is where we find new "(" after closing ")",the real position is always the position in the original string (all string not when divided into new parts
                        if i=="(":
                            left.append(realpos)
                            curr=pos  # current position, changes in each new string (new part, in case the main strign was more than one part)
                            templeft.append(realpos)
                        elif i==")":
                            right.append(realpos)
                            tempright.append(realpos)
                            nextl=pos

                        if nextl<curr and nextl>0:  # here means that we have more one part, by finding "(" after ")"
                            endofcurrentpart=str(left[-1]-1)
                            startofnextpart=str(right[-1]+1)
                            left.pop()
                            templeft.pop()
                            realpos-=1
                            st=st[pos:]  # make the next part as new string
                            partslist.append(int(endofcurrentpart))
                            partslist.append(int(startofnextpart))
                            cut=realpos  # cut at the real position of thr whole original string
                            pos=0  # position=zero for the new string (new part)
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

                if (len(lstlstr) != len(lstlstl)) or ((")" not in hairpin) and ('(' not in hairpin)):  # checking the parts
                    flag=1

                else:
                    for i in range(len(lstlstr)):
                        if len(lstlstr[i])!=len(lstlstl[i]):  # check if the foldings not equal
                            flag=1

                if stat=="old":
                    x=spos
                    y=epos
                elif stat=="new":
                    x=newspos
                    y=newepos
                partslstlen=len(partslist)
                i=0
                # if the parts are overlapping, like: ..(((.(((((...)))))..((((..)))).)))...,,,in this case we take only the parts inside and ignore the outer foldings, in this section we correct the parts

                if warnlength==1:
                    log.debug("here warn")
                    log.debug(partslist)
                    # for i in range(0,partslstlen):
                    while i <=len(partslist)-1:
                        if i==0:
                            warnleft=[]
                            warnright=[]
                            warnhairpin=hairpin[0:(int(partslist[0])+1)]
                            log.debug(warnhairpin)
                            for k in range(0,len(warnhairpin)):
                                if hairpin[k]=="(":
                                    warnleft.append(k)
                                if hairpin[k]==")":
                                    warnright.append(k)
                            if len(warnleft)>len(warnright):
                                log.debug("R")
                                newsindex=len(warnleft)-1-len(warnright)
                                partslist.insert(0,-1)
                                partslist.insert(1,warnleft[newsindex]+1)
                                if len(partslist)>4:
                                    i=4
                                elif len(partslist)<=4:
                                    i=3
                                log.debug(partslist)
                                log.debug(i)
                            if len(warnleft)<len(warnright): # which is impossible to happen
                                log.debug("L")
                                diff=len(warnright)-len(warnleft)
                                neweindex=warnright[len(warnright)-diff]-1
                                log.debug(neweindex)
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
                            log.debug("end")
                            log.debug(partslist)
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

                            if len(warnleft)<len(warnright): # which is impossible to happen
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
                            log.debug("mid")
                            log.debug(partslist)
                            for k in range(int(partslist[i-1]), int(partslist[i])+1):
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

                            if len(warnleft)<len(warnright) and len(partslist)>4 and (int(len(partslist))-i)>2:  # which is impossible to happen
                                diff=len(warnright)-len(warnleft)
                                neweindex=warnright[len(warnright)-diff]-1
                                partslist[i]=neweindex
                                i=i+2
                            elif (len(warnleft)<len(warnright) and len(partslist)>4 and (int(len(partslist))-i)==2) or (len(warnleft)<len(warnright) and len(partslist)==4):  # which is impossible to happen
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
                # at this section the the parts list is corrected ( if needed )
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
                        if i==0 and partslist[0]!=-1:  # if the partslist starting from position zero, i.e. the first part starts at position zero
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
                            elif (((x in range (0,int(partslist[0])+1)) and (y not in list(range(0,int(partslist[0])+1)))) or ((x not in list(range(0,int(partslist[0])+1))) and (y in range (0,int(partslist[0])+1))))and stat=="old" :#if x or y out of the  hairpin
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
                            elif (((x in range (0,int(partslist[0])+1)) and (y not in list(range(0,int(partslist[0])+1)))) or ((x not in list(range(0,int(partslist[0])+1))) and (y in range (0,int(partslist[0])+1))))and stat=="new" :#if x or y out of the  hairpin
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
                            elif (((x not in list(range(int(partslist[i]),len(hairpin)))) and (y in range(int(partslist[i]),len(hairpin)))) or ((y not in list(range(int(partslist[i]),len(hairpin)))) and (x in range(int(partslist[i]),len(hairpin))))) and stat=="old":
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
                            elif (((x not in list(range(int(partslist[i]),len(hairpin)))) and (y in range(int(partslist[i]),len(hairpin)))) or ((y not in list(range(int(partslist[i]),len(hairpin)))) and (x in range(int(partslist[i]),len(hairpin))))) and stat=="new":
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
                            elif ((x not in list(range(int(partslist[i-1]),int(partslist[i])+1)) and (y in range(int(partslist[i-1]),int(partslist[i])+1))) or ((y not in list(range(int(partslist[i-1]),int(partslist[i])+1))) and (x in range(int(partslist[i-1]),int(partslist[i])+1)))) and stat=="old":
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
                            elif ((x not in list(range(int(partslist[i-1]),int(partslist[i])+1)) and (y in range(int(partslist[i-1]),int(partslist[i])+1))) or ((y not in list(range(int(partslist[i-1]),int(partslist[i])+1))) and (x in range(int(partslist[i-1]),int(partslist[i])+1)))) and stat=="new":
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

                    # all if here and elif here, are not broken when len>0 then we compare the loop
                    if len(lstcompstart)>0 and len(lstcompinside)==0 and len(lstcompend)==0:  #start only
                        finalcomp=lstcompstart
                        posi=finalcomp[2]
                        finalcomp.append(0)
                        finalcomp.append(int(partslist[0]))
                    elif len(lstcompstart)==0 and len(lstcompinside)>0 and len(lstcompend)==0:  #middle only
                        finalcomp=lstcompinside
                        posi=finalcomp[2]  # the index in partslist, that stores the real end position of a given sub-stem

                        finalcomp.append(int(partslist[posi-1]))
                        finalcomp.append(int(partslist[posi]))
                    elif len(lstcompstart)==0 and len(lstcompinside)==0 and len(lstcompend)>0:  #end only
                        finalcomp=lstcompend
                        posi=finalcomp[2]
                        finalcomp.append(int(partslist[posi]))
                        finalcomp.append(len(hairpin)-1)
                    elif len(lstcompstart)>0 and len(lstcompinside)>0:  #start + middle
                        if lstcompstart[1]==lstcompinside[1] and lstcompstart[3]>=lstcompinside[3]:  # if both true or false loops we compare the number of foldings
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
                        #CAVH
                        for i in currhairpin:
                            if i != ")": #CAVH To detect ')' at the beginning
                                index = currhairpin.index(i)
                                currhairpin=currhairpin[index:]
                                break
                        oldlstlstr=[]
                        oldlstlstl=[]
                        for nuc in range(0,len(currhairpin)):
                            if currhairpin[nuc]==")":
                                oldlstlstr.append(nuc+finalcomp[4])
                            elif currhairpin[nuc]=="(":
                                oldlstlstl.append(nuc+finalcomp[4])

                    if len(finalcomp)>0 and stat=="new":
                        currhairpin=hairpin[finalcomp[4]:finalcomp[5]+1]
                        for i in currhairpin:
                            if i != ")":  # To detect ')' at the beginning
                                index = currhairpin.index(i)
                                currhairpin=currhairpin[index:]
                                break
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

                    if (((x in range(loopstart,loopend+1)) and (y not in list(range(loopstart,loopend+1)))) or ((x not in list(range(loopstart,loopend+1))) and (y in range(loopstart,loopend+1)))) and stat=="old":
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

                    elif (((x in range(loopstart,loopend+1)) and (y not in list(range(loopstart,loopend+1)))) or ((x not in list(range(loopstart,loopend+1))) and (y in range(loopstart,loopend+1)))) and stat=="new":
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

                    log.debug(["old: oldbroken= ", oldbroken, "/oldloop= ", oldloop, "start= ", oldhairpstart, "end= ", oldhairpend, "Counts= ", oldncounts,  "old score= ", oldscore, "old parts= ", oldparts,  "old energy",  oldscore])
                    log.debug(["new: newbroken= ", newbroken, "/newloop= ", newloop, "start= ", newhairpstart, "end= ", newhairpend, "Counts= ", newncounts,  "new score= ", newscore, "new parts= ", newparts,  "new energy",  newscore])
                    temcorrectsplit=precid.split()
                    precId=temcorrectsplit[0].strip()
                    precides=temcorrectsplit[1].strip()
                    if (oldbroken==False and newbroken==False):  # To save the old when are possibly flipped into new, so no need to add old when new is broken!
                        log.debug("save old good")
                        mirstarcorr, mirstarsposcorr, mirstareposcorr, miroriencorr=getmirstar(spos, epos, matseq, oldlstlstl, oldlstlstr, oldprecseq, oldhairpstart, oldhairpend)

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
                if done:  # Check energies for both
                    if oldscore > -10 and newscore > -10:
                        log.debug("Both sequence energies are really high for a miRNA")
                        listofboth.append(precid)
                        listofboth.append(oldprecseq)
                        checkenergy = False
                    else:  #At least one of the sequences has < -10 MFE
                        checkenergy = True

                if done and checkenergy:
                    if oldbroken and newbroken:
                        log.debug("remove old and no new")
                        listofboth.append(precid)
                        listofboth.append(oldprecseq)
                    elif (not oldbroken and not oldloop and oldparts==1) and ((oldscore<=newscore) or (newbroken)):
                        log.debug("old is Better1")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("new is Better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)

                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("remove old and no new")
                        listofboth.append(precid)
                        listofboth.append(oldprecseq)
                    elif (oldscore<=newscore and not oldbroken and oldloop and not newbroken and not newloop and newparts==1 and newncounts> oldncounts and newscore<=-10.00):
                        log.debug("new is Better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("new is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old is better2")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("new is better with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("remove old and no news")
                        listofboth.append(precid)
                        listofboth.append(oldprecseq)

                    elif (newscore<oldscore and not newbroken and newloop and not oldbroken and not oldloop):
                        log.debug("old is better3")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old with warnings")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old is better4")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("new is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old is better5")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old is better6")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old is Better7")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("new is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("new is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("remove old and no new")
                        listofboth.append(precid)
                        listofboth.append(oldprecseq)

                    elif newscore<oldscore and (newbroken or newloop) and newparts>=1 and not oldbroken and not oldloop and oldparts>=1:
                        log.debug("old is better")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("new is better here")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("new warning better here")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(newspos,newepos,matseq,newlstlstl,newlstlstr,newseq,newhairpstart,newhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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
                        log.debug("old with warning")
                        mirstar,mirstarspos,mirstarepos,mirorien=getmirstar(spos,epos,matseq,oldlstlstl,oldlstlstr,oldprecseq,oldhairpstart,oldhairpend)
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug([mirstar,mirstarspos,mirstarepos])
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

                    elif (newbroken==False and newncounts>0 and newscore<=-10.00):
                        splitgoodprec=precid.split()
                        listgoodnew.append(str(splitgoodprec[0]).strip())
                        listgoodnew.append(str(splitgoodprec[1]).strip())
                        log.debug(["here look test:",newspos,newepos,matseq,newlstlstl,newlstlstr,newseq1,newhairpstart,newhairpend])
                        mirstar, mirstarspos, mirstarepos, mirorien=getmirstar(newspos, newepos, matseq, newlstlstl, newlstlstr, newseq1, newhairpstart, newhairpend)
                        log.debug(["test 1515", mirstar, mirstarspos, mirstarepos, mirorien])
                        if mirstarspos!=-1 and mirstarepos!=-1:
                            log.debug(["also here",mirstar,mirstarspos,mirstarepos])
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
                            if mirorien=='5p' and ycut>10:
                                newseq1=newseq1[:mirstarepos+11]
                            elif mirorien=='5p' and ycut<=10:
                                newseq1=newseq1
                            elif mirorien=='3p' and xcut>10:
                                newseq1=newseq1[mirstarspos-10:]
                            elif mirorien=='3p' and xcut<=10:
                                newseq1=newseq1
                            listgoodnew.append(newseq1)  # [precID[0](x-mir-y),precid(MI),matstardesc,mirseq,mirstar,start..end,newseq]
                        else:
                            del listgoodnew[-2:]

        return oldlstlstr, oldlstlstl, oldparts, finaloldcomp, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listgoodnew

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def getmirstar(spos,epos,mature,lstl,lstr,precursor,hairpstart,hairpend):
    logid = scriptname+'.getmirstar: '
    mirstarepos=0
    limit_mirmature = 19
    try:
        log.debug("get mirstar here 18")
        mirflag=False
        if epos>=lstr[0]:
            orien="3p"
        elif epos<lstr[0]:
            orien="5p"

        if ((orien == '5p' and (spos > lstl[-1] or epos>=lstr[0])) or (orien == '3p' and (spos < lstl[-1] or epos < lstr[0])) or epos<=spos):
             log.debug("-1 wrong")
             mirflag=True
             return "", -1, -1, 'p'

        if orien=="5p":
            tempr=lstr
            templ=lstl
            rev=tempr[::-1]
            log.debug("get mirstar here 19")
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
                        sind=templ.index(i)  # first ( after spos
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
                            sind=templ.index(i)  # after the end of the mat
                            sind1=templ[(sind-1)]  # before the end of the mat
                            diff=epos-sind1
                            mirstarspos=rev[sind-1]-diff+2
                            break
            if mirstarspos<=epos:  # To avoid overlapping between mir and mir*, this would happen in case the mir is in the loop
                mirstarepos=mirstarepos+(epos-mirstarspos+1)
                mirstarspos=mirstarspos+(epos-mirstarspos+1)

            if mirstarepos>len(precursor)-1:  # To avoid mir* end position going outside the precursor
                mirstarepos=len(precursor)-1  # TODO: what if after reduction mirstar is too short?
            mirstarspos = comp_5p(lstl, lstr, spos, mirstarspos, precursor, 1)  # Comparing 5'ends of mir and mir*
            mirstar=precursor[mirstarspos:mirstarepos+1]
            mirstar= mirstar.replace("T","U")  # here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
            log.debug(["get mirstar here 20",mirstar,mirstarspos,orien])
            if len(mirstar) < limit_mirmature:
                log.debug("no predicted mir*")
                return "",-1,-1,'p'
            else:
                mirflag=True
                return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))

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
                        sind=tempr.index(i)  # first ( after spos
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
                            sind=tempr.index(i)  # after the end of the mat
                            sind1=tempr[sind-1]  # before the end of the mat
                            diff=epos-sind1
                            mirstarspos=rev[sind-1]-diff+2
                            break

            if mirstarspos<=0:
                mirstarspos=0

            if mirstarepos>=spos:  # To avoid overlapping between mir and mir*, this would happen in case the mir is in the loop
                mirstarspos=mirstarspos-(mirstarepos-spos+1)
                mirstarepos=mirstarepos-(mirstarepos-spos+1)

            mirstarspos = comp_5p(lstr, lstl, spos, mirstarspos, precursor, 1)  # Comparing 5'ends of mir and mir*
            mirstar=precursor[mirstarspos:mirstarepos+1]
            mirstar=mirstar.replace("T","U")  # here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
            log.debug(["get mirstar here 24",mirstar,mirstarspos,orien])
            if len(mirstar) < limit_mirmature:
                log.debug("no predicted mir*")
                return "", -1, -1, 'p'
            else:
                mirflag=True
            return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))

        if not mirflag:
            log.debug("no predicted mir*")
            return "", -1, -1, 'p'

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def getmirstarbak(spos, epos, mature, lstl, lstr, precursor, hairpstart, hairpend):
    logid = scriptname+'.getmirstarbak: '
    try:
        log.debug("get mirstar here 18")
        mirflag=False
        if epos>=lstr[0]:
            orien="3p"
        elif epos<lstr[0]:
            orien="5p"

        if (spos in range(hairpstart,lstr[0]) and (epos>=lstr[0] or epos<spos)) or (spos in range(lstr[0],hairpend+1) and (epos<spos or epos>hairpend)):
             log.debug("-1 wrong")
             mirflag=True
             return "",-1,-1,'p'

        if orien=="5p":
            tempr=lstr
            templ=lstl
            rev=tempr[::-1]
            log.debug("get mirstar here 19")
            for i in templ:
                if i==int(spos):
                    sind=templ.index(i)  #list index(in the list of foldings), of the nucleotide representing the first fold, after the first nucleotide of the mature
                    sx=templ[sind]  #get the item at sind, which is the position in the precursor, of the first folding nucleotide after the first nucleotide of mature
                    sposstar=rev[sind]  #start position of the mature*
                    mirstarspos=sposstar-len(mature)+1
                    mirstarepos=sposstar
                    mirstar=precursor[mirstarspos+2:mirstarepos+3]
                    if mirstarspos<=epos:  #To avoid overlapping between mir and mir*, this would happen in case the mir is in the loop
                        mirstarepos=mirstarepos+(epos-mirstarspos+1)
                        mirstarspos=mirstarspos+(epos-mirstarspos+1)

                    if mirstarepos>len(precursor)-1:  #To avoid mir* end position going outside the precursor
                        mirstarepos=len(precursor)-1

                    mirstarspos = comp_5p(lstl, lstr, spos, mirstarspos+2, precursor, 1)  # Comparing 5'ends of mir and mir*
                    mirstar=precursor[mirstarspos:mirstarepos+1]
                    mirstar= mirstar.replace("T","U")  # here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
                    log.debug(["get mirstar here 20",mirstar,mirstarspos,orien])
                    mirflag=True
                    return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))
                    break

                if i>int(spos):
                    log.debug("get mirstar here 21")
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
                    log.debug(["get mirstar here 22",mirstar,mirstarspos,orien])
                    mirflag=True
                    return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))
                    break

        elif orien=="3p":
            tempr=lstr
            templ=lstl
            rev=templ[::-1]
            log.debug("get mirstar here 23")
            log.debug(rev)
            for i in tempr:
                log.debug("i "+str(i))
                if i==int(spos):
                    log.debug("i=spos "+str(spos))
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
                    log.debug(["get mirstar here 24",mirstar,mirstarspos,orien])
                    mirflag=True
                    return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))
                    break

                if i>int(spos):
                    log.debug("i>spos")
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

                    log.debug("mir*")
                    log.debug([mirstarspos,mirstarepos])
                    mirstar=precursor[mirstarspos:mirstarepos+1]#
                    mirstar=mirstar.replace("T","U")#here it is minus because we are in the 3p arm, the sposstar is actually the last nucleotide in the mir* which is the firt one folding to mir
                    log.debug(["get mirstar here 25",mirstar,mirstarspos,orien])
                    mirflag=True
                    return (str(mirstar),int(mirstarspos),int(mirstarepos),str(orien))
                    break

        if not mirflag:
            log.debug("no predicted mir*")
            return "",-1,-1,'p'

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def predict(align,matId,newmatID,matfile,filename,precdescrip,mapfile,directory,listremovedbroken,listremovedscore,listremovedcomposition):
    logid = scriptname+'.predict: '
    try:
        if directory[-1]!="/":# to work in both cases, user puts / at the end of the directory or not
            directory=str(directory)+"/"
        stockfile=alignTostock(align)
        alignment = AlignIO.read(stockfile, "stockholm")
        listrecords=[]
        for record in alignment:
            if matId == record.id:
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
                # the window is taking substring from the aligned mature line with szie of  real mature, and move nuc by nuc on the aligned prec to see how many nucleotides matching
                # First evaluate if there is not a high number of Ns in the sequence
                (proportionN, totalN, totalseq) = countNSeq(originalseq)
                if proportionN <= 40 and proportionN > -1: #Number based on RNAfold rules, have to be N content <= 40 from the complete seq
                    for i in range (0,lastindex+1):#moving along the precursor, and stop at the last index, where after that will be out of range
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
                        log.debug(["broken",famsplit[1]])
                        if famsplit[1] not in listremovedbroken:
                            listremovedbroken.append(famsplit[1].strip())
                else:
                    log.debug(["Sequence lot of Ns",precdescrip.split()[1]])
                    if precdescrip.split()[1] not in listremovedcomposition:
                        listremovedcomposition.append(precdescrip.split()[1].strip())
        return listremovedbroken,listremovedscore,listremovedcomposition

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def checknomat(precfile,mapfile,matfile,directory,precfilename,listremovedbroken,listremovedscore,listremovedN,nomats,listnomat):
    logid = scriptname+'.checknomat: '
    try:
        flagnomatexists=False
        directory=os.path.abspath(directory)+"/"
        matIDs=[]
        precnomat=[]
        countnomat=0#count number of precs without mature

        for record in SeqIO.parse(openfile(precfile), 'fasta'):
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
            log.debug("here is True")
            for i in matIDs:
                mtf = openfile(matfile)
                for record in SeqIO.parse(mtf, 'fasta'):
                    if i.strip() in record.description:
                        with open(directory+'tempmat.fa','a') as tempmaturefile:
                            tempmaturefile.write(">"+record.description+"\n"+str(record.seq)+"\n")

            for prec in SeqIO.parse(openfile(directory+'nomat-'+precfilename+'.fa'),'fasta'):
                listofmat=[]
                listofscore=[]
                splitprecdisc=prec.description.split()#because the clustal reads only id, so it will only write in thr result the id without MI00....
                prectempid=splitprecdisc[0]+"-"+splitprecdisc[1]
                listnomat.append((splitprecdisc[1]).strip())
                for mat in SeqIO.parse(openfile(directory+'tempmat.fa'),'fasta'):
                    splitmatdisc=mat.description.split()
                    mattempid=splitmatdisc[0]+"-"+splitmatdisc[1]#because the clustal reads only id, so it will only write in thr result the id without MIMAT...
                    with open(directory+'temptoalign.fa','w') as temptoalign:
                        temptoalign.write(">"+prectempid+"\n"+str(prec.seq)+"\n"+">"+mattempid+"\n"+str(mat.seq)+"\n")
                    infile=directory+'temptoalign.fa'
                    outfile=directory+'temptoalign.aln'
                    cline = ClustalwCommandline("clustalw", infile=infile, outfile=outfile)

                    stdout,stder=cline()
                    with open(directory+'temptoalign.txt','w') as fi:
                        fi.write(stdout)

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
                log.debug(["best for ",prectempsplit[1]," is ",bestMatID, "of ", listofscore, listofmat])

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

                    listremovedbroken,listremovedscore,listremovedN=predict(tempoutfilepredict,matId,newmatID,matfile,filename,prec.description,mapfile,directory,listremovedbroken,listremovedscore,listremovedN)
                    stktemptopredict=tempoutfilepredict+".stk"
                    os.remove(str(stktemptopredict))

                elif maxscore<21:
                    prectempsplit=prec.description.split()
                    log.debug(["if not score",prectempsplit[1]])
                    if prectempsplit[1] not in listremovedscore:
                        listremovedscore.append(prectempsplit[1].strip())

        elif countnomat>0 and len(matIDs)==0:
            flagnomatexists=True
            nomats=-1
            return flagnomatexists,nomats,listremovedbroken,listremovedscore,listremovedN,listnomat
        elif countnomat<=0:
            log.debug("do the normal procedure")
            #flip here
            nomats=0
            flagnomatexists=False
            log.debug("all precursors has mature, at least one")
        return flagnomatexists,nomats,listremovedbroken,listremovedscore,listremovedN,listnomat
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def readfoldpredict(foldfile,predictedspos,predictedepos):#,predictedspos,predictedepos):
    logid = scriptname+'.readfoldpredict: '
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
                partslist=[]
                parts=0
                pos=0
                realpos=0
                left=[]
                right=[]
                lstlstr=[]#list of list for right ")"
                lstlstl=[]#list of list for left "("
                curr=0
                next=0
                cut=0
                flag=0
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
                            if (((x in range(0,loopstart)) and (y in range(0,loopstart))) or ((x in range(loopend+1,int(partslist[0])+1)) and (y in range(loopend+1,int(partslist[0])+1)))):
                                log.debug("a")
                                flagvalid=1
                                finalpredictedspos=predictedspos
                                finalpredictedepos=predictedepos

                            elif (x in range(0,loopstart)) and (y in range(loopstart,loopend+1)):
                                flagvalid=1
                                shiftby=y-(loopstart-1)
                                if x<=shiftby:
                                    log.debug("b")
                                    x=0
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    log.debug("c")
                                    x=x-shiftby
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif (x in range(loopstart,loopend+1)) and (y in range(loopend+1,int(partslist[0])+1)):
                                flagvalid=1
                                shiftby=(loopend+1)-x
                                if (y+shiftby)>=(len(hairpin)-1):
                                    log.debug("d")
                                    y=len(hairpin)-1
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                    break
                                else:
                                    log.debug("e")
                                    y=y+shiftby
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))):
                                log.debug("f")
                                finalpredictedspos=-1
                                finalpredictedepos=-1
                            elif (((x in range (0,int(partslist[0])+1)) and (y not in list(range(0,int(partslist[0])+1)))) or ((x not in list(range(0,int(partslist[0])+1))) and (y in range (0,int(partslist[0])+1)))):#if x or y out of the  hairpin
                                log.debug("g")
                                finalpredictedspos=-1
                                finalpredictedepos=-1
                            elif ((x in range(0,loopstart)) and (y in range(loopend+1,int(partslist[0])+1))) or ((x in range(loopend+1,int(partslist[0])+1)) and ((y in range(0,loopstart)))):
                                log.debug("h")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                        if i==len(partslist)-1 and flagvalid!=1:
                            loopstart=hairpin.rfind('(',int(partslist[i])+1)+1
                            loopend=hairpin.find(')',int(partslist[i])+1)-1

                            if (((x in range(int(partslist[i]),loopstart)) and (y in range(int(partslist[i]),loopstart) )) or ((x in range(loopend+1,len(hairpin))) and (y in range(loopend+1,len(hairpin))))):
                                log.debug("1")
                                flagvalid=1
                                finalpredictedspos=predictedspos
                                finalpredictedepos=predictedepos
                            elif (x in range(int(partslist[i]),loopstart)) and (y in range(loopstart,loopend+1)):
                                flagvalid=1
                                shiftby=y-(loopstart-1)
                                if x<=shiftby:
                                    log.debug("2")
                                    x=0
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    log.debug("3")
                                    x=x-shiftby
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif  (x in range(loopstart,loopend+1)) and (y in range(loopend+1,len(hairpin))):
                                shiftby=(loopend+1)-x
                                if (y+shiftby)>=(len(hairpin)-1):
                                    log.debug("4")
                                    flagvalid=1
                                    y=len(hairpin)-1
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    log.debug("5")
                                    flagvalid=1
                                    y=y+shiftby
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))):
                                log.debug("6")
                                finalpredictedspos=-1
                                finalpredictedepos=-1
                            elif (((x not in list(range(int(partslist[i]),len(hairpin)))) and (y in range(int(partslist[i]),len(hairpin)))) or ((y not in list(range(int(partslist[i]),len(hairpin)))) and (x in range(int(partslist[i]),len(hairpin))))):
                                log.debug([x,y])
                                log.debug("7")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                            elif (((x in range(int(partslist[i]),loopstart)) and (y in range(loopend+1,len(hairpin)))) or ((x in range(loopend+1,len(hairpin))) and (y in range(int(partslist[i]),loopstart)))):
                                log.debug("8")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                        if i!=0 and i!=len(partslist)-1 and i%2==0 and flagvalid!=1:
                            loopstart=hairpin.rfind('(',int(partslist[i-1]),int(partslist[i])+1)+1
                            loopend=hairpin.find(')',int(partslist[i-1]),int(partslist[i])+1)-1
                            if (((x in range(int(partslist[i-1]),loopstart)) and (y in range(int(partslist[i-1]),loopstart))) or ((x in range(loopend+1,int(partslist[i])+1)) and (y in range(loopend+1,int(partslist[i])+1)))):
                                log.debug("9")
                                flagvalid=1
                                finalpredictedspos=predictedspos
                                finalpredictedepos=predictedepos

                            elif (x in range(int(partslist[i-1]),loopstart)) and (y in range(loopstart,loopend+1)):
                                flagvalid=1
                                shiftby=y-(loopstart-1)
                                if x<=shiftby:
                                    log.debug("10")
                                    x=0
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    log.debug("11")
                                    x=x-shiftby
                                    y=y-shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif  (y in range(loopend+1,int(partslist[i])+1)) and (x in range(loopstart,loopend+1)):
                                flagvalid=1
                                shiftby=(loopend+1)-x
                                if (y+shiftby)>=(len(hairpin)-1):
                                    log.debug("12")
                                    y=len(hairpin)-1
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y
                                else:
                                    log.debug("13")
                                    y=y+shiftby
                                    x=x+shiftby
                                    finalpredictedspos=x
                                    finalpredictedepos=y

                            elif ((x in range(loopstart,loopend+1)) and (y in range(loopstart,loopend+1))):
                                log.debug("14")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                            elif ((x not in list(range(int(partslist[i-1]),int(partslist[i])+1)) and (y in range(int(partslist[i-1]),int(partslist[i])+1))) or ((y not in list(range(int(partslist[i-1]),int(partslist[i])+1))) and (x in range(int(partslist[i-1]),int(partslist[i])+1)))):
                                log.debug("15")
                                finalpredictedspos=-1
                                finalpredictedepos=-1

                            elif (((x in range(int(partslist[i-1]),loopstart)) and (y in range(loopend+1,int(partslist[i])+1))) or ((x in range(loopend+1,int(partslist[i])+1)) and (y in range(int(partslist[i-1]),loopstart)))):
                                log.debug("16")
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
        log.debug([finalpredictedspos,(finalpredictedepos+1)])
        flagfinalparts=0
        if (finalpredictedspos!=-1 and finalpredictedepos!=-1 and len(partslist)>1):
            for i in range(0,len(partslist)):
                if i==0:
                    if (finalpredictedspos in range(0,int(partslist[0])+1)) and (finalpredictedepos  in range(0,int(partslist[0])+1)):
                        flagfinalparts=1
                        break

                if i!=0 and i!=len(partslist)-1 and i%2==0:
                    if  (finalpredictedspos in range(int(partslist[i-1]),int(partslist[i])+1)) and (finalpredictedepos in range(int(partslist[i-1]),int(partslist[i])+1)):
                        log.debug("here 2")
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

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def correct(corid,flanking,countcorrected,countcorrectedTonew,listofnew,listofnewloop,listoldstatus,templong,listmisalignedcorr,listcorrected,listcorrectedori,listgoodnew):
    logid = scriptname+'.correct: '
    try:
        log.debug("def correct")
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
        log.debug(listidsnewnewloop)

        if corid in listoldstatus:
            corrind=int(listoldstatus.index(corid))#the x-mir-y
            corriddes=str(listoldstatus[corrind+1])#the MI00001 id
            if corid in listidsnewnewloop and corriddes in templong:
                indexlongmat=int(templong.index(corriddes.strip()))
                longseq=str(templong[indexlongmat+1]).replace('T','U')
                log.debug(["long seq",longseq,templong[indexlongmat]])
                corrind=int(listoldstatus.index(corid))
                correctseq=str(listoldstatus[corrind+2])
                correctmir=str(listoldstatus[corrind+3])
                correctmirstar=str(listoldstatus[corrind+4])
                orien=str(listoldstatus[corrind+5])
                (coor1, coor2, correctfinalseq) = find_positions(longseq, correctseq,
                                                 correctmir, correctmirstar,
                                                 flanking)

                if coor2<coor1:
                    tempcor=correctmir[:]
                    correctmir=correctmirstar[:]
                    correctmirstar=tempcor[:]

                start_mir = correctfinalseq.find(correctmir)
                end_mir = (start_mir + len(correctmir)) - 1
                correctfinalseqnomir = correctfinalseq[end_mir + 1:]
                index_to_update = correctfinalseq.find(correctfinalseqnomir)
                listmisalignedcorr.append(corriddes)
                listmisalignedcorr.append(correctfinalseq)
                listmisalignedcorr.append(int(start_mir))
                listmisalignedcorr.append(int(end_mir))
                start_mirstar_temp = correctfinalseqnomir.find(correctmirstar)
                star_mirstar = start_mirstar_temp + index_to_update
                end_mirstar = (start_mirstar + len(correctmirstar)) - 1
                listmisalignedcorr.append(int(star_mirstar))
                listmisalignedcorr.append(int(end_mirstar))
                listmisalignedcorr.append(str(orien))
                countcorrected=countcorrected+1
                listcorrected.append(corriddes.strip())
                log.debug(["corrected seq",correctfinalseq])

        if corid.strip() in listgoodnew and corid.strip() not in listidsnewnewloop and len(listgoodnew)>2:
            corrind=int(listgoodnew.index(corid.strip()))#x-mir-y
            corriddes=listgoodnew[corrind+1]#MI00
            if corriddes.strip() in templong:
                indexlongmat=int(templong.index(corriddes.strip()))
                longseq=str(templong[indexlongmat+1]).replace('T','U')
                log.debug(["long seq",longseq,templong[indexlongmat]])
                correctseq=str(listgoodnew[corrind+6])
                correctmir=str(listgoodnew[corrind+3])
                correctmirstar=str(listgoodnew[corrind+4])
                orien=str(listgoodnew[corrind+6])
                (coor1, coor2, correctfinalseq) = find_positions(longseq, correctseq,
                                                 correctmir, correctmirstar,
                                                 flanking)

                if coor2<coor1:
                    tempcor=correctmir[:]
                    correctmir=correctmirstar[:]
                    correctmirstar=tempcor[:]

                start_mir = correctfinalseq.find(correctmir)
                end_mir = (start_mir + len(correctmir)) - 1
                correctfinalseqnomir = correctfinalseq[end_mir + 1:]
                index_to_update = correctfinalseq.find(correctfinalseqnomir)
                listmisalignedcorr.append(corriddes)
                listmisalignedcorr.append(correctfinalseq)
                listmisalignedcorr.append(int(start_mir))
                listmisalignedcorr.append(int(end_mir))
                start_mirstar_temp = correctfinalseqnomir.find(correctmirstar)
                star_mirstar = start_mirstar_temp + index_to_update
                end_mirstar = (start_mirstar + len(correctmirstar)) - 1
                listmisalignedcorr.append(int(star_mirstar))
                listmisalignedcorr.append(int(end_mirstar))
                listmisalignedcorr.append(str(orien))
                countcorrectedTonew=int(countcorrectedTonew)+1
                listcorrectedori.append(corriddes.strip())
                log.debug(["corrected seq",correctfinalseq])

        return int(countcorrected),int(countcorrectedTonew),listmisalignedcorr,listcorrected,listcorrectedori

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def index_genomes(genomes_file):
    logid = scriptname+'.index_genomes: '
    try:
        listofgenomes = []
        with open(genomes_file) as files:
            for genome in files:
                genome = str(genome.rstrip())
                # Blast db files
                # file_names_list = [genome+".ndb", genome+".nhr", genome+".nin",
                #                   genome+".not", genome+".nsq", genome+".ntf", genome+".nto"]
                file_names_list = [genome+".nin",genome+".nhr",genome+".nsq"]
                if all(list(map(os.path.isfile,file_names_list))):
                    log.debug('Using detected blast index from '+str(genome))
                    listofgenomes.append(genome)
                    continue
                else:
                    log.debug('Building blast index on '+str(genome))
                    cline = NcbimakeblastdbCommandline(dbtype="nucl",input_file=genome)
                    cline()
                    listofgenomes.append(genome)
        return listofgenomes

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def sublist(queue, configurer, level, filename, args):
    logid = scriptname+'.sublist: '
    configurer(queue, level)
    log.debug(logid+'Starting to process '+str(filename))
    try:
        filesdir=str(args.famdir)#dir for families
        command="list"
        mapfile=str(args.mapping)#mapping of mir to mirfam
        matfile=str(args.mature)#mature sequences
        matrdir=args.maturedir#directory for mature files
        genomes_file=args.genomes #file that contains all genomes
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
        listremovedN=[]
        listofmirstar=[]
        listnomat=[]
        list2mat=[]
        list2matcoor=[]
        list1matcoor=[]
        listmatcoor=[]
        listnogenomes=[]#list of IDs that no existing genome to search, for their species
        listnotingenome=[]#list of IDs that are not found in the genome searched, for their species
        listremovenoloop = [] #Collect candidates without loop after curation
        listhighmfe = [] #Collect all candidates with high MFE after all curation
        listlonghairpin = [] #Collect all candidates with a final reported hairpin
        listbadmature = [] #Collect all candidates that failed mature prediction
        nomats=0
        templong=[]
        userflanking=0
        listmisalignedcorr=[]
        listgoodnew=[]
        listcorrected=[]
        listcorrectedori=[]

        log.debug(logid+filename)
        numberoffamilies+=1
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
        listnomatremoved=[]
        flagnomatexists=False

        ## prepare output directory for current family
        filename = str( filename).strip()
        outdir = str( args.outdir) + filename + ".out/"
        makeoutdir( outdir, args.force)
        filen = filesdir + filename + ".fa"

        OldShanon=0
        NewShanon=0
        tempcountsucnomat=0
        infile=""
        outfile=""
        infile=filen[:]
        outfile=outdir+filename+"-tempshan.aln"
        clustaline = ClustalwCommandline("clustalw2", infile=infile, outfile=outfile)
        stdoutshan,stdershan=clustaline()
        alignTostock(outfile)
        OldShanon=CalShanon(outfile+'.stk')
        log.debug(logid+str(["OldShanon",OldShanon]))
        os.remove(outfile)
        os.remove(filesdir+filename+".dnd")
        os.remove(outfile+".stk")
        infile=""
        outfile=""
        userflanking=int(args.extension)
        indexed_genomes = index_genomes(genomes_file) # Create blastn databases

        with openfile(filen) as fl:
            for rec in SeqIO.parse(fl,'fasta'):
                pidsplit=rec.description.split()
                pid=str(pidsplit[1])
                mat2seq=str(rec.seq)
                coorflag=0
                smat=""
                emat=""
                with openfile(mapfile) as mf:
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

        log.debug(logid+str(["list 2 mat is here",list2mat]))

        if os.path.isfile(outdir+'nomat-'+filename.strip()+'.fa'):#before calling checknomat in Submit
            log.debug(logid+"The file "+filename.strip()+" already processed, will be done again")
            os.remove(outdir+'nomat-'+filename+'.fa')

        flagnomatexists,nomats,listremovedbroken,listremovedscore,listremovedN,listnomat=checknomat(filen,mapfile,matfile,outdir,filename,listremovedbroken,listremovedscore,listremovedN,nomats,listnomat)

        if flagnomatexists and nomats!=-1:
            log.debug(logid+"flagnomatexists"+str(flagnomatexists)+';'+str(nomats))
            if os.path.isfile(outdir+filename.strip()+"-new.fa"):
                os.remove(outdir+filename.strip()+"-new.fa")

            listnomatremoved=listremovedbroken+listremovedscore+listremovedN
            log.debug(logid+str(["all removed",listnomatremoved]))

            if len(listnomatremoved)>0:
                log.debug(logid+"listnomatremoved"+str(listnomatremoved))
                fl = openfile(filen)
                for record in SeqIO.parse(fl, 'fasta'):
                    tempdes=record.description
                    tempdeslst=record.description.split()
                    if tempdeslst[1].strip() not in listremovedscore and tempdeslst[1].strip() not in listremovedbroken and tempdeslst[1].strip() not in listremovedN:# and tempdeslst[1].strip() not in lst2mat:
                        newprecfile=open(outdir+filename+"-new.fa",'a')
                        newprecfile.write(">"+str(tempdes)+"\n"+str(record.seq)+"\n")
                        newprecfile.close()

            for f in [outdir+"tempfold.fa", outdir+"tempmat.fa", outdir+"temptoalign.aln",outdir+"temptoalign.dnd", outdir+"temptoalign.fa", outdir+"temptoalign.txt", outdir+"temptofold.fa", outdir+"temptopredict.aln", outdir+"temptopredict.dnd", outdir+"temptopredict.fa"]:
                if os.path.isfile(f):
                    os.remove(f)

            if len(listnomatremoved)>0:
                listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listnogenomes, listnotingenome, templong, listgoodnew=flip(filename.strip(), outdir+filename+"-new.fa", outdir, mapfile, matfile, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, list2mat, listnogenomes, listnotingenome, templong, listgoodnew, indexed_genomes, args)
            else:
                listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listnogenomes, listnotingenome, templong, listgoodnew=flip(filename.strip(), filen, outdir, mapfile, matfile, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, list2mat, listnogenomes, listnotingenome, templong, listgoodnew, indexed_genomes, args)

        elif flagnomatexists and nomats==-1:
            log.debug(logid+"flagnomatexists"+str(flagnomatexists)+';'+str(nomats))
            with open(outdir+filename.strip()+"-summ.txt","a") as summaryfile:
                summaryfile.write("no matures for the sequences and no related matures in the mapping file\n")

        elif not flagnomatexists:
            log.debug(logid+str(["this flip",flagnomatexists]))
            listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, listnogenomes, listnotingenome, templong, listgoodnew=flip(filename.strip(), filen, outdir, mapfile, matfile, listofnew, listofnewloop, listoldstatus, listofoldloop, listofold, listofboth, listofmirstar, listnomat, list2mat, listnogenomes, listnotingenome, templong, listgoodnew, indexed_genomes, args)#filename: filename/family, filen: the file itself(with the directory)

        log.debug(logid+"listofnew: "+str(listofnew))
        if os.path.isfile(outdir+filename+"-new.fa"):
            os.remove(outdir+filename+"-new.fa")

        if len(listofnew)>0:
            for i in range(0,len(listofnew)):
                if i%2==0:
                    log.debug(logid+"new"+str([">"+listofnew[i],listofnew[i+1]]))
                    with open(outdir+filename+"-res.fa","a") as familyfileres:
                        familyfileres.write(">"+str(listofnew[i])+"\n"+str(listofnew[i+1].replace('T','U'))+"\n")

        if len(listofnewloop)>0:
            for i in range(0,len(listofnewloop)):
                if i%2==0:
                    log.debug(logid+"new loop"+str([">"+listofnewloop[i],listofnewloop[i+1]]))
                    with open(outdir+filename+"-res.fa","a") as familyfileres:
                        familyfileres.write(">"+str(listofnewloop[i])+"\n"+str(listofnewloop[i+1].replace('T','U'))+"\n")

        if len(listofold)>0:
            for i in range(0,len(listofold)):
                if i%2==0:
                    log.debug(logid+"old"+str([">"+listofold[i],listofold[i+1]]))
                    with open(outdir+filename+"-res.fa","a") as familyfileres:
                        familyfileres.write(">"+str(listofold[i])+"\n"+str(listofold[i+1].replace('T','U'))+"\n")

        if len(listofoldloop)>0:
            for i in range(0,len(listofoldloop)):
                if i%2==0:
                    log.debug(logid+"old loop"+str([">"+listofoldloop[i],listofoldloop[i+1]]))
                    with open(outdir+filename+"-res.fa","a") as familyfileres:
                        familyfileres.write(">"+str(listofoldloop[i])+"\n"+str(listofoldloop[i+1].replace('T','U'))+"\n")

        listnomatbroken=listremovedbroken
        listnomatscore=listremovedscore
        listnomatN = listremovedN

        if len(list2mat)>0:
            log.debug(logid+"list2mat: "+str(list2mat))
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
                    mspos=None
                    mepos=None
                    mspos=pseq.find(fmatseq)
                    mepos=pseq.rfind(ematseq)  # rfind returns the last index of the match

                    xcutseq=len(pseq[:mspos])#used in case the seq not found in the genome
                    ycutseq=len(pseq[mepos+1:])

                    if precID==pid:
                        # On miRBase some mature were assigned for the same miRNAs.
                        if mspos == -1 or mepos == -1:
                            log.error(logid+'Referred mir or mir* from '+pid+' did not fit into the precursor, maybe some reference is repeated on mapping file?')
                            sys.exit()

                        long2matseq=""
                        long2matseq, listnogenomes, listnotingenome=getindex2mat(pseq.replace("U", "T"), specie, precID, precDes, listnogenomes, listnotingenome, args)  # from here on we have the reverse complement if seq on minus strand
                        log.debug(logid+'long2matseq: '+str(long2matseq))

                        if long2matseq!="":
                            (fspos, fepos, cutpseq) = find_positions(long2matseq, pseq.replace("U", "T"),
                                                                             fmatseq.replace("U", "T"), ematseq.replace("U", "T"),
                                                                             userflanking)

                            with open(outdir+filename+"-res.fa","a") as familyfileres:
                                familyfileres.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")
                            with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                                familyfileresfinal.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")

                        elif long2matseq=="":
                            if xcutseq>=xcut and xcut>=0 and xcut<=50:
                                fspos=mspos-xcut
                            elif xcutseq<xcut and xcut>=0 and xcut<=50:#if requested flankings more than available for this seq, because not found in genome,then take all available
                                fspos=0
                            elif xcutseq>10 and (xcut>50 or xcut<0):
                                fspos=mspos-10
                            elif xcutseq<=10 and (xcut>50 or xcut<0):
                                fspos=0

                            if ycutseq>=ycut and ycut>=0 and ycut<=50:
                                fepos=mepos+ycut+1
                            elif ycutseq<ycut and ycut>=0 and ycut<=50:
                                fepos=len(pseq)
                            elif ycutseq>10 and (ycut>50 or ycut<0):
                                fepos=mepos+11
                            elif ycutseq<=10 and (ycut>50 or ycut<0):
                                fepos=len(pseq)

                            cutpseq=pseq[fspos:fepos+len(ematseq)]
                            with open(outdir+filename+"-res.fa","a") as familyfileres:
                                familyfileres.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")
                            with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                                familyfileresfinal.write(">"+str(record.description)+"\n"+str(cutpseq).replace('T','U')+"\n")


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
            log.debug(["res",resprecid,mat2seq,mat1seq])

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
                    curmatseq = None
                    curmatstar = None
                    # Here use find because both matures were reported at the
                    # Beginning and we are iterating over the complete mature file
                    for record in SeqIO.parse(mtf, 'fasta'):
                        sequence=str(record.seq)
                        if firstmat.strip() == record.description.split(" ")[1]: # specific label
                            first = 'Found'
                            curmatseq = sequence

                        if lastmat.strip() == record.description.split(" ")[1]: # specific label
                            second = 'Found'
                            curmatstar = sequence
                    (startmat, startmatstar, finalseq) = find_positions(mat1seq, mat1seq,curmatseq, curmatstar,userflanking)
                    log.debug(logid+str([curmatseq,startmat, startmatstar, finalseq]))
                    if startmat == None or startmatstar == None:
                        log.debug(logid+'Not possible to locate the mapping referred mir or mir* on the mature file for '+resprecid+' with '+mat2seq)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        finalseq=curmatseq
                        list2matcoor.append("NULL")#matstarID/here no mir* found
                        list2matcoor.append("NULL")#matstarID/here no mir* found
                        list2matcoor.append(startmat)
                        list2matcoor.append(endmat)
                        list2matcoor.append(startmatstar)
                        list2matcoor.append(endmatstar)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        listbadmature.append(resprecid)
                        break

                    endmat = startmat+len(curmatseq)-1
                    endmatstar = (startmatstar+len(curmatstar))-1
                    log.debug(logid+str(["heres new",mat2seq,finalseq,startmat,endmat,startmatstar,endmatstar,curmatseq,curmatstar]))
                    if first and second and startmat != None and startmatstar != None:
                        coorflag=1

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
                    startmat=0
                    endmat=0
                    startmatstar=0
                    endmatstar=0
                    list1matcoor.append(resprecid)
                    finalseq="A"
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

                    # Here genome did not exists, then it is reeplaced by the same precursor seq to obtain
                    # the mature coordinates
                    (coortemp1, coortemp2, finalseq) = find_positions(mat1seq, mat1seq,
                                                 curmatseq, curmatstar,
                                                 userflanking)

                    if coortemp2<coortemp1:
                        tempseqex=curmatseq[:]
                        curmatseq=curmatstar[:]
                        curmatstar=tempseqex[:]

                    log.debug(logid+str(["no long",coortemp1,coortemp2,curmatseq,curmatstar]))
                    startmat=coortemp1
                    endmat=(startmat+len(curmatseq))-1

                    if star:
                        startmatstar=coortemp2
                        endmatstar=startmatstar+int(len(curmatstar)-1)
                        log.debug(logid+str(['coor star not',startmatstar,endmatstar,finalseq]))#,startmatstarlong,endmatstarlong)

                    if star and nstar and startmatstar!=-1 and endmatstar!=-1 and startmatstar>endmat:
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                            familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break
                    elif star and nstar and startmatstar!=-1 and endmatstar!=-1 and startmatstar<endmat:#to put them in order
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
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
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
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
                    curmatseq=None
                    curmatstar=None
                    startmat=0
                    endmat=0
                    startmatstar=0
                    endmatstar=0
                    longseq=""
                    startmat=0
                    finalseq="A"
                    startmatstar=0
                    endmatstar=0
                    list1matcoor.append(resprecid)
                    indexlongmat=int(templong.index(resprecid.strip()))
                    longseq=str(templong[indexlongmat+1]).replace('T','U')

                    with openfile(matfile) as mtf:
                        for reco in SeqIO.parse(mtf, 'fasta'):
                            splitreco=reco.description.split()
                            if curmatID == splitreco[1]:
                                curmatseq=str(reco.seq)
                                nstar=True
                                break

                    log.debug(logid+'Curmatseq: '+str(curmatseq))

                    with openfile(outdir+filename.strip()+"-mirstar.fa") as msfa:
                        for starrec in SeqIO.parse(msfa,'fasta'):
                            curmatsplit=(starrec.description).split()
                            curmatsplit1=(curmatsplit[1]).split('-')
                            starrecID=curmatsplit1[0]
                            if (curmatID).strip()==(starrecID).strip():
                                curmatstar=str(starrec.seq)
                                star=True
                                break

                    if not curmatstar:
                        log.error(logid+'Not possible to define curmatstar for '+ resprecid +' in '+str(matfile)+' and '+str(outdir+filename.strip()+"-mirstar.fa"))
                        sys.exit()

                    log.debug(["coor1temp",longseq,curmatseq,curmatstar,resprecid])
                    (coortemp1, coortemp2, finalseq) = find_positions(longseq, mat2seq,
                                                 curmatseq, curmatstar,
                                                 userflanking)

                    if coortemp1 == None or coortemp2 == None:
                        log.debug(logid+'Not possible to locate miR or miR* in ' + resprecid + " " + curmatseq)
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        finalseq=curmatseq
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
                        listbadmature.append(resprecid)
                        break

                    if coortemp2<coortemp1:
                        tempseqex=curmatseq[:]
                        curmatseq=curmatstar[:]
                        curmatstar=tempseqex[:]

                    #Look for mir in final seq:
                    startmirfinal = int(finalseq.find(curmatseq))
                    endmirfinal = int(startmirfinal + len(curmatseq) - 1)
                    subsetnomir = finalseq[endmirfinal:]
                    indexupdate = finalseq.find(subsetnomir)
                    #Look mirstar in final seq
                    startmirstarfinal = subsetnomir.find(curmatstar)
                    endmirstarfinal = int(startmirstarfinal + len(curmatstar) - 1)
                    #update coordinates
                    startmat=startmirfinal
                    endmat=endmirfinal
                    startmatstar = startmirstarfinal + indexupdate
                    endmatstar = endmirstarfinal + indexupdate
                    log.debug(['coor',startmat, endmat, startmatstar,endmatstar])
                    # Previous evaluations:
                    (structure, mfe, loopsize, length_precursor) = evaluate_final_hairpin(finalseq, startmat, endmat, startmatstar, endmatstar, resprecid)

                    if mfe > -10:
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
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
                        listhighmfe.append(resprecid)
                        break
                    if length_precursor > 200:
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
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
                        listlonghairpin.append(resprecid)
                        break
                    if loopsize <= 1:
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
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
                        listremovenoloop.append(resprecid)
                        break

                    if star and nstar and startmatstar!=-1 and endmatstar!=-1 and startmatstar>endmat and loopsize > 1:
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                            familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    elif star and nstar and startmatstar!=-1 and endmatstar!=-1 and startmatstar<endmat and loopsize > 1:#to put them in order
                        list1matcoor.append(curmatsplit[1].strip())#matstarID
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
                            familyfileresfinal.write(">"+resprecdesc+"\n"+finalseq+"\n")
                        startmat=0
                        endmat=0
                        startmatstar=0
                        endmatstar=0
                        break

                    # Here is not possible to calculate the loop
                    elif (not star and nstar) or startmatstar==-1 or endmatstar==-1:
                        startmat=int(mat1seq.find(curmatseq))
                        endmat=(startmat+len(curmatseq))-1
                        finalseq=mat1seq
                        startmstar=0
                        endmatstar=len(mat1seq)-1
                        list1matcoor.append(curmatID.strip())#mat original
                        list1matcoor.append("NULL")#matstarID/here no mir* found
                        list1matcoor.append(startmat)
                        list1matcoor.append(endmat)
                        list1matcoor.append(startmatstar)
                        list1matcoor.append(endmatstar)
                        with open(outdir+filename.strip()+"-Final.fasta","a") as familyfileresfinal:
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

        log.debug([len(list2matcoor), len(list1matcoor), list2matcoor, list1matcoor])
        listmatcoor=list2matcoor+list1matcoor

        mi=0
        mk=0
        log.debug(len(listmatcoor))
        log.debug(len(listmatcoor)/7)

        while mi<=int(len(listmatcoor))-7:
            mk=mk+1
            r=0

            for n in range(mk+1,int(len(listmatcoor)/7)+1):
                with open(outdir+filename.strip()+"-Final.anc","a") as anchorcoorfile:
                    if listmatcoor[mi+1] == "NULL" or listmatcoor[mi+3] == -1 or listmatcoor[mi+10+r] == -1 or listmatcoor[mi+3] == 0 or listmatcoor[mi+10+r] == 0:
                        continue
                    else:
                        anchorcoorfile.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+3]+1)+" "+str(listmatcoor[mi+10+r]+1)+" "+str(22)+" "+str(1)+"\n")
                    if listmatcoor[mi+1] == "NULL" or listmatcoor[mi+5] == -1 or listmatcoor[mi+12+r] == -1 or listmatcoor[mi+5] == 0 or listmatcoor[mi+12+r] == 0:
                        continue
                    else:
                        anchorcoorfile.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+5]+1)+" "+str(listmatcoor[mi+12+r]+1)+" "+str(22)+" "+str(1)+"\n")
                r=r+7

            mi=mi+7

        log.debug(["here listmatcoor 1 ",listmatcoor])
        maxidesc=0

        finalstk=open(outdir+filename.strip()+'.stk','a')
        finalstk.write('# STOCKHOLM 1.0\n')
        if matrdir:
            fs=os.environ["DIALIGN2_DIR"]=matrdir
            log.debug(matrdir)
        f1=os.popen("dialign2-2 -n -anc -fa "+outdir+filename.strip()+'-Final.fasta')
        log.debug(f1)
        f1.close()

        if os.path.isfile(outdir+filename.strip()+'-Final.fa'):
            doalifold(outdir+filename.strip()+"-Final.fa",outdir)
            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-Final.fa'),'fasta'):
                if len(rec.description.strip())>maxidesc:
                    maxidesc=len(rec.description.strip())

                if maxidesc<len('#=GC SS_cons'):
                    maxidesc=len('#=GC SS_cons')

            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-Final.fa'),'fasta'):
                finalstk.write(str(rec.description.strip())+" "*(maxidesc-len(rec.description.strip())+2)+str(rec.seq)+"\n")

            struct=getstructure(outdir+'alifoldtemp.txt')
            os.remove(outdir+'alifoldtemp.txt')
            os.remove(outdir+filename.strip()+"-Final.ali")
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
            log.debug([sthalfsum,ndhalfsum,totalstnucnum,totalndnucnum])
            # Here calculated the average of nt on left side of align (stnucavg) and on the right side (ndnucavg) on the alignment
            stnucavg=totalstnucnum/numofseqs
            ndnucavg=totalndnucnum/numofseqs
            alignment=AlignIO.read(outdir+filename.strip()+'.stk',"stockholm")
            # Here iterate again but with calculated average scores
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
                    log.debug([record.id,'st'])
                    log.debug([stnucnum,ndnucnum,nucnum])
                    log.debug(["list new good",corid,listgoodnew])
                    countcorrected,countcorrectedTonew,listmisalignedcorr,listcorrected,listcorrectedori=correct(corid.strip(),userflanking,countcorrected,countcorrectedTonew,listofnew,listofnewloop,listoldstatus,templong,listmisalignedcorr,listcorrected,listcorrectedori,listgoodnew)

                if (ndnucnum>(ndnucavg*1.5) and (stnucnum<(stnucavg/2) or stnucnum==0)) or (ndnucnum>0.7*nucnum and (stnucnum<(stnucavg/2) or stnucnum==0)):#-(ndnucavg*0.1)):
                    log.debug([record.id,'nd'])
                    log.debug([stnucnum,ndnucnum,nucnum])
                    log.debug(["list new good",corid,listgoodnew])
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
            log.debug(["here listmatcoor 2 ",listmatcoor])

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

            log.debug(["here listmatcoor 3 ",listmatcoor])
            mi=0
            mk=0
            while mi<=int(len(listmatcoor))-7:
                mk=mk+1
                r=0
                for n in range(mk+1,int(len(listmatcoor)/7)+1):
                    with open(outdir+filename.strip()+"-corrected.anc","a") as anchorcoorfilecorrected:
                        anchorcoorfilecorrected.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+3])+" "+str(listmatcoor[mi+10+r])+" "+str(22)+" "+str(1)+"\n")
                        anchorcoorfilecorrected.write(str(mk)+" "+str(n)+" "+str(listmatcoor[mi+5])+" "+str(listmatcoor[mi+12+r])+" "+str(22)+" "+str(1)+"\n")
                    r=r+7

                mi=mi+7

            maxidesc=0
            finalstkcorrected=open(outdir+filename.strip()+'corrected.stk','a')
            finalstkcorrected.write('# STOCKHOLM 1.0\n')
            if matrdir:
                fe=os.environ["DIALIGN2_DIR"]=matrdir
            f11=os.popen("dialign2-2 -n -anc -fa  "+outdir+filename.strip()+'-corrected.fasta')
            f11.close()

            doalifold(outdir+filename.strip()+"-corrected.fa",outdir)

            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-corrected.fa'),'fasta'):
                if len(rec.description.strip())>maxidesc:
                    maxidesc=len(rec.description.strip())

                if maxidesc<len('#=GC SS_cons'):
                    maxidesc=len('#=GC SS_cons')

            for rec in SeqIO.parse(openfile(outdir+filename.strip()+'-corrected.fa'),'fasta'):
                finalstkcorrected.write(str(rec.description.strip())+" "*(maxidesc-len(rec.description.strip())+2)+str(rec.seq)+"\n")

            struct=getstructure(outdir+'alifoldtemp.txt')
            os.remove(outdir+'alifoldtemp.txt')
            os.remove(outdir+filename.strip()+"-corrected.fa")
            os.remove(outdir+filename.strip()+"-corrected.ali")
            finalstkcorrected.write('#=GC SS_cons'+" "*(maxidesc-len('#=GC SS_cons')+2)+str(struct))
            finalstkcorrected.close()
            NewShanon=CalShanon(outdir+filename.strip()+'corrected.stk')

        else:
            if os.path.isfile(outdir+filename.strip()+'-Final.fa'):
                os.remove(outdir+filename.strip()+"-Final.fa")
                NewShanon=CalShanon(outdir+filename.strip()+'.stk')
                log.debug(["stk file studied is: "+outdir+filename.strip()+'.stk'])
                log.debug(["final.fa exists",NewShanon])

            else:
                NewShanon=OldShanon
                log.debug(["final.fa DO NOT exists",NewShanon])

        log.debug(["new shanon",NewShanon])
        finalcoor=open(outdir+filename+"-FinalCoor.txt","a")
        log.debug(["here listmatcoor 4 ",len(listmatcoor),listmatcoor])

        for l in range(0,len(listmatcoor),7):
            finalcoor.write(str(listmatcoor[l]).strip()+" "+str(listmatcoor[l+1]).strip()+" "+str(listmatcoor[l+2]).strip()+" "+str(listmatcoor[l+3]).strip()+" "+str(listmatcoor[l+4]).strip()+" "+str(listmatcoor[l+5]).strip()+" "+str(listmatcoor[l+6]).strip()+"\n")

        finalcoor.close()

        with open(outdir+filename.strip()+"-summ.txt","a") as summaryfile:
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
                    if listnomat[k] not in listnomatbroken and listnomat[k] not in listnomatscore and listnomat[k] not in listnomatN:
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
                summaryfile.write("--->All precursors have annotated matures\n")

            summaryfile.write("\n")
            summaryfile.write("---------------------------Original precursors totally removed ----------------------------\n")

            if len(listofboth)>0:
                for i in range(0,len(listofboth)):
                    if i%2==0:
                        tempsplit=listofboth[i].split()
                        summaryfile.write(str(tempsplit[1].strip())+"\n")
            elif len(listnomatN)>0:
                for j in range(0, len(listnomatN)):
                    summaryfile.write(str(listnomatN[j].strip()) + "\n")
            elif len(listremovenoloop)>0:
                for k in range(0, len(listremovenoloop)):
                    summaryfile.write(str(listremovenoloop[k].strip()) + "\n")
            elif len(listhighmfe)>0:
                for k in range(0, len(listhighmfe)):
                    summaryfile.write(str(listhighmfe[k].strip()) + "\n")
            elif len(listbadmature)>0:
                for k in range(0, len(listbadmature)):
                    summaryfile.write(str(listbadmature[k].strip()) + "\n")
            elif len(listlonghairpin)>0:
                for k in range(0, len(listlonghairpin)):
                    summaryfile.write(str(listlonghairpin[k].strip()) + "\n")
            else:
                summaryfile.write("---> NO Original precursors totally removed ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("#Precursors without loop regions (<= 1 nt)---------------------------\n")
            if len(listremovenoloop)>0:
                for j in range(0, len(listremovenoloop)):
                    summaryfile.write(str(listremovenoloop[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported a valid loop region ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("# Precursors with high N content---------------------------\n")
            if len(listnomatN)>0:
                for j in range(0, len(listnomatN)):
                    summaryfile.write(str(listnomatN[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported valid nucleotide composition with N content > 60 % ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("# Precursors with high MFE---------------------------\n")
            if len(listhighmfe)>0:
                for j in range(0, len(listhighmfe)):
                    summaryfile.write(str(listhighmfe[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported a valid MFE range ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("# Precursors with unsuccessful mature prediction---------------------------\n")
            if len(listbadmature)>0:
                for j in range(0, len(listbadmature)):
                    summaryfile.write(str(listbadmature[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported canonical mature positions ----------------------------\n")

            summaryfile.write("\n")
            summaryfile.write("# Precursors with reported precursor size > 200 nt---------------------------\n")
            if len(listlonghairpin)>0:
                for j in range(0, len(listlonghairpin)):
                    summaryfile.write(str(listlonghairpin[j].strip()) + "\n")
            else:
                summaryfile.write("---> All precursors reported a valid size ----------------------------\n")

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
        Removed=(len(listofboth)/2 + len(listnomatN) + len(listremovenoloop) + len(listhighmfe) + len(listbadmature) + len(listlonghairpin))
        TotalnumberofSequences=Processed+Removed
        #no mats/predicted-no prediction
        predicted=tempcountsucnomat
        noprediction=len(listnomat)-tempcountsucnomat-len(listnomatN)-len(listremovenoloop)-len(listhighmfe)-len(listbadmature)-len(listlonghairpin)

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

        with open(outdir+filename.strip()+"-summ.txt","a") as summaryfile:
            summaryfile.write("---------------------------Results In Numbers----------------------------\n")
            summaryfile.write("*Number of remained precursors= "+str(int((len(listofold)/2)+(len(list2mat)/3)))+"\n")
            summaryfile.write("*Number of remained precursors with bad positioned matures= "+str(int((len(listofoldloop)/2)))+"\n")
            summaryfile.write("*Number of flipped(changed) precursors= "+str(int((len(listofnew)/2)))+"\n")
            summaryfile.write("*Number of flipped(changed) precursors with bad positioned matures= "+str(int((len(listofnewloop)/2)))+"\n")
            summaryfile.write("*Number of removed precursors= "+str(int(Removed))+"\n")
            summaryfile.write("*Number of precursors without a given matures= "+str(int(len(listnomat)))+"\n")
            summaryfile.write("*Number of precursors with successfully predicted matures= "+str(int(len(listnomat)-len(listremovedbroken)-len(listremovedscore)-len(listremovedN)-len(listremovenoloop)-len(listhighmfe)-len(listbadmature)-len(listlonghairpin)))+"\n")
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
            if listnomat[k] not in listnomatbroken and listnomat[k] not in listnomatscore and listnomat[k] not in listnomatN and listnomat[k] not in listremovenoloop and listnomat[k] not in listhighmfe and listnomat[k] not in listbadmature and listnomat[k] not in listlonghairpin:
                listsuccpred.append(str(listnomat[k].strip()))

        listremovedjson=[]
        for i in range(0,len(listofboth)):
            if i%2==0:
                tempsplit=listofboth[i].split()
                listremovedjson.append(str(tempsplit[1].strip()))
        for i in range(0,len(listnomatN)):
            listremovedjson.append(listnomatN[i].strip())
        for j in range(0, len(listremovenoloop)):
            listremovedjson.append(listremovenoloop[j].strip())
        for k in range(0, len(listhighmfe)):
            listremovedjson.append(listhighmfe[k].strip())
        for k in range(0, len(listbadmature)):
            listremovedjson.append(listbadmature[k].strip())
        for k in range(0, len(listlonghairpin)):
            listremovedjson.append(listlonghairpin[k].strip())
        data = {
            "Processed":int(Processed),
            "Precursors totally removed":{
                "Number":int(len(listremovedjson)),
                "IDs":listremovedjson,
            },
            "Remained precursors":int((len(listofold)/2)+(len(list2mat)/3)-len(listnomatN)-len(listremovenoloop)-len(listhighmfe)-len(listbadmature)-len(listlonghairpin)),
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
            "Precursors with high N content":{
                "Number": int(len(listnomatN)),
                "IDs": listnomatN,
            },
            "Precursors without loop":{
                "Number": int(len(listremovenoloop)),
                "IDs": listremovenoloop,
            },
            "Precursors with high MFE":{
                "Number": int(len(listhighmfe)),
                "IDs": listhighmfe,
            },
            "Precursors with unsuccessful mature prediction":{
                "Number": int(len(listbadmature)),
                "IDs": listbadmature,
            },
            "Precursors with long precursors":{
                "Number": int(len(listlonghairpin)),
                "IDs": listlonghairpin,
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
        del listremovedN[:]
        del listnomat[:]
        del list2mat[:]
        del listmatcoor[:]
        del listnogenomes[:]#list of IDs that no existing genome to search, for their species
        del listnotingenome[:]

        familyfileres.close()
        log.debug("done")
        log.debug([listofold,listofnew,listofnewloop,listofoldloop,listremovedbroken,listremovedscore,listremovedN,listofmirstar])
        log.debug(list2mat)
        if os.path.isfile(outdir+filename.strip()+"-res.fa"):
            os.remove(outdir+filename.strip()+"-res.fa")

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def parseargs():
    parser = argparse.ArgumentParser(description='MIRfix automatically curates miRNA datasets by improving alignments of their precursors, the consistency of the annotation of mature miR and miR* sequence, and the phylogenetic coverage. MIRfix produces alignments that are comparable across families and sets the stage for improved homology search as well as quantitative analyses.')
    parser.add_argument("-j", "--cores", type=int, default=1, help='Number of parallel processes to run')
    parser.add_argument("-f", "--families", type=str, required=True, help='Path to list of families to work on')
    parser.add_argument("-i", "--famdir", type=str, required=True, help='Directory where family files are located')
    parser.add_argument("-g", "--genomes", type=str, required=True, help='Genome FASTA files to parse')
    parser.add_argument("-m", "--mapping", type=str, required=True, help='Mapping between precursor and families')
    parser.add_argument("-a", "--mature", type=str, required=True, help='FASTA files containing mature sequences')
    parser.add_argument("-d", "--maturedir", type=str, default='', help='Directory of matures')
    parser.add_argument("-o", "--outdir", type=str, default='', help='Directory for output')
    parser.add_argument("--force", action='store_true', help='Force MIRfix to overwrite existing output directories')
    parser.add_argument("-e", "--extension", type=int, default=10, help='Extension of nucleotides for precursor cutting')
    parser.add_argument("-l", "--logdir", type=str, default='LOGS', help='Directory to write logfiles to')
    parser.add_argument("--loglevel", type=str, default='WARNING', choices=['WARNING','ERROR','INFO','DEBUG','CRITICAL'], help="Set log level")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def main(args):
    try:
        #  Logging configuration
        logdir = args.logdir
        logfile = str.join(os.sep,[os.path.abspath(logdir),scriptname+'.log'])
        makelogdir(logdir)
        makelogfile(logfile)
        #  Multiprocessing with spawn
        nthreads = args.cores or 2
        multiprocessing.set_start_method('spawn')  # set multiprocessing start method to safe spawn
        pool = multiprocessing.Pool(processes=nthreads, maxtasksperchild=1)
        queue = multiprocessing.Manager().Queue(-1)
        listener = multiprocessing.Process(target=listener_process, args=(queue, listener_configurer, logfile, args.loglevel))
        listener.start()
        worker_configurer(queue, args.loglevel)
        log.info(logid+'Running '+scriptname+' on '+str(args.cores)+' cores.')
        log.info(logid+'CLI: '+sys.argv[0]+' '+'{}'.format(' '.join( [shlex.quote(s) for s in sys.argv[1:]] )))
        lfams = []
        with openfile(args.families) as filelist:
            for line in filelist:
                line = line.strip()
                if line.endswith( ".fa"):
                    line = line[:-3]
                ## check if there are already output directories
                ## to prevent errors
                if os.path.exists( args.outdir + line + ".out/") and not args.force:
                    sys.exit( '\nError: At least one output directory already exists! Won\'t override!\nTo override the original output folder, please specify the \'--force\' option.\n')
                lfams.append( line)

        for fam in lfams:
            pool.apply_async(sublist, args=(queue, worker_configurer, args.loglevel, fam, args))
        pool.close()
        pool.join()
        queue.put_nowait(None)
        listener.join()
        sys.exit()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

##############################MAIN##############################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args = parseargs()

        # find_executable('clustalw2') or sys.exit('Please install clustalw2 to run this')
        find_executable('dialign2-2') or sys.exit('Please install dialign2-2 to run this')

        main(args)
        sys.exit()

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))
