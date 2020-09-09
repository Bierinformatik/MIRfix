import os, sys, inspect
from lib.logger import *
##other modules
import numpy as np
import heapq
from operator import itemgetter
from natsort import natsorted, ns
import traceback as tb
import sys
import re
import pprint
from io import StringIO
import os
import gzip
import math
from collections import defaultdict
##import Bio modules
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

scriptname = os.path.basename(inspect.stack()[-1].filename).replace('.py','')

# Code:All subs from here on
#Files
def openfile(f):
    logid = scriptname+'.openfile: '
    try:
        return open(f,'r') if not '.gz' in f[-4:] else gzip.open(f,'rt')
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def getfilename(dirfile):#to get the name of the file without directory or extension, in principal was used to get famname
    logid = scriptname+'.getfilename: '
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
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

#foldings
def dofold(listnewold,oldid,precseq,newid,newseq):
    logid = scriptname+'.dofold: '
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
        log.debug(["folding old is here:",oldstruct,oldscore])
        log.debug(["folding new is here:",newstruct,newscore])
        listnewold.extend((oldid,str(precseq),oldstruct,oldscore,newid,str(newseq),newstruct,newscore))
        return listnewold

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def foldnomat(inputfasta,outputfasta):#fold the new and the old sequences, using temp file every time I get the new sequence from the original
    logid = scriptname+'.foldnomat: '
    try:
        f=os.popen("RNAfold -d3 --noPS --noLP <"+inputfasta)
        fi=f.read()
        wr=open(outputfasta,"w")
        f.close()
        wr.write(fi)
        wr.close
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def doalifold(alnfile,outdir):
    logid = scriptname+'.doalifold: '
    try:
        f=os.popen("RNAalifold --noPS -r "+alnfile) # -r = RIBOSUM scoring
        alifoldtemp=outdir+'alifoldtemp.txt'
        foldrestemp=open(alifoldtemp,'w')
        fi=f.read()
        foldrestemp.write(fi)
        foldrestemp.close()
        f.close()
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def getstructure(alifoldtemp):
    logid = scriptname+'.getstructure: '
    try:
        for line in (openfile(alifoldtemp)):
            item=line.split()
            hairpin = ''
            if len(item)>0 and 'A' not in item[0] and 'C' not in item[0] and 'G' not in item[0] and 'T' not in item[0] and 'a' not in item[0] and 'c' not in item[0] and 'g' not in item[0] and 't' not in item[0] and 'N' not in item[0] and 'n' not in item[0] and ('.' in item[0] or '(' in item[0] or ')' in item[0]):
                hairpin=item[0]
                return str(hairpin)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def cal_ent(s):
    logid = scriptname+'.cal_ent: '
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
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def CalShanon(stkfile):
    logid = scriptname+'.CalShanon: '
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
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

#comparison
def comp_5p(ll, lr, sm, ss, precursor, run):
    logid = scriptname+'.comp_5p: '
    try:
        log.debug(['comp5p', ll, lr, sm, ss, precursor, run])
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
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


#format
def alignTostock(align):
    logid = scriptname+'.alignTostock: '
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
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


### CAVH functions

def count_repetitions(longseq,pattern):
    results=0
    len_pattern = len(pattern)
    for i in range(len(longseq)):
        if longseq[i:i+len_pattern] == pattern:
            results += 1
            print(longseq[i:i+len_pattern],i,i+len_pattern)
    return results

def search_pattern(longseq, pattern):
    start = longseq.find(pattern)
    end = int((start + len(pattern)) -1)
    return (start, end)


def generate_array(longseq, pattern, number, mode):
    pattern_len = len(pattern)
    array = []
    short_sequence = None
    update_index = None
    for j in range(0, number):  # (0, n-1)
        if j == 0:  # Search on the complete sequence
            start = longseq.find(pattern)
            end = int((start + pattern_len) - 1)
            assembly_data = [start, end, mode]
            array.append(assembly_data)
            short_sequence = longseq[end + 1:]
            update_index = longseq.find(short_sequence)
        else:
            start_temp = short_sequence.find(pattern)
            if start_temp != -1:
                end_temp = int((start_temp + pattern_len) - 1)
                start = int(start_temp + update_index)
                end = int(end_temp + update_index)
                assembly_data = [start, end, mode]
                array.append(assembly_data)
                short_sequence = short_sequence[end_temp + 1:]
                update_index = longseq.find(short_sequence)
    return array


def detect_overlapp_in_array(complete_array):
    columns = len(complete_array)
    for i in range(0, columns -1):
        if complete_array[i+1][0] in range(complete_array[i][0],
           complete_array[i][1] + 1):
            return 1
    return 0


def detect_pairs(complete_array, longseq, mat1seq, curmatseq, curmatstar):
    columns = len(complete_array)  # Complete number columns
    # Detect overlapping
    overlapping = detect_overlapp_in_array(complete_array)
    pair = defaultdict(dict)
    number_pair = 0
    if overlapping == 0:
        for i in range(0, columns-1):  # In array =columns -1,range not reach end#
            if complete_array[i][-1] != complete_array[i+1][-1]:
                pair[number_pair] = [complete_array[i], complete_array[i+1]]
                number_pair = number_pair + 1
    else:  # At least one pair overlapped with another one
        substring = None
        update_index = None
        for j in range(0, columns -1):  # (0, n-1)
            if j == 0:  # Search on the complete sequence
                identity_first = complete_array[j][-1]  # M or S
                if identity_first == "M":
                    start_m = longseq.find(curmatseq)  # Start M
                    end_m = (start_m + len(curmatseq)) - 1  # End M
                    substring = longseq[end_m + 1:]  # Restrict seq no M
                    start_temp_s = substring.find(curmatstar)  # Detect start S
                    update_index = longseq.find(substring)  # Get index
                    start_s = start_temp_s + update_index  # Update index S
                    end_s = (start_s + len(curmatstar)) -1
                    element1 = [start_m, end_m, "M"]
                    element2 = [start_s, end_s, "S"]
                    pair[number_pair] = [element1, element2]
                    number_pair = number_pair + 1
                else:  # When S is first, search M
                    start_s = longseq.find(curmatstar)  # Start S
                    end_s = (start_s + len(curmatstar))-1  # End S
                    substring = longseq[end_s + 1:]  # Restrict seq no S
                    update_index = longseq.find(substring)
                    start_temp_m = substring.find(curmatseq)  # Search M
                    start_m = start_temp_m + update_index  # Start M
                    end_m = (start_m + len(curmatseq)) - 1  # End M
                    element1 = [start_s, end_s, "S"]
                    element2 = [start_m, end_m, "M"]
                    pair[number_pair] = [element1, element2]
                    number_pair = number_pair + 1
            else:  # This is not the first element
                if complete_array[j][0] in range(complete_array[j-1][0],
                                                 complete_array[j-1][1]):
                    continue
                else:
                    identity_first = complete_array[j][-1]  # M or S
                    if identity_first == "M":
                        start_temp_m = substring.find(curmatseq)  # Start M in sub
                        end_temp_m = (start_temp_m + len(curmatseq)) - 1
                        start_m = start_temp_m + update_index  # Updated M Start
                        end_m = (start_m + len(curmatseq)) - 1  # Updated M end
                        substring = substring[end_temp_m + 1:]   # String no start M
                        update_index = longseq.find(substring)
                        start_temp_s = substring.find(curmatstar)
                        if (start_temp_s == -1):  # Not possible to map S
                            continue
                        start_s = start_temp_s + update_index  # Updated start S
                        end_s = (start_s + len(curmatstar)) -1  # Updated end S
                        element1 = [start_m, end_m, "M"]
                        element2 = [start_s, end_s, "S"]
                        pair[number_pair] = [element1, element2]
                        number_pair = number_pair + 1
                    else:  # When S is first, search M
                        start_temp_s = substring.find(curmatstar)  # Start M in sub
                        end_temp_s = (start_temp_s + len(curmatstar)) - 1
                        start_s = start_temp_s + update_index  # Updated S Start
                        end_s = (start_s + len(curmatstar)) - 1  # Updated S end
                        substring = substring[end_temp_s + 1:]   # String no start S
                        update_index = longseq.find(substring)
                        start_temp_m = substring.find(curmatseq)
                        start_m = start_temp_m + update_index  # Update start M
                        if (start_temp_m == -1):  # Not possible to map M
                            continue
                        end_m = (start_s + len(curmatseq)) -1  # Updated end M
                        element1 = [start_s, end_s, "S"]
                        element2 = [start_m, end_m, "M"]
                        pair[number_pair] = [element1, element2]
                        number_pair = number_pair + 1
    return (pair, number_pair)


def define_best_pair(pairs, distances_loop, startprecursor,
                     endprecursor, precursorlen, flanking, longseq):
    all_distances_sorted = sorted(distances_loop, key=itemgetter(-1))
    for i in range(len(all_distances_sorted)):
        if (all_distances_sorted[i][1][0][0] in range(startprecursor, endprecursor + 1)
                and all_distances_sorted[i][1][0][1] in range(startprecursor, endprecursor + 1)
                and all_distances_sorted[i][1][1][0] in range(startprecursor, endprecursor + 1)
                and all_distances_sorted[i][1][1][1] in range(startprecursor, endprecursor + 1)):
            coord1 = all_distances_sorted[i][1][0][0]
            coord2 = all_distances_sorted[i][1][1][0]
            finalcoord2 = all_distances_sorted[i][1][1][1]
            startfinalseq = coord1-flanking
            endfinalseq = finalcoord2+flanking
            correctfinalseq = longseq[startfinalseq:endfinalseq+1]
            break
    return (coord1, coord2, correctfinalseq)


def find_positions(longseq, mat1seq, curmatseq, curmatstar, flanking):
    # Features sequences
    precursorlen = len(mat1seq)
    # Locate
    startprecursor = longseq.find(mat1seq)
    endprecursor = int((startprecursor + precursorlen)-1)
    # Search the number of substrings along the long sequence
    mir = count_repetitions(longseq,curmatseq)
    mirstar = count_repetitions(longseq,curmatstar)
    mir_array = generate_array(longseq, curmatseq, mir, "M")
    mirstar_array = generate_array(longseq, curmatstar, mirstar, "S")
    # All array
    all_array = mir_array + mirstar_array
    all_array_sorted = sorted(all_array, key=itemgetter(0))
    (pairs, number_pairs) = detect_pairs(all_array_sorted, longseq,
                                         mat1seq, curmatseq, curmatstar)
    distances_loop = []
    for pair in pairs:
        loopsize = abs(pairs[pair][0][1] - pairs[pair][1][0]) - 1
        distances_loop.append([pair, pairs[pair], loopsize])
    (coord1, coord2, correctfinalseq) = define_best_pair(pairs, distances_loop,
                                                         startprecursor,
                                                         endprecursor,
                                                         precursorlen,
                                                         flanking,
                                                         longseq)
    return (coord1, coord2, correctfinalseq)

#Utils
def removekey(d, key):
    logid = scriptname+'.removekey: '
    try:
        r = dict(d)
        del r[key]
        return r
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def getlowest_list(a, n):
    logid = scriptname+'.getlowest_list: '
    try:
        if n > len(a) - 1:
            b = len(a) - 1
        else:
            b = n
        if len(a) > 0 and n > 0:
            return list(np.partition(a, b)[:n])
        else:
            return list(None for i in range(n))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def gethighest_list(a, n):
    logid = scriptname+'.gethighest_list: '
    try:
        if len(a)-n < 0:
            b = len(a)-1
        else:
            b = len(a)-n
        if len(a) > 0 and n > 0:
            return list(np.partition(a, b)[-n:])
        else:
            return list(None for i in range(n))
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def getlowest_dict(a, n):
    logid = scriptname+'.getlowest_dict: '
    try:
        if n > len(a):
            b = len(a)
        else:
            b = n
        if len(a) > 0:
            return dict(heapq.nsmallest(b,a.items(), key=itemgetter(1)))
        else:
            return dict({i:None for i in range(n)})
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def gethighest_dict(a, n):
    logid = scriptname+'.gethighest_dict: '
    try:
        if n > len(a):
            b = len(a)
        else:
            b = n
        if len(a) > 0:
            return dict(heapq.nlargest(b,a.items(), key=itemgetter(1)))
        else:
            return dict({i:None for i in range(n)})
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def convertcol(entry):
    logid = scriptname+'.convertcol: '
    try:
        if isinvalid(entry):
#       if entry is None or entry == 'NA' or entry == 'nan' or entry is np.nan:
            return np.nan
        else:
            return float(entry)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


def isvalid(x=None):
    logid = scriptname+'.isvalid: '
    try:
        if x:
            if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
                return False
            else:
                return True
        else:
            return False
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def isinvalid(x=None):
    logid = scriptname+'.isinvalid: '
    try:
        if x:
            if x in ('None', 'nan', 'none', 'NA', 'NAN') or x is None or x is np.nan:
                return True
            else:
                return False
        else:
            return True
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def makeoutdir(outdir):
    logid = scriptname+'.makeoutdir: '
    try:
        if not os.path.isabs(outdir):
            outdir =  os.path.abspath(outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        return outdir
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def parseseq(sequence):
    logid = scriptname+'.parseseq: '
    try:
        if (isinstance(sequence, StringIO)):
            seq = sequence

        elif ( isinstance(sequence, str) and sequence == 'random' ):
            rand = "\n".join(createrandseq(length, gc, number, alphabet))
            seq = StringIO(rand)
            o = gzip.open('Random.fa.gz','wb')
            o.write(bytes(rand,encoding='UTF-8'))
            o.close()

        elif (isinstance(sequence, str) and os.path.isfile(sequence)):
            if '.gz' in sequence :
                seq = gzip.open(sequence,'rt')
            else:
                seq = open(sequence,'rt')
        else:
            header = ">Seq1:default:nochrom:(.)"
            s = sequence
            seq = StringIO("{header}\n{s}".format(header=header, s=s))

        return seq

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def npprint(a, o=None):#, format_string ='{0:.2f}'):
    logid = scriptname+'.npprint: '
    try:
        out = ''
        it = np.nditer(a, flags=['f_index'])
        while not it.finished:
            out += "%d\t%0.7f" % (it.index+1,it[0])+"\n"
            it.iternext()
        if o:
            o.write(bytes(out,encoding='UTF-8'))
        else:
            print(out)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def idfromfa(id):
    logid = scriptname+'.idfromfa: '
    goi, chrom, strand = [None, None, None]
    try:
        goi, chrom = id.split(':')[::2]
        strand = str(id.split(':')[3].split('(')[1][0])
    except:
        log.error(logid+'Fasta header is not in expected format, you will loose information on strand and chromosome')
        goi = id
        chrom, strand = ['na','na']

    if goi and chrom and strand:
        return [str(goi), str(chrom), str(strand)]
    else:
        log.error(logid+'Could not assign any value from fasta header, please check your fasta files')
        sys.exit('Could not assign any value from fasta header, please check your fasta files')

def print_globaldicts():
    logid = scriptname+'.print_globaldicts: '
    try:
        for name, value in globals().copy().items():
            print(name, value)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))


def print_globallists():
    logid = scriptname+'.print_globallists: '
    try:
        for name, value in globals().deepcopy().items():
            print(name, value)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        log.error(logid+''.join(tbe.format()))

#def bpp_callback(v, v_size, i, maxsize, what, data):
#   if what & RNA.PROBS_WINDOW_BPP:
#       data['bpp'].extend([{'i': i, 'j': j, 'p': p} for j, p in enumerate(v) if (p is not None)])# and (p >= 0.01)])
#
#def up_callback(v, v_size, i, maxsize, what, data):
#   if what & RNA.PROBS_WINDOW_UP:
#       #    data['up'].extend([{ 'i': i, 'up': v}])
#       data['up'].extend([v])
#
# Collection.py ends here
