### Collection.py ---
##
## Filename: Collection.py
## Description:
## Author: Joerg Fallmann
## Maintainer:
## Created: Thu Sep  6 09:02:18 2018 (+0200)
## Version:
## Package-Requires: ()
## Last-Updated: Tue Nov  5 15:56:40 2019 (+0100)
##           By: Joerg Fallmann
##     Update #: 252
## URL:
## Doc URL:
## Keywords:
## Compatibility:
##
######################################################################
##
### Commentary:
###import os, sys, inspect
# # realpath() will make your script run, even if you symlink it :)
# cmd_folder = os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) ))
# if cmd_folder not in sys.path:
#     sys.path.insert(0, cmd_folder)
#
# # Use this if you want to include modules from a subfolder
# cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath( inspect.getfile( inspect.currentframe() )) )),"lib")
# if cmd_subfolder not in sys.path:
#     sys.path.insert(0, cmd_subfolder)
#
# # Info:
# # cmd_folder = os.path.dirname(os.path.abspath(__file__)) # DO NOT USE __file__ !!!
# # __file__ fails if the script is called in different ways on Windows.
# # __file__ fails if someone does os.chdir() before.
# # sys.argv[0] also fails, because it doesn't not always contains the path.
##
##
######################################################################
##
### Change Log:
##
##
######################################################################
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.
##
######################################################################
##
### Code:
### IMPORTS
import os, sys, inspect

##load own modules
from lib.logger import makelogdir, setup_multiprocess_logger
# Create log dir
makelogdir('logs')
# Define loggers
scriptname=os.path.basename(__file__)
#streamlog = setup_multiprocess_logger(name='', log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')
clog = setup_multiprocess_logger(name=scriptname, log_file='logs/'+scriptname, logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level='WARNING')
from lib.Collection import *

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
#Biopython stuff
from Bio import SeqIO
from Bio.Seq import Seq

# Code:All subs from here on
def create_kmers(choices, length):
    logid = scriptname+'.create_kmers: '
    try:
        #     choices=['A','T','G','C']
        bases = list(choices)
        k = length

        return [''.join(p) for p in itertools.product(bases, repeat=k)]

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def randseq(alphabet, length):
    logid = scriptname+'.randseq: '
    try:
        l=''
        for i in range(length):
            l += str(''.join(choice(alphabet)))
            return str(l)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def weightedrandseq(alphabet, probs , length):
    logid = scriptname+'.weightedrandseq: '
    try:
        seq=[]
        for i in alphabet:
            weight=int(next(probs))
            for a in range(weight):
                seq+=i

        if len(seq) < length:
            for i in range(length-len(seq)):
                seq += str(''.join(choice(alphabet)))

        shuffle(seq)
        return seq
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

def toarray(file, ulim=None):
    logid = scriptname+'.toarray: '
    try:
        if not ulim:
            ulim = 1
        x = np.loadtxt(str(file), usecols = (ulim), delimiter = '\t', unpack = True, converters = {ulim: lambda s: convertcol(s.decode("utf-8"))})
        return x
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

def parse_annotation_bed(bed, annotated=None):
    logid = scriptname+'.parse_annotation_bed: '
    anno = defaultdict(list)
    if os.path.isfile(os.path.abspath(bed)):
        if '.gz' in bed:
            f = gzip.open(bed,'rt')
        else:
            f = open(bed,'rt')
    else:
        f = bed
    try:
        for line in f:
            entries = line.rstrip().split('\t')
            goi = entries[3]
            if annotated:
                start = int(entries[10])
                end   = int(entries[11])-1
            else:
                start = int(entries[1])
                end   = int(entries[2])-1
            anno[str(goi)].append('-'.join([str(start),str(end)]))
        return anno
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def readConstraintsFromBed(bed, linewise=None):
    logid = scriptname+'.readConstraintsFromBed: '
    cons = defaultdict(list)
    try:
        for line in bed:
            entries = line.rstrip().split('\t')
            start = int(entries[1])+1
            end = entries[2]
            goi = entries[3]
            if linewise:
                cons['lw'].append('-'.join([str(start),str(end)]))
            else:
                cons[str(goi)].append('-'.join([str(start),str(end)]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def readPairedConstraintsFromBed(bed, linewise=None):
    logid = scriptname+'.readPairedConstraintsFromBed: '
    cons = defaultdict(list)
    try:
        for line in bed:
            entries = line.rstrip().split('\t')
            if len(entries) % 2:
                raise Exception('Unbalanced paired bed, please make sure the paired bed consists of equal number of fields for both constraint entries')
            else:
                second = int((len(entries)/2)+1)
                clog.debug(logid+'For line '+str(entries)+' start of second bed part is field '+str(second))
            if int(entries[1]) > -1 and int(entries[second]) > -1:
                start_one = int(entries[1])+1
                end_one = entries[2]
                goi = entries[3]
                start_two = int(entries[second])+1
                end_two = int(entries[second+1])
                if linewise:
                    cons['lw'].append('-'.join([str(start_one), str(end_one), str(start_two), str(end_two)]))
                else:
                    cons[str(goi)].append('-'.join([str(start_one), str(end_one), str(start_two), str(end_two)]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def readConstraintsFromCSV(csv, linewise=None):
    logid = scriptname+'.readConstraintsCSV: '
    cons = defaultdict(
        lambda: defaultdict(list)
    )

    try:
        for line in csv:
            entries = split(',',line.rstrip())
            if linewise:
                cons['def'].append('-'.join([str(entries[1]),str(entries[2])]))
            else:
                cons[entries[3]].append('-'.join([str(entries[1]),str(entries[2])]))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def readConstraintsFromGeneric(generic, linewise=None):
    logid = scriptname+'.readConstraintsFromGeneric: '
    cons = defaultdict(
        lambda: defaultdict(list)
    )

    try:
        for line in csv:
            entries = re.split(r'[ ,|;"]+', line.rstrip())
            if len(entries > 2):
                if linewise:
                    cons['lw'].append('-'.join([str(entries[1]),str(entries[2])]))
                else:
                    cons[entries[0]].append('-'.join([str(entries[1]),str(entries[2])]))
            else:
                if linewise:
                    cons['lw'].append('-'.join([str(entries[1]),str(entries[2])]))
                else:
                    cons['generic'].append('-'.join(str(entries[1:2])))
        return cons
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

def plot_data(fa, raw, consu, consp, const, xs, cons, saveas, outdir):
    logid = scriptname+'.plot_data: '
    anime = []
    #define xs for constraint line
    consl = []
    try:
        for x in const:
            consl.append(1.25)
        width = 16/100*len(fa.seq)
        height = 9
        fig = plt.figure(figsize=(width,height),dpi=80)
        ax1 = fig.add_subplot(111)

        ax2 = ax1.twiny()
        #   line, = ax.plot([], [], lw=2)
        plt.title("Blue-- = Unconstraint, Green-. = Unpaired, Red = Paired, Gray = Constraint",y=1.075)
        ax1.set_ylabel('Prob unpaired')
        ax1.set_xlabel('Nucleotides')
        #   plt.xticks(range(0,len(fa.seq)+1),(' '+fa.seq),size='small')
        #add lines to plot
        ax1.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-')
        ax1.set_xlim(0,len(fa.seq)+1)
        ax2.set_xlim(ax1.get_xlim())
        ax1.set_xticks(range(0,len(fa.seq)+1))
        ax1.set_xticklabels((' '+fa.seq), ha="right")
        #   ax2.set_xlabel(r"Modified x-axis: $1/(1+X)$")
        ax2.set_xticks(range(1,len(fa.seq)+1))
        ax2.set_xticklabels(range(1,len(fa.seq)+1), rotation=45, ha="right")
        # We change the fontsize of minor ticks label
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax1.tick_params(axis='both', which='minor', labelsize=4)
        ax2.tick_params(axis='both', which='major', labelsize=5)
        ax2.tick_params(axis='both', which='minor', labelsize=3)
        goi, chrom = fa.id.split(':')[::2]
        strand = str(fa.id.split(':')[3].split('(')[1][0])
        fig.savefig('StruCons_'+goi+'_'+cons+'.'+saveas)
        plt.close()
        #   anime.append(plt.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-'))
        #   return anime
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def plot_temp(fa, raw, temp, xs, saveas, outdir):
    logid = scriptname+'.plot_temp: '
    try:
        anime = []
        #define xs for constraint line
        width = 16/100*len(fa.seq)
        height = 9
        fig = plt.figure(figsize=(width,height),dpi=80)
        ax1 = fig.add_subplot(111)
        plt.title("Blue-- = "+temp+" degree",y=1.075)
        ax1.set_ylabel('Prob unpaired')
        ax1.set_xlabel('Nucleotides')
        #add lines to plot
        ax1.plot(xs, raw, 'b-')
        ax1.set_xlim(0,len(fa.seq)+1)
        ax1.set_xticks(range(0,len(fa.seq)+1))
        ax1.set_xticklabels((' '+fa.seq))
        # We change the fontsize of minor ticks label
        ax1.tick_params(axis='both', which='major', labelsize=8)
        ax1.tick_params(axis='both', which='minor', labelsize=4)
        goi, chrom = fa.id.split(':')[::2]
        strand = str(fa.id.split(':')[3].split('(')[1][0])
        fig.savefig('TempCons_'+goi+'_'+temp+'.'+saveas)
        plt.close()
        #   anime.append(plt.plot(xs, raw, 'b-', xs, consu, 'g-', xs, consp, 'r-', const, consl, 'k-'))
        #   return anime
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def calc_gibbs(fc):
    logid = scriptname+'.calc_gibbs: '
    try:
        return fc.pf()[1]
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        clog.error(logid+''.join(tbe.format()))

def get_bppm(tmp, start, end):
    logid = scriptname+'.get_bppm: '
    bppm = []
    try:
#        clog.debug('\t'.join(map(str,[tmp,start,end])))
        if start < 0 or end > len(tmp):
            clog.warning(logid+'start of constraint '+str(start)+' end of constraint '+str(end)+' while length of bpp matrix '+str(len(tmp))+'! Skipping!')
            return None

        for item in tmp:
            for i in range(int(start),int(end)+1):
                if item[i] > 0.0:
                    bppm.append(str.join('\t',[str(tmp.index(item)), str(i), str(item[i])]))
        return bppm
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def get_ddg(file):
    logid = scriptname+'.parseseq: '
    try:
        ret = defaultdict()
        if (isinstance(file, str) and os.path.isfile(file)):
            if '.gz' in file :
                res = gzip.open(file,'rt')
            else:
                res = open(file,'rt')

            for line in res:
                clog.debug(logid+line)
                if 'Condition' in line[0:15]:
                    continue
                else:
                    cond, gibbs, dg, nrg, cons = line.rstrip().split('\t')
                    if not str(cons) in ret:
                        ret[str(cons)] = defaultdict()
                    ret[str(cons)][cond] = float(dg)
        return ret

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def calc_ddg(ddgs):
    logid = scriptname+'.calc_ddg: '

    try:
        clog.debug(logid+str(ddgs))
        ddg = ddgs['constraint_unpaired']+ddgs['secondconstraint_unpaired']-ddgs['bothconstraint_unpaired']-ddgs['unconstraint']
        """Yi-Hsuan Lin, Ralf Bundschuh, RNA structure generates natural cooperativity between single-stranded RNA binding proteins targeting 5' and 3'UTRs, Nucleic Acids Research, Volume 43, Issue 2, 30 January 2015, Pages 1160-1169, https://doi.org/10.1093/nar/gku1320"""

        return ddg

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        clog.error(logid+''.join(tbe.format()))

def calc_bpp(bppm):
    logid = scriptname+'.calc_bpp: '
    bpp = 0.0;
    try:
        for entry in bppm:
            base, mate, prob = map(float,entry.split('\t'))
            bpp += prob
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        clog.error(logid+''.join(tbe.format()))

    return bpp

def calc_nrg(bpp):
    logid = scriptname+'.calc_nrg: '
    #set kT for nrg2prob and vice versa calcs
    kT = 0.61632077549999997
    nrg = 0.0;
    try:
        if bpp > 0.0:
            nrg = -1 * kT * math.log(bpp)
        return nrg
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def print_region_up(data, seqlength=None, region=None):
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    #   pp.pprint(data)
    logid = scriptname+'.print_region_up: '
    try:
        if data:
            ups=''
            x = int(region)
            for i in range(int(seqlength)):
                if isinvalid(data[i][x]):
                    data[i][x] = np.nan
                else:
                    data[i][x] = round(data[i][x],7)
                ups+=str(i+1)+"\t"+str(data[i][x])+"\n"
            return ups
        else:
            clog.error(logid+'No up data to print')
            return ups

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def print_up(data=None, seqlength=None, region=None):
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    #   pp.pprint(data)
    logid = scriptname+'.print_up: '
    try:
        if data:
            ups=''
            for i in range(int(seqlength)):
                for x in range(1,region+1):
                    if isinvalid(data[i][x]):
                        data[i][x] = np.nan
                    else:
                        data[i][x] = round(data[i][x],7)
                ups+=str(i+1)+"\t"+"\t".join(map(str,data[i][1:region+1]))+"\n"
            return ups
        else:
            clog.error(logid+'No up data to print')
            return ups
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def up_to_array(data=None, region=None, seqlength=None):
    #   pp = pprint.PrettyPrinter(indent=4)#use with pp.pprint(datastructure)
    #   pp.pprint(data[165553:165588])
    logid = scriptname+'.up_to_array: '
    try:
        if data:
            entries=[]
            if not seqlength:
                seqlength = len(data)
            if not region:
                region = slice(1,len(data[0]))
            for i in range(seqlength):
                entries.append([])
                for e in range(len(data[i])):
                    if isinvalid(data[i][e]):
                        data[i][e] = np.nan
                    else:
                        data[i][e] = round(data[i][e],8)
                entries[i].extend(data[i][region])
            return np.array(entries)
        else:
            clog.error(logid+'No up data to print')
            return np.array()
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))

def printdiff(a, o=None):
    logid = scriptname+'.printdiff: '
    try:
        np.save(o, a)
        #np.savetxt(o, a, delimiter='\t', encoding='bytes')
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def read_precalc_plfold(data, name, seq):
    logid = scriptname+'.read_precalc_plfold: '
    try:
        for i in range(len(seq)):
            data.append([])
            data[i] = []
        with gzip.open(name,'rt') as o:
            for line in o:
                cells = line.rstrip().split('\t')
                data[int(cells[0])-1].append([])
                data[int(cells[0])-1][0] = None
                for a in range(1,len(cells)):
                    data[int(cells[0])-1].append([])
                    data[int(cells[0])-1][a] = float(cells[a])
        return data
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        clog.error(logid+''.join(tbe.format()))

def pl_to_array(name, ulim, fmt='npy'):
    logid = scriptname+'.pl_to_array: '
    try:
        clog.debug('\t'.join([logid,name]))
        if fmt == 'txt':
            return np.array(np.loadtxt(name, usecols=ulim, unpack=True, delimiter='\t', encoding='bytes'))
        elif fmt == 'npy':
            return np.array(np.load(name)[:,ulim])
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
        exc_type, exc_value, exc_tb,
        )
        clog.error(logid+' '+name+': '.join(tbe.format()))

def idfromfa(id):
    logid = scriptname+'.idfromfa: '
    goi, chrom, strand = [None, None, None]
    try:
        goi, chrom = id.split(':')[::2]
        strand = str(id.split(':')[3].split('(')[1][0])
    except:
        clog.error(logid+'Fasta header is not in expected format, you will loose information on strand and chromosome')
        goi = id
        chrom, strand = ['na','na']

    if goi and chrom and strand:
        return [str(goi), str(chrom), str(strand)]
    else:
        clog.error(logid+'Could not assign any value from fasta header, please check your fasta files')
        sys.exit('Could not assign any value from fasta header, please check your fasta files')

def constrain_paired(fc, start, end):
    logid = scriptname+'.constrain_paired: '
    try:
        for x in range(start, end+1):
            fc.hc_add_bp_nonspecific(x,0) #0 means without direction  ( $ d < 0 $: pairs upstream, $ d > 0 $: pairs downstream, $ d == 0 $: no direction)
        return fc
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        clog.error(logid+''.join(tbe.format()))

def constrain_unpaired(fc, start, end):
    logid = scriptname+'.constrain_unpaired: '
    try:
        for x in range(start, end+1):
            fc.hc_add_up(x)
        return fc
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
            )
        clog.error(logid+''.join(tbe.format()))

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
        clog.error(logid+''.join(tbe.format()))


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
        clog.error(logid+''.join(tbe.format()))

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
