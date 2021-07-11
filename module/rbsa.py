#!/usr/bin/env python
# coding: utf-8

import os
from collections import defaultdict
import argparse
import sys


def nucmer_align(ref,scf):
    with open("nucmer.sh","w") as f:
        f.write('''
nucmer --mum -c 200 -l 100 {:s} {:s}  -p scf
delta-filter -i 85 -l 100 -1 scf.delta > scf.flt.delta
show-coords -cdlqoTH scf.flt.delta > scf.flt.coords
mummerplot --fat --large --png -p {:s} scf.flt.delta
        '''.format(ref,scf,scf))
    os.system('sh nucmer.sh')

def anchor2pos():
    pos = open("ScafInChr.pos","w")
    with open("scf.flt.coords","r") as f:
        for each in f:
            line = each.strip().split()
            if int(line[3]) > int(line[2]):
                pos.write(line[13]+"\t"+line[14]+"\t"+line[14]+"_"+line[2]+"_"+line[3]+"\t"+line[2]+"\t"+line[3]+"\t+\t"\
                + str(int(line[3])-int(line[2]))+"\t"+line[0]+"\t"+line[1]+"\n")
            else:
                 pos.write(line[13]+"\t"+line[14]+"\t"+line[14]+"_"+line[3]+"_"+line[2]+"\t"+line[3]+"\t"+line[2]+"\t-\t"\
                + str(int(line[2])-int(line[3]))+"\t"+line[0]+"\t"+line[1]+"\n")
    pos.close()


def mcscan_anchor2pos(A,B,A_B_anchor):
     A_bed = {}
     B_bed = {}
     with open(A,"r") as f:
          for each  in f:
             line = each.strip().split()
             A_bed[line[3]] = [line[3],line[0],line[1],line[2]]
        
     with open(B,"r") as f:
         for each  in f:
             line = each.strip().split()
             B_bed[line[3]] = [line[3],line[0],line[1],line[2]]
     pos = open("ScafInChr.pos","w")
     with open(A_B_anchor,"r") as f:
         for each in f:
             line = each.strip().split()
             i = 0
             while i < int(line[4]):
                 i += 1
                 pos.write(A_bed[line[0]][1]+'\t'+ B_bed[line[2]][1] + "\t" + \
                 A_bed[line[0]][0]+"_"+ A_bed[line[1]][0] + "_" + B_bed[line[2]][0] + "_" + B_bed[line[3]][0] + "\t" + \
                 B_bed[line[2]][2] + "\t" + B_bed[line[3]][3] + "\t" + line[5] + "\t" + \
                 str(int(A_bed[line[1]][3])-int(A_bed[line[0]][2])) + "\t" + A_bed[line[0]][2] + "\t" + A_bed[line[1]][3] + "\n")
     pos.close()




def conflict_scaf():
    #Identify scaffold in conflict
    scafInChr = {}
    exclude = {}
    with open('scaf.id','w') as w:
        with open('ScafInChr.pos','r')  as f:
            for each in f:
                line = each.strip().split()
                if line[1] not in  scafInChr:
                    scafInChr[line[1]] = line[0]
                else:
                    if scafInChr[line[1]] == line[0]:
                        continue
                    else:
                        if line[1] not in exclude:
                            w.write(line[1]+"\n")
                            exclude[line[1]] = 0         

def rm_marker():
    
    #read scaf id
    scaffold = {}
    with open('scaf.id','r') as f:
        for each in f:
            line = each.strip().split()
            scaffold[line[0]] = 0
    scfNum = {}
    scfOnChr = {}
    scfOnChrNum = defaultdict(dict)
    
    #Count the number of scaffold on each chromosome
    with open('ScafInChr.pos','r')  as f:
        for each in f:
            line = each.strip().split()
            if line[1]  in scaffold:
                if line[1] not in scfOnChrNum or line[0] not in scfOnChrNum[line[1]]:
                    scfOnChrNum[line[1]][line[0]] = 0
                else:
                    scfOnChrNum[line[1]][line[0]] += 1
    
    # assign Scaffold to chromosomes
    for scf in scfOnChrNum:
        for chr in scfOnChrNum[scf]:
            if scf not in scfNum:
                scfNum[scf] = scfOnChrNum[scf][chr]
                scfOnChr[scf] = chr
            elif scfOnChrNum[scf][chr] > scfNum[scf]:
                scfNum[scf] = scfOnChrNum[scf][chr]
                scfOnChr[scf] = chr
    
    #output result
    with open('ScafInChr.pos.rm.pos','w') as pos:
        with open('rm_marker.lst','w') as rm:
            with open('ScafInChr.pos','r') as f:
                for each in f:
                    each = each.strip()
                    line = each.split()
                    if line[1] not in scaffold:
                        pos.write(each+"\n")
                    else:
                        if line[0] in scfOnChr[line[1]] and line[0] != '':
                            pos.write(each+'\n')
                        else:
                            rm.write(each+'\n')
                                     
    

def ordering_orientation():
    #Determine the orientation and position of scaffolds on chromosome
    ori_plus = defaultdict(int)
    ori_min = defaultdict(int)
    scaf_len = {}
    scaf_cm = {}
    scaf_chr = {}
    with open("ScafInChr.pos.rm.pos","r") as f:
        for  each in f:
            line = each.strip().split()
            if line[5] == '+':
                ori_plus[line[1]] += 1
            else:
                ori_min[line[1]] += 1
            
            if line[1] not in scaf_len or scaf_len[line[1]] < line[6]:
                scaf_len[line[1]] = line[6]
                scaf_cm[line[1]] = line[7]+"\t"+line[8]
                scaf_chr[line[1]] = line[0]
                
    #output
    with open("CHR.scaf.chrpos","w") as f:
        for key in scaf_chr:
            if ori_plus[key] >= ori_min[key]:
                ori = '1'
            else:
                ori = '-1'
            f.write(scaf_chr[key]+"\t"+key+"\t"+scaf_cm[key]+"\t"+ori+"\n")
    
    # order scaffolds according its position on the chromosome
    os.system('cat CHR.scaf.chrpos  |sort -k1,1 -k 3,3g > CHR.scaf.chrpos.sort')
                
            

def unplaced_scaf(scf):
    #Find the unplaced scaffold's ids
    scafOnChr = {}
    with open("CHR.scaf.chrpos.sort","r") as f:
        for each in f:
            line = each.strip().split()
            scafOnChr[line[1]] = line[0]
    
    with open("non-chr.lst","w") as out:
        with open(scf,"r") as f:
            for each in f:
                line = each.strip().split()
                if  line  and line[0].find(">") != -1:
                    line[0] = line[0].lstrip(">")
                    if line[0] not in scafOnChr:
                        out.write("0\t"+line[0]+"\n")

def reverse(seq):
    #Take the reverse complementary sequence of the sequence
    seq = ''.join(reversed(seq))
    trans = str.maketrans('atcgATCG','tagcTAGC')
    return seq.translate(trans)
    

def anchor(scaffold):
    chrom = defaultdict(list)
    scf = defaultdict(str)
    #Store scaffold and orientation in the order of chromosomes
    with open("CHR.scaf.chrpos.sort","r") as f:
        for  each in f:
            line = each.strip().split()
            chrom[line[0]].append([line[1],line[-1]])
    
    #Storing scaffold sequences
    with open(scaffold,"r") as f:
        for each in f:
            line = each.strip().split()
            if line and line[0].find(">") == 0:
                line[0] = line[0].lstrip(">")
                name = line[0]
                scf[name] = ''
            elif line:
                scf[name]+=line[0]
    #Export the agp file and the chromosome file            
    sequence = ''
    strand = ''
    with open("Chr.ScafCut.fa","w") as fa:
        with open("ScafInChr.lst.agp","w") as f:
            
            for key in chrom.keys():
                sequence = ''
                strand = ''
                i = -1
                for scaffold,tag in chrom[key]:
                    i+=1
                    j = str(i)
                    if sequence != '':
                        start = len(sequence) + 1
                        end = str(start + 100)
                        start = str(start)
                        sequence += 'N'*100
                        f.write(key+"\t"+start+"\t"+end+"\t"+j+"\tU\t100\tcontig\tno\tna\n")

                    start = len(sequence) + 1
                    if tag == '-1':
                        sequence += reverse(scf[scaffold])
                        strand = '-'
                    else:
                        sequence += scf[scaffold]
                        strand = '+'
                    end = len(sequence)
                    length = str(end-start+1)
                    start = str(start)
                    end = str(end)
                    i+=1
                    j = str(i)
                    f.write(key+"\t"+start+"\t"+end+"\t"+j+"\tW\t"+scaffold+"\t1\t"+length+"\t"+strand+"\n")

                fa.write(">"+key+"\n"+sequence+"\n")
           
            #Output unanchored scaffold sequences 
            with open("non-chr.lst","r") as f:
                for each in f:
                    line =each.strip().split()
                    if line[1] in scf:
                        fa.write(">"+line[1]+"\n"+scf[line[1]]+"\n")
            unplace = open("unplace.agp","w")
            with open("non-chr.lst","r") as f:
                for each in f:
                    line =each.strip().split()
                    if line[1] in scf:
                        unplace.write("unpl_"+line[1]+"\t1\t"+str(len(scf[line[1]]))+"\t1\tW\t"+line[1]+"\t1\t"+str(len(scf[line[1]]))+"\t+\n")
            unplace.close()

if __name__ == '__main__':
    parser =  argparse.ArgumentParser(usage = """
    rbsa.py  --type nucmer --ref REF --scf SCF
    or rbsa.py  --type mcscan --A_bed A_BED --B_bed B_BED --anchor ANCHOR --scf SCF
    """,description = "Using reference genomes  to anchor scaffolds")
    parser.add_argument('--type',type=str,help="nucmer or mcscan anchor",choices = ["nucmer","mcscan"])
    parser.add_argument('--ref',type=str,help="This is reference genome sequences which has  anchored to the chromosomal level")
    parser.add_argument('--scf',type=str,help="This is the scaffold sequences ",required = True )
    parser.add_argument('--A_bed',type=str,help="species_a bed file")
    parser.add_argument('--B_bed',type=str,help="species_b bed file")
    parser.add_argument('--anchor',type=str,help="anchor.simple files of species_a and species_b")
    arg = vars(parser.parse_args())
    if arg["type"] == "nucmer":
        nucmer_align(arg["ref"],arg["scf"])
        anchor2pos()
        conflict_scaf()
        rm_marker()
        ordering_orientation()
        unplaced_scaf(arg["scf"])
        anchor(arg["scf"])
        nucmer_align(arg["ref"],"Chr.ScafCut.fa")
    else:
        mcscan_anchor2pos(arg["A_bed"],arg["B_bed"],arg["anchor"])
        conflict_scaf()
        rm_marker()
        ordering_orientation()
        unplaced_scaf(arg["scf"])
        anchor(arg["scf"])
        nucmer_align(arg["ref"],"Chr.ScafCut.fa")
        

    











