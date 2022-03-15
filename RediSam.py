import argparse
import os
import math

parser = argparse.ArgumentParser(description="RediSam")
# input
parser.add_argument("-Ir", "--reads", help="Input redi reads")
parser.add_argument("-Ip", "--posreads", help="Input redi pos reads")
parser.add_argument("-B", "--bam", help="Input_bam")
parser.add_argument("-Bb", "--bam_bai", help="Input bam index")
parser.add_argument("-R", "--reference", help="Input reference index")
parser.add_argument("-Rf", "--reference_fai", help="Input reference index")


# output
parser.add_argument("-Or", "--output_reads_psl", help="output reads.psl", default='reads.psl')
parser.add_argument("-Ob", "--output_badreads", help="output badreads", default='badreads.txt')
parser.add_argument("-Oe", "--output_bed", help="output bed", default='bed')
parser.add_argument("-O1", "--output_1", help="output 1", default='bed.bam')
parser.add_argument("-O2", "--output_2", help="output 2", default='bed_ns.bam')
parser.add_argument("-O3", "--output_3", help="output 3", default='bed_ns_fx.bam')
parser.add_argument("-O4", "--output_4", help="output 4", default='bed_ns_fx_st.bam')
parser.add_argument("-O5", "--output_5", help="output 5", default='dedup.bam')
parser.add_argument("-O6", "--output_6", help="output 6", default='dedup.bam.bai')



args = parser.parse_known_args()[0]
print(args)


##INPUT
READS = args.reads
POSREADS = args.posreads
BAM = args.bam
BAMbai = args.bam_bai
REFERENCE = args.reference
REFERENCEfai = args.reference_fai

##OUTPUT
READSpsl = args.output_reads_psl
BADREADS = args.output_badreads
BED = args.output_bed
OT1 = args.output_1
OT2 = args.output_2
OT3 = args.output_3
OT4 = args.output_4
OT5 = args.output_5
OT6 = args.output_6


with open(READSpsl,'w') as OUT:
    a=0
with open(BADREADS,'w') as OUT:
    a=0
with open(BED, 'w') as OUT:
    a=0
with open(OT1, 'w') as OUT:
    a=0
with open(OT2, 'w') as OUT:
    a=0
with open(OT3, 'w') as OUT:
    a=0
with open(OT4, 'w') as OUT:
    a=0
with open(OT5, 'w') as OUT:
    a=0
with open(OT6, 'w') as OUT:
    a=0

# pblat -t=dna -q=rna -stepSize=5 -repMatch=2253 -min- Score=20 -minIdentity=0genome.fa outReads reads.psl
cmd = 'pblat -t=dna -q=rna -stepSize=5 -repMatch=2253 -minScore=20 -minIdentity=0 %s %s %s' %(REFERENCE, READS, READSpsl)
os.system(cmd)


#### for blat
def getPS ( line ) :
    pid = (100.0 - (pslCalcMilliBad(line) * 0.1))
    score = pslScore(line)
    # print "The percentage:",pid
    # print "Score:",score
    return pid, score
def pslScore ( cols ) :
    sizeMul = 1
    return sizeMul * (int(cols[0]) + (int(cols[2]))) - sizeMul * int(cols[1]) - int(cols[4]) - int(cols[6])
def round ( number ) :
    return int(number + .5);
def pslCalcMilliBad ( cols ) :
    sizeMul = 1
    # cols[0]  matches
    # cols[1]  misMatches
    # cols[2]  repMaches
    # cols[4]  qNumInsert
    # cols[6]  tNumInsert
    # cols[11] qStart
    # cols[12] qEnd
    # cols[15] tStart
    # cols[16] tEnd
    qAliSize = sizeMul * (int(cols[12]) - int(cols[11]))
    tAliSize = int(cols[16]) - int(cols[15])
    # I want the minimum of qAliSize and tAliSize
    if qAliSize < tAliSize :
        aliSize = qAliSize  # ? $aliSize = $qAliSize : $aliSize =  $tAliSize;
    else :
        aliSize = tAliSize
    # return 0 is AliSize == 0
    if aliSize <= 0 : return 0
    # size diff
    sizeDiff = qAliSize - tAliSize
    if sizeDiff < 0 : sizeDiff = 0
    # insert Factor
    insertFactor = int(cols[4])
    # $insertFactor += $cols[6];
    milliBad = (1000 * (int(cols[1]) * sizeMul + insertFactor + round(3 * math.log(1 + sizeDiff)))) / (
                sizeMul * (int(cols[0]) + int(cols[2]) + int(cols[1])))
    return milliBad
def com ( num, list ) :
    for i in list :
        if i[0] <= num <= i[1] : return 1
    return 0
def min95 ( val, score ) :
    # if val < (score*95.0)/100: return 1
    if val < (score * 0.95) : return 1
    return 0
def readLines ( lines ) :
    res = []
    for line in lines :
        pidd, score = getPS(line)
        # print pidd,score
        sp = [int(x) for x in (line[18].strip(',')).split(',')]
        tstarts = [int(x) for x in (line[20].strip(',')).split(',')]
        ex = [(tstarts[x] + 1, tstarts[x] + sp[x]) for x in range(len(sp))]
        nl = [line[9], score, str(int(line[11]) + 1), line[12], str(line[10]), pidd, line[13], line[8],
              int(line[15]) + 1, int(line[16]), ex, int(line[0])]
        res.append((int(line[0]), nl))  # score
    # if d.has_key(line[9]): d[line[9]].append((score,nl))
    # else: d[line[9]]=[(score,nl)]
    return res
def comp ( ri, hits ) :
    g, ng = 0, 0
    hits.sort()
    hits.reverse()
    if len(hits) == 1 :  # unique hit with editing candidate position included
        if hits[0][1][6] == ri[2] and com(ri[1], hits[0][1][10]) :
            g += 1  # float(hits[0][1][5])>=90.0
        else :
            ng += 1
    elif len(hits) > 1 :  # multiple hits
        if hits[0][1][6] == ri[2] and min95(hits[1][0],
                                            hits[0][0]) :  # if second best score less than 95% of first best score
            if com(ri[1], hits[0][1][10]) :
                g += 1  # if first best hit include editing position
            else :
                ng += 1
        else :
            ng += 1
    if g > ng : return 1
    return 0
def readPSL ( infile, outfile ) :
    f = open(infile)
    o = open(outfile, 'w')
    name, lines, xx = '', [], 0
    while 1 :
        line = f.readline()
        if not line :
            if name == '' : break
            nn = name.split('$')
            oread = (name, int(nn[2]), nn[1])
            bread = readLines(lines)
            badr = ''
            if len(bread) == 0 :
                badr = name
            else :
                if not comp(oread, bread) : badr = name
            if badr != '' :
                # o.write(name[:-2]+' '+name[-1]+'\n')
                o.write(name.split('_')[0] + ' ' + name.split('$')[0][-1] + '\n')
                xx += 1
            break
        if line.strip() == '' : continue
        if line.startswith('psL') : continue
        if (line.strip()).startswith('match') : continue
        if line.startswith('-') : continue
        l = (line.strip()).split('\t')
        if l[9] != name :
            if len(lines) != 0 :
                nn = name.split('$')
                # (rname,pileupcolumn.pos+1,chr)
                oread = (name, int(nn[2]), nn[1])  # dread[name]
                bread = readLines(lines)
                badr = ''
                if len(bread) == 0 :
                    badr = name
                else :
                    if not comp(oread, bread) : badr = name
                if badr != '' :
                    # o.write(name[:-2]+' '+name[-1]+'\n')
                    o.write(name.split('_')[0] + ' ' + name.split('$')[0][-1] + '\n')
                    xx += 1
            lines = [l]
            name = l[9]
        else :
            lines.append(l)
    f.close()
    o.close()
    return xx
def readgf ( infile ) :
    f = open(infile)
    for i in f :
        if 'Server ready for queries!' in i :
            f.close()
            return 1
    f.close()
    return 0
def parse ( line ) :
    l = (line.strip()).split('\t')
    cc = (int(l[3]), int(l[4]))
    return cc
# readPsl.py reads.psl badreads.txt
readPSL(READSpsl,BADREADS)
# sort -k1,1 -k2,2n -k3,3n first/DnaRna_51144481/outPosReads_51144481 | mergeBed > bed
cmd = 'sort -k1,1 -k2,2n -k3,3n %s | mergeBed > %s' %(POSREADS, BED)
os.system(cmd)

# samtools view -@ 4 -L bed -h -b../../../Alignment/SRR1258218_Aligned.sortedByCoord.out.bam > SRR1258218_bed.bam
cmd = 'samtools view -L %s -h -b %s > %s' %(BED, BAM, OT1)
os.system(cmd)
cmd = 'rm %s' %(BED)
os.system(cmd)
# samtools sort -@ 4 -n SRR1258218_bed.bam -o SRR1258218_bed_ns.bam
cmd = 'samtools sort -n %s -o %s' %(OT1, OT2)
os.system(cmd)
cmd = 'rm %s' %(OT1)
os.system(cmd)
# samtools fixmate -@ 4 -m SRR1258218_bed_ns.bam SRR1258218_bed_ns_fx.bam
cmd = 'samtools fixmate -@ 4 -m %s %s' %(OT2, OT3)
os.system(cmd)
cmd = 'rm %s' %(OT2)
os.system(cmd)
# samtools sort -@ 4 SRR1258218_bed_ns_fx.bam -o SRR1258218_bed_ns_fx_st.bam
cmd = 'samtools sort %s -o %s' %(OT3, OT4)
os.system(cmd)
cmd = 'rm %s' %(OT3)
os.system(cmd)
# samtools markdup -r -@ 4 SRR1258218_bed_ns_fx_st.bam SRR1258218_bed_dedup.bam
cmd = 'samtools markdup -r %s %s' %(OT4, OT5)
os.system(cmd)
cmd = 'rm %s' %(OT4)
os.system(cmd)
# samtools index SRR1258218_bed_dedup.bam
cmd = 'samtools index %s' %(OT5)
os.system(cmd)