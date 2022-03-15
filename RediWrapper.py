import argparse
import pysam

parser = argparse.ArgumentParser(description="RediWrapper")
# input
parser.add_argument("-I", "--redi", help="Input redi")
parser.add_argument("-R", "--rmsk", help="Input rmsk")
parser.add_argument("-Rt", "--rmsk_tbi", help="Input tbi rmsk")
parser.add_argument("-S", "--snp", help="Input snp")
parser.add_argument("-St", "--snp_tbi", help="Input snp tbi")
parser.add_argument("-A", "--atlas", help="Input atlas")
parser.add_argument("-At", "--atlas_tbi", help="Input atlas tbi")
# output
parser.add_argument("-Ot2", "--output_pos_txt", help="output pos file txt", default='pos.txt')
parser.add_argument('-Ot1', '--output_posalu_txt', help='output posalu file txt', default='posalu.txt')
parser.add_argument("-Og2", "--output_pos_gff", help="output pos file gff", default='pos.gff')
parser.add_argument('-Og1', '--output_posalu_gff', help='output posalu file gff', default='posalu.gff')
parser.add_argument("-Ogz2", "--output_pos_gff_gz", help="output pos file gff.gz", default='pos.gff.gz')
parser.add_argument('-Ogz1', '--output_posalu_gff_gz', help='output posalu gff.gz', default='posalu.gff.gz')
parser.add_argument("-Ogzt2", "--output_pos_gff_gz_tbi", help="output pos file gff.gz.tbi", default='pos.gff.gz.tbi')
parser.add_argument('-Ogzt1', '--output_posalu_gff_gz_tbi', help='output posalu file gff.gz.tbi', default='posalu.gff.gz.tbi')
parser.add_argument("-KN", "--output_known", help="output known editing", default='KnownEditing')

# OTHER
parser.add_argument("-H", "--header", help="header", action= "store_true")


args = parser.parse_known_args()[0]
print(args)

###INPUT
REDI = args.redi
RMSK = args.rmsk
SNP = args.snp
ATLAS = args.atlas
RMSKtbi = args.rmsk_tbi
SNPtbi = args.snp_tbi
ATLAStbi = args.atlas_tbi



###OUTPUT
POStxt = args.output_pos_txt
POSALUtxt = args.output_posalu_txt
POSgff = args.output_pos_gff
POSALUgff = args.output_posalu_gff
POSgffgz = args.output_pos_gff_gz
POSALUgffgz = args.output_posalu_gff_gz
POSgffgztbi = args.output_pos_gff_gz_tbi
POSALUgffgztbi = args.output_posalu_gff_gz_tbi
KN = args.output_known

###OTHER
HEAD = args.header

with open(POStxt, 'w') as OUT1 :
    a = 0
with open(POSALUtxt, 'w') as OUT1 :
    a=0
with open(POSgff, 'w') as OUT1 :
    a=0
with open(POSALUgff, 'w') as OUT1 :
    a=0
with open(POSgffgz, 'w') as OUT1 :
    a=0
with open(POSALUgffgz, 'w') as OUT1 :
    a=0
with open(POSgffgztbi, 'w') as OUT1 :
    a=0
with open(POSALUgffgztbi, 'w') as OUT1 :
    a=0
with open(KN, 'w') as OUT1 :
    a=0


def Set_Chr_Nr_ (Chr):
    """ Sort by chromosome """
    if Chr:
        New = Chr[3:]
        if New == 'X': New = 23
        elif New == 'Y': New = 24
        else: New = int(New)
    else:
        New = 0
    return New

def parse ( res ) :
    d = {'+' : {}, '-' : {}}
    anns = '+'
    for i in res :
        if i[3] == '+' :
            if d['+'].has_key(i[1]) :
                if i[0] not in d['+'][i[1]][0] : d['+'][i[1]][0] = d['+'][i[1]][0] + ',' + i[0]
                if i[2] + '-' + i[0] not in d['+'][i[1]][1] : d['+'][i[1]][1] = d['+'][i[1]][1] + ',' + i[2] + '-' + i[
                    0]
            else :
                d['+'][i[1]] = [i[0], i[2] + '-' + i[0]]
        elif i[3] == '-' :
            if d['-'].has_key(i[1]) :
                if i[0] not in d['-'][i[1]][0] : d['-'][i[1]][0] = d['-'][i[1]][0] + ',' + i[0]
                if i[2] + '-' + i[0] not in d['-'][i[1]][1] : d['-'][i[1]][1] = d['-'][i[1]][1] + ',' + i[2] + '-' + i[
                    0]
            else :
                d['-'][i[1]] = [i[0], i[2] + '-' + i[0]]
    gip = '$'.join(d['+'].keys())
    featp = '$'.join([d['+'][x][0] for x in d['+'].keys()])
    tip = '$'.join([d['+'][x][1] for x in d['+'].keys()])
    gim = '$'.join(d['-'].keys())
    featm = '$'.join([d['-'][x][0] for x in d['-'].keys()])
    tim = '$'.join([d['-'][x][1] for x in d['-'].keys()])
    p = [featp, gip, tip]
    m = [featm, gim, tim]
    pm = [(featp + '&' + featm).strip('&'), (gip + '&' + gim).strip('&'), (tip + '&' + tim).strip('&')]
    if len(d['+']) == 0 and len(d['-']) != 0 : anns = '-'
    if len(d['+']) == 0 : p = ['-', '-', '-']
    if len(d['-']) == 0 : m = ['-', '-', '-']
    if len(d['+']) == 0 and len(d['-']) == 0 :
        pm = ['-', '-', '-']
        anns = '+-'
    if len(d['+']) != 0 and len(d['-']) != 0 : anns = '+-'
    return (p, m, pm, anns)


RMSKtabix = pysam.Tabixfile(RMSK)
RMSKcontig = RMSKtabix.contigs
SNPtabix = pysam.Tabixfile(SNP)
SNPcontig = SNPtabix.contigs
ATLAStabix = pysam.Tabixfile(ATLAS)
ATLAScontig = ATLAStabix.contigs

with open (REDI,'r') as INP1:
    if HEAD:
        header = INP1.readline()
    else:
        for line in INP1:
            lines = line.rstrip('\n').split('\t')
            CHR = lines[0]
            if CHR == 'chrM' : continue
            POS = int(lines[1])
            ### First filter on variant positions (10)
            VAR = lines[7]
            if VAR == '-':
                continue
            REF = lines[2]
            COV = int(lines[4])
            FREQ = float(lines[8])
            ACGT = eval(lines[6])
            STRAND = lines[3]
            nACGT = ['A', 'C', 'G', 'T']
            iACGT = [0, 1, 2, 3]
            ### filter SNP
            snp = []
            if CHR in SNPcontig :
                snp = [(kk.feature, kk.gene_id, kk.transcript_id, kk.strand) for kk in
                        SNPtabix.fetch(reference=CHR, start=POS - 1, end=POS, parser=pysam.asGTF())]
            if len(snp) > 0 : continue
            snp = ['-', '-']
            ### ALU - nonALU, nonREP
            CTRLalu = []
            if CHR in RMSKcontig :
                CTRLalu = [(kk.feature, kk.gene_id, kk.transcript_id, kk.strand) for kk in
                        RMSKtabix.fetch(reference=CHR, start=POS - 1, end=POS, parser=pysam.asGTF())]
            if len(CTRLalu) > 0 :
                ann = parse(CTRLalu)
                CTRLalu_res = ann[2]
            else :
                CTRLalu_res = ['-', '-']
            ### filter ALU
            if CTRLalu_res[1][0:3] == 'Alu':
                ###SEL1
                if COV < 5 : continue
                if FREQ == 0.0: continue
                SUP = 1
                CTRLv = 0
                for nn in iACGT :
                    if nACGT[nn] == REF : continue
                    if ACGT[nn] >= SUP :
                        CTRLv += 1
                if CTRLv == 0 : continue
                ### ED
                if CHR in ATLAScontig :
                    ED = [(kk.feature, kk.gene_id, kk.transcript_id, kk.strand) for kk in
                            ATLAStabix.fetch(reference=CHR, start=POS -1 , end=POS, parser=pysam.asGTF())]
                if len(ED) == 0 :
                    CTRLed = ['-']
                else:
                    CTRLed = ['ed']
                lines.extend(CTRLalu_res[:2])
                lines.extend(snp)
                lines.extend(CTRLed)
                nline = '\t'.join(lines) + '\n'
                if CTRLed[0] == '-':
                    with open(POSALUtxt, 'a') as OUT:
                        OUT.write(nline)
                    if STRAND == '0' :
                        strand = '-'
                    else :
                        strand = '+'
                    CHR_POS = '{0}-{1}'.format(CHR, POS)
                    gffLine = [CHR, 'reditoolTable', 'pos', str(POS), str(POS), '.', strand, '.', CHR_POS]
                    nline = '\t'.join(gffLine) + '\n'
                    with open(POSALUgff, 'a') as OUT:
                        OUT.write(nline)
                else:
                    with open(KN,'a') as OUT:
                        OUT.write(nline)
            else:
                ###SEL2
                if COV < 10 : continue
                if FREQ < 0.1: continue
                SUP = 3
                CTRLv = 0
                for nn in iACGT :
                    if nACGT[nn] == REF : continue
                    if ACGT[nn] >= SUP :
                        CTRLv += 1
                if CTRLv == 0 : continue
                ### ED
                if CHR in ATLAScontig :
                    ED = [(kk.feature, kk.gene_id, kk.transcript_id, kk.strand) for kk in
                            ATLAStabix.fetch(reference=CHR, start=POS -1 , end=POS, parser=pysam.asGTF())]
                if len(ED) == 0 :
                    CTRLed = ['-']
                else :
                    CTRLed = ['ed']
                cSINE = CTRLalu_res[0]
                if cSINE == '-' :
                    lines.extend(CTRLalu_res[:2])
                    lines.extend(snp)
                    lines.extend(CTRLed)
                    nline = '\t'.join(lines) + '\n'
                    if CTRLed[0] == '-' :
                        with open(POStxt, 'a') as OUT :
                            OUT.write(nline)
                        if STRAND == '0' :
                            strand = '-'
                        else :
                            strand = '+'
                        CHR_POS = '{0}-{1}'.format(CHR, POS)
                        gffLine = [CHR, 'reditoolTable', 'pos', str(POS), str(POS), '.', strand, '.', CHR_POS]
                        nline = '\t'.join(gffLine) + '\n'
                        with open(POSgff, 'a') as OUT :
                            OUT.write(nline)
                    else:
                        with open(KN, 'a') as OUT :
                            OUT.write(nline)

                else:
                    if cSINE == 'Simple_repeat' : continue
                    if cSINE == 'Low_complexity' : continue
                    lines.extend(CTRLalu_res[:2])
                    lines.extend(snp)
                    lines.extend(CTRLed)
                    nline = '\t'.join(lines) + '\n'
                    if CTRLed[0] == '-' :
                        with open(POStxt, 'a') as OUT :
                            OUT.write(nline)
                        if STRAND == '0' :
                            strand = '-'
                        else :
                            strand = '+'
                        CHR_POS = '{0}-{1}'.format(CHR, POS)
                        gffLine = [CHR, 'reditoolTable', 'pos', str(POS), str(POS), '.', strand, '.', CHR_POS]
                        nline = '\t'.join(gffLine) + '\n'
                        with open(POSgff, 'a') as OUT :
                            OUT.write(nline)
                    else:
                        with open(KN, 'a') as OUT :
                            OUT.write(nline)


pysam.tabix_index(POSgff, preset='gff', force=True)
pysam.tabix_index(POSALUgff, preset='gff', force=True)
