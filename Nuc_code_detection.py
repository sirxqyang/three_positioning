#!/usr/bin/python

from Genome import whole_genome
from Genome import chrlength

def Nucleosome_code(positionfile, spename):
    # Positionfile is BED format file
    genome = whole_genome(spename)
    seq = genome.getnucleosomecode()
    aa = [0]*251
    tt = [0]*251
    total_aa = 0
    total_tt = 0
    aattat = [0] * 251
    total = 0
    for line in open(positionfile).xreadlines():
        line = line.strip().split('\t')
        if line[0] not in chrlength[spename].keys():
            continue
        if line[5] == '+':
            chr = line[0]
            position = int(line[1]) - 1
            begin = position - 50
            if begin < 0:
                continue
            end = position + 201
            if end > chrlength[spename][chr]:
                continue
            total += 1
            subseq = seq[chr][begin:end]
            for i in range(251):
                aattat[i] += int(subseq[i])
        elif line[5] == '-':
            chr = line[0]
            position = int(line[2]) - 1
            begin = position - 200
            if begin < 0:
                continue
            end = position + 51
            if end > chrlength[spename][chr]:
                continue
            total += 1
            subseq = ''
            subseqold = seq[chr][begin:end]
            for i in reversed(subseqold):
                subseq += i
            for i in range(251):
                aattat[i] += int(subseq[i])
    for i in range(251):
        aattat[i] = aattat[i] / float(total)
    return aattat

if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    
    usage = 'usage: python Nuc_code_detection.py options positionfile'
    parser = OptionParser(usage)
    
    parser.add_option("-s", "--spename", dest="spename",
                      default="yeast_SGD", help="species name" )
    parser.add_option("-n", "--name", dest="name",
                      default="NUC", help="name of the data, e.g. ACF_1_pile")
        
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_help()
        sys.exit()
    aattat = Nucleosome_code(args[0], options.spename)
    for i in range(251):
        aattat[i] = str(aattat[i])
    print 'Code' + ' <- c(' + ', '.join(aattat) + ')'
    #thing =
    #f = open('thing', 'w')
    #f.write=(thing)
