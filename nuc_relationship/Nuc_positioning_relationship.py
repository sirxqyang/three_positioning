#!/usr/bin/python

from SeqTag import SeqTag
from Genome import chrlength

def Positioning_relationship(tagfile, positionfile, ofile, spename):
    # Both tagfile and Positionfile are BED format files
    
    tagplus = SeqTag(tagfile, spename, extension = 1, shift = 0, strand = '+')
    tagminus = SeqTag(tagfile, spename, extension = 1, shift = 0, strand = '-')
    
    same_strand_value = [0] * 2001
    oppo_strand_value = [0] * 2001
    for line in open(positionfile).xreadlines():
        line = line.strip().split('\t')
        if line[0] not in chrlength[spename].keys():
            continue
        if line[5] == '+':
            chr = line[0]
            position = int(line[1]) - 1
            begin = max((position - 1000), 0)
            end = min(position + 1001, chrlength[spename][chr])
            for i in range(begin, end):
                same_strand_value[1000 + i - position] += tagplus.values[chr][i]
                oppo_strand_value[1000 + i - position] += tagminus.values[chr][i]
        elif line[5] == '-':
            chr = line[0]
            position = int(line[2]) - 1
            begin = max((position - 1000), 0)
            end = min(position + 1001, chrlength[spename][chr])
            for i in range(begin, end):
                same_strand_value[1000 - i + position] += tagminus.values[chr][i]
                oppo_strand_value[1000 - i + position] += tagplus.values[chr][i]
    
    output = open(ofile,'w')
    
    for i in range(2001):
        same_strand_value[i] = str(same_strand_value[i])
        oppo_strand_value[i] = str(oppo_strand_value[i])
    
    output.write('Same_strand' + ' <- c(' + ', '.join(same_strand_value) + ')' + '\n')
    output.write('Oppo_strand' + ' <- c(' + ', '.join(oppo_strand_value) + ')' + '\n')
    
    output.close()
    #return (same_strand_value, oppo_strand_value)

if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    usage = 'usage: python Nuc_positioning_relationship.py options tagfile positionfile'
    parser = OptionParser(usage)
    
    parser.add_option("-s", "--spename", dest="spename",
                      default="zv9", help="species name" )
    parser.add_option("-n", "--name", dest="name",
                      default="NUC", help="name of the data, e.g. ACF_1_pile")
    parser.add_option("-o", "--output", dest="output", help="output file")
    
    (options, args) = parser.parse_args()
    if len(args) < 2:
        parser.print_help()
        sys.exit()

    Positioning_relationship(args[0], args[1], options.output, options.spename)
    #for i in range(2001):
     #   same_strand_value[i] = str(same_strand_value[i])
      #  oppo_strand_value[i] = str(oppo_strand_value[i])
    
    #print 'Same_strand_' + options.name + ' <- c(' + ', '.join(same_strand_value) + ')'
    #print 'Oppo_strand_' + options.name + ' <- c(' + ', '.join(oppo_strand_value) + ')'
