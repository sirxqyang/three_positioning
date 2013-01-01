#!/usr/bin/python

def get_uniq_reads(inputbed, outputbed):

    reads = {}
    f = open(inputbed, 'r')
    f_uniq = open(outputbed,'w')

    for line in f:

        line_parser = line.strip().split()
        chr = line_parser[0]
        start = line_parser[1]
        end = line_parser[2]
        strand = line_parser[5]

        key = chr+','+start+','+end+','+strand

        reads[key] = 0

    for key in reads.keys():
       key_parser = key.split(',')
       f_uniq.write( key_parser[0] + '\t' + key_parser[1] + '\t' + key_parser[2] + '\t' + '0'+'\t'+'0'+'\t'+key_parser[3]+'\n')

    f_uniq.close()

if __name__=='__main__':
    import sys

    from optparse import OptionParser
    usage = 'usage: python get_uniq_reads.py -t inputbed -o outputbed'
    parser = OptionParser(usage)

    parser.add_option("-t", dest="infile", help="input raw bed" )
    parser.add_option("-o", dest="ofile", help="output bed of uniq location reads")

    (options, args) = parser.parse_args()

    if not options.infile or not options.ofile:   #only required argument
        parser.print_help()
        sys.exit(1)
    
    get_uniq_reads(options.infile, options.ofile)

