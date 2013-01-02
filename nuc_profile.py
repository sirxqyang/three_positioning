chrlength = {}
chrlength['Yeast'] = {'chr1' : 230208,  'chr2' : 813178,  'chr3' : 316617,  'chr4' : 1531918,
'chr5' : 576869,  'chr6' : 270148,  'chr7' : 1090947, 'chr8' : 562643,
'chr9' : 439885,  'chr10' : 745745,  'chr11' : 666454,  'chr12' : 1078175,
'chr13' : 924429,  'chr14' : 784333,  'chr15' : 1091289, 'chr16' : 948062,
'chrM' : 85779}

chrlength['zv9'] = { 'chr1':60348388,'chr2':60300536,'chr3':63268876,'chr4':62094675,
'chr5':75682077,'chr6':59938731,'chr7':77276063,'chr8':56184765,'chr9':58232459,'chr10':46591166,'chr11':46661319,
'chr12':50697278,'chr13':54093808,'chr14':53733891,'chr15':47442429,'chr16':58780683,'chr17':53984731,'chr18':49877488,
'chr19':50254551,'chr20':55952140,'chr21':44544065,'chr22':42261000,'chr23':46386876,'chr24':43947580,'chr25':38499472 }

chrlength['hg19'] = { 'chr1':249250621,'chr2':243199373,'chr3':198022430,'chr4':191154276,
'chr5':180915260,'chr6':171115067,'chr7':159138663,'chr8':146364022,'chr9':141213431,'chr10':135534747,'chr11':135006516,
'chr12':133851895,'chr13':115169878,'chr14':107349540,'chr15':102531392,'chr16':90354753,'chr17':81195210,'chr18':78077248,
'chr19':59128983,'chr20':63025520,'chr21':48129895,'chr22':51304566,'chrX':155270560,'chrY':59373565,'chrM':16571 }


def tag2hash(bedfile, spename, extension = 73, shift = 36):
    values = {}
    for chrname in chrlength[spename].keys():
        values[chrname] = [0] * chrlength[spename][chrname]
        
    for line in open(bedfile).xreadlines():
        line = line.strip().split()
        if line[0] not in chrlength[spename].keys():
            continue
        if line[5] == '+':           # Tag in plus strand
            b = int(line[1]) + shift
            e = b + extension
            for k in xrange(max(b, 1), min(e, chrlength[spename][line[0]] + 1)): # the step is 1
                values[line[0]][k - 1] += 1
        elif line[5] == '-':        # Tag in minus strand
            e = int(line[2]) - shift
            b = e - extension
            for k in xrange(max(b, 1), min(e, chrlength[spename][line[0]] + 1)): # the step of wig file is 1
                values[line[0]][k - 1] += 1
        else:
            continue
    return values


def tag_to_wig(hash_values, path_name, spename, step = 1):
    namelist = chrlength[spename].keys()
    namelist.sort()
    for name in namelist:
        fileout = open(name +'_' + path_name + '.wig', 'w')
        print >>fileout, "track type=wiggle_0 name=" + name
        print >>fileout, "fixedStep  chrom=" + name + "  start=1  step=" + str(step)
        for k in xrange(0, len(hash_values[name]), step):
            print >>fileout, hash_values[name][k]
        fileout.close()


def nucleosome_tss_profile(hash_values, genefile, spename, upstream = 1000, downstream = 2000):
    f_gene = open(genefile,'r')
    total_tss = 0
    length = upstream + downstream + 1
    all_gene = [0]*length
    
    for line in f_gene:
        line_gene = line.strip().split()
        if line_gene[0] not in chrlength[spename].keys():
           continue

        if line_gene[5] == '+':
            begin = int(line_gene[1])
            for k in range(length):
                if (begin - upstream + k) < 0:
                    pass
                elif (begin - upstream + k) > len(values[line_gene[0]]) - 1:
                    pass
                else:
                    all_gene[k] += values[line_gene[0]][begin-upstream +k]
            total_tss += 1
        elif line_gene[5] == '-':
            begin = int(line_gene[2])
            for k in range(length):
                if (begin + upstream - k) < 0:
                    pass
                elif (begin + upstream - k) > len(values[line_gene[0]]) -1:
                    pass
                else:
                    all_gene[k] += values[line_gene[0]][begin + upstream - k]
            total_tss += 1
        else:
            continue
    meta_gene = [float(k)/total_tss for k in all_gene]
    f_gene.close()
    
    return (range(-1*upstream, downstream + 1), meta_gene)


def nucleosome_tts_profile(hash_values, genefile, spename, upstream = 2000, downstream = 1000):
    f_gene = open(genefile,'r')
    total_tts = 0
    length = upstream + downstream + 1
    all_gene = [0]*length
    
    for line in f_gene:
        line_gene = line.strip().split()
        if line_gene[0] not in chrlength[spename].keys():
            continue

        if line_gene[5] == '+':
            begin = int(line_gene[2])
            for k in range(length):
                if (begin - upstream + k) < 0:
                    pass
                elif (begin - upstream + k) > len(values[line_gene[0]]) - 1:
                    pass
                else:
                    all_gene[k] += values[line_gene[0]][begin-upstream +k]
            total_tts += 1
        elif line_gene[5] == '-':
            begin = int(line_gene[1])
            for k in range(length):
                if (begin + upstream - k) < 0:
                    pass
                elif (begin + upstream - k) > len(values[line_gene[0]]) -1:
                    pass
                else:
                    all_gene[k] += values[line_gene[0]][begin + upstream - k]
            total_tts += 1
        else:
            continue
    meta_gene = [float(k)/total_tts for k in all_gene]
    f_gene.seek(0)
    
    return (range(-1*upstream, downstream + 1), meta_gene)


def nucleosome_tss_profile_each_gene(hash_values, genefile, ofile, spename, upstream = 600, downstream = 1000):
    f_gene = open(genefile,'r')
    total_tss = 0
    length = upstream + downstream + 1
    each_gene = [0]*length
    
    ofile = open(ofile, 'w')
    for line in f_gene:
        line_gene = line.strip().split()
        if line_gene[0] not in chrlength[spename].keys():
            continue
        if line_gene[5] == '+':
            begin = int(line_gene[1])
            for k in range(length):
                if (begin - upstream + k) < 0:
                    pass
                elif (begin - upstream + k) > len(values[line_gene[0]]) - 1:
                    pass
                else:
                    each_gene[k] = values[line_gene[0]][begin-upstream +k]
            ofile.write(line_gene[3] + '\t')
            for l in xrange(0, 1601, 10):
                if l < 1600:
                    ofile.write(str(each_gene[l]) + '\t')
                elif l == 1600:
                    ofile.write(str(each_gene[l]) + '\n')
                else:
                    pass
                each_gene = [0]*length

        elif line_gene[5] == '-':
            begin = int(line_gene[2])
            for k in range(length):
                if (begin + upstream - k) < 0:
                    pass
                elif (begin + upstream - k) > len(values[line_gene[0]]) -1:
                    pass
                else:
                    each_gene[k] += values[line_gene[0]][begin + upstream - k]
            ofile.write(line_gene[3] + '\t')
            for l in xrange(0, 1601, 10):
                if l < 1600:
                    ofile.write(str(each_gene[l]) + '\t')
                elif l == 1600:
                    ofile.write(str(each_gene[l]) + '\n')
                else:
                    pass
            each_gene = [0]*length

        else:
            continue
    f_gene.close()
    ofile.close()


def nucleosome_tts_profile_each_gene(hash_values, genefile, ofile, spename, upstream = 1000, downstream = 600):
    f_gene = open(genefile,'r')
    total_tts = 0
    length = upstream + downstream + 1
    each_gene = [0]*length

    ofile = open(ofile, 'w')
    for line in f_gene:
        line_gene = line.strip().split()
        if line_gene[0] not in chrlength[spename].keys():
            continue

        if line_gene[5] == '+':
            begin = int(line_gene[2])
            for k in range(length):
                if (begin - upstream + k) < 0:
                    pass
                elif (begin - upstream + k) > len(values[line_gene[0]]) - 1:
                    pass
                else:
                    each_gene[k] = values[line_gene[0]][begin-upstream +k]
            ofile.write(line_gene[3] + '\t')
            for l in xrange(0, 1601, 10):
                if l < 1600:
                    ofile.write(str(each_gene[l]) + '\t')
                elif l == 1600:
                    ofile.write(str(each_gene[l]) + '\n')
                else:
                    pass
            each_gene = [0]*length

        elif line_gene[5] == '-':
            begin = int(line_gene[1])
            for k in range(length):
                if (begin + upstream - k) < 0:
                    pass
                elif (begin + upstream - k) > len(values[line_gene[0]]) -1:
                    pass
                else:
                    each_gene[k] += values[line_gene[0]][begin + upstream - k]
            ofile.write(line_gene[3] + '\t')
            for l in xrange(0, 1601, 10):
                if l < 1600:
                    ofile.write(str(each_gene[l]) + '\t')
                elif l == 1600:
                    ofile.write(str(each_gene[l]) + '\n')
                else:
                    pass
            each_gene = [0]*length

        else:
            continue
    f_gene.close()
    ofile.close()


if __name__ == '__main__':
    import sys
    from optparse import OptionParser
    
    usage = "usage: python nucleosome_profile.py -i input_bed -w output_wig --tss output_tss --tts output_tts -g gene"
    parser = OptionParser(usage)
    
    parser.add_option("-i",dest="input",help="input Bed file")
    parser.add_option("-w",dest="wig",help="output wiggle directory")
    parser.add_option("-e",dest="extension",default=73,type="int",help="extension")
    parser.add_option("-t",dest="shift",default=36,type="int",help="shift")
    parser.add_option("--tss",dest="tss",help="output tss profile")
    parser.add_option("--tts",dest="tts",help="output tts profile")
    parser.add_option("--tsseach",dest="tsseach",help="output tss each gene profile")
    parser.add_option("--ttseach",dest="ttseach",help="output tts each gene profile")
    parser.add_option("-s", "--spename", dest="spename",default="zv9", type="string", help="species name" )
    parser.add_option("-g",dest="gene",default='/mnt/Storage/home/fuk/Zebrafish/gene/zebrafish_refseq', help="refseq gene list")
    parser.add_option("--name",dest="name",default='NA',help="Run name")
    (options,args) = parser.parse_args()
    
    if not options.input or not options.gene:
        parser.print_help()
        sys.exit()
    
    values = tag2hash(options.input, spename=options.spename, extension=int(options.extension), shift=int(options.shift))
    #if options.wig[-1] == '/':
    #    wigdir = options.wig[:-1]
    #else:
    #    wigdir = options.wig
    #tag_to_wig(values, wigdir, step = 1)
    
    #if not options.tss or not options.tts or not options.tsseach or not options.ttseach:
    #    sys.exit()    
    
    #tag_to_wig(values, options.name, options.spename, step = 10)
    
    if options.tss:
        (coordinate_tss, meta_gene_tss) = nucleosome_tss_profile(values,options.gene, options.spename)
        f_metagene_tss = open(options.tss,'w')
        for i in range(0, len(coordinate_tss), 10):
            f_metagene_tss.write(str(coordinate_tss[i]) + '\t' + str(meta_gene_tss[i]) + '\n')
        f_metagene_tss.close()
    if options.tts:
        (coordinate_tts, meta_gene_tts) = nucleosome_tts_profile(values,options.gene, options.spename)
        f_metagene_tts = open(options.tts,'w')
        for i in range(0, len(coordinate_tts), 10):
            f_metagene_tts.write(str(coordinate_tts[i]) + '\t' + str(meta_gene_tts[i]) + '\n')
        f_metagene_tts.close()
    
    if options.tsseach:
        nucleosome_tss_profile_each_gene(values,options.gene,options.tsseach,options.spename)
    if options.ttseach:
        nucleosome_tts_profile_each_gene(values,options.gene,options.ttseach,options.spename)


