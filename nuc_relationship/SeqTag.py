#!/usr/bin/python

# Yong Zhang @ DFCI

from Genome import chrlength
import sys
import math
#from rpy import r

class SeqTag(object):
    """ The class SeqTag(bedfile, spename, extension = 73, shift = 36)
        tag_to_wig: from bed format file to wig file
        normalize_subtraction:  two wig files normalization and subtraction
        wig_genomic_feature: plot average tag density around certain genomic regions
    """
    def __init__(self, bedfile, spename, extension = 73, shift = 36, strand = 'both'):
        """ strand: 'both', '+' or '-' """
        self.spename = spename
        self.values = {}
        for chrname in chrlength[spename].keys():
            self.values[chrname] = [0] * chrlength[spename][chrname]
            
        for line in open(bedfile).xreadlines():
            line = line.strip().split()
            if line[0] not in chrlength[spename].keys():
                continue
            try:
                if line[5] == '+' and strand != '-':           # Tag in plus strand
                    b = int(line[1]) + shift
                    e = b + extension
                    for k in xrange(max(b, 0), min(e, chrlength[spename][line[0]])):
                        self.values[line[0]][k - 1] += 1
                            
                elif line[5] == '-' and strand != '+':        # Tag in minus strand
                    e = int(line[2]) - shift
                    b = e - extension
                    for k in xrange(max(b, 0), min(e, chrlength[spename][line[0]])):
                        self.values[line[0]][k - 1] += 1
                else:
                    continue
            except:
                print >> sys.stderr, 'Tag position file error: ', sys.exc_info()[0], sys.exc_info()[1]
                sys.exit()
                
    def tag_to_wig(self, wigfile):
        fileout = open(wigfile, 'w')
        namelist = chrlength[self.spename].keys()
        namelist.sort()
        for name in namelist:
            print >>fileout, "track type=wiggle_0 name=" + name
            print >>fileout, "fixedStep  chrom=" + name + "  start=0  step=1"
            for k in xrange(len(self.values[name])):
                print >>fileout, self.values[name][k]
        fileout.close()
    
    def return_hash(self):
        return self.values
    
    def normalize_subtraction(self, wigfile, othervalue, ratio1, ratio2):
        fileout = open(wigfile, 'w')
        namelist = chrlength[self.spename].keys()
        namelist.sort()
        for name in namelist:
            print >>fileout, "track type=wiggle_0 name=" + name
            print >>fileout, "fixedStep  chrom=" + name + "  start=0  step=1"
            for k in xrange(len(self.values[name])):
                print >>fileout, (self.values[name][k] * ratio1) - othervalue[name][k] * ratio2
        fileout.close()
    
    def enriched_loci(self, cutoff, outputfile, strand = '+'):
        fileout = open(outputfile, 'w')
        namelist = chrlength[self.spename].keys()
        namelist.sort()
        for name in namelist:
            for k in xrange(len(self.values[name])):
                if self.values[name][k] >= cutoff:
                    print >>fileout, name + '\t' + str(k + 1) + '\t' + strand + '\t' + str(self.values[name][k])
        fileout.close()
        
        clusters = []
        chrold = ''
        posold = 0
        for line in open(outputfile).xreadlines():
            line = line.strip().split()
            if line[0] != chrold:
                chrold = line[0]
                clusters.append([int(line[1])])
                posold = int(line[1])
                continue
            if int(line[1]) - posold > 50:
                clusters.append([int(line[1])])
                posold = int(line[1])
                continue
            clusters[-1].append(int(line[1]))
            posold = int(line[1])
        
        value = [0] * 150
        number = 0
        if strand == '+':
            for k in clusters:
                if len(k) > 2:
                    number += 1
                    for i in xrange(1, len(k)):
                        if (k[i] - k[0]) < 150:
                            value[k[i] - k[0]] += 1
        else:
            for k in clusters:
                if len(k) > 2:
                    number += 1
                    for i in xrange(len(k) - 1):
                        if (k[-1] - k[i]) < 150:
                            value[k[-1] - k[i]] += 1
        for k in xrange(len(value)):
            value[k] = str(value[k] / float(number))
        
        return value
        
    def enriched_loci_center_pattern(self, cutoff, outputfile, strand = '+'):
        fileout = open(outputfile, 'w')
        namelist = chrlength[self.spename].keys()
        namelist.sort()
        for name in namelist:
            for k in xrange(len(self.values[name])):
                if self.values[name][k] >= cutoff:
                    print >>fileout, name + '\t' + str(k + 1) + '\t' + strand + '\t' + str(self.values[name][k])
        fileout.close()
        
        clusters = []
        chrold = ''
        posold = 0
        for line in open(outputfile).xreadlines():
            line = line.strip().split()
            if line[0] != chrold:
                chrold = line[0]
                clusters.append([[int(line[1]), int(line[3])]])
                posold = int(line[1])
                continue
            if int(line[1]) - posold > 50:
                clusters.append([[int(line[1]), int(line[3])]])
                posold = int(line[1])
                continue
            clusters[-1].append([int(line[1]), int(line[3])])
            posold = int(line[1])
        
        value = [0] * 75
        number = 0
        for k in clusters:
            if len(k) > 3:# and (k[-1][0] - k[0][0])>= 10:
                number += 1
                pos = []
                h = 0
                p = 0
                for i in range(len(k)):
                    pos.append(k[i][0])
                    if k[i][1] >= k[h][1]:
                        h = i
                        p = k[h][0] 
                for i in xrange(len(pos)):
                    if abs(pos[i] - p) == 0:
                        continue
                    if abs(pos[i] - p) <= 75:
                            value[abs(pos[i] - p) - 1] += 1
                            
        for k in xrange(len(value)):
            value[k] = str(value[k] / float(number))
            
        print outputfile, number
        return value
        
    def wig_genomic_feature(self, genomic_bedfile, width = 1501, strand = 'both', remove = 0.0, plot = 'F'):
        """ plot average tag density around certain genomic regions 
        strand: 'both' or '+'
        """
        halfwidth = width / 2
        length = halfwidth * 2 + 1
        wigs_5end = [0] * length
        wigs_3end = [0] * length
                
        genelist_5end = {}
        genelist_3end = {}
        for lineold in open(genomic_bedfile).xreadlines():
            line = lineold.strip().split()
            if not self.values.has_key(line[0]):
                continue
            if len(line) < 6 or line[5] == '+':
                begin = int(line[1])
                end = int(line[2])
                genelist_5end[lineold.strip()] = 0
                genelist_3end[lineold.strip()] = 0
                
                for k in range(length):
                    if (begin - halfwidth + k) < 0:
                        pass
                    elif (begin - halfwidth + k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        genelist_5end[lineold.strip()] += self.values[line[0]][begin - halfwidth + k]
                        
                    if (end - halfwidth + k) < 0:
                        pass
                    elif (end - halfwidth + k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        genelist_3end[lineold.strip()] += self.values[line[0]][end - halfwidth + k]
                        
            elif strand == 'both':
                begin = int(line[2])
                end = int(line[1])
                genelist_5end[lineold.strip()] = 0
                genelist_3end[lineold.strip()] = 0
                
                for k in range(length):
                    if (begin + halfwidth - k) < 0:
                        pass
                    elif (begin + halfwidth - k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        genelist_5end[lineold.strip()] += self.values[line[0]][begin + halfwidth - k]
                        
                    if (end + halfwidth - k) < 0:
                        pass
                    elif (end + halfwidth - k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        genelist_3end[lineold.strip()] += self.values[line[0]][end + halfwidth - k]
                        
        if remove <= 0.0:
            genelist5 = genelist_5end.keys()
            genelist3 = genelist_3end.keys()
        elif remove < 1.0:
            a = [ (value, key) for key, value in genelist_5end.iteritems() ]
            a.sort()
            total = len(genelist_5end.keys())
            genelist5 = []
            for i in range(int(total * remove / 2), total - int(total * remove / 2)):
                genelist5.append(a[i][1])
                
            a = [ (value, key) for key, value in genelist_3end.iteritems() ]
            a.sort()
            total = len(genelist_3end.keys())
            genelist3 = []
            for i in range(int(total * remove / 2), total - int(total * remove / 2)):
                genelist3.append(a[i][1])
        else:
            sys.exit('Wrong remove value!')
        
        # 5 end
        for line in genelist5:
            line = line.strip().split()
            if len(line) < 6 or line[5] == '+':
                begin = int(line[1])
                end = int(line[2])
                for k in range(length):
                    if (begin - halfwidth + k) < 0:
                        pass
                    elif (begin - halfwidth + k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        wigs_5end[k] += self.values[line[0]][begin - halfwidth + k]
                        
            elif strand == 'both':
                begin = int(line[2])
                end = int(line[1])
                for k in range(length):
                    if (begin + halfwidth - k) < 0:
                        pass
                    elif (begin + halfwidth - k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        wigs_5end[k] += self.values[line[0]][begin + halfwidth - k]
        
        # 3 end
        for line in genelist3:
            line = line.strip().split()
            if len(line) < 6 or line[5] == '+':
                begin = int(line[1])
                end = int(line[2])
                for k in range(length):
                    if (end - halfwidth + k) < 0:
                        pass
                    elif (end - halfwidth + k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        wigs_3end[k] += self.values[line[0]][end - halfwidth + k]
                        
            elif strand == 'both':
                begin = int(line[2])
                end = int(line[1])
                for k in range(length):
                    if (end + halfwidth - k) < 0:
                        pass
                    elif (end + halfwidth - k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        wigs_3end[k] += self.values[line[0]][end + halfwidth - k]
  
        for k in range(length):
            wigs_5end[k] = wigs_5end[k] / float(len(genelist5))
            wigs_3end[k] = wigs_3end[k] / float(len(genelist3))
        
        #fileout = open(wigfile, 'w')
        #print >> fileout, '# 5\' end centered profile'
        #for k in range(length):
        #    print >> fileout, str(k - halfwidth) + '\t' + str(wigs_5end[k])
        #print >> fileout, '# 3\' end centered profile'
        #for k in range(length):
        #    print >> fileout, str(k - halfwidth) + '\t' + str(wigs_3end[k])
        #fileout.close()
        return (wigs_5end, wigs_3end)
        #if plot == 'T':
        #     r.pdf(wigfile + '.pdf', width = 5, height = 8.5)
        #     r.par(mfrow=(2, 1))
        #     r.plot(range(-1 * halfwidth, halfwidth + 1), wigs_5end, col = 'red', type = 'l', main = '5\' End Tag Profile', xlab = 'Locations relative to 5\' end (bp)', ylab = 'average tag density')
        #     r.plot(range(-1 * halfwidth, halfwidth + 1), wigs_3end, col = 'red', type = 'l', main = '3\' End Tag Profile', xlab = 'Locations relative to 3\' end (bp)', ylab = 'average tag density')
        #     r.dev_off()

    def wig_genomic_feature_normalize(self, genomic_bedfile, width = 1501, referencetag = 3.0, strand = 'both', remove = 0.0, plot = 'F'):
        """ plot average tag density around certain genomic regions 
        strand: 'both' or '+'
        """
        halfwidth = width / 2
        length = halfwidth * 2 + 1
        wigs_5end = [0] * length
        wigs_3end = [0] * length
                
        genelist_5end = {}
        genelist_3end = {}
        for lineold in open(genomic_bedfile).xreadlines():
            line = lineold.strip().split()
            if not self.values.has_key(line[0]):
                continue
            if len(line) < 6 or line[5] == '+':
                begin = int(line[1])
                end = int(line[2])
                genelist_5end[lineold.strip()] = 0
                genelist_3end[lineold.strip()] = 0
                
                for k in range(length):
                    if (begin - halfwidth + k) < 0:
                        pass
                    elif (begin - halfwidth + k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        genelist_5end[lineold.strip()] += self.values[line[0]][begin - halfwidth + k]
                        
                    if (end - halfwidth + k) < 0:
                        pass
                    elif (end - halfwidth + k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        genelist_3end[lineold.strip()] += self.values[line[0]][end - halfwidth + k]
                        
            elif strand == 'both':
                begin = int(line[2])
                end = int(line[1])
                genelist_5end[lineold.strip()] = 0
                genelist_3end[lineold.strip()] = 0
                
                for k in range(length):
                    if (begin + halfwidth - k) < 0:
                        pass
                    elif (begin + halfwidth - k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        genelist_5end[lineold.strip()] += self.values[line[0]][begin + halfwidth - k]
                        
                    if (end + halfwidth - k) < 0:
                        pass
                    elif (end + halfwidth - k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        genelist_3end[lineold.strip()] += self.values[line[0]][end + halfwidth - k]
                        
        if remove <= 0.0:
            genelist5 = genelist_5end.keys()
            genelist3 = genelist_3end.keys()
        elif remove < 1.0:
            a = [ (value, key) for key, value in genelist_5end.iteritems() ]
            a.sort()
            total = len(genelist_5end.keys())
            genelist5 = []
            for i in range(int(total * remove / 2), total - int(total * remove / 2)):
                genelist5.append(a[i][1])
                
            a = [ (value, key) for key, value in genelist_3end.iteritems() ]
            a.sort()
            total = len(genelist_3end.keys())
            genelist3 = []
            for i in range(int(total * remove / 2), total - int(total * remove / 2)):
                genelist3.append(a[i][1])
        else:
            sys.exit('Wrong remove value!')
        
        # 5 end
        for lineold in genelist5:
            line = lineold.strip().split()
            if len(line) < 6 or line[5] == '+':
                begin = int(line[1])
                end = int(line[2])
                ref = max(referencetag, genelist_5end[lineold] / float(width))
                for k in range(length):
                    if (begin - halfwidth + k) < 0:
                        pass
                    elif (begin - halfwidth + k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        wigs_5end[k] += ((self.values[line[0]][begin - halfwidth + k]) / float(ref))
                        
            elif strand == 'both':
                begin = int(line[2])
                end = int(line[1])
                ref = max(referencetag, genelist_5end[lineold] / float(width))
                for k in range(length):
                    if (begin + halfwidth - k) < 0:
                        pass
                    elif (begin + halfwidth - k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        wigs_5end[k] += ((self.values[line[0]][begin + halfwidth - k]) / float(ref))
        
        # 3 end
        for line in genelist3:
            line = line.strip().split()
            if len(line) < 6 or line[5] == '+':
                begin = int(line[1])
                end = int(line[2])
                ref = max(referencetag, genelist_3end[lineold] / float(width))
                for k in range(length):
                    if (end - halfwidth + k) < 0:
                        pass
                    elif (end - halfwidth + k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        wigs_3end[k] += ((self.values[line[0]][end - halfwidth + k]) / float(ref))
                        
            elif strand == 'both':
                begin = int(line[2])
                end = int(line[1])
                ref = max(referencetag, genelist_3end[lineold] / float(width))
                for k in range(length):
                    if (end + halfwidth - k) < 0:
                        pass
                    elif (end + halfwidth - k) > len(self.values[line[0]]) - 1:
                        pass
                    else:
                        wigs_3end[k] += ((self.values[line[0]][end + halfwidth - k]) / float(ref))
  
        for k in range(length):
            wigs_5end[k] = wigs_5end[k] / float(len(genelist5))
            wigs_3end[k] = wigs_3end[k] / float(len(genelist3))
        
        #if plot == 'T':
        #     r.pdf(wigfile + '.pdf', width = 5, height = 8.5)
        #     r.par(mfrow=(2, 1))
        #     r.plot(range(-1 * halfwidth, halfwidth + 1), wigs_5end, col = 'red', type = 'l', main = '5\' End Tag Profile', xlab = 'Locations relative to 5\' end (bp)', ylab = 'average tag density')
        #     r.plot(range(-1 * halfwidth, halfwidth + 1), wigs_3end, col = 'red', type = 'l', main = '3\' End Tag Profile', xlab = 'Locations relative to 3\' end (bp)', ylab = 'average tag density')
        #     r.dev_off()
        return (wigs_5end, wigs_3end)

    def gene_promoter_profile(self, bedfile, outfile):
        total = 0
        len = 0
        for chrname in chrlength[self.spename].keys():
            len += chrlength[self.spename][chrname]
            total += sum(self.values[chrname])
        
        average = float(total) / len
        
        f1 = open(outfile + '_TSS.xls', 'w')
        names = ['GeneID']
        for i in range(-450, 251, 10):
            names.append('TSS_' + str(i))
        print >>f1, '\t'.join(names)
        
        f2 = open(outfile + '_TTS.xls', 'w')
        names = ['GeneID']
        for i in range(-250, 451, 10):
            names.append('TTS_' + str(i))
        print >>f2, '\t'.join(names)
        
        # TSS & TTS
        for line in open(bedfile):
            if line[:3] != 'chr':
                continue
            line = line.strip().split()
            if line[5] == '+':
                TSSbegin = int(line[1]) - 450
                TSSend = int(line[1]) + 251
                TTSbegin = int(line[2]) - 250
                TTSend = int(line[2]) + 451             
                step = 10
            else:
                TSSbegin = int(line[2]) + 450
                TSSend = int(line[2]) - 251
                TTSbegin = int(line[1]) + 250
                TTSend = int(line[1]) - 451              
                step = -10             
                
            if min(TSSbegin, TSSend, TTSbegin, TTSend) < 0 or max(TSSbegin, TSSend, TTSbegin, TTSend) >= chrlength[self.spename][line[0]]:
                continue
            TSSave = float(sum(self.values[line[0]][min(TSSbegin, TSSend): max(TSSbegin, TSSend)])) / 700
            TTSave = float(sum(self.values[line[0]][min(TTSbegin, TTSend): max(TTSbegin, TTSend)])) / 700
            
            TSSvalue = [line[3]]
            for i in range(TSSbegin, TSSend, step):
                if self.values[line[0]][i] >=  max(average, TSSave):
                    TSSvalue.append(str(math.log(min(5, float(self.values[line[0]][i]) / max(average, TSSave)))/math.log(2.0)))
                else:
                    TSSvalue.append(str(math.log(max(0.2, float(self.values[line[0]][i]) / max(average, TSSave)))/math.log(2.0)))
            print >>f1, '\t'.join(TSSvalue)
            
            TTSvalue = [line[3]]
            for i in range(TTSbegin, TTSend, step):
                if self.values[line[0]][i] >=  max(average, TTSave):
                    TTSvalue.append(str(math.log(min(5, float(self.values[line[0]][i]) / max(average, TTSave)))/math.log(2.0)))
                else:
                    TTSvalue.append(str(math.log(max(0.2, float(self.values[line[0]][i]) / max(average, TTSave)))/math.log(2.0)))                
            print >>f2, '\t'.join(TTSvalue)
        f1.close()
        f2.close()

    def gene_promoter_profile_2(self, bedfile, outfile):
        total = 0
        len = 0
        for chrname in chrlength[self.spename].keys():
            len += chrlength[self.spename][chrname]
            total += sum(self.values[chrname])
        
        average = float(total) / len
        
        f1 = open(outfile + '_TSS.xls', 'w')
        names = ['GeneID']
        for i in range(-450, 251, 10):
            names.append('TSS_' + str(i))
        print >>f1, '\t'.join(names)
        
        f2 = open(outfile + '_TTS.xls', 'w')
        names = ['GeneID']
        for i in range(-250, 451, 10):
            names.append('TTS_' + str(i))
        print >>f2, '\t'.join(names)
        
        # TSS & TTS
        for line in open(bedfile):
            if line[:3] != 'chr':
                continue
            line = line.strip().split()
            if line[5] == '+':
                TSSbegin = int(line[1]) - 450
                TSSend = int(line[1]) + 251
                TTSbegin = int(line[2]) - 250
                TTSend = int(line[2]) + 451             
                step = 10
            else:
                TSSbegin = int(line[2]) + 450
                TSSend = int(line[2]) - 251
                TTSbegin = int(line[1]) + 250
                TTSend = int(line[1]) - 451              
                step = -10             
                
            if min(TSSbegin, TSSend, TTSbegin, TTSend) < 0 or max(TSSbegin, TSSend, TTSbegin, TTSend) >= chrlength[self.spename][line[0]]:
                continue
            TSSave = float(sum(self.values[line[0]][min(TSSbegin, TSSend): max(TSSbegin, TSSend)])) / 700
            TTSave = float(sum(self.values[line[0]][min(TTSbegin, TTSend): max(TTSbegin, TTSend)])) / 700
            
            TSSvalue = [line[3]]
            for i in range(TSSbegin, TSSend, step):
                TSSvalue.append(str(float(self.values[line[0]][i]) / max(average, TSSave)))
            print >>f1, '\t'.join(TSSvalue)
            
            TTSvalue = [line[3]]
            for i in range(TTSbegin, TTSend, step):
                TTSvalue.append(str(float(self.values[line[0]][i]) / max(average, TTSave)))               
            print >>f2, '\t'.join(TTSvalue)
        f1.close()
        f2.close()

def test_tag_to_wig():
    # E. coli
    names = ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']
    for i in names:
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/UCSC_version/BEDfile/' + i + '_Ecoli.bed', 'ecoli', extension = 10, shift = 68)
        tag.tag_to_wig('/home/zy/collaboration/Nuc_Struhl/Data/UCSC_version/WIGfile_10bp/' + i + '_Ecoli.wig')
    
    # Yease UCSC format
    for i in names:
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/UCSC_version/BEDfile/' + i + '_yeast.bed', 'yeast_UCSC', extension = 10, shift = 68)
        tag.tag_to_wig('/home/zy/collaboration/Nuc_Struhl/Data/UCSC_version/WIGfile_10bp/' + i + '_Yeast.wig')
        
    # Yease SGD format
    for i in names:
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/SGD_version/BEDfile/' + i + '_yeast.bed', 'yeast_SGD', extension = 10, shift = 68)
        tag.tag_to_wig('/home/zy/collaboration/Nuc_Struhl/Data/SGD_version/WIGfile_10bp/' + i + '_Yeast.wig')

def wig_genomic_feature_four_sample_plot(genename, extension = 73, shift = 36, strand = 'both', end = 5, version = 'SGD', width = 1501, multiple = 1, remove = 0.0):
    samples = ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']
    reftags = [0.3 * extension, 0.27 * extension, 0.1 * extension, 0.11 * extension]
    #ratio = [0.93, 1.04, 2.83, 2.55]
    
    mylist = {}
    
    for i in range(len(samples)):
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/' + version + '_version/BEDfile/' + samples[i] + '_yeast.bed', 'yeast_' + version, extension = extension, shift = shift, strand = strand)
        if strand == 'both':
            (end5, end3) = tag.wig_genomic_feature_normalize('/home/zy/collaboration/Nuc_Struhl/Data/Genome/Annotation/TSS_cluster/SGD_gene_' + version + genename, width = width, referencetag = reftags[i], remove = remove)
        else:
            (end5, end3) = tag.wig_genomic_feature_normalize('/home/zy/collaboration/Nuc_Struhl/Data/Genome/Annotation/TSS_cluster/SGD_gene_' + version + genename, width = width, referencetag = reftags[i], strand = '+', remove = remove)
        mylist[samples[i]] = []
        if end == 5:
            for k in range(len(end5)):
                mylist[samples[i]].append(str(end5[k]))
                #mylist[samples[i]].append(str(multiple * end5[k] * ratio[i] / extension))
        else:
            for k in range(len(end3)):
                mylist[samples[i]].append(str(end3[k]))
                #mylist[samples[i]].append(str(multiple * end3[k] * ratio[i] / extension))
    return mylist

def wig_genomic_feature_plot(Rscript, version = 'SGD', remove = 0.0, cluster = 0):
    fo = open(Rscript, 'a')
    #values = {'_5_end_remove_MNase_biased.bed' : 'All 4,918 Genes 5 end (Removed MNase biased regions)', '_5_end_convergent_remove_MNase_biased.bed' : '1,696 5 End Convergent Genes (Removed MNase biased regions)', '_3_end_remove_MNase_biased.bed' : 'All 5,200 Genes 3 end (Removed MNase biased regions)', '_3_end_convergent_remove_MNase_biased.bed' : '1,036 3 End Convergent Genes (Removed MNase biased regions)'}
    values = {'_C' + str(cluster) + '.bed' : 'All Genes in Cluster ' + str(cluster), '_C' + str(cluster) + '_5_end_convergent.bed' : '5 End Convergent Genes in Cluster ' + str(cluster)}
    
    x = []
    begin = -750
    for k in range(1501):
        x.append(str(begin))
        begin += 1
    print >>fo, 'x <- c(' + ', '.join(x) + ')\n'
    
    # 5' end 10 bp extension
    print >>fo, "pdf('Yeast_5_end_10bp_" + version + ".pdf', width = 10, height = 4.25)"
    print >>fo, "par(mfrow=c(1, 2))"
    
    for gene in values.keys(): #['_5_end_remove_MNase_biased.bed', '_5_end_convergent_remove_MNase_biased.bed']:
        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 10, shift = 68, strand = 'both', end = 5, version = version, width = 1501, multiple = 1, remove = remove)
        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0, 1), main = '5 End Tag Profile (" + values[gene] + ")', xlab = 'Locations relative to 5 end (bp)', ylab = 'average tag profile')"
        k = 4
        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
            k += 1
        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
        
#    for gene in ['.bed', '_5end_convergent.bed']:
#        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 10, shift = -5, strand = '+', end = 5, version = version, width = 1501, multiple = 1, remove = remove)
#        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
#        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0.1, 0.25), main = '5 End Plus Tag Profile (" + values[gene] + ", plus strand genes)', xlab = 'Locations relative to 5 end (bp)', ylab = 'average plus tag density')"
#        k = 4
#        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
#            k += 1
#        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
        
#    for gene in ['.bed', '_5end_convergent.bed']:
#        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 10, shift = -5, strand = '-', end = 5, version = version, width = 1501, multiple = 1, remove = remove)
#        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
#        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0.1, 0.25), main = '5 End Minus Tag Profile (" + values[gene] + ", plus strand genes)', xlab = 'Locations relative to 5 end (bp)', ylab = 'average minus tag density')"
#        k = 4
#        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
#            k += 1
#        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
    
    print >>fo, 'dev.off()'
    
    # 5' end 73 bp extension
    print >>fo, "pdf('Yeast_5_end_73bp_" + version + ".pdf', width = 10, height = 4.25)"
    print >>fo, "par(mfrow=c(1, 2))"
    
    for gene in values.keys(): #['_5_end_remove_MNase_biased.bed', '_5_end_convergent_remove_MNase_biased.bed']:
        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 73, shift = 36, strand = 'both', end = 5, version = version, width = 1501, multiple = 1, remove = remove)
        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0, 1), main = '5 End Tag Profile (" + values[gene] + ")', xlab = 'Locations relative to 5 end (bp)', ylab = 'average tag profile')"
        k = 4
        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
            k += 1
        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
        
#    for gene in ['.bed', '_5end_convergent.bed']:
#        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 73, shift = -36, strand = '+', end = 5, version = version, width = 1501, multiple = 1, remove = remove)
#        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
#        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0.1, 0.25), main = '5 End Plus Tag Profile (" + values[gene] + ", plus strand genes)', xlab = 'Locations relative to 5 end (bp)', ylab = 'average plus tag density')"
#        k = 4
#        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
#            k += 1
#        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
        
#    for gene in ['.bed', '_5end_convergent.bed']:
#        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 73, shift = -36, strand = '-', end = 5, version = version, width = 1501, multiple = 1, remove = remove)
#        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
#        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0.1, 0.25), main = '5 End Minus Tag Profile (" + values[gene] + ", plus strand genes)', xlab = 'Locations relative to 5 end (bp)', ylab = 'average minus tag density')"
#        k = 4
#        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
#            k += 1
#        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
    
    print >>fo, 'dev.off()'

    """     
    # 3' end 10 bp extension
    print >>fo, "pdf('Yeast_3_end_10bp_" + version + ".pdf', width = 10, height = 4.25)"
    print >>fo, "par(mfrow=c(1, 2))"

    for gene in ['_3_end_remove_MNase_biased.bed', '_3_end_convergent_remove_MNase_biased.bed']:
        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 10, shift = 68, strand = 'both', end = 3, version = version, width = 1501, multiple = 1, remove = remove)
        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0, 1), main = '3 End Tag Profile (" + values[gene] + ")', xlab = 'Locations relative to 3 end (bp)', ylab = 'average tag profile')"
        k = 4
        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
            k += 1
        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
        
#    for gene in ['.bed', '_3end_convergent.bed']:
#        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 10, shift = -5, strand = '+', end = 3, version = version, width = 1501, multiple = 1, remove = remove)
#        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
#        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0.1, 0.25), main = '3 End Plus Tag Profile (" + values[gene] + ", plus strand genes)', xlab = 'Locations relative to 3 end (bp)', ylab = 'average plus tag density')"
#        k = 4
#        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
#            k += 1
#        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
        
#    for gene in ['.bed', '_5end_convergent.bed']:
#        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 10, shift = -5, strand = '-', end = 3, version = version, width = 1501, multiple = 1, remove = remove)
#        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
#        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0.1, 0.25), main = '3 End Minus Tag Profile (" + values[gene] + ", plus strand genes)', xlab = 'Locations relative to 3 end (bp)', ylab = 'average minus tag density')"
#        k = 4
#        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
#            k += 1
#        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
    
    print >>fo, 'dev.off()'
    
    # 5' end 73 bp extension
    print >>fo, "pdf('Yeast_3_end_73bp_" + version + ".pdf', width = 10, height = 4.25)"
    print >>fo, "par(mfrow=c(1, 2))"
    
    for gene in ['_3_end_remove_MNase_biased.bed', '_3_end_convergent_remove_MNase_biased.bed']:
        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 73, shift = 36, strand = 'both', end = 3, version = version, width = 1501, multiple = 1, remove = remove)
        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0, 1), main = '3 End Tag Profile (" + values[gene] + ")', xlab = 'Locations relative to 3 end (bp)', ylab = 'average tag profile')"
        k = 4
        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
            k += 1
        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
        
#    for gene in ['.bed', '_3end_convergent.bed']:
#        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 73, shift = -36, strand = '+', end = 3, version = version, width = 1501, multiple = 1, remove = remove)
#        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
#        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0.1, 0.25), main = '3 End Plus Tag Profile (" + values[gene] + ", plus strand genes)', xlab = 'Locations relative to 3 end (bp)', ylab = 'average plus tag density')"
#        k = 4
#        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
#            k += 1
#        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
        
#    for gene in ['.bed', '_3end_convergent.bed']:
#        mylist = wig_genomic_feature_four_sample_plot(gene, extension = 73, shift = -36, strand = '-', end = 3, version = version, width = 1501, multiple = 1, remove = remove)
#        for n in ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
#        print >>fo, "plot(x, ACF1, col = 2, type = 'l', ylim = c(0.1, 0.25), main = '3 End Minus Tag Profile (" + values[gene] + ", plus strand genes)', xlab = 'Locations relative to 3 end (bp)', ylab = 'average minus tag density')"
#        k = 4
#        for sample in ['SALT', 'C_MNase', 'C_SONIC']:
#            print >>fo, "lines(x, " + sample + ", col = " + str(k) + ", type = 'l')"
#            k += 1
#        print >>fo, "legend('topright',  c('" + "', '".join(['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']) + "'), col = c(2,4,5,6), lty = 1)"
    """    
    print >>fo, 'dev.off()'    
    fo.close()
    

def wig_genomic_feature(Rfile):
    names = ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']
    name2 = ['ACF', 'SALT', 'MNase_ctl', 'Sonication_ctl']
    genes = ['', '_5end_convergent', '_5end_shared_others_5end', '_5end_shared_others_3end', '_3end_convergent', '_3end_shared_others_5end', '_3end_shared_others_3end']
    ratio = [0.93, 1.04, 2.83, 2.55]
    
    # Yeast SGD
    fo = open(Rfile, 'a')
    
    print >>fo, "## Yeast SGD"
    
    for i in range(len(names)):
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/SGD_version/BEDfile/' + names[i] + '_yeast.bed', 'yeast_SGD', extension = 10, shift = 68)
        for j in genes:
            end5 = []
            end3 = []
            (end5, end3) = tag.wig_genomic_feature('/home/zy/collaboration/Nuc_Struhl/Data/Genome/Annotation/SGD_gene_SGD' + j + '.bed')
            for k in range(len(end5)):
                end5[k] = str(end5[k] * ratio[i])
                end3[k] = str(end3[k] * ratio[i])
            print >>fo, names[i] + '_end5_' + j + ' <- c(' + ', '.join(end5) + ')'
            print >>fo, names[i] + '_end3_' + j + ' <- c(' + ', '.join(end3) + ')'
    
    x = []
    begin = int(-1 * (len(end5) / 2))
    for k in range(len(end5)):
        x.append(str(begin))
        begin += 1
    
    print >>fo, 'x <- c(' + ', '.join(x) + ')\n' 
    print >>fo, "pdf('Yeast_5_end_SGD.pdf', width = 10, height = 8.5)"
    print >>fo, "par(mfrow=c(2, 2))"
    gene5end = ['', '_5end_convergent', '_5end_shared_others_5end', '_5end_shared_others_3end']
    mains = ['All genes (6,604)', 'Convergent genes (2,258)', 'Sharing other gene 5 end (1,895)', 'Sharing other gene 3 end (2,395)']
    for i in range(len(gene5end)):
        print >>fo, "plot(x, ACF1_end5_" + gene5end[i] +", col = 2, type = 'l', main = '5 End Tag Profile (" + mains[i] + ")', xlab = 'Locations relative to 5 end (bp)', ylab = 'average tag density')"
        for k in range(1, 4):
            print >>fo, "lines(x, " + names[k] + "_end5_" + gene5end[i]  + ", col = " + str(k + 3) + ", type = 'l')"
        print >>fo, "legend('topright',  c('" + "', '".join(name2) + "'), col = c(2,4,5,6), lty = 1)"
    print >>fo, "dev.off()"
    
    print >>fo, "pdf('Yeast_3_end_SGD.pdf', width = 10, height = 8.5)"
    print >>fo, "par(mfrow=c(2, 2))"
    gene5end = ['', '_3end_convergent', '_3end_shared_others_3end', '_3end_shared_others_5end']
    mains = ['All genes (6,604)', 'Convergent genes (1,317)', 'Sharing other gene 3 end (2,970)', 'Sharing other gene 5 end (2,451)']
    for i in range(len(gene5end)):
        print >>fo, "plot(x, ACF1_end3_" + gene5end[i] +", col = 2, type = 'l', main = '3 End Tag Profile (" + mains[i] + ")', xlab = 'Locations relative to 3 end (bp)', ylab = 'average tag density')"
        for k in range(1, 4):
            print >>fo, "lines(x, " + names[k] + "_end3_" + gene5end[i]  + ", col = " + str(k + 3) + ", type = 'l')"
        print >>fo, "legend('topright',  c('" + "', '".join(name2) + "'), col = c(2,4,5,6), lty = 1)"
    print >>fo, "dev.off()"
    
    sys.exit()
    
    # Yeast UCSC
    
    print >>fo, "## Yeast UCSC"
    for i in range(len(names)):
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/UCSC_version/BEDfile/' + names[i] + '_yeast.bed', 'yeast_UCSC', extension = 10, shift = 68)
        for j in genes:
            end5 = []
            end3 = []
            (end5, end3) = tag.wig_genomic_feature('/home/zy/collaboration/Nuc_Struhl/Data/Genome/Annotation/SGD_gene_UCSC' + j + '.bed')
            for k in range(len(end5)):
                end5[k] = str(end5[k] * ratio[i])
                end3[k] = str(end3[k] * ratio[i])
            print >>fo, names[i] + '_end5_' + j + ' <- c(' + ', '.join(end5) + ')'
            print >>fo, names[i] + '_end3_' + j + ' <- c(' + ', '.join(end3) + ')'
    
    x = []
    begin = int(-1 * (len(end5) / 2))
    for k in range(len(end5)):
        x.append(str(begin))
        begin += 1
    
    print >>fo, 'x <- c(' + ', '.join(x) + ')\n' 
    print >>fo, "pdf('Yeast_5_end_UCSC.pdf', width = 10, height = 8.5)"
    print >>fo, "par(mfrow=c(2, 2))"
    gene5end = ['', '_5end_convergent', '_5end_shared_others_5end', '_5end_shared_others_3end']
    mains = ['All genes (5,768)', 'Convergent genes (2,559)', 'Sharing other gene 5 end (1,419)', 'Sharing other gene 3 end (1,854)']
    for i in range(len(gene5end)):
        print >>fo, "plot(x, ACF1_end5_" + gene5end[i] +", col = 2, type = 'l', main = '5 End Tag Profile (" + mains[i] + ")', xlab = 'Locations relative to 5 end (bp)', ylab = 'average tag density')"
        for k in range(1, 4):
            print >>fo, "lines(x, " + names[k] + "_end5_" + gene5end[i]  + ", col = " + str(k + 3) + ", type = 'l')"
        print >>fo, "legend('topright',  c('" + "', '".join(name2) + "'), col = c(2,4,5,6), lty = 1)"
    print >>fo, "dev.off()"
    
    print >>fo, "pdf('Yeast_3_end_UCSC.pdf', width = 10, height = 8.5)"
    print >>fo, "par(mfrow=c(2, 2))"
    gene5end = ['', '_3end_convergent', '_3end_shared_others_3end', '_3end_shared_others_5end']
    mains = ['All genes (5,768)', 'Convergent genes (1,404)', 'Sharing other gene 3 end (2,579)', 'Sharing other gene 5 end (1,810)']
    for i in range(len(gene5end)):
        print >>fo, "plot(x, ACF1_end3_" + gene5end[i] +", col = 2, type = 'l', main = '3 End Tag Profile (" + mains[i] + ")', xlab = 'Locations relative to 3 end (bp)', ylab = 'average tag density')"
        for k in range(1, 4):
            print >>fo, "lines(x, " + names[k] + "_end3_" + gene5end[i]  + ", col = " + str(k + 3) + ", type = 'l')"
        print >>fo, "legend('topright',  c('" + "', '".join(name2) + "'), col = c(2,4,5,6), lty = 1)"
    print >>fo, "dev.off()"
    fo.close()
        
def test():
    tagc1 = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/BEDfile/C_MNase_Ecoli.bed', 'ecoli')
    othervalue1 = tagc1.return_hash()
    tagc2 = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/BEDfile/C_SONIC_Ecoli.bed', 'ecoli')
    othervalue2 = tagc2.return_hash()
    tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/BEDfile/ACF1_Ecoli.bed', 'ecoli')
    tag.normarlize('/home/zy/collaboration/Nuc_Struhl/Data/WIGfile_normalized/ACF1_Ecoli_MNase.wig', othervalue1, 1.24, 0.61)
    tag.normarlize('/home/zy/collaboration/Nuc_Struhl/Data/WIGfile_normalized/ACF1_Ecoli_SONIC.wig', othervalue2, 1.24, 1.09)
    tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/BEDfile/SALT_Ecoli.bed', 'ecoli')
    tag.normarlize('/home/zy/collaboration/Nuc_Struhl/Data/WIGfile_normalized/SALT_Ecoli_MNase.wig', othervalue1, 3.86, 0.61)
    tag.normarlize('/home/zy/collaboration/Nuc_Struhl/Data/WIGfile_normalized/SALT_Ecoli_SONIC.wig', othervalue2, 3.86, 1.09)
    
    tagc1 = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/BEDfile/C_MNase_yeast.bed', 'yeast')
    othervalue1 = tagc1.return_hash()
    tagc2 = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/BEDfile/C_SONIC_yeast.bed', 'yeast')
    othervalue2 = tagc2.return_hash()
    tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/BEDfile/ACF1_yeast.bed', 'yeast')
    tag.normarlize('/home/zy/collaboration/Nuc_Struhl/Data/WIGfile_normalized/ACF1_yeast_MNase.wig', othervalue1, 0.93, 2.83)
    tag.normarlize('/home/zy/collaboration/Nuc_Struhl/Data/WIGfile_normalized/ACF1_yeast_SONIC.wig', othervalue2, 0.93, 2.55)
    tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/BEDfile/SALT_yeast.bed', 'yeast')
    tag.normarlize('/home/zy/collaboration/Nuc_Struhl/Data/WIGfile_normalized/SALT_yeast_MNase.wig', othervalue1, 1.04, 2.83)
    tag.normarlize('/home/zy/collaboration/Nuc_Struhl/Data/WIGfile_normalized/SALT_yeast_SONIC.wig', othervalue2, 1.04, 2.55)

def application_enriched_loci():
    names = ['ACF1', 'SALT', 'C_MNase', 'C_SONIC']
    cutoffs = [3, 3, 3, 3]
    fo = open('/home/zy/collaboration/Nuc_Struhl/Analysis/Enriched_loci/Yeast_enriched_loci_center_pattern_SGD.R', 'w')
    print >>fo, "pdf('Yeast_enriched_loci_center_patternSGD.pdf', width = 10, height = 8.5)"
    print >>fo, "par(mfrow=c(2, 2))"
    print >>fo, "x <- seq(1, 75)"
    for i in range(4):
        name = names[i]
        cutoff = cutoffs[i]
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/SGD_version/BEDfile/' + name + '_yeast.bed', 'yeast_SGD', extension = 1, shift = 0, strand = '+')
        value = tag.enriched_loci_center_pattern(cutoff, '/home/zy/collaboration/Nuc_Struhl/Analysis/Enriched_loci/Yeast_enriched_loci_' + name + '_+_SGD.txt', strand = '+')
        print >>fo, "value1 <-c(" + ', '.join(value) + ")"
        #print >>fo, "plot(x, value1, col = 2, type = 'l', main = '" + name + "', ylim =c(0, 0.3))"
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/SGD_version/BEDfile/' + name + '_yeast.bed', 'yeast_SGD', extension = 1, shift = 0, strand = '-')
        value = tag.enriched_loci_center_pattern(cutoff, '/home/zy/collaboration/Nuc_Struhl/Analysis/Enriched_loci/Yeast_enriched_loci_' + name + '_-_SGD.txt', strand = '-')
        print >>fo, "value2 <-c(" + ', '.join(value) + ")"
        print >>fo, "value <- (value1 + value2) / 2"
        print >>fo, "plot(x, value, col = 2, type = 'l', main = '" + name + "', ylim =c(0, 0.2), xlab = 'Distance to local summit locus (bp)', ylab = '% cluster having tag enriched loci')"
        print >>fo, "abline(v = c(10, 20, 30, 40, 50, 60, 70), col = 'grey', lty = 3)"
        #print >>fo, "lines(x, value2, col = 4, type = 'l')"
        #print >>fo, "legend('topright', c('+', '-'), col = c(2,4), lty = 1)"
    print >>fo, "dev.off()"
    fo.close()


def wig_genomic_feature_plot_plosbio(Rscript, version = 'SGD', remove = 0.0, cluster = 0):
    fo = open(Rscript, 'a')
    values = {'_C' + str(cluster) + '.bed' : 'All Genes in Cluster ' + str(cluster), '_C' + str(cluster) + '_5_end_convergent.bed' : '5 End Convergent Genes in Cluster ' + str(cluster)}
    
    x = []
    begin = -750
    for k in range(1501):
        x.append(str(begin))
        begin += 1
    print >>fo, 'x <- c(' + ', '.join(x) + ')\n'
    
    # 5' end 10 bp extension
    print >>fo, "pdf('Yeast_5_end_in_vivo_10bp_" + version + "_cluster.pdf', width = 10, height = 17)"
    print >>fo, "par(mfrow=c(4, 2))"
    
    for gene in values.keys(): 
        mylist = wig_genomic_feature_two_sample_plot_plosbio(gene, extension = 10, shift = 68, strand = 'both', end = 5, version = version, width = 1501, multiple = 1, remove = remove)
        for n in ['normal', 'heatshock']:
            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
        print >>fo, "plot(x, normal, col = 2, type = 'l', ylim = c(0, 1), main = '" + values[gene] + " ()', xlab = 'Locations relative to 5 end (bp)', ylab = 'average tag profile')"
        print >>fo, "lines(x, heatshock, col = 4, type = 'l')"
        print >>fo, "legend('topright',  c('" + "', '".join(['normal', 'heatshock']) + "'), col = c(2,4), lty = 1)"
        
    print >>fo, 'dev.off()'
    
    # 5' end 73 bp extension
    print >>fo, "pdf('Yeast_5_end_in_vivo_73bp_" + version + "_cluster.pdf', width = 10, height = 17)"
    print >>fo, "par(mfrow=c(4, 2))"
    
    for gene in values.keys(): 
        mylist = wig_genomic_feature_two_sample_plot_plosbio(gene, extension = 73, shift = 36, strand = 'both', end = 5, version = version, width = 1501, multiple = 1, remove = remove)
        for n in ['normal', 'heatshock']:
            print >>fo, n + ' <- c(' + ', '.join(mylist[n]) + ')'
        print >>fo, "plot(x, normal, col = 2, type = 'l', ylim = c(0, 1), main = '" + values[gene] + " ()', xlab = 'Locations relative to 5 end (bp)', ylab = 'average tag profile')"
        print >>fo, "lines(x, heatshock, col = 4, type = 'l')"
        print >>fo, "legend('topright',  c('" + "', '".join(['normal', 'heatshock']) + "'), col = c(2,4), lty = 1)"
        
    print >>fo, 'dev.off()'    
    fo.close()
    
def wig_genomic_feature_two_sample_plot_plosbio(genename, extension = 73, shift = 36, strand = 'both', end = 5, version = 'SGD', width = 1501, multiple = 1, remove = 0.0):
    samples = ['normal', 'heatshock']
    reftags = [0.042 * extension, 0.086 * extension]
    #ratio = [0.93, 1.04, 2.83, 2.55]
    
    mylist = {}
    
    for i in range(len(samples)):
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/in_vivo/PlosBio_SGD/BEDfile/PlosBio_' + samples[i] + '.bed', 'yeast_' + version, extension = extension, shift = shift, strand = strand)
        if strand == 'both':
            (end5, end3) = tag.wig_genomic_feature_normalize('/home/zy/collaboration/Nuc_Struhl/Data/Genome/Annotation/TSS_cluster/SGD_gene_' + version + genename, width = width, referencetag = reftags[i], remove = remove)
        else:
            (end5, end3) = tag.wig_genomic_feature_normalize('/home/zy/collaboration/Nuc_Struhl/Data/Genome/Annotation/TSS_cluster/SGD_gene_' + version + genename, width = width, referencetag = reftags[i], strand = '+', remove = remove)
        mylist[samples[i]] = []
        if end == 5:
            for k in range(len(end5)):
                mylist[samples[i]].append(str(end5[k]))
                #mylist[samples[i]].append(str(multiple * end5[k] * ratio[i] / extension))
        else:
            for k in range(len(end3)):
                mylist[samples[i]].append(str(end3[k]))
                #mylist[samples[i]].append(str(multiple * end3[k] * ratio[i] / extension))
    return mylist

def plosbio():
    for name in ['heatshock', 'normal']:
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/in_vivo/PlosBio_SGD/BEDfile/PlosBio_' + name + '.bed', 'yeast_SGD', extension = 73, shift = 36)
        tag.tag_to_wig('/home/zy/collaboration/Nuc_Struhl/Data/in_vivo/PlosBio_SGD/WIGfile_73bp/PlosBio_' + name + '.wig')
        tag = SeqTag('/home/zy/collaboration/Nuc_Struhl/Data/in_vivo/PlosBio_SGD/BEDfile/PlosBio_' + name + '.bed', 'yeast_SGD', extension = 10, shift = 68)
        tag.tag_to_wig('/home/zy/collaboration/Nuc_Struhl/Data/in_vivo/PlosBio_SGD/WIGfile_10bp/PlosBio_' + name + '.wig')
    for i in range(4):
        wig_genomic_feature_plot_plosbio('/home/zy/collaboration/Nuc_Struhl/Analysis/Gene_features/profile/Gene_5_end_in_vivo_feature_SGD_080604.R', version = 'SGD', cluster = i)
def NFR_xls():
    import sys
    tagfile = sys.argv[1]
    genefile = sys.argv[2]
    output = sys.argv[3]
    tag = SeqTag(tagfile, 'yeast_SGD', extension = 73, shift = 36)
    tag.gene_promoter_profile(genefile, output)

def NFR_xls_2():
    import sys
    tagfile = sys.argv[1]
    genefile = sys.argv[2]
    output = sys.argv[3]
    tag = SeqTag(tagfile, 'yeast_SGD', extension = 73, shift = 36)
    tag.gene_promoter_profile_2(genefile, output)
    
if __name__ == '__main__':
    #wig_genomic_feature_plot('/home/zy/collaboration/Nuc_Struhl/Analysis/Gene_features/Gene_5_3_end_feature_SGD_080509.R', version = 'SGD')
    #for i in range(4):
    #    wig_genomic_feature_plot('/home/zy/collaboration/Nuc_Struhl/Analysis/Gene_features/profile/Gene_5_end_feature_SGD_080603.R', version = 'SGD', cluster = i)
    #application_enriched_loci()
    #plosbio()
    NFR_xls_2()
