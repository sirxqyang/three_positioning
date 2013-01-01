#!/usr/bin/python

# Yong Zhang @ DFCI

import sys, time
#from numpy import *

chrlength = {}

chrlength['ecoli'] = {'U00096' : 4639675}
chrlength['yeast_UCSC'] = {'chr1' : 230208,  'chr2' : 813136,  'chr3' : 316613,  'chr4' : 1531914,
                      'chr5' : 576869,  'chr6' : 270148,  'chr7' : 1090944, 'chr8' : 562639,
                      'chr9' : 439885,  'chr10' : 745446,  'chr11' : 666445,  'chr12' : 1078173,
                      'chr13' : 924430,  'chr14' : 784328,  'chr15' : 1091285, 'chr16' : 948060,
                      'chrM' : 85779}

chrlength['yeast_SGD'] = {'chr01' : 230208,  'chr02' : 813178,  'chr03' : 316617,  'chr04' : 1531918,
                      'chr05' : 576869,  'chr06' : 270148,  'chr07' : 1090947, 'chr08' : 562643,
                      'chr09' : 439885,  'chr10' : 745745,  'chr11' : 666454,  'chr12' : 1078175,
                      'chr13' : 924429,  'chr14' : 784333,  'chr15' : 1091289, 'chr16' : 948062,
                      'chrmt' : 85779}

chrlength['zv9_test'] = { 'chr1':60348388}

chrlength['zv9'] = { 'chr1':60348388,'chr2':60300536,'chr3':63268876,'chr4':62094675,
'chr5':75682077,'chr6':59938731,'chr7':77276063,'chr8':56184765,'chr9':58232459,'chr10':46591166,'chr11':46661319,
'chr12':50697278,'chr13':54093808,'chr14':53733891,'chr15':47442429,'chr16':58780683,'chr17':53984731,'chr18':49877488,
'chr19':50254551,'chr20':55952140,'chr21':44544065,'chr22':42261000,'chr23':46386876,'chr24':43947580,'chr25':38499472}

chrlength['zv9_part1'] = { 'chr1':60348388,'chr2':60300536,'chr3':63268876,'chr4':62094675,
'chr5':75682077,'chr6':59938731,'chr7':77276063,'chr8':56184765,'chr9':58232459,'chr10':46591166,'chr11':46661319,
'chr12':50697278,'chr13':54093808}

chrlength['zv9_part2'] = { 'chr14':53733891,'chr15':47442429,'chr16':58780683,'chr17':53984731,'chr18':49877488,
'chr19':50254551,'chr20':55952140,'chr21':44544065,'chr22':42261000,'chr23':46386876,'chr24':43947580,'chr25':38499472}

chrlength['hg19'] = { 'chr1':249250621,'chr2':243199373,'chr3':198022430,'chr4':191154276,
'chr5':180915260,'chr6':171115067,'chr7':159138663,'chr8':146364022,'chr9':141213431,'chr10':135534747,'chr11':135006516,
'chr12':133851895,'chr13':115169878,'chr14':107349540,'chr15':102531392,'chr16':90354753,'chr17':81195210,'chr18':78077248,
'chr19':59128983,'chr20':63025520,'chr21':48129895,'chr22':51304566,'chrX':155270560,'chrY':59373565,'chrM':16571 }


import time
time.asctime()
class Genome(object):
    
    def __init__(self, filename):
        temp = []
        for line in open(filename).xreadlines():
            if line[0] == '>':
                continue
            temp.append(line.strip())
        self.seq = ''.join(temp)
    def length(self):
        return len(self.seq)
    
    def ATcontent(self, length, width = 146):
        at = []
        for k in range(length - width):
            subseq = self.seq[k : k + width]
            numat = subseq.count('A') + subseq.count('T') + subseq.count('a') + subseq.count('t')
            numcg = subseq.count('C') + subseq.count('G') + subseq.count('c') + subseq.count('g')
            if (numat + numcg) >= 0.8 * width:
                at.append(str(k + width / 2) + '\t' + str(int(1000 * numat / (numat + numcg)) / 10.0))
        return at

lastfix = {'ecoli' : '.fna', 'yeast_SGD' : '.fsa', 'yeast_UCSC' : '.fa','zv9':'.fa','zv9_test':'.fa', 'hg19':'.fa'}

dirname = '/mnt/Storage/home/yangxq/project/Nucleosome/annotation/human/chrfa'

class whole_genome(object):
    def __init__(self, spename):
        if not chrlength.has_key(spename):
            print 'Wrong species name:'
            sys.exit()
        self.seq = {}
        for chrname in chrlength[spename].keys():
            temp = []
            for line in open(dirname + '/' + chrname + lastfix[spename]).xreadlines():
                if line[0] == '>':
                    continue
                temp.append(line.strip().upper())
            self.seq[chrname] = ''.join(temp)
            
    def getseq(self):
        return self.seq
    def getnucleosomecode(self):
        code = {}
        codenew = {}
        for chrname in self.seq.keys():
            code[chrname] = []
            for i in range(len(self.seq[chrname]) - 1):
                pa = self.seq[chrname][i:i+2]
                if pa == 'AT' or pa == 'AA' or pa == 'TA' or pa == 'TT':
                    code[chrname].append('1')
                else:
                    code[chrname].append('0')
            code[chrname].append('0')
            codenew[chrname] = ''.join(code)
        return code
        
   

def test():
#    dirname = '/home/zy/collaboration/Nuc_Struhl/Data/Genome/ecoli/'
#    fo = open('/home/zy/collaboration/Nuc_Struhl/Analysis/ATcontent/Ecoli_AT.wig', 'w')
#    for name in chrlength['ecoli'].keys():
#        genome = Genome(dirname + name + '.fna')
#        at = genome.ATcontent(chrlength['ecoli'][name])
#        print >>fo, "track type=wiggle_0"
#        print >>fo, "variableStep  chrom=" + name
#        for k in xrange(len(at)):
#            print >>fo, at[k]
#    fo.close()
    
    dirname = '/mnt/Storage/home/yangxq/project/Nucleosome/annotation/human/chrfa'
    fo = open('/mnt/Storage/home/yangxq/project/Nucleosome/GSM945580/run/Human_1k.wig', 'w')
    for name in chrlength['hg19'].keys():
        genome = Genome(dirname + name + '.fa')
        length = genome.length()
        print name + ' : ' + str(length)
        at = genome.ATcontent(length, width= 1000)
        print >>fo, "track type=wiggle_0"
        print >>fo, "variableStep  chrom=" + name
        for k in xrange(len(at)):
            print >>fo, at[k]
    fo.close()

def test_convergent_genes():
    gene = Gene('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD.bed')
    gene.convergent_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_5end_convergent.bed', end = 5)
    gene.convergent_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_3end_convergent.bed', end = 3)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_5end_shared_others_5end.bed', end = 5, end_others = 5)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_5end_shared_others_3end.bed', end = 5, end_others = 3)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_3end_shared_others_5end.bed', end = 3, end_others = 5)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_SGD_3end_shared_others_3end.bed', end = 3, end_others = 3)
    
    gene = Gene('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC.bed')
    gene.convergent_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_5end_convergent.bed', end = 5)
    gene.convergent_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_3end_convergent.bed', end = 3)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_5end_shared_others_5end.bed', end = 5, end_others = 5)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_5end_shared_others_3end.bed', end = 5, end_others = 3)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_3end_shared_others_5end.bed', end = 3, end_others = 5)
    gene.share_ends_genes('/Users/yzhang/Bioinfo/Nuc_Struhl/Genes_features/Annotation/SGD_gene_UCSC_3end_shared_others_3end.bed', end = 3, end_others = 3)
    
if __name__ == '__main__':
    #test_convergent_genes()
    test()
    
