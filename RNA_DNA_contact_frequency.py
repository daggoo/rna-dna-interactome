import ipyparallel
import collections
import os
import zmq as pyz
import gc
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import matplotlib.ticker as ticker
from matplotlib.ticker import Locator
from collections import Counter
import re
import itertools
import matplotlib as mpl
import time
mpl.rcParams['agg.path.chunksize'] = 30000

import gzip

class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in xrange(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))


#os.system("ipcluster stop --profile=natalia &")
#os.system("ipcluster start --profile=natalia -n 4 &")
#os.system("ipcluster start --profile=natalia -n 2 &")

number_of_cpu = 8

path = "/home/ubuntu/RNAdata/Dekker/"

k562_rna = "/home/ubuntu/RNAdata/K562.annot.full.tab"
k562_dna = "/home/ubuntu/RNAdata/Dekker/K562_all_contacts.pairs.txt.gz"
introns_file = "/home/ubuntu/RNAdata/Dekker/gencode.v27.gffutils.intron.processed.bed"
exons_file = "/home/ubuntu/RNAdata/Dekker/gencode.v27.gffutils.exon.processed.bed"
genes_file = "/home/ubuntu/RNAdata/Dekker/gencode.v27.gffutils.genes.processed.bed"

k562_stat = "/home/ubuntu/RNAdata/Dekker/hg19_chrms_length"

chrms = []
for j in range (0, 5):    
    tmp = []
    for i in range (1, 5):
        tmp.append('chr' + str(j * 4 + i))
    chrms.append(tmp)

tmp = []
tmp.append('chr21')
tmp.append('chr22')
tmp.append('chrX')
tmp.append('chrY')
chrms.append(tmp)

chrms_flat = ["chr" + str(i) for i in range(1,23)]
chrms_flat.append("chrX")
chrms_flat.append("chrY")

print "Chrms flat list"
print chrms_flat

exons = {}
exons_by_chr_by_bin = {}

exons_strands = []
last_exons = {}

introns = {}
introns_by_chr_by_bin = {}

introns_strands = []
genes = {}
genes_length = {}
genes_by_chr_by_bin = {}
max_gene_length = 0
min_gene_length = 1000000000

bin_size = 0
bin_size_to_split_annotations = 1000000
binId = lambda x: int(x/bin_size)

chrms_length = {}
chrms_ordered_by_length = []
#
#
def load_length(fl_in):
    global chrms_length, chrms_ordered_by_length
    with open(fl_in, 'r') as tab_fl:
        for i, l in enumerate(tab_fl):
            ln = l.rstrip().split("\t")
            if len(ln)==2:
                #print ln[0] + "\t" + ln[-2]
                chrms_length[ln[0]] = int(ln[1])
                chrms_ordered_by_length.append(ln[0])
#   
#    
def load_coordinates(fl_in, df, mode, df_by_bin):
    df_tmp = {}
    count_strange = 0
    with open(fl_in, 'r') as fl_raw:
            for i, l in enumerate(fl_raw):                
                ln = l.rstrip().split("\t")                
                
                if ln[8]!='protein_coding':
                    continue
                
                if len(ln)<9:
                    print "Warning: load_coordinates short line " + str(ln)
                    continue
                
                chr = ln[0]
                start = int(ln[1])
                end = int(ln[2])
                strand = ln[5]
                id = ln[7]
                bn = int(start/bin_size_to_split_annotations)
                
                length = end - start
                tul = [start, end, strand, length]  
                
                if start > end:
                    count_strange+=1
                    
                if mode=="exons":
                    tul2 = [start, end, strand, length, chr]   
                    if id not in last_exons:
                        last_exons[id] = tul2
                    else:
                        if start>=last_exons[id][0]:
                            last_exons[id] = tul2                         
                
                
                if chr not in df_tmp:
                    df_tmp[chr] = []
                    df_by_bin[chr] = {}
                
                if bn not in df_by_bin[chr]:
                    df_by_bin[chr][bn] = []
                        
                df_tmp[chr].append(tul)
                df_by_bin[chr][bn].append(tul)
                
                                 
    for chr in sorted(df_tmp.keys()):
        df[chr] = np.array(df_tmp[chr])
        print chr ,
        print len(df_tmp[chr]) ,
        print "| " ,
    
    print "\ncount_strange = " ,
    print count_strange                    
#
#
def load_gene_data (fl_in, df):
    global max_gene_length
    global min_gene_length
    df_tmp = {}
    
    count_pos_strand = 0
    count_neg_strand = 0
    with open(fl_in, 'r') as fl_raw:
            for i, l in enumerate(fl_raw):                
                ln = l.rstrip().split("\t")  
                if ln[8]!='protein_coding':
                    continue                
                if len(ln)<9:
                    print "Warning: load_coordinates short line " + str(ln)
                    continue
                
                id = ln[7]
                chr = ln[0]
                start = int(ln[1])
                end = int(ln[2])
                strand = ln[5]
                length = abs(end - start)
                
                '''
                if binId(length) - (binId(end) - binId(start)) !=-1:
                    print id ,
                    print binId(length) - (binId(end) - binId(start))
                    print start ,
                    print end , 
                    print length
                '''    
                
                if length > max_gene_length:
                    max_gene_length = length
                if length < min_gene_length:
                    min_gene_length = length 
                    
                if id not in genes_length:
                    if strand=="+":
                        genes_length[id] = [start, end, strand, length, chr]          
                        bn = int(start/bin_size_to_split_annotations)
                    elif strand=="-":
                        genes_length[id]=[start, end, strand, length, chr]          
                        bn = int(start/bin_size_to_split_annotations)                        
                        #genes_length[id] = [end, start, strand, length, chr] 
                        #bn = int(end/bin_size_to_split_annotations)
                else: 
                    continue
                    #print l.rstrip()
                    #print id ,
                    #print genes_length[id]
                    
                if strand=="+":
                    tul=[start, end, strand, length, id, chr] 
                    count_pos_strand+=1
                elif strand=="-":
                    #tul=[end, start, strand, length, id, chr]
                    tul = [start, end, strand, length, id, chr] 
                    count_neg_strand+=1
                else:
                    print "Warning: load_gene_data unknown strand " + l
                    break
                
                if chr not in df_tmp:
                    df_tmp[chr] = []    
                    genes_by_chr_by_bin[chr] = {}
                
                if bn not in genes_by_chr_by_bin[chr]:
                    genes_by_chr_by_bin[chr][bn] = []
                
                df_tmp[chr].append(tul)
                genes_by_chr_by_bin[chr][bn].append(tul)
                '''
                if id=="AC073333.1":
                    print chr
                    print bn
                    print tul
                    print genes_by_chr_by_bin[chr]
                    print "Strange"
                    print genes_by_chr_by_bin[chr][bn]
                '''
    print "Number of genes per chromosome"                
    for chr in sorted(df_tmp.keys()):
        df[chr] = np.array(df_tmp[chr])
        print chr ,
        print len(df_tmp[chr]) ,     
        print "| " ,
    print "\n count_pos_strand = " ,
    print count_pos_strand
    
    print "\n count_neg_strand = " ,
    print count_neg_strand    
#
#
def split_mrnas_by_chr_of_rna_part(fl_in, fl_out_prefix):
    list_of_chmrs = []
    chrm_files = {}
    for chr in chrms_flat:
        chrm_files[chr] = open(fl_out_prefix + "_" + chr + ".tab", 'w+')
    with open(fl_in, 'r') as fl_raw:         
        for i, l in enumerate(fl_raw):  
            if 'protein_coding'  in l:
                ln = l.rstrip().split("\t")    
                chrm_files[ln[0]].write(l)
    for chr in chrms_flat:
        chrm_files[chr].close()      
#
#
def split_dna_contacts_by_chr_of_first_dna_part(fl_in, fl_out_prefix):
    list_of_chmrs = []
    chrm_files = {}
    add_chrms = []
    for chr in chrms_flat:
        chrm_files[chr] = gzip.open(fl_out_prefix + "_" + chr + ".tab", 'wt')
    with gzip.open(fl_in, 'r') as fl_raw:         
        for i, l in enumerate(fl_raw):
            ln = l.rstrip().split("\t") 
            if ln[1] in chrms_flat:
                chrm_files[ln[1]].write(l)
            else:
                if ln[1] not in add_chrms: 
                    add_chrms.append(ln[1])
    print "Additional chromosomes"
    print add_chrms
    for chr in chrms_flat:
        chrm_files[chr].close() 
#
#
#def core_parse_chrms_data_vs_annotated_coordinates(lines, fl_out, mode, df_starts, df_ends, df_strands):    
def core_parse_chrms_data_vs_annotated_coordinates( lines, mode, df_starts, df_ends, df_strands):  
    result = []
    #import numpy as np
    if mode=='RNA':
        for line in lines:
            if line is None:
                continue
            ln = line.rstrip().split("\t")
            coord1 = int(ln[1])
            bn = int(coord1/bin_size_to_split_annotations)
            chr = ln[0]            
            rna_strand = ln[8]   
            if numpy.any((df_starts[bn] <= coord1) & (df_ends[bn] >= coord1) & (df_strands[bn] == rna_strand)):
                result.append(line) 
            else:
                result.append("\n")                 
            
    elif mode=='DNA':   
        for line in lines:
            if line is None:
                continue        
            ln = line.rstrip().split("\t") 
            coord1 = int(ln[2])
            bn = int(coord1/bin_size_to_split_annotations)
            chr = ln[1]  
            if numpy.any((df_starts[bn] <= coord1) & (df_ends[bn] >= coord1)):
                result.append(line) 
            else:
                result.append("\n")            
    else:
        return "Warning: core_parse_chrms_data_vs_annotated_coordinates unknown mode\n"   
    
    return result
                            
#
#
def parse_chrms_data_vs_annotated_coordinates(path, df, df_name, mode, number_of_cpu):
          
    filenames_rna = [] 
    filenames_dna = []
    
    for i in range(len(chrms_ordered_by_length)):
        if chrms_ordered_by_length[i]=="chrY":
            continue
        filenames_rna.append("k562_rna_" + chrms_ordered_by_length[i]+".tab")
        filenames_dna.append("k562_dna_" + chrms_ordered_by_length[i]+".tab")
        
  
    client = ipyparallel.Client(profile='natalia') 
    dview = client[:]
    print "clients ids ",
    print client.ids
    dview.block = False  
        
    results = [None]*number_of_cpu
    if mode=="DNA":
        for fl in filenames_dna:     
            line_counter=0
            with open (os.path.join(path , fl), 'r') as inLines,\
                 open (os.path.join(path , fl.split(".tab")[0] + "_" + df_name + ".tab"), 'w+') as out_file:            
                chr=fl.split("_")[2].split(".tab")[0]
                if chr not in chrms_flat:
                    print "Warning: parse_chrms_data_vs_annotated_coordinates unknow chromosome ",
                    print chr
                    return         
                print chr ,
                print fl 
                
                df_starts = np.array(df[chr][:][:,0],  dtype='int')
                df_ends = np.array(df[chr][:][:,1], dtype='int')
                #df_strands=np.array(df[chr][:,2], dtype='str')     
                df_strands = []    
             
                with dview.sync_imports():
                    import numpy
                t = time.time()    
                
                #for k, l in enumerate (inLines):
                for lines in itertools.izip_longest(*[inLines]*((1000*number_of_cpu)*number_of_cpu)):
                    line_counter+=1000*number_of_cpu*number_of_cpu
                    
                    if (line_counter%(1000*number_of_cpu*1000))==0:
                        t3 = time.time() - t
                        print line_counter ,
                        print('%f secs (multicore)' % t3)
                        t = time.time() 
                    for engine in range(number_of_cpu):
                        #print engine , 
                        dview.targets = [engine]
                        line1 = engine * 1000 * number_of_cpu
                        line3 = engine * 1000 * number_of_cpu + 1000 * number_of_cpu - 1
                        
                        results[engine] = dview.apply(core_parse_chrms_data_vs_annotated_coordinates, 
                                                 lines[line1:line3],
                                                 "DNA",
                                                 df_starts, 
                                                 df_ends, 
                                                 df_strands)
                                
                    
                    #print '\nWaiting dnas...\n'
                    dview.wait(results)                            
                    for engine in range(number_of_cpu):
                        #if results[engine][0]!="":
                        for line in results[engine][0]:
                            out_file.write(line)
                print 'Number of lines appr',
                print line_counter
           # client.close()
            gc.collect()     
    elif mode=="RNA":
        for fl in filenames_rna:     
            line_counter=0
            with open (os.path.join(path , fl), 'r') as inLines,\
                 open (os.path.join(path , fl.split(".tab")[0] + "_" + df_name + ".tab"), 'w+') as out_file:            
                chr=fl.split("_")[2].split(".tab")[0]
                if chr not in chrms_flat:
                    print "Warning: parse_chrms_data_vs_annotated_coordinates unknow chromosome ",
                    print chr
                    return         
                print chr ,
                print fl 
                
                df_starts=np.array(df[chr][:][:,0],  dtype='int')
                df_ends=np.array(df[chr][:][:,1], dtype='int')
                df_strands=np.array(df[chr][:][:,2], dtype='str')     
             
                with dview.sync_imports():
                    import numpy
                t = time.time()    
                
                for lines in itertools.izip_longest(*[inLines]*((1000*number_of_cpu)*number_of_cpu)):
                    line_counter+=1000*number_of_cpu*number_of_cpu
                    
                    if (line_counter%(1000*number_of_cpu*1000))==0:
                        t3 = time.time() - t
                        print line_counter ,
                        print('%f secs (multicore)' % t3)
                        t = time.time() 
                    for engine in range(number_of_cpu):
                        dview.targets = [engine]
                        line1 = engine * 1000 * number_of_cpu
                        line3 = engine * 1000 * number_of_cpu + 1000 * number_of_cpu - 1
                        
                        results[engine] = dview.apply(core_parse_chrms_data_vs_annotated_coordinates, 
                                                 lines[line1:line3],
                                                 "RNA",
                                                 df_starts, 
                                                 df_ends, 
                                                 df_strands)

                    dview.wait(results)                            
                    for engine in range(number_of_cpu):
                        for line in results[engine][0]:
                            out_file.write(line)
                print 'Number of lines appr',
                print line_counter
            gc.collect()             
        
                
#
#
def filter_out_last_exon(path):
    files =["k562_dna_all_genes_exons.tab",  "k562_rna_all_genes_exons.tab"]
    
    for fl in files:
        with open(fl, 'r') as inLines, \
             open (fl.split(".tab")[0]+"_wo_last_exon.tab", 'w+') as output:

            for i, l in enumerate(inLines):
                ln=l.rstrip().split("\t")
                if "rna" in fl:
                    id=ln[10]
                    start=int(ln[1])
                    chr=ln[0]
                    
                elif "dna" in fl:
                    id=ln[8]
                    start=int(ln[2])
                    chr=ln[1]
                    
                if id not in genes_length:
                    if "_" in id:
                        tmp_id=id.split("_")[0]
                        if tmp_id in genes_length:
                            id=tmp_id
                        else:
                            #print id
                            continue                    
                if (start >= last_exons[id][0]) and (start<=last_exons[id][1]):
                    pass
                else:
                    output.write(l)
                
#
#
def core_annotate_DNA_contacts_left_part(lines, genes_by_bin):      
    result=[]    
    pos=2            
    for line in lines:
        if line in ['\n', '\r\n']:
            continue              
        if line is None:
            continue        
        ln=line.rstrip().split("\t")
        bn=int(int(ln[2])/bin_size_to_split_annotations)
        
        '''
        for bn in genes_by_bin:
            for tmp in genes_by_bin[bn]:
                if tmp[5]!=ln[1]:
                    return [line.rstrip() + "\t" + str(tmp)]
        '''
        if bn not in genes_by_bin:
            continue
        for tmp in genes_by_bin[bn]:    
            if tmp[5]!=ln[1]:
                return ["Warning:" + line.rstrip() + "\t" + str(tmp)]
            
                
                #
            if int(ln[2])>=int(tmp[0]) and int(ln[2])<=int(tmp[1]) and tmp[2]==ln[5]:  
                result.append(line.rstrip() + "\t" + tmp[4] + "\n")   
        
    return result
#
#
def annotate_DNA_contacts_left_part (path, df_name, number_of_cpu):
    print "Annotating " + df_name 
    filenames_dna=[]
    results=[None]*number_of_cpu
    
    for i in range(len(chrms_ordered_by_length)):
        if chrms_ordered_by_length[i] in ["chrY"]:
            continue
        filenames_dna.append("k562_dna_" + chrms_ordered_by_length[i]+".tab")
        
    client = ipyparallel.Client(profile='natalia') 
    dview = client[:]
    print "clients ids ",
    print client.ids
    dview.block = False 
    
        
    for fl in filenames_dna:  
        with open (os.path.join(path , fl.split(".tab")[0] + "_" + df_name + ".tab"), 'r') as inLines, \
             open (os.path.join(path , fl.split(".tab")[0] + "_" + df_name + "_annotated_genes.tab"), 'w+') as out_file:   
              
            chr=fl.split("_")[2].split(".tab")[0]
            print chr ,
            print os.path.join(path , fl.split(".tab")[0] + "_" + df_name + ".tab") ,
            print os.path.join(path , fl.split(".tab")[0] + "_" + df_name + "_annotated_genes.tab")
            
            #dview['genes_by_bin']=genes_by_chr_by_bin[chr]
            
            #print genes_by_chr_by_bin[chr]
            
            dview['bin_size_to_split_annotations']=bin_size_to_split_annotations
            
            print 'computing'
            if chr not in chrms_flat:
                print "Warning: parse_chrms_data_vs_annotated_coordinates unknown chromosome ",
                print chr
                return         
            
            line_counter=0
            t = time.time() 
            for lines in itertools.izip_longest(*[inLines]*((100000*number_of_cpu)*number_of_cpu)):
                line_counter+=100000*number_of_cpu*number_of_cpu
                
                
                if (line_counter%(100000*number_of_cpu*100))==0:
                    t3 = time.time() - t
                    print line_counter ,
                    print('%f secs (multicore)' % t3)
                    t = time.time() 
                   
                for engine in range(number_of_cpu):
                    dview.targets = [engine]
                    line1=engine*100000*number_of_cpu
                    line3=engine*100000*number_of_cpu+100000*number_of_cpu-1                            
                    results[engine] = dview.apply(core_annotate_DNA_contacts_left_part, 
                                             lines[line1:line3],
                                             genes_by_chr_by_bin[chr])
                  
                dview.wait(results)                            
                for engine in range(number_of_cpu):  
                    for line in results[engine][0]: 
                        if "Warning" in line:
                            print line
                            #time.sleep(5)
                            break
                        out_file.write(line)      
                
    client.close()                        
#
#
def core_parse_gene_data_vs_annotated_coordinates(lines, genes_length, mode, flank):      
    #import numpy
    result=[]
    if mode=='DNA':   
        for line in lines:
            if line in ['\n', '\r\n']:
                continue              
            if line is None:
                continue        
            ln=line.rstrip().split("\t")
            if ln[1] != ln[3] :
                continue 
                
            #df_starts=numpy.array(genes[ln[1]][:,0],  dtype='int')
            #df_ends=numpy.array(genes[ln[1]][:,1],  dtype='int')
            #ids=numpy.array(genes[ln[1]][:,4],  dtype='str')
            '''
            b=(df_starts <= int(ln[2])) & (df_ends >= int(ln[2])) & (df_starts <= int(ln[4])) & (df_ends >= int(ln[4]))          
            for i in range(len(b)):
                if b[i]==1:
                    result.append(line.rstrip() + "\t" + ids[i] + "\n") 
            '''
            
            '''        
            if numpy.any((df_starts <= int(ln[2])) & (df_ends >= int(ln[2])) & (df_starts <= int(ln[4])) & (df_ends >= int(ln[4]))):
                result.append(line.rstrip() + "\t" + tmp[4]) 
            '''

            if int(ln[2])>=int(genes_length[ln[8]][0]) - flank  and int(ln[2])<=int(genes_length[ln[8]][1]) + flank  and int(ln[4])>=int(genes_length[ln[8]][0]) - flank   and int(ln[4])<=int(genes_length[ln[8]][1])  + flank:                    
                result.append(line)  
                    
            
    else:
        return ["Warning: core_parse_chrms_data_vs_annotated_coordinates unknown mode\n"]       
    return result
#
#
def parse_gene_data_vs_annotated_coordinates (path, df_name, mode, number_of_cpu, flank):
    print "Parsing " + df_name + "_" + mode
    filenames_rna=[] 
    filenames_dna=[]
    results=[None]*number_of_cpu
    for i in range(len(chrms_ordered_by_length)):
        if chrms_ordered_by_length[i]=="chrY":
            continue
        filenames_rna.append("k562_rna_" + chrms_ordered_by_length[i]+".tab")
        filenames_dna.append("k562_dna_" + chrms_ordered_by_length[i]+".tab")
        
        
    if mode=="DNA_exon_intron":
        client = ipyparallel.Client(profile='natalia') 
        dview = client[:]
        print "clients ids ",
        print client.ids
        dview.block = False         
        with open (os.path.join(path , "k562_dna_all_genes" + "_" + df_name + ".tab"), 'w+') as out_file:
            for fl in filenames_dna:     

                with open (os.path.join(path , fl.split(".tab")[0] + "_" + df_name + "_annotated_genes.tab"), 'r') as inLines:   
                    
                    chr=fl.split("_")[2].split(".tab")[0]
                    print chr

                    line_counter=0
                    if chr not in chrms_flat:
                        print "Warning: parse_chrms_data_vs_annotated_coordinates unknown chromosome ",
                        print chr
                        return         
                    
                    print chr ,
                    print os.path.join(path , fl.split(".tab")[0] + "_" + df_name + ".tab") 
                    #print os.system("cat " + os.path.join(path , fl.split(".tab")[0] + "_" + df_name + ".tab") + " | sed '/^\s*$/d' | wc -l" )

                    t = time.time() 
                    for lines in itertools.izip_longest(*[inLines]*((100000*number_of_cpu)*number_of_cpu)):
                        line_counter+=100000*number_of_cpu*number_of_cpu
    
                        if (line_counter%(100000*number_of_cpu*100))==0:
                            t3 = time.time() - t
                            print line_counter ,
                            print('%f secs (multicore)' % t3)
                            t = time.time() 
                            
                        for engine in range(number_of_cpu):
                            dview.targets = [engine]
                            line1=engine*100000*number_of_cpu
                            line3=engine*100000*number_of_cpu+100000*number_of_cpu-1                            
                            results[engine] = dview.apply(core_parse_gene_data_vs_annotated_coordinates, 
                                                     lines[line1:line3],
                                                     genes_length,
                                                     #df_starts,
                                                     #df_ends,
                                                     #ids,
                                                     "DNA",
                                                     flank)
                            
                        dview.wait(results)                            
                        for engine in range(number_of_cpu):                            
                            for line in results[engine][0]:
                                out_file.write(line)   

    elif mode=="RNA_exon_intron":
        count_strange=0
        with open (os.path.join(path , "k562_rna_all_genes" + "_" + df_name + "_" + str(flank) + ".tab"), 'w+') as out_file:
            for fl in filenames_rna:     
                
                with open (os.path.join(path , fl.split(".tab")[0] + "_" + df_name + ".tab"), 'r') as inLines:            
                    chr=fl.split("_")[2].split(".tab")[0]
                    if chr not in chrms_flat:
                        print "Warning: parse_chrms_data_vs_annotated_coordinates unknown chromosome ",
                        print chr
                        return         
                    print chr ,
                    print os.path.join(path , fl.split(".tab")[0] + "_" + df_name + ".tab")
                    
                    for i, line in enumerate (inLines): 
                        if line in ['\n', '\r\n']:
                            continue                        
                        ln=line.rstrip().split("\t")
                        id=ln[10] 
                        if id not in genes_length:
                            if "_" in id:
                                tmp_id=id.split("_")[0]
                                if tmp_id in genes_length:
                                    id=tmp_id
                                else:
                                    print id
                                    continue 
                        if int(ln[12]) > int(ln[13]):
                            count_strange+=1
                        if int(ln[1]) > int (ln[2]):
                            count_strange+=1
                        if int(ln[6]) > int (ln[7]):
                            count_strange+=1
                                
                        if (ln[11] == genes_length[id][4]) and (int(ln[12]) >= (int(genes_length[id][0]) - flank)) and (int(ln[12]) <= (int(genes_length[id][1]) + flank)):
                            out_file.write(line)
        print "\ncount_strange = " ,
        print count_strange                
    elif mode=="DNA_all":
        with open (os.path.join(path , "k562_dna_all_genes" + "_" + "all_regions" + ".tab"), 'w+') as out_file:
            with open (k562_dna, 'r') as inLines:  
                for i, line in enumerate (inLines): 
                    ln=line.rstrip().split("\t")
                    if (ln[1] != ln[3]) :
                        continue
                    for tmp in genes[ln[1]]:
                        if int(ln[2]) > tmp[0] and int(ln[2]) < tmp[1] and int(ln[4])  > tmp[0] and int(ln[4]) < tmp[1]:
                            out_file.write(line)
                            break
    elif mode=="RNA_all":
        with open (os.path.join(path , "k562_rna_all_genes" + "_" + "all_regions" + ".tab"), 'w+') as out_file:
            with open (k562_rna, 'r') as inLines:  
                for i, line in enumerate (inLines): 
                    ln=line.rstrip().split("\t")
                    if (ln[0] != ln[11]) or ln[9]!='protein_coding':
                        continue                    
                    if (int(ln[12]) > int(ln[6])) and (int(ln[12]) < int(ln[13])):
                        out_file.write(line)    
#
#
def parse_gene_data_vs_length (path, df, df_name, mode, portion):
    tmp=[sublist[3] for chr in genes.keys() for sublist in genes[chr]]
    tmp_genes=[]
    result=[]
    if portion>0:
        threshold=sorted(tmp[len(tmp)/portion], reverse=True)
        tmp_genes=[sublist[4] for chr in genes.keys() for sublist in genes[chr] if sublist[3] > threshold] 
    elif portion<0:
        threshold=sorted(tmp[len(tmp)/(-1*portion)], reverse=False)
        tmp_genes=[sublist[4] for chr in genes.keys() for sublist in genes[chr] if sublist[3] < threshold] 
        
    with open (os.path.join(path , "k562_" + mode + "_all_genes" + "_" + df_name + ".tab"), 'r') as inLines:
        for i, l in enumerate(inLines):
            ln=l.rstrip().split("\t")
            if ln[10] in tmp_genes:
                result.append(l)
    return result
    
#
#
def make_intragene_calculation (mode,  title, df_name, fl_in, fl_out_intra, fl_out_intra_by_gene, 
                                fl_out_raw_intra, statistics , chrms_needed, table_file_path):
    print fl_in
    
    statistics.write(fl_in + "\n")
    paths_pics="pics/"
    values_intra=open(fl_out_intra, 'w+')
    values_intra_by_gene=open(fl_out_intra_by_gene, 'w+')
    values_raw_intra_by_gene=open(fl_out_raw_intra, 'w+')
    
    table_file = open (table_file_path, 'w+')
    
    if mode=="RNA":
        bin_1_pos=1
        chr_1_pos=0
        bin_2_pos=12        
        chr_2_pos=11
        id_pos=10
        
    elif mode=="DNA":
        bin_1_pos=2
        chr_1_pos=1
        bin_2_pos=4        
        chr_2_pos=3
        id_pos=8
    else:
        return
    
    genes_lengths_sorted=sorted([int(x[3]) for x in genes_length.values()], reverse=True)
    #print genes_lengths_sorted[0]
    #print genes_lengths_sorted[-1]
    print "Max gene length " ,
    print max_gene_length
    print "Min gene length " ,
    print min_gene_length    
    print "# of genes " ,
    print len(genes_length)
    x_ax=range(-1* binId(genes_lengths_sorted[0]) -1 , binId(genes_lengths_sorted[0])  + 2 , 1)
    #x_ax=range(0, binId(genes_lengths_sorted[0])  + 1 , 1)
    #print x_ax
    dist_intra={}
    dist_by_bins_intra={}
    bins_intra={}
   
    if mode=="RNA":
        fl_raw=open(fl_in, 'r')
    elif mode=="DNA":
        fl_raw=open(fl_in, 'r')
        
    count=0
    norm_value={}
    
    count_empty_line=0
    count_error_line=0
    count_to_check=0
    norm_coef=0
    count_change_sign=0
    count_pos_strand=0
    count_neg_strand=0
    
    d_range=range(-2500, 2500)
    
    for i, l in enumerate(fl_raw): 
        if l in ['\n', '\r\n']:
            count_empty_line+=1
            continue
        if "Warning" in l:
            count_error_line+=1
            continue
        ln=l.rstrip().split("\t")
        if (mode=='RNA' and len(ln)<15) or (mode=='DNA' and len(ln)<8 ):
            print "Warning: make_intrachromosomal_calculation short line "  +str(ln)
            continue
        if ln[chr_2_pos] not in chrms_needed:
            continue
        
        id=ln[id_pos]
        
        if id not in genes_length:
            if "_" in id:
                tmp_id=id.split("_")[0]
                if tmp_id in genes_length:
                    id=tmp_id
                else:
                    print id
                    continue
                
        if  not (int(ln[bin_1_pos])>=int(genes_length[id][0])  and int(ln[bin_1_pos])<=int(genes_length[id][1])) :
            continue                
        if  not (int(ln[bin_2_pos])>=int(genes_length[id][0])  and int(ln[bin_2_pos])<=int(genes_length[id][1])) :
            continue                
        
        binRNA=binId(int(ln[bin_1_pos]))
        binDNA=binId(int(ln[bin_2_pos]))      
        
        d=binDNA-binRNA
        
        #if d==0 or (d ) > 100 or (d) < -100:
        if  not (d in d_range) :
        #if d==0 or not (d in d_range) :
            norm_coef+=1
            continue
        if mode=='DNA' and binDNA<binRNA:
            count_to_check+=1
           
        chr1=ln[chr_1_pos]
        chr2=ln[chr_2_pos]
        
        
        
        #d=abs(binDNA-binRNA)
        
        strand=genes_length[id][2]
        
        
        if genes_length[id][2]=="-":
            count_neg_strand+=1
            #continue
            
        if genes_length[id][2]=="+":
            count_pos_strand+=1
            #continue
            
        if mode=='DNA' and strand=='-':
            strand='+'
            
        norm_coef+=1
        
        if id in norm_value:
            norm_value[id]=norm_value[id]+1
        else:
            norm_value[id]=1
            
        if strand=="-":    
        #if mode=='RNA' and genes_length[id][2]=="-":
        #if mode!='DNA' and genes_length[id][2]=="-":
            d = (-1) * d
            count_change_sign+=1
        
        '''
        #Changed
        if mode=='DNA' and genes_length[id][2]=="-":
            tmp_value=binRNA
            binRNA=binDNA
            binDNA=tmp_value
            #d = (-1) * d
            #count_change_sign+=1
        '''        
        if id not in dist_intra:
            dist_intra[id]=[{}, None]            
        
        geneLengthInBins = max(binId(int(genes_length[id][3])), binId(genes_length[id][1]) - binId(genes_length[id][0]))
        
        if strand == "+":
            relativeBinRNA = binRNA - binId(genes_length[id][0])
            relativeBinDNA = binDNA - binId(genes_length[id][0])
            bn_range=range (max(-1 * relativeBinRNA, d_range[0]), min(geneLengthInBins  - relativeBinRNA + 1, d_range[-1] + 1))
          
        elif strand == "-":
            relativeBinRNA = binId(genes_length[id][1]) - binRNA
            relativeBinDNA = binId(genes_length[id][1]) - binDNA
            bn_range=range (max(-1 * relativeBinRNA, d_range[0]), min(geneLengthInBins - relativeBinRNA + 1, d_range[-1] + 1))
            
            if relativeBinRNA > geneLengthInBins:
                print id ,
                print relativeBinRNA ,
                print binRNA ,
                print binId(genes_length[id][0]) ,
                print binId(genes_length[id][1]) ,
                print binId(int(genes_length[id][3])) ,
                print geneLengthInBins
                print "\n"
                
        else:
            print "\nError: Unknown strand, exiting ..\n"
            break
        
        if relativeBinRNA not in dist_intra[id][0]:
            dist_intra[id][0][relativeBinRNA]={}
            dist_intra[id][1]=bn_range

            for t in d_range:
                if t in bn_range:
                    dist_intra[id][0][relativeBinRNA][t]=0
                #else:
                #    dist_intra[id][0][relativeBinRNA][t]=float('nan')
                
  
        if d not in bn_range:
            print id
            print genes_length[id]
            print geneLengthInBins
            print relativeBinRNA            
            print d
            print bn_range            
            print l
            
        if np.isnan(dist_intra[id][0][relativeBinRNA][d]) or (d not in d_range ):# is None:           
            print "Error with d and d_range interevals ! exiting .."
            print d
            print dist_intra[id][0][relativeBinRNA]
            print binId(int(genes_length[id][3])) + 2
            print d_range[-1]
            break            
        else:
            dist_intra[id][0][relativeBinRNA][d]=dist_intra[id][0][relativeBinRNA][d] + 1             
        
    fl_raw.close() 
    
    y_ax_intra=[]      
    number_of_bins=0
    number_of_non_zero_pos=0   
    number_of_non_one_pos=0 
    signal_by_d={}    
    norm_value_by_d={}
    empty_line=0
    number_of_empty_lines=0    
    
    for d in range (-2600, 2602):
        signal_by_d[d]=0
        norm_value_by_d[d]=0
        
     
    for id in sorted(dist_intra.keys()):
        number_of_bins+=len(dist_intra[id][0])        
        for bin in dist_intra[id][0]:            
            table_file.write(id + "," 
                             + str(max (binId(int(genes_length[id][3])), 
                                        binId(genes_length[id][1]) - binId(genes_length[id][0]))) 
                             + "," + str(bin) + "," + genes_length[id][2] 
                             + "," + str(bin_size))            
            #for d in dist_intra[id][0][bin]:     
            for d in d_range:   
                if d in dist_intra[id][0][bin]:
                    table_file.write("\t" + str(dist_intra[id][0][bin][d]))
                    if np.isnan(dist_intra[id][0][bin][d]) or dist_intra[id][0][bin][d]==0:
                        pass
                    else:
                        empty_line=1
                    if not np.isnan(dist_intra[id][0][bin][d]) :
                        if dist_intra[id][0][bin][d]!=0:
                            number_of_non_zero_pos+=1  
                            if dist_intra[id][0][bin][d]==1:
                                number_of_non_one_pos+=1                            
                        signal_by_d[d]+=dist_intra[id][0][bin][d]
                        norm_value_by_d[d]+=1  
                else:
                    table_file.write("\t" + "nan")
                    
                                 
            table_file.write ("\n")
            if empty_line==0:
                number_of_empty_lines+=1
            empty_line=0
    
    '''
    for d in range(-150, 151):
        print signal_by_d[d] ,
    print "\n"
    
    for d in range(-150, 151):
        print norm_value_by_d[d] ,
    '''
    
    print "\nNumber of ids :" ,
    print len(dist_intra.keys()) 
    
    print "number_of_bins :" ,
    print number_of_bins 
    
    print "number_of_non_zero_pos" ,
    print number_of_non_zero_pos     
    
    print "number_of_non_one_pos" ,
    print number_of_non_one_pos       
    
    print "number_of_empty_lines" ,
    print number_of_empty_lines    
    
    print "count_change_sign = ",
    print count_change_sign
    
    print "count_pos_strand = ",
    print count_pos_strand
    
    print "count_neg_strand = ",
    print count_neg_strand
    
    print "number of lines in the file = ",
    print i
    
    print "count_to_check = ",
    print count_to_check
    
    print "count_empty_line = ",
    print count_empty_line
    statistics.write(str(count_empty_line) + "\n")
    
    print "count_error_line = ",
    print count_error_line
    statistics.write(str(count_error_line) + "\n")
    
    
    count_protruding_x=0
    #norm_coef=1
    for x in x_ax:       
        y_ax_intra_tmp=[]
        
        for gene in dist_intra.keys():  
            if (x > 0 and x > binId(int(genes_length[gene][3]))+1) or (x < 0 and x < -1 * binId(int(genes_length[gene][3]))-1):
                continue
            else:
                count_protruding_x+=1
            for bn in dist_intra[gene][0]:
                if x in dist_intra[gene][0][bn]: 
                    y_ax_intra_tmp.append(dist_intra[gene][0][bn][x])               
                
                else:
                    y_ax_intra_tmp.append(float('nan'))
                    #pass
        '''
        if np.isnan(y_ax_intra_tmp).all():
            print " x: all nans" ,
            print x
        if len(y_ax_intra_tmp)==0:
            print " x: len y_tmp==0" ,
            print x            
        '''
        if len(y_ax_intra_tmp)!=0 and (not np.isnan(y_ax_intra_tmp).all()) and np.nanmean(y_ax_intra_tmp)!=0:
            y_ax_intra.append(np.nanmean(y_ax_intra_tmp)/(norm_coef * 1.0))
            values_intra.write(str(x) + "\t" + str(np.nanmean(y_ax_intra_tmp)/(norm_coef * 1.0)) + "\t" + str(norm_coef) + "\t")
            for y in y_ax_intra_tmp:
                values_intra.write(str(y) + ";")
            values_intra.write("\n")
        else:
            y_ax_intra.append(float('nan'))
            values_intra.write(str(x) + "\t" + str(float('nan')) + "\t" + str(norm_coef) + "\t;" + "\n")  
        
    x_ax_np=np.array(x_ax).astype(np.integer)
    y_ax_intra_np=np.array(y_ax_intra).astype(np.double)
    
    print "count_protruding_x = " ,
    print count_protruding_x ,
    print "\n"
    fig, ax = plt.subplots(figsize=(6, 6))    
    
    ax.plot(x_ax, y_ax_intra, linewidth=1, label='Intragenes ' )
    #ax.plot(x_ax_np[y_ax_inter_mask], y_ax_intra_np[y_ax_inter_mask], linewidth=1, label='Intrachromosomal, ' + chr)
    

    tmp_title=title + " " +  "_" + df_name
    
    ax.set_title(tmp_title, fontsize=14)    
    ax.legend(loc='lower center', prop={'size': 6})
    
    #plt.savefig(path + paths_pics + "chrm_territory_" + "_" +df_name+"_"+mode+"_DNA_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    
    ax.set_yscale("log", nonposy='clip')
    #ax.set_ylim(0.0000001, 0.01)
    #plt.savefig(path + paths_pics + "chrm_territory_" + "_" +df_name+"_"+mode+"_DNA_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    
    #ax.set_xscale("log", nonposx='clip')
    #ax.set_xscale("symlog", linthreshy=1)
    ax.set_xscale("symlog")
    plt.savefig(path + paths_pics + "intragene"+"_" + df_name + "_" + mode + "_DNA_ylog_xlog" + title.split("intragene")[-1] + "_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
   
    
    plt.close()    
    values_intra.close()
    values_intra_by_gene.close()
    values_raw_intra_by_gene.close()
    table_file.close()
#
#
def make_intragene_calculation_relative_bin (mode,  title, df_name, fl_in, fl_out_intra, fl_out_intra_by_gene, fl_out_raw_intra, statistics , chrms_needed):
    
    print fl_in
    
    statistics.write(fl_in + "\n")
    paths_pics="pics/"
    values_intra=open(fl_out_intra, 'w+')
    values_intra_by_gene=open(fl_out_intra_by_gene, 'w+')
    values_raw_intra_by_gene=open(fl_out_raw_intra, 'w+')

    if mode=="RNA":
        bin_1_pos=1
        chr_1_pos=0
        bin_2_pos=12        
        chr_2_pos=11
        id_pos=10
        
    elif mode=="DNA":
        bin_1_pos=2
        chr_1_pos=1
        bin_2_pos=4        
        chr_2_pos=3
        id_pos=8
    else:
        return
    
    genes_lengths_sorted=sorted([int(x[3]) for x in genes_length.values()], reverse=True)

    print "Max gene length " ,
    print max_gene_length
    print "Min gene length " ,
    print min_gene_length    
    print "# of genes " ,
    print len(genes_length)
    max_dist=7
    dist_range=range(-7, 8)
    
    x_ax=dist_range
    #x_ax=range(-1* binId(genes_lengths_sorted[0]) -1 , binId(genes_lengths_sorted[0])  + 2 , 1)
    #x_ax=range(0, binId(genes_lengths_sorted[0])  + 1 , 1)
    
    dist_intra={}
    dist_by_bins_intra={}
    bins_intra={}
   
    if mode=="RNA":
        fl_raw=open(fl_in, 'r')
    elif mode=="DNA":
        fl_raw=open(fl_in, 'r')
        
    count=0
    norm_value={}
    
    count_empty_line=0
    count_error_line=0
    count_to_check=0
    norm_coef=0
    count_change_sign=0
    count_pos_strand=0
    count_neg_strand=0
    
    for i, l in enumerate(fl_raw): 
        if l in ['\n', '\r\n']:
            count_empty_line+=1
            continue
        if "Warning" in l:
            count_error_line+=1
            continue
        ln=l.rstrip().split("\t")
        if (mode=='RNA' and len(ln)<15) or (mode=='DNA' and len(ln)<8 ):
            print "Warning: make_intrachromosomal_calculation short line "  +str(ln)
            continue
        if ln[chr_2_pos] not in chrms_needed:
            continue
        
        id=ln[id_pos]
        
        if id=="STPG4_3276" or id=="SDHD_3854":
            continue

        
        if id not in genes_length:
            if "_" in id:
                tmp_id=id.split("_")[0]
                if tmp_id in genes_length:
                    id=tmp_id
                else:
                    print id
                    continue
        
        global bin_size
        bin_size=1+(genes_length[id][3])/6
        
        
        binStart=binId(int(genes_length[id][0]))
            
        binRNA=binId(int(ln[bin_1_pos]) - int(genes_length[id][0]))
        binDNA=binId(int(ln[bin_2_pos]) - int(genes_length[id][0]))      
        
        if genes_length[id][2]=="-":
            binRNA = 5 - binRNA
            binDNA = 5 - binDNA
        
        d=binDNA-binRNA
        
        if mode=='DNA' and binDNA<binRNA:
            count_to_check+=1
           
        chr1=ln[chr_1_pos]
        chr2=ln[chr_2_pos]
        
        norm_coef+=1
         
        if genes_length[id][2]=="-":
            count_neg_strand+=1
            
        if genes_length[id][2]=="+":
            count_pos_strand+=1
        
        '''   
        if genes_length[id][2]=="-":
            d = (-1) * d
            count_change_sign+=1
        '''
        
        if abs(d) > max_dist or abs(binRNA) > 5 or abs(binDNA) > 5 or binRNA*binDNA<0 :
            print bin_size
            print binRNA
            print binDNA
            print ln[bin_1_pos]
            print ln[bin_2_pos]
            print l
            print genes_length[id]
            print ""
            norm_coef-=1
            continue
        
        if id not in dist_intra:
            dist_intra[id]=[{}, {}]    
        
            for bn in range(0,6):  
                dist_intra[id][0][bn]={}
                for t in dist_range:
                    dist_intra[id][0][bn][t]=0   
                    
        if binRNA not in  dist_intra[id][0] or d not in dist_intra[id][0][binRNA]:
            print "binRNA " ,
            print binRNA
            print "d " ,
            print d
            print l
            
        dist_intra[id][0][binRNA][d]=dist_intra[id][0][binRNA][d] + 1
        
        '''
        if binRNA not in dist_intra[id][1][d]:
            dist_intra[id][1][d].append(binRNA)
        '''

    fl_raw.close()          
    y_ax_intra=[] 
    
    print "count_change_sign = ",
    print count_change_sign
    
    print "count_pos_strand = ",
    print count_pos_strand
    
    print "count_neg_strand = ",
    print count_neg_strand
    
    print "number of lines in the file = ",
    print i
    
    print "count_to_check = ",
    print count_to_check
    
    print "count_empty_line = ",
    print count_empty_line
    statistics.write(str(count_empty_line) + "\n")
    
    print "count_error_line = ",
    print count_error_line
    statistics.write(str(count_error_line) + "\n")
    
    
    count_protruding_x=0
    
    
    y_ax_sliding={}
    y_ax_fixed={}
    y_ax_sum_fixed=[]
    y_ax_sum_sliding=[]
    
    
    def total_signal_for_gene(gene_data):
        tmp=0
        for bn in gene_data:
            for x in gene_data[bn]:
                tmp+=gene_data[bn][x]
        return float(tmp)
            
        
    for bn in range(0, 6):
        y_ax_sliding[bn]=[]
        y_ax_fixed[bn]=[]
        
    for x in x_ax:  
        sum_fixed=[]
        #for bn in [0, 5]:
        for bn in range(0, 6):
            y_ax_intra_tmp=[]
            for gene in dist_intra.keys():  
                if (x > 0 and x > binId(int(genes_length[gene][3])) + 1) or (x < 0 and x < -1 * binId(int(genes_length[gene][3])) - 1):
                    count_protruding_x+=1
                    continue            
                
                if x in dist_intra[gene][0][bn]:
                    if total_signal_for_gene(dist_intra[gene][0])==0:
                        print dist_intra[gene][0]
                    #y_ax_intra_tmp.append(dist_intra[gene][0][bn][x])   
                    y_ax_intra_tmp.append(dist_intra[gene][0][bn][x]/total_signal_for_gene(dist_intra[gene][0]))
            
            #y_ax_fixed[bn].append(np.sum(y_ax_intra_tmp)/(norm_coef * 1.0))
            if np.mean(y_ax_intra_tmp)/(1.0)!=0:
                y_ax_fixed[bn].append(np.mean(y_ax_intra_tmp)/(1.0))
                
            else:
                y_ax_fixed[bn].append(float('nan'))
            sum_fixed.extend(y_ax_intra_tmp)
            
        y_ax_sum_fixed.append(np.mean(sum_fixed))    
        
    for bn in range(0, 6):
    #for bn in [0, 5]:
        for j in range(len(y_ax_fixed[bn])):
            if (j - bn) < -5:
                y_ax_sliding[bn].append( float('nan'))
                continue

            y_ax_sliding[bn].append( y_ax_fixed[bn][j - bn])
        
        
    x_ax_np=np.array(x_ax).astype(np.integer)
    y_ax_intra_np=np.array(y_ax_intra).astype(np.double)
    
    print "count_protruding_x = " ,
    print count_protruding_x ,
    print "\n"
    
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 12))    
    
    print "y_ax_fixed[0]"
    print y_ax_fixed[0]
    print "y_ax_sliding[0]"
    print y_ax_sliding[0]
    
    print "y_ax_fixed[1]"
    print y_ax_fixed[1]
    print "y_ax_sliding[1]"
    print y_ax_sliding[1]
    
    print "y_ax_fixed[2]"
    print y_ax_fixed[2]
    print "y_ax_sliding[2]"
    print y_ax_sliding[2]
    
    for bn in range(0, 6):   
    #for bn in [0, 5]:
        ax[0,0].plot(x_ax, y_ax_fixed[bn], linewidth=1.5, label='fixed bin = ' + str(bn) )
        ax[0,0].set_yscale("log", nonposy='clip')
    ax[0,1].plot(x_ax, y_ax_sum_fixed, linewidth=1.5, label='sum fixed bn' )
    ax[0,1].set_yscale("log", nonposy='clip') 
    for bn in range(0, 6):   
    #for bn in [0, 5]:        
        ax[1,0].plot(x_ax, y_ax_sliding[bn], linewidth=1.5, label='sliding bin = ' + str(bn) )
        ax[1,0].set_yscale("log", nonposy='clip')
    #ax[2,2].plot(x_ax, y_ax_sum_fixed, linewidth=1, label='sum fixed bn' )

    plt.savefig(path + paths_pics + "intragene"+"_" + df_name + "_" + mode + "_DNA_ylog_xlog" + title.split("intragene")[-1] + "_rel_bin.png", dpi=300, figsize=(12, 12)) 
       
    plt.close()    
    values_intra.close()
    values_intra_by_gene.close()
    values_raw_intra_by_gene.close()
#
#
'''

for gene in dist_by_bins_intra.keys():            
            if x in dist_by_bins_intra[gene]:   
                values_raw_intra_by_gene.write(str(x) + "\t" + str(gene) + "\t")                
                for bn in dist_by_bins_intra[gene][x]:
                    y_ax_intra_tmp.append(dist_by_bins_intra[gene][x][bn])
                    values_raw_intra_by_gene.write(str(dist_by_bins_intra[gene][x][bn]) + ";")                
                values_raw_intra_by_gene.write("\t" + str(norm_value[gene]) + "\n") 
                
            else:
                values_raw_intra_by_gene.write(str(x) + "\t" + str(gene) + "\t" + str(0) + ";\t" + str(norm_value[gene]) + "\n")       
                
'''


'''
            if (length + 1 - x)==0:
                    print gene ,
                    print "x=",
                    print x ,
                    print "dist_intra[gene][x]=",
                    print dist_intra[gene][x] ,
                    print "length=",
                    print length
                    print "norm_value[gene]=",
                    print norm_value[gene] , 
                    print "binId(genes_length[gene][3])=",
                    print binId(genes_length[gene][3]) ,
                    print 'genes_length[gene]',
                    print genes_length[gene]   
                    print " "
'''
#
#
def make_combined_graph_genes(fl1_intra_in, fl2_intra_in, fl3_intra_in, fl4_intra_in, name):
    path_pics="pics/"
            
    fl1=open(fl1_intra_in, 'r')
    fl2=open(fl2_intra_in, 'r')
    fl3=open(fl3_intra_in, 'r')
    fl4=open(fl4_intra_in, 'r')
    
    fl1=list(enumerate(fl1))
    fl2=list(enumerate(fl2))
    fl3=list(enumerate(fl3))
    fl4=list(enumerate(fl4))
    
    
    x_ax=[]
    y_ax1=[]
    y_ax2=[]
    y_ax3=[]
    y_ax4=[]
    
    y_ax1_ratio=[] 
    y_ax2_ratio=[]  
    y_ax3_ratio=[]  
    y_ax4_ratio=[]
    
    
    for i in range(len(fl1)):
        x_ax.append(int(fl1[i][1].rstrip().split("\t")[0]))
        
        y_ax1.append(float(fl1[i][1].rstrip().split("\t")[1])) 
        y_ax2.append(float(fl2[i][1].rstrip().split("\t")[1])) 
        y_ax3.append(float(fl3[i][1].rstrip().split("\t")[1])) 
        y_ax4.append(float(fl4[i][1].rstrip().split("\t")[1])) 
        
         
    
    fig, ax = plt.subplots(figsize=(6, 6))    
    
    l=(len(x_ax) -1)/2 
    
    ax.plot(x_ax, y_ax1,  linewidth=0.7, color='blue', linestyle=":", label='Intragene, ' + fl1_intra_in.split("_values.tab")[0])    
    ax.plot(x_ax[(l):(l+2)], y_ax1[(l):(l+2)], linewidth=0.7, color='blue', linestyle=":", label='_nolegend_')    
    
    ax.plot(x_ax, y_ax2,  linewidth=0.7, color='orange', linestyle="-.", label='Intragene, ' + fl2_intra_in.split("_values.tab")[0])    
    ax.plot(x_ax[(l):(l+2)], y_ax2[(l):(l+2)], linewidth=0.7, color='orange', linestyle="-.", label='_nolegend_')
    
    ax.plot(x_ax, y_ax3,  linewidth=0.7, color='yellow', linestyle="--", label='Intragene, ' + fl3_intra_in.split("_values.tab")[0])
    ax.plot(x_ax[(l-2):(l+2)], y_ax3[(l-2):(l+2)], linewidth=0.7, color='yellow', linestyle="--", label='_nolegend_')
    
    ax.plot(x_ax, y_ax4,  linewidth=0.7, color='red', linestyle="-", label='Intragene, ' + fl4_intra_in.split("_values.tab")[0])
    ax.plot(x_ax[(l-2):(l+2)], y_ax4[(l-2):(l+2)], linewidth=0.7, color='red', linestyle="-", label='_nolegend_')
    
    
    ax.legend(loc='lower center', prop={'size': 6}, frameon=False)
    
    #plt.savefig(path_pics+"chrm_territory_"+"_combined_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6))
    ax.set_yscale("log")
    #ax.set_ylim(0.000000001, 0.01)
    #plt.savefig(path_pics+"chrm_territory_"+"_combined_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    #ax.set_xscale("symlog", nonposx='clip', linthreshy=1)
    ax.set_xscale("symlog")
    plt.savefig(path_pics + "intragene" + "_combined_ylog_xlog" + name.split("intragene")[-1] + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close()       

    
    x_ax_ratio = range(1, x_ax[-1])  
    #print x_ax_ratio
    for i in x_ax_ratio:
        '''
        print i , 
        print y_ax2[i+l] ,
        print y_ax2[l-i] ,
        print y_ax1[i+l] ,
        print y_ax1[l-i] ,
        print l ,
        print i+l ,
        print l-i
        '''
        if not np.isnan(y_ax1[i+l]) and not np.isnan(y_ax1[l-i]):
            y_ax1_ratio.append(y_ax1[i+l]/y_ax1[l-i])
        else:
            y_ax1_ratio.append(float('nan'))
            
        if not np.isnan(y_ax2[i+l]) and not np.isnan(y_ax2[l-i]):
            y_ax2_ratio.append(y_ax2[i+l]/y_ax2[l-i])                
        else:
            y_ax2_ratio.append(float('nan'))        
        
        if not np.isnan(y_ax3[i+l]) and not np.isnan(y_ax3[l-i]):
            y_ax3_ratio.append(y_ax3[i+l]/y_ax3[l-i])   
        else:
            y_ax3_ratio.append(float('nan'))            
        
        if not np.isnan(y_ax4[i+l]) and not np.isnan(y_ax4[l-i]):
            y_ax4_ratio.append(y_ax4[i+l]/y_ax4[l-i])     
        else:
            y_ax4_ratio.append(float('nan'))            
            
    
    fig2, ax2 = plt.subplots(figsize=(6, 6))    
    #plt.axes().yaxis.set_minor_locator(ticker.AutoMinorLocator())
    #plt.axes().xaxis.set_minor_locator(ticker.AutoMinorLocator())
    #ax2.tick_params(axis='x',which='minor',bottom=False)
    #ax2.tick_params(axis='y',which='minor',bottom=False)
    
    #ax2.axis([1, 10, 1, 1000])

    
    ax2.plot(x_ax_ratio, y_ax1_ratio,  
            linewidth=0.5, 
            color='blue', 
            linestyle=":", 
            label='Intragene, ' + fl1_intra_in.split("_values")[0]) 
    
    print x_ax_ratio[0:4]
    
    print y_ax1_ratio[0:4]
    '''
    ax2.plot(x_ax_ratio[0:4], y_ax1_ratio[0:4] ,  'ro-', 
            linewidth=0.5, 
            color='blue')
            #linestyle="ro-", 
            #label='Intragene, ' + fl1_intra_in.split("_values")[0]) 
    '''
    
    ax2.plot(x_ax_ratio, y_ax2_ratio, 
            linewidth=0.5, 
            color='orange', 
            linestyle="-.", 
            label='Intragene, ' + fl2_intra_in.split("_values")[0])    
        
    ax2.plot(x_ax_ratio, y_ax3_ratio, 
            linewidth=0.5, 
            color='yellow', 
            linestyle="--", 
            label='Intragene, ' + fl3_intra_in.split("_values")[0])
    
    ax2.plot(x_ax_ratio, y_ax4_ratio, 
            linewidth=0.5, 
            color='red', 
            linestyle="-", 
            label='Intragene, ' + fl4_intra_in.split("_values")[0])    
    
    ax2.legend(loc='upper center', prop={'size': 6}, frameon=False)
    
    #plt.savefig(path_pics+"chrm_territory_"+"_combined_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6))
    #ax2.set_yscale("symlog")
    #ax.set_ylim(1, 1e+1)
    yaxis = plt.gca().yaxis
    #yaxis.set_minor_locator(MinorSymLogLocator(1e-1))    
    yaxis.set_minor_locator(ticker.AutoMinorLocator())    
    ax2.set_xscale("symlog")
    xaxis = plt.gca().xaxis
    xaxis.set_minor_locator(MinorSymLogLocator(1e-1))
    
    #plt.savefig(path_pics+"chrm_territory_"+"_combined_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    #ax.set_xscale("symlog", nonposx='clip', linthreshy=1)
    

    plt.savefig(path_pics + "intragene" + "_combined_ratio" + name.split("intragene")[-1] + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close()        
        
#
def make_combined_graph_genes_wo_last_exon( fl2_intra_in, fl4_intra_in, name):
    path_pics="pics/" 
    
    print fl4_intra_in
    print fl2_intra_in
    fl2=open(fl2_intra_in, 'r')
    fl4=open(fl4_intra_in, 'r')  
    
    fl2=list(enumerate(fl2))
    fl4=list(enumerate(fl4))   

    x_ax=[]
    y_ax2=[]
    y_ax4=[]    
    
    def return_conf_interval(vals): 
        vals=np.array(vals).astype(np.double)
        if len(vals)==0 or np.count_nonzero(~np.isnan(vals))==0:
            return [None, None]
        std=np.nanstd(vals, ddof=1)        
        mn=np.nanmean(vals)        
        length=np.count_nonzero(~np.isnan(vals))
        return [mn - 1.96 * (std/np.sqrt(length)) , mn + 1.96 * (std/np.sqrt(length))]
    
    def return_bootstrap_values_for_mean(vals, number_of_iterations):
        mean_estimates=[]
        vals=np.array(vals).astype(np.double)
        for _ in range(number_of_iterations): 
            if len(vals)==0 or np.count_nonzero(~np.isnan(vals)) == 0:
                mean_estimates.append(None)
                continue
            re_sample_idx = np.random.randint(0, len(vals), len(vals))
            mean_estimates.append(np.nanmean(vals[re_sample_idx]))        
        return mean_estimates
    
    
    BS_i=200
    
    y_ax2_CI_high=[]    
    y_ax2_CI_low=[]    
    y_ax4_CI_high=[]
    y_ax4_CI_low=[]
    y_ax2_BS={}
    y_ax4_BS={}
    
    for i in range(BS_i):
        y_ax2_BS[i]=[]
        y_ax4_BS[i]=[]
    
    for i in range(len(fl2)):
        x=int(fl2[i][1].rstrip().split("\t")[0])
        
        x_ax.append(int(fl2[i][1].rstrip().split("\t")[0]))
        
        y_ax2.append(float(fl2[i][1].rstrip().split("\t")[1])) 
        y_ax4.append(float(fl4[i][1].rstrip().split("\t")[1])) 
        
        #---------------------------------------------------------------            
        norm_value=float(fl2[i][1].rstrip().split("\t")[2])
        flat_counts=[]        

        
        for val in fl2[i][1].rstrip().split("\t")[3].split(";"):
            try:
                if val!='':                           
                    flat_counts.append(float(val)/(1.0 * norm_value))
            except ValueError:
                print val            
        if len(flat_counts)>0:
            vals = np.array(flat_counts)            
            #print "CI and bootstrapping DNA-DNA ... "
            CI=return_conf_interval(vals)            
            y_ax2_CI_high.append(CI[1])            
            y_ax2_CI_low.append(CI[0])
            bs_vals=return_bootstrap_values_for_mean(flat_counts, BS_i)
            for j in range(BS_i):            
                y_ax2_BS[j].append(bs_vals[j])
        else:

            y_ax2_CI_high.append(None)            
            y_ax2_CI_low.append(None)
            for j in range(BS_i):            
                y_ax2_BS[j].append(None) 
                
        #---------------------------------------------------------------            
        norm_value=float(fl4[i][1].rstrip().split("\t")[2])
        flat_counts=[]
        
        
        for val in fl4[i][1].rstrip().split("\t")[3].split(";"):
            try:
                if val!='':                           
                    flat_counts.append(float(val)/(1.0 * norm_value))
            except ValueError:
                print val                    
        if len(flat_counts)>0:
            vals = np.array(flat_counts)
            #print "CI and bootstrapping RNA-DNA ... "
            CI=return_conf_interval(vals)
            y_ax4_CI_high.append(CI[1])            
            y_ax4_CI_low.append(CI[0])
            bs_vals=return_bootstrap_values_for_mean(flat_counts, BS_i)
            for j in range(BS_i):            
                y_ax4_BS[j].append(bs_vals[j])
        else:
            
            y_ax4_CI_high.append(None)            
            y_ax4_CI_low.append(None)
            for j in range(BS_i):            
                y_ax4_BS[j].append(None)        
        
        
        '''
        if int(fl2[i][1].rstrip().split("\t")[0]) ==0 :
            print float(fl2[i][1].rstrip().split("\t")[1])
            print float(fl4[i][1].rstrip().split("\t")[1])
        '''
    
    fig, ax = plt.subplots(figsize=(6, 6))        
        
    ax.plot(x_ax, y_ax2,  linewidth=0.5, color='orange', label='Intragene, ' + fl2_intra_in.split("_values.tab")[0])    
    ax.plot(x_ax, y_ax4,  linewidth=0.5, color='red', label='Intragene, ' + fl4_intra_in.split("_values.tab")[0])
    
    l=(len(x_ax) -1)/2 
    ax.plot(x_ax[(l):(l+2)], y_ax2[(l):(l+2)], linewidth=0.5, color='orange', label='_nolegend_')
    ax.plot(x_ax[(l-2):(l+2)], y_ax4[(l-2):(l+2)], linewidth=0.5, color='red', label='_nolegend_')
        
    ax.legend(loc='lower center', prop={'size': 6}, frameon=False)
    
    #plt.savefig(path_pics+"chrm_territory_"+"_combined_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6))
    #ax.set_yscale("log", nonposy='clip')
    ax.set_yscale("log")
    #ax.set_ylim(0.000000001, 0.01)
    #plt.savefig(path_pics+"chrm_territory_"+"_combined_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    ax.set_xscale("symlog", linthreshy=0.1)
    plt.savefig(path_pics + "intragene" + "_combined_ylog_xlog" + name.split("intragene")[-1] + "_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close()   
    
    
    #-------------------------------------------------------BS-------------------------------------------------        
    l=(len(x_ax) -1)/2 
    fig, ax = plt.subplots(figsize=(6, 6)) 
    for i in range(BS_i):                
        ax.plot(x_ax, y_ax4_BS[i], linewidth=0.1, alpha=0.2, color='grey', linestyle="-")
        ax.plot(x_ax, y_ax2_BS[i], linewidth=0.1, alpha=0.2, color='grey', linestyle="-")
        #ax.plot(x_ax[(l-6):(l+6)], y_ax2_BS[i][(l-6):(l+6)], linewidth=0.1, alpha=0.4, color='green', linestyle="-", label='_nolegend_')
        
    ax.plot(x_ax, y_ax2,  linewidth=0.5, color='orange', label='Intragene, ' + fl2_intra_in.split("_values.tab")[0])    
    ax.plot(x_ax, y_ax4,  linewidth=0.5, color='red', label='Intragene, ' + fl4_intra_in.split("_values.tab")[0])
    
    
    #ax.plot(x_ax[(l):(l+2)], y_ax2[(l):(l+2)], linewidth=0.5, color='orange', label='_nolegend_')
    #ax.plot(x_ax[(l-2):(l+2)], y_ax4[(l-2):(l+2)], linewidth=0.5, color='red', label='_nolegend_')
    
    
    ax.legend(loc='lower center', prop={'size': 6}, frameon=False)
    ax.set_yscale("log")
    ax.set_xscale("symlog", linthreshy=0.1)
    plt.savefig(path_pics + "intragene" + "_combined_ylog_xlog" + name.split("intragene")[-1] + "_BS_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close()    
    
    #-------------------------------------------------------CI-------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 6))
    
    
    x_ax_np=np.array(x_ax).astype(np.integer)
    
    y_ax4_CI_low_np=np.array(y_ax4_CI_low).astype(np.double)
    y_ax4_CI_high_np=np.array(y_ax4_CI_high).astype(np.double)
    y_ax4_mask=np.isfinite(np.log(y_ax4_CI_high_np)) &  np.isfinite(np.log(y_ax4))
    
    y_ax4_CI_low_masked=np.ma.masked_where( ~y_ax4_mask, y_ax4_CI_low_np)
    y_ax4_CI_high_masked=np.ma.masked_where(~y_ax4_mask, y_ax4_CI_high_np)    
    
    
    print np.log(y_ax4)[0:100]
    print y_ax4_mask[0:100]
    print y_ax4_CI_low_np[0:100]
    
    '''
    print y_ax4[len(x_ax)/2 -100: len(x_ax)/2 + 100]
    print y_ax4_CI_low[len(x_ax)/2 -100: len(x_ax)/2 + 100]
    print y_ax4_CI_low_np[len(x_ax)/2 -100: len(x_ax)/2 + 100]
    print y_ax4_mask[len(x_ax)/2 -100: len(x_ax)/2 + 100]
    print y_ax4_CI_low_masked[len(x_ax)/2 -100: len(x_ax)/2 + 100]
    '''
    
    ax.plot(x_ax_np, y_ax4_CI_low_masked, linewidth=0.1, alpha=0.4, color='grey')
    ax.plot(x_ax_np, y_ax4_CI_high_masked, linewidth=0.1, alpha=0.4, color='grey')
    ax.fill_between(x_ax, y_ax4_CI_high_masked, y_ax4_CI_low_masked, color='grey', alpha=0.4)
        
    y_ax2_CI_low_np=np.array(y_ax2_CI_low).astype(np.double)
    y_ax2_CI_high_np=np.array(y_ax2_CI_high).astype(np.double)
    y_ax2_mask=np.isfinite(y_ax2_CI_high_np)  
    
    y_ax2_CI_low_masked=np.ma.masked_where(~y_ax2_mask, y_ax2_CI_low_np)
    y_ax2_CI_high_masked=np.ma.masked_where(~y_ax2_mask, y_ax2_CI_high_np)  
    
    ax.plot(x_ax_np, y_ax2_CI_high_masked, linewidth=0.1, alpha=0.4, color='grey', linestyle="-")
    ax.plot(x_ax_np, y_ax2_CI_low_masked, linewidth=0.1, alpha=0.4, color='grey', linestyle="-")
    ax.fill_between(x_ax, y_ax2_CI_high_masked, y_ax2_CI_low_masked, color='grey', alpha=0.4)
    '''
    ax.plot(x_ax[(l-6):(l+6)], y_ax2_CI_high_masked[(l-6):(l+6)], linewidth=0.1, alpha=0.4, color='blue', linestyle="-")
    ax.plot(x_ax[(l-6):(l+6)], y_ax2_CI_low_masked[(l-6):(l+6)], linewidth=0.1, alpha=0.4, color='blue', linestyle="-")
    ax.fill_between(x_ax[(l-6):(l+6)], y_ax2_CI_high_masked[(l-6):(l+6)], y_ax2_CI_low_masked[(l-6):(l+6)], color='green', alpha=0.4)
    '''
    
    ax.plot(x_ax, y_ax2,  linewidth=0.5, color='orange', label='Intragene, ' + fl2_intra_in.split("_values.tab")[0])    
    ax.plot(x_ax, y_ax4,  linewidth=0.5, color='red', label='Intragene, ' + fl4_intra_in.split("_values.tab")[0])
    
    
    l=(len(x_ax)-1)/2 
    
    '''
    ax.plot(x_ax[(l):(l+2)], y_ax2[(l):(l+2)], linewidth=0.5, color='orange', label='_nolegend_')
    ax.plot(x_ax[(l-2):(l+2)], y_ax4[(l-2):(l+2)], linewidth=0.5, color='red', label='_nolegend_')
    '''
    
    ax.legend(loc='lower center', prop={'size': 6}, frameon=False)
    ax.set_yscale("log")
    ax.set_xscale("symlog", linthreshy=0.1)
    plt.savefig(path_pics + "intragene" + "_combined_ylog_xlog" + name.split("intragene")[-1] + "_CI_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close()   
    
#
#
def make_intrachromosomal_and_interchromosomal_calculation (mode, chr, title, df_name, fl_in, fl_out_intra, fl_out_raw_intra, fl_out_inter, fl_out_raw_inter, statistics , chrms_needed):
    print fl_in
    
    statistics.write(fl_in + "\n")
    paths_pics="pics/"
    values_intra=open(fl_out_intra, 'w+')
    values_raw_intra=open(fl_out_raw_intra, 'w+')
    values_inter=open(fl_out_inter, 'w+')
    value_raw_inter=open(fl_out_raw_inter, 'w+')
    
    '''
    if chr=='chr1':
        if mode=='RNA':
            norm_value=4271640
        if mode=='DNA':
            norm_value=118745965
    else:
        return
    '''
    if mode=="RNA":
        bin_1_pos=1
        bin_2_pos=12        
        chr_2_pos=11
        
    elif mode=="DNA":
        bin_1_pos=2
        bin_2_pos=4        
        chr_2_pos=3
    else:
        return
    
    x_ax=range(0, binId(chrms_length[chr])  + 1 ,1)
    dist_intra={}
    dist_by_bins_intra={}
    bins_intra={}
    dist_inter={}
    dist_by_bins_inter={}
    bins_inter={}
    
    if mode=="RNA":
        fl_raw=open(fl_in, 'r')
    elif mode=="DNA":
        #fl_raw=gzip.open(fl_in, 'r')
        fl_raw=open(fl_in, 'r')
        
    count=0
    norm_value=0
    
    count_empty_line=0
    count_error_line=0
    
    for i, l in enumerate(fl_raw): 
        if l in ['\n', '\r\n']:
            count_empty_line+=1
            continue
        if "Warning" in l:
            count_error_line+=1
            continue
        ln=l.rstrip().split("\t")
        if (mode=='RNA' and len(ln)<15) or (mode=='DNA' and len(ln)<8 ):
            print "Warning: make_intrachromosomal_calculation " + str(ln)
            continue
        if ln[chr_2_pos] not in chrms_needed:
            continue
        
        norm_value+=1
        binRNA=binId(int(ln[bin_1_pos]))
        binDNA=binId(int(ln[bin_2_pos]))   
                    
        if count >0:
            print ln
            print ln[bin_1_pos] ,
            print ln[bin_2_pos] ,
            print ln[chr_2_pos]
            count-=1
        
        #if mode=="DNA" and ln[chr_2_pos]!=chr and ln[chr_2_pos]!='chrX' and ln[chr_2_pos]!='chrY':
            #print l
            
        if ln[chr_2_pos] not in dist_inter: 
            dist_inter[ln[chr_2_pos]]=1
        else:
            dist_inter[ln[chr_2_pos]]+=1
        

          
        if ln[chr_2_pos]==chr:                
            d=abs(binDNA-binRNA)
            if d not in dist_intra:
                
                for t in range(0, binId(chrms_length[chr])  + 1):
                    dist_intra[d]=0
                    
                dist_intra[d]=1
                dist_by_bins_intra[d]={}
                dist_by_bins_intra[d][binRNA]=1
            else:
                dist_intra[d]+=1
                if binRNA not in dist_by_bins_intra[d]:
                    dist_by_bins_intra[d][binRNA]=1
                else:
                    dist_by_bins_intra[d][binRNA]=dist_by_bins_intra[d][binRNA] + 1

                    
    print sorted(dist_inter.keys())
    
    fl_raw.close()          
    y_ax_intra=[] 
    y_ax_inter_chr1=[] 
    y_ax_inter_chr1_A=[]
    y_ax_inter_chr10=[]
    y_ax_inter_chr21=[]
    
    num_bins_genome=0
    sum_for_chr=0
    
    for ch in chrms_length:
        if ch!=chr:
            if ch in dist_inter:
                sum_for_chr+=dist_inter[ch]
        num_bins_genome+=binId(chrms_length[ch]) + 1 
        
    for ch in sorted(chrms_length.keys()):
        #print ch , 
        #print binId(chrms_length[ch]) ,
        if ch in dist_inter:
            #print dist_inter[ch]
            value_raw_inter.write(ch +"\t" + str(binId(chrms_length[ch]) + 1) + "\t" + str(dist_inter[ch]) + "\t" + str(norm_value) +"\n")            
        else:
            #print "0"
            value_raw_inter.write(ch +"\t" + str(binId(chrms_length[ch]) + 1) + "\t" + str(0) + "\t" + str(norm_value) +"\n")
            
        
    print "norm_value (== non empty lines) = ",
    print  norm_value
    statistics.write(str(norm_value) + "\n")
    
    print "count_empty_line = ",
    print count_empty_line
    statistics.write(str(count_empty_line) + "\n")
    
    print "count_error_line = ",
    print count_error_line
    statistics.write(str(count_error_line) + "\n")
    
    value_for_chr10=0
    value_for_chr21=0
    
    if 'chr10' in dist_inter:
        value_for_chr10=dist_inter['chr10']
    else:
        value_for_chr10=0
        
    if 'chr21' in dist_inter:
        value_for_chr21=dist_inter['chr21']
    else:
        value_for_chr21=0    
        
    for x in x_ax:
        norm_bins=(num_bins_genome - binId(chrms_length[chr]) -1) * (binId(chrms_length[chr]) + 1)
        
        y_ax_inter_chr1_A.append((sum_for_chr * 1.0)/ (norm_bins * norm_value))
        values_inter.write(str(x) + "\t" + str((sum_for_chr * 1.0)/ (norm_bins * norm_value))  + "\n")
        
        y_ax_inter_chr10.append((value_for_chr10 * 1.0)/ ( (binId(chrms_length[chr]) + 1) * (binId(chrms_length['chr10']) + 1) * norm_value) )
        y_ax_inter_chr21.append((value_for_chr21 * 1.0)/ ( (binId(chrms_length[chr]) + 1) * (binId(chrms_length['chr21']) + 1) * norm_value) )                
        
        # 10, 21 old norm coefficients!
        #(num_bins_genome - binId(chrms_length[chr]) -1) * (binId(chrms_length[chr]) + 1)
        #y_ax_inter_chr10.append((value_for_chr10 * 1.0)/ ( (binId(chrms_length[chr]) + 1) * (binId(chrms_length['chr10']) + 1) * norm_value) )
        #y_ax_inter_chr21.append((value_for_chr21 * 1.0)/ ( (binId(chrms_length[chr]) + 1) * (binId(chrms_length['chr21']) + 1) * norm_value) )        
        #y_ax_inter_chr1.append((sum_for_chr+ dist_inter[chr]* 1.0)/ ( num_bins_genome * binId(chrms_length[chr]) *norm_value))
       
        if x in dist_intra:
            #sm=np.mean(dist_by_bins_intra[x].values())
            sm = (dist_intra[x] * 1.0) / (binId(chrms_length[chr]) + 1 - x)
            y_ax_intra.append((sm * 1.0) / (norm_value )  )
            values_intra.write(str(x) + "\t" + str((sm * 1.0) / (norm_value) ) + "\n")
            values_raw_intra.write(str(x) + "\t")
            
            for bn in dist_by_bins_intra[x]:
                values_raw_intra.write(str(dist_by_bins_intra[x][bn]) + ";")
            
            values_raw_intra.write("\t" +  str(norm_value) + "\n")
            
            #y_ax_intra.append((dist_intra[x]*1.0) / (norm_value * (binId(chrms_length[chr]) + 1 - x)) )             
            #values_intra.write(str(x) + "\t" + str((dist_intra[x]*1.0) / (norm_value * (binId(chrms_length[chr]) + 1 - x)) ) + "\n")
            #values_raw_intra.write(str(x) + "\t" + str(dist_intra[x]) + "\t" +  str(norm_value) + "\t" +  str(binId(chrms_length[chr]) + 1 - x) + "\n")
        else:
            y_ax_intra.append(None)
            values_intra.write(str(x) + "\t" + str(0) + "\n")
            values_raw_intra.write(str(x) + "\t" + str(0) + ";\t" +  str(norm_value) + "\n")
    
    x_ax_np=np.array(x_ax).astype(np.integer)
    y_ax_intra_np=np.array(y_ax_intra).astype(np.double)
    y_ax_inter_mask=np.isfinite(y_ax_intra_np)
    
    fig, ax = plt.subplots(figsize=(6, 6))    
    
    
    #ax.plot(x_ax, y_ax_intra, linewidth=1, label='Intrachromosomal, ' + chr)
    ax.plot(x_ax_np[y_ax_inter_mask], y_ax_intra_np[y_ax_inter_mask], linewidth=1, label='Intrachromosomal, ' + chr)
    
    #ax.plot(x_ax, y_ax_inter_chr1, '--', linewidth=0.5, label='Interchromosomal, ' + chr +'_inc chr1')
    #ax.plot(x_ax, y_ax_inter_chr1_A, '--', linewidth=0.5, label='Interchromosomal, ' + chr +'_ex chr1')    
    ax.plot(x_ax, y_ax_inter_chr1_A, '--', linewidth=0.5, label='Interchromosomal, ' + chr +'')    
    ax.plot(x_ax, y_ax_inter_chr10, '--',linewidth=0.5, label='Interchromosomal, ' + chr +'-chr10')
    ax.plot(x_ax, y_ax_inter_chr21, '--', linewidth=0.5, label='Interchromosomal, ' + chr +'-chr21')
    
    tmp_title=title + " " + chr + "_" + df_name
    
    ax.set_title(tmp_title, fontsize=14)    
    ax.legend(loc='upper right', prop={'size': 6})
    
    #plt.savefig(path + paths_pics + "chrm_territory_"+chr+"_" +df_name+"_"+mode+"_DNA_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    
    ax.set_yscale("log", nonposy='clip')
    ax.set_ylim(0.0000001, 0.01)
    #plt.savefig(path + paths_pics + "chrm_territory_"+chr+"_" +df_name+"_"+mode+"_DNA_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    
    #ax.set_xscale("log", nonposx='clip')
    ax.set_xscale("symlog", linthreshy=1)
    #plt.savefig(path + paths_pics + "chrm_territory_"+chr+"_" +df_name+"_"+mode+"_DNA_ylog_xlog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    
    
    plt.close()    
    values_intra.close()
    values_raw_intra.close()
    values_inter.close()
    value_raw_inter.close()    
#
#
def make_combined_graph(fl1_intra_in, fl2_intra_in, fl3_intra_in, fl4_intra_in, fl1_inter_in, fl2_inter_in, fl3_inter_in, fl4_inter_in, name, chr):
    path_pics="pics/"
    fl1=open(fl1_intra_in, 'r')
    fl2=open(fl2_intra_in, 'r')
    fl3=open(fl3_intra_in, 'r')
    fl4=open(fl4_intra_in, 'r')
    
    fl1=list(enumerate(fl1))
    fl2=list(enumerate(fl2))
    fl3=list(enumerate(fl3))
    fl4=list(enumerate(fl4))
    
    
    fl1_inter=open(fl1_inter_in, 'r')
    fl2_inter=open(fl2_inter_in, 'r')
    fl3_inter=open(fl3_inter_in, 'r')
    fl4_inter=open(fl4_inter_in, 'r')  
    
    fl1_inter=list(enumerate(fl1_inter))
    fl2_inter=list(enumerate(fl2_inter))
    fl3_inter=list(enumerate(fl3_inter))
    fl4_inter=list(enumerate(fl4_inter))
    
    x_ax=[]
    y_ax1=[]
    y_ax2=[]
    y_ax3=[]
    y_ax4=[]
    
    y_inter_ax1=[]
    y_inter_ax2=[]
    y_inter_ax3=[]
    y_inter_ax4=[]    
    
    for i in range(len(fl1)):
        x_ax.append(int(fl1[i][1].rstrip().split("\t")[0]))
        
        y_ax1.append(float(fl1[i][1].rstrip().split("\t")[1])) 
        y_ax2.append(float(fl2[i][1].rstrip().split("\t")[1])) 
        y_ax3.append(float(fl3[i][1].rstrip().split("\t")[1])) 
        y_ax4.append(float(fl4[i][1].rstrip().split("\t")[1])) 
        
        y_inter_ax1.append(float(fl1_inter[i][1].rstrip().split("\t")[1])) 
        y_inter_ax2.append(float(fl2_inter[i][1].rstrip().split("\t")[1])) 
        y_inter_ax3.append(float(fl3_inter[i][1].rstrip().split("\t")[1])) 
        y_inter_ax4.append(float(fl4_inter[i][1].rstrip().split("\t")[1]))         
    
    fig, ax = plt.subplots(figsize=(6, 6))    
        
    ax.plot(x_ax, y_ax1,  linewidth=0.5, color='blue', label='Intrachromosomal, ' + fl1_intra_in.split("_values.tab")[0])
    ax.plot(x_ax, y_ax2,  linewidth=0.5, color='orange', label='Intrachromosomal, ' + fl2_intra_in.split("_values.tab")[0])    
    ax.plot(x_ax, y_ax3,  linewidth=0.5, color='green', label='Intrachromosomal, ' + fl3_intra_in.split("_values.tab")[0])
    ax.plot(x_ax, y_ax4,  linewidth=0.5, color='red', label='Intrachromosomal, ' + fl4_intra_in.split("_values.tab")[0])
    
    ax.plot(x_ax, y_inter_ax1,  '--',  linewidth=0.5, color='blue', label='Intrachromosomal, ' + fl1_inter_in.split("_values_inter.tab")[0])
    ax.plot(x_ax, y_inter_ax2,  '--', linewidth=0.5,  color='orange', label='Intrachromosomal, ' + fl2_inter_in.split("_values_inter.tab")[0])    
    ax.plot(x_ax, y_inter_ax3,  '--',  linewidth=0.5, color='green', label='Intrachromosomal, ' + fl3_inter_in.split("_values_inter.tab")[0])
    ax.plot(x_ax, y_inter_ax4,  '--',  linewidth=0.5, color='red', label='Intrachromosomal, ' + fl4_inter_in.split("_values_inter.tab")[0])    
    
    ax.legend(loc='upper right', prop={'size': 6},frameon=False)
    
    #plt.savefig(path_pics+"chrm_territory_" + chr + "_combined_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6))
    ax.set_yscale("log", nonposy='clip')
    ax.set_ylim(0.0000001, 0.01)
    #plt.savefig(path_pics+"chrm_territory_" + chr + "_combined_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    ax.set_xscale("symlog")
    plt.savefig(path_pics+"chrm_territory_" + chr + "_combined_ylog_xlog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close()    

#
#
def make_combined_graph_for_whole_genome(path, path_pics, chrms):
    names={}    
    ax1=["dna", "rna"]
    ax2=["introns", "exons"]   
    ax3= ["inter", "intra"]
    
    #--------------------------------Init dictionary-------------------------------- 
    for (a,b,c) in list(itertools.product(ax1, ax2, ax3)):
        names[a+b+c]={}
        names[a+b+c]["file"]={}
        for chr in chrms:
            if c=="intra":
                names[a+b+c]["file"][chr]=path + "k562_" + a + "_" + chr + "_" + b + "_raw_values_" + c + "_" +str(bin_size)+".tab"
            if c=="inter":
                names[a+b+c]["file"][chr]=path + "k562_" + a + "_" + chr + "_" + b + "_raw_values_by_chr_" + c + "_" +str(bin_size)+".tab"
        names[a+b+c]["values_from_file"]={}
        names[a+b+c]["values_to_plot"]=[]
    
    #--------------------------------Reading into dictionary-------------------------------- 
    tmp_x_ax=[]
    for chr in chrms:
        for (a,b,c) in list(itertools.product(ax1, ax2, ax3)):
            names[a+b+c]["values_from_file"][chr]={}
            with open(names[a+b+c]["file"][chr], 'r+') as inLines:
                for i, l in enumerate(inLines):
                    ln=l.rstrip().split("\t")
                    if c=="intra":                        
                        names[a+b+c]["values_from_file"][chr][int(ln[0])]= [ln[1], int(ln[2])]
                        tmp_x_ax.append(int(ln[0]))
                    else:
                        #print l
                        names[a+b+c]["values_from_file"][chr][ln[0]] = [int(ln[1]), float(ln[2]), float(ln[3])]
    
    #--------------------------------X axis--------------------------------  
    tmp_set = set(tmp_x_ax)    
    x_ax=np.array(sorted(list(tmp_set))).astype(np.integer)
    print "x_ax " ,
    print x_ax
    
    #--------------------------------INTRA-------------------------------- 
    
    def return_conf_interval(vals):    
        std=np.std(vals, ddof=1)
        mn=np.mean(vals)        
        return [mn - 1.96 * (std/np.sqrt(len(vals))) , mn + 1.96 * (std/np.sqrt(len(vals)))]
    
    def return_bootstrap_values_for_mean(vals, number_of_iterations):
        mean_estimates=[]
        vals=np.array(vals)
        for _ in range(number_of_iterations):            
            re_sample_idx = np.random.randint(0, len(vals), len(vals))
            mean_estimates.append(np.mean(vals[re_sample_idx]))        
        return mean_estimates
       
    def return_bootsrap_values_for_sampling_whole_table_intra(table, number_of_iterations, chr_list):        
        curves={}       
        dist_range={}  
        dist_new_table={}  
        tmp_numerator={}
        tmp_denominator2={}          
        tmp_denominator={}
        for chr in chr_list:
            dist_range[chr] = binId(chrms_length[chr]) + 1              
        for i in range(number_of_iterations): 
            dist_new_table[i]={}
            tmp_numerator[i]=0
            tmp_denominator2[i]=0 
            tmp_denominator[i]=0
            curves[i]=[]
            for chr in chr_list:
                dist_new_table[i][chr]={}          
        
        for x in x_ax: 
            for i in range(number_of_iterations): 
                tmp_numerator[i]=0
                tmp_denominator2[i]=0
                tmp_denominator[i]=0
            flat_counts=[]

            for chr in table:  
                #flat_counts=[]
                if x not in table[chr]:
                    continue       
                for i in range(number_of_iterations):
                    #tmp_denominator[i]+= float(table[chr][x][1])
                    tmp_denominator2[i]+= float(binId(chrms_length[chr]) + 1 - x)   
                for val in table[chr][x][0].split(";"):
                    try:
                        if val!='':                           
                            flat_counts.append(float(val)/(1.0 * float(table[chr][x][1])))                            
                    except ValueError:
                        print val 
                
            vals = np.array(flat_counts)                
            for i in range(number_of_iterations): 
                re_sample_idx=np.random.randint(0, len(vals), len(vals))
                sum_counts=sum(vals[re_sample_idx])
                #tmp_numerator[i] += ( 1.0 * sum_counts ) / (1.0 * float(table[chr][x][1]))
                tmp_numerator[i] = ( 1.0 * sum_counts ) 
                
            for i in range(number_of_iterations):
                if tmp_numerator[i]!=0:
                    #curves[i].append(tmp_numerator[i]/(tmp_denominator[i] * tmp_denominator2[i]))
                    curves[i].append(tmp_numerator[i]/( tmp_denominator2[i]))
                else:
                    curves[i].append(None) 
                        
                        
        '''        
                new_d_list=np.random.randint(0, dist_range[chr], count)                
                for i in range(number_of_iterations):
                    for d in new_d_list:
                        if d not in dist_new_table[i][chr]:
                            dist_new_table[i][chr][d]=0
                        dist_new_table[i][chr][d]+=1    
                        tmp_numerator[i]+= ( 1.0 )                
                    tmp_numerator[i] = ( 1.0 * tmp_numerator[i] ) / (1.0 * float(table[chr][x][1]))
                    tmp_denominator2[i]+= float(binId(chrms_length[chr]) + 1 - x)
            
            for i in range(number_of_iterations):        
                if tmp_numerator!=0:
                    curves[i].append(tmp_numerator[i]/tmp_denominator2[i])
                else:
                    curves[i].append(None)                         
                        
        '''                 
        '''
        vals = np.array(flat_counts)        
        for i in range(number_of_iterations): 
            re_sample_idx=np.random.randint(0, len(vals), len(vals))
            for d in vals[re_sample_idx]:
                if d not in dist_new_table[i][chr]:
                    dist_new_table[i][chr][d]=0
                dist_new_table[i][chr][d]+=1
        '''    
                        
        '''            
        for i in range(number_of_iterations):                
            curves[i]=[]
            for x in x_ax:
                tmp_numerator=0
                tmp_denominator2=0
                for chr in table:
                    if x not in dist_new_table[i][chr]:                        
                        continue
                    tmp_numerator += ( 1.0 * dist_new_table[i][chr][x]) / (1.0 * float(table[chr][x][1]))
                    tmp_denominator2 += float(binId(chrms_length[chr]) + 1 - x)
                if tmp_numerator!=0:
                    curves[i].append(tmp_numerator/tmp_denominator2)
                else:
                    curves[i].append(None)    
         '''           
                    
        return curves
    
    
      
    names_statistics={}
    bootsrtap_i=10
    for (a,b,c) in list(itertools.product(ax1, ax2, ["intra"])):
        names_statistics[a+b+c]={}
        names_statistics[a+b+c]["values_to_plot_conf_int_low"]=[]
        names_statistics[a+b+c]["values_to_plot_conf_int_high"]=[]
        names_statistics[a+b+c]["values_to_plot_bootstrap"]={}
        names_statistics[a+b+c]["bootstrap_curves"]={}
    stat_lines={}
    for x in x_ax:
        for (a,b,c) in list(itertools.product(ax1, ax2, ["intra"])):    
            tmp_numerator=[]
            tmp_divided=[]
            tmp_denominator=0    
            tmp_denominator2=0 
            
            for chr in chrms:
                if chr not in stat_lines:
                    stat_lines[chr]=""
                count_is_one=0
                all_counts=[]
                if chr not in names[a+b+c]["values_from_file"]:
                    continue
                if x not in names[a+b+c]["values_from_file"][chr]:
                    continue
                for val in names[a+b+c]["values_from_file"][chr][x][0].split(";"):
                    try:
                        if val!='':
                            t=int(val)
                            all_counts.append(t)
                            if t==1:
                                count_is_one+=1
                            #tmp_numerator.append((t * 1.0) / (1.0 * (binId(chrms_length[chr]) + 1 - x)))
                            tmp_numerator.append((t * 1.0) / (1.0 * float(names[a+b+c]["values_from_file"][chr][x][1])))
                    except ValueError:
                        print val
                stat_lines[chr] = stat_lines[chr] + " " + str(count_is_one) + "-" + str(len(all_counts))
                #tmp_numerator.extend(map(int, [val for val in names[a+b+c]["values_from_file"][chr][x][0].split(";") if val!="" and int(val) > 0 ] ))
                tmp_denominator2+=float(binId(chrms_length[chr]) + 1 - x)
                #tmp_denominator+=float(names[a+b+c]["values_from_file"][chr][x][1])
            if len(tmp_numerator)!=0: 
                tmp_divided = [(val * 1.0) / (1.0 * tmp_denominator2) for val in tmp_numerator]
                names[a+b+c]["values_to_plot"].append(sum(tmp_numerator)/(tmp_denominator2))
                #names[a+b+c]["values_to_plot"].append(np.mean(tmp_numerator)/tmp_denominator)
                #names[a+b+c]["values_to_plot"].append(np.mean(tmp_divided))
                #ci=return_conf_interval(tmp_divided)
                #names_statistics[a+b+c]["values_to_plot_conf_int_low"].append(ci[0])
                #names_statistics[a+b+c]["values_to_plot_conf_int_high"].append(ci[1])
                #names_statistics[a+b+c]["values_to_plot_bootstrap"][x]=return_bootstrap_values_for_mean(tmp_divided, 200)
            else:
                names[a+b+c]["values_to_plot"].append(None)    
                names_statistics[a+b+c]["values_to_plot_conf_int_low"].append(None)
                names_statistics[a+b+c]["values_to_plot_conf_int_high"].append(None)   
                names_statistics[a+b+c]["values_to_plot_bootstrap"][x]=None
    '''
    for chr in sorted(chrms):
        print chr ,
        print stat_lines[chr]
    '''    
    print "Bootstrapping ... "
    for (a,b,c) in list(itertools.product(ax1, ax2, ["intra"])): 
        print (a,b,c)
        chr_list=names[a+b+c]["values_from_file"].keys()
        names_statistics[a+b+c]["bootstrap_curves"]=return_bootsrap_values_for_sampling_whole_table_intra(names[a+b+c]["values_from_file"], bootsrtap_i, chr_list)
        
    dna_introns_intra_y_ax_np=np.array(names["dnaintronsintra"]["values_to_plot"]).astype(np.double)
    dna_introns_intra_y_ax_mask=np.isfinite(np.log(dna_introns_intra_y_ax_np))
        
    dna_exons_intra_y_ax_np=np.array(names["dnaexonsintra"]["values_to_plot"]).astype(np.double)
    dna_exons_intra_y_ax_mask=np.isfinite(np.log(dna_exons_intra_y_ax_np))
    
    rna_introns_intra_y_ax_np=np.array(names["rnaintronsintra"]["values_to_plot"]).astype(np.double)
    rna_introns_intra_y_ax_mask=np.isfinite(np.log(rna_introns_intra_y_ax_np))
    
    rna_exons_intra_y_ax_np=np.array(names["rnaexonsintra"]["values_to_plot"]).astype(np.double)
    rna_exons_intra_y_ax_mask=np.isfinite(np.log(rna_exons_intra_y_ax_np))
        
    print rna_exons_intra_y_ax_np[-500:]
    
    #--------------------------------INTER--------------------------------
    
    num_bins_genome=0
    for ch in chrms_length:
        num_bins_genome+=binId(chrms_length[ch])         
    for (a,b,c) in list(itertools.product(ax1, ax2, ["inter"])): 
        tmp_numerator=0
        tmp_denominator=0
        tmp_denominator2=0        
        for chr in chrms:
            if chr not in names[a+b+c]["values_from_file"]:
                continue
            for chr_2 in names[a+b+c]["values_from_file"][chr]:            
                if chr!=chr_2:
                    #tmp_numerator+=float(names[a+b+c]["values_from_file"][chr][chr_2][1])
                    tmp_numerator+=float(names[a+b+c]["values_from_file"][chr][chr_2][1])/ (1.0 * float(names[a+b+c]["values_from_file"][chr][chr][2]))
                else:                
                    chr_length=binId(chrms_length[chr]) 
                    #tmp_numerator+=float(names[a+b+c]["values_from_file"][chr][chr_2][1])
                    #tmp_denominator2+=float(names[a+b+c]["values_from_file"][chr][chr][2])
            tmp_denominator+=float((num_bins_genome - chr_length) * chr_length)         
        for x in x_ax:        
            #names[a+b+c]["values_to_plot"].append(tmp_numerator/(1.0 * tmp_denominator2 * tmp_denominator )) 
            names[a+b+c]["values_to_plot"].append(tmp_numerator/(1.0 * tmp_denominator )) 
            
            
    dna_introns_inter_y_ax=np.array(names["dnaintronsinter"]["values_to_plot"]).astype(np.double)
    dna_exons_inter_y_ax=np.array(names["dnaexonsinter"]["values_to_plot"]).astype(np.double)
    rna_introns_inter_y_ax=np.array(names["rnaintronsinter"]["values_to_plot"]).astype(np.double)
    rna_exons_inter_y_ax=np.array(names["rnaexonsinter"]["values_to_plot"]).astype(np.double)       
    
    
    for (a,b,c) in list(itertools.product(ax1, ax2, ["inter"])):
        names_statistics[a+b+c]={}
        names_statistics[a+b+c]["bootstrap_curves"]={}
    
    def return_bootsrap_values_for_sampling_whole_table_inter(table, number_of_iterations, chr_list):
        curves={}       
 
        tmp_numerator={}
        tmp_denominator=0
           
        for i in range(number_of_iterations): 
            tmp_numerator[i]=0
            curves[i]=[]

        flat_counts=[]
        
        for chr in table:
            if chr not in table[chr]:
                continue
            for chr_2 in table[chr]:            
                if chr!=chr_2:
                    flat_counts.append(float(table[chr][chr_2][1]) / 
                                       (1.0 * float(table[chr][chr][2])))
                else:                
                    chr_length=binId(chrms_length[chr]) 
                    
            tmp_denominator+=float((num_bins_genome - chr_length) * chr_length)        
            
        vals = np.array(flat_counts)                
        for i in range(number_of_iterations): 
            re_sample_idx=np.random.randint(0, len(vals), len(vals))
            sum_counts=sum(vals[re_sample_idx])
            tmp_numerator[i] = ( 1.0 * sum_counts ) 
    
        for x in x_ax: 
            for i in range(number_of_iterations):
                if tmp_numerator[i]!=0:
                    curves[i].append(tmp_numerator[i]/( tmp_denominator))
                else:
                    curves[i].append(None) 
        return curves            
                

    
    print "Bootstrapping ... "
    for (a,b,c) in list(itertools.product(ax1, ax2, ["inter"])): 
        print (a,b,c)
        chr_list=names[a+b+c]["values_from_file"].keys()
        names_statistics[a+b+c]["bootstrap_curves"]=return_bootsrap_values_for_sampling_whole_table_inter(names[a+b+c]["values_from_file"], bootsrtap_i, chr_list)
    
    #--------------------------------Normal picture--------------------------------
    
    fig, ax = plt.subplots(figsize=(6, 6))    
        
    ax.plot(x_ax[dna_introns_intra_y_ax_mask], dna_introns_intra_y_ax_np[dna_introns_intra_y_ax_mask], 
            linewidth=0.5, color='blue', label='Intrachromosomal, dna_introns_intra')
    ax.plot(x_ax[dna_exons_intra_y_ax_mask], dna_exons_intra_y_ax_np[dna_exons_intra_y_ax_mask], 
            linewidth=0.5, color='orange', label='Intrachromosomal, dna_exons_intra')
    ax.plot(x_ax[rna_introns_intra_y_ax_mask], rna_introns_intra_y_ax_np[rna_introns_intra_y_ax_mask], 
            linewidth=0.5, color='green', label='Intrachromosomal, rna_introns_intra')
    ax.plot(x_ax[rna_exons_intra_y_ax_mask], rna_exons_intra_y_ax_np[rna_exons_intra_y_ax_mask], 
            linewidth=0.5, color='red', label='Intrachromosomal, rna_exons_intra')   
   
    ax.plot(x_ax, dna_introns_inter_y_ax, 
            '--', linewidth=0.5, color='blue', label='Intrachromosomal, dna_introns_inter')
    ax.plot(x_ax, dna_exons_inter_y_ax, 
            '--', linewidth=0.5, color='orange', label='Intrachromosomal, dna_exons_inter')
    ax.plot(x_ax, rna_introns_inter_y_ax, 
            '--', linewidth=0.5, color='green', label='Intrachromosomal, rna_introns_inter')
    ax.plot(x_ax, rna_exons_inter_y_ax, 
            '--', linewidth=0.5, color='red', label='Intrachromosomal, rna_exons_inter')
    
    ax.legend(loc = 'upper right', prop={'size': 6}, frameon=False)
    
    #plt.savefig(path_pics+"whole_genome_combined_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6))
    ax.set_yscale("log", nonposy='clip')
    #ax.set_ylim(0.000000001, 0.01)
    #plt.savefig(path_pics+"whole_genome_combined_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    ax.set_xscale("symlog")
    plt.savefig(path_pics+"whole_genome_combined_ylog_xlog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close()     
    
    #--------------------------------Normal picture with CI area--------------------------------
    '''
    fig, ax = plt.subplots(figsize=(6, 6))    
        
    ax.plot(x_ax[dna_introns_intra_y_ax_mask], dna_introns_intra_y_ax_np[dna_introns_intra_y_ax_mask], 
            linewidth=0.5, color='blue', label='Intrachromosomal, dna_introns_intra')
    ax.plot(x_ax[dna_exons_intra_y_ax_mask], dna_exons_intra_y_ax_np[dna_exons_intra_y_ax_mask], 
            linewidth=0.5, color='orange', label='Intrachromosomal, dna_exons_intra')
    ax.plot(x_ax[rna_introns_intra_y_ax_mask], rna_introns_intra_y_ax_np[rna_introns_intra_y_ax_mask], 
            linewidth=0.5, color='green', label='Intrachromosomal, rna_introns_intra')
    ax.plot(x_ax[rna_exons_intra_y_ax_mask], rna_exons_intra_y_ax_np[rna_exons_intra_y_ax_mask], 
            linewidth=0.5, color='red', label='Intrachromosomal, rna_exons_intra')   
   
    ax.plot(x_ax, names_statistics["dnaintronsintra"]["values_to_plot_conf_int_low"], 
            linewidth=0.1, alpha=0.2, color='blue')
    ax.plot(x_ax, names_statistics["dnaintronsintra"]["values_to_plot_conf_int_high"], 
            linewidth=0.1, alpha=0.2, color='blue')    
    ax.fill_between(x_ax, names_statistics["dnaintronsintra"]["values_to_plot_conf_int_high"], 
                    names_statistics["dnaintronsintra"]["values_to_plot_conf_int_low"], 
                    color='blue', alpha=0.2)
    
    ax.plot(x_ax, names_statistics["dnaexonsintra"]["values_to_plot_conf_int_low"], 
            linewidth=0.1, alpha=0.2, color='orange')
    ax.plot(x_ax[dna_introns_intra_y_ax_mask], names_statistics["dnaexonsintra"]["values_to_plot_conf_int_high"], 
            linewidth=0.1, alpha=0.2, color='orange')    
    ax.fill_between(x, names_statistics["dnaexonsintra"]["values_to_plot_conf_int_high"], 
                    names_statistics["dnaexonsintra"]["values_to_plot_conf_int_low"], 
                    color='orange', alpha=0.2)
    
    
    ax.plot(x_ax[dna_introns_intra_y_ax_mask], names_statistics["rnaintronsintra"]["values_to_plot_conf_int_low"], 
            linewidth=0.1, alpha=0.2, color='green')
    ax.plot(x_ax[dna_introns_intra_y_ax_mask], names_statistics["rnaintronsintra"]["values_to_plot_conf_int_high"], 
            linewidth=0.1, alpha=0.2, color='green')    
    ax.fill_between(x, names_statistics["rnaintronsintra"]["values_to_plot_conf_int_high"], 
                    names_statistics["rnaintronsintra"]["values_to_plot_conf_int_low"], 
                    color='green', alpha=0.2)    


    ax.plot(x_ax[dna_introns_intra_y_ax_mask], names_statistics["rnaexonsintra"]["values_to_plot_conf_int_low"], 
            linewidth=0.1, alpha=0.2, color='red')
    ax.plot(x_ax[dna_introns_intra_y_ax_mask], names_statistics["rnaexonsintra"]["values_to_plot_conf_int_high"], 
            linewidth=0.1, alpha=0.2, color='red')    
    ax.fill_between(x, names_statistics["rnaexonsintra"]["values_to_plot_conf_int_high"], 
                    names_statistics["rnaexonsintra"]["values_to_plot_conf_int_low"], 
                    color='red', alpha=0.2)
    
    
    ax.plot(x_ax, dna_introns_inter_y_ax, 
            '--', linewidth=0.5, color='blue', label='Intrachromosomal, dna_introns_inter')
    ax.plot(x_ax, dna_exons_inter_y_ax, 
            '--', linewidth=0.5, color='orange', label='Intrachromosomal, dna_exons_inter')
    ax.plot(x_ax, rna_introns_inter_y_ax, 
            '--', linewidth=0.5, color='green', label='Intrachromosomal, rna_introns_inter')
    ax.plot(x_ax, rna_exons_inter_y_ax, 
            '--', linewidth=0.5, color='red', label='Intrachromosomal, rna_exons_inter')
    
    ax.legend(loc = 'upper right', prop={'size': 6}, frameon=False)
    
    #plt.savefig(path_pics+"whole_genome_combined_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6))
    ax.set_yscale("log", nonposy='clip')
    #ax.set_ylim(0.000000001, 0.01)
    #plt.savefig(path_pics+"whole_genome_combined_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    ax.set_xscale("symlog")
    plt.savefig(path_pics+"whole_genome_combined_ylog_xlog_CI" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close() 
    
    '''
    #--------------------------------Normal picture with bootsrtapped curves--------------------------------
    
    fig, ax = plt.subplots(figsize=(6, 6))    
    
    
    
    for i in range(bootsrtap_i):                
        ax.plot(x_ax, names_statistics["dnaintronsintra"]["bootstrap_curves"][i], 
            linewidth=0.1, alpha=0.1, color='grey')
        ax.plot(x_ax, names_statistics["dnaexonsintra"]["bootstrap_curves"][i], 
                linewidth=0.1, alpha=0.1, color='grey')
        ax.plot(x_ax, names_statistics["rnaintronsintra"]["bootstrap_curves"][i], 
                linewidth=0.1, alpha=0.1, color='grey')
        ax.plot(x_ax, names_statistics["rnaexonsintra"]["bootstrap_curves"][i], 
                linewidth=0.1, alpha=0.1, color='grey')     
        
        
    ax.plot(x_ax[dna_introns_intra_y_ax_mask], dna_introns_intra_y_ax_np[dna_introns_intra_y_ax_mask], 
            linewidth=0.5, color='blue', label='Intrachromosomal, dna_introns_intra')
    ax.plot(x_ax[dna_exons_intra_y_ax_mask], dna_exons_intra_y_ax_np[dna_exons_intra_y_ax_mask], 
            linewidth=0.5, color='orange', label='Intrachromosomal, dna_exons_intra')
    ax.plot(x_ax[rna_introns_intra_y_ax_mask], rna_introns_intra_y_ax_np[rna_introns_intra_y_ax_mask], 
            linewidth=0.5, color='green', label='Intrachromosomal, rna_introns_intra')
    ax.plot(x_ax[rna_exons_intra_y_ax_mask], rna_exons_intra_y_ax_np[rna_exons_intra_y_ax_mask], 
            linewidth=0.5, color='red', label='Intrachromosomal, rna_exons_intra')   
    
    
    for i in range(bootsrtap_i):                
        ax.plot(x_ax, names_statistics["dnaintronsinter"]["bootstrap_curves"][i], 
            linewidth=0.1, alpha=0.1, color='grey')
        ax.plot(x_ax, names_statistics["dnaexonsinter"]["bootstrap_curves"][i], 
                linewidth=0.1, alpha=0.1, color='grey')
        ax.plot(x_ax, names_statistics["rnaintronsinter"]["bootstrap_curves"][i], 
                linewidth=0.1, alpha=0.1, color='grey')
        ax.plot(x_ax, names_statistics["rnaexonsinter"]["bootstrap_curves"][i], 
                linewidth=0.1, alpha=0.1, color='grey')     
        
 
    
    ax.plot(x_ax, dna_introns_inter_y_ax, 
            '--', linewidth=0.5, color='blue', label='Intrachromosomal, dna_introns_inter')
    ax.plot(x_ax, dna_exons_inter_y_ax, 
            '--', linewidth=0.5, color='orange', label='Intrachromosomal, dna_exons_inter')
    ax.plot(x_ax, rna_introns_inter_y_ax, 
            '--', linewidth=0.5, color='green', label='Intrachromosomal, rna_introns_inter')
    ax.plot(x_ax, rna_exons_inter_y_ax, 
            '--', linewidth=0.5, color='red', label='Intrachromosomal, rna_exons_inter')
    
    ax.legend(loc = 'upper right', prop={'size': 6}, frameon=False)
    
    #plt.savefig(path_pics+"whole_genome_combined_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6))
    ax.set_yscale("log", nonposy='clip')
    #ax.set_ylim(0.000000001, 0.01)
    #plt.savefig(path_pics+"whole_genome_combined_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    ax.set_xscale("symlog")
    plt.savefig(path_pics+"whole_genome_combined_ylog_xlog_BS_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close()       
    
    
    '''
    #--------------------------------Normal picture with bootsrtapped dots--------------------------------
    
    fig, ax = plt.subplots(figsize=(6, 6))    
        
    ax.plot(x_ax[dna_introns_intra_y_ax_mask], dna_introns_intra_y_ax_np[dna_introns_intra_y_ax_mask], 
            linewidth=0.5, color='blue', label='Intrachromosomal, dna_introns_intra')
    ax.plot(x_ax[dna_exons_intra_y_ax_mask], dna_exons_intra_y_ax_np[dna_exons_intra_y_ax_mask], 
            linewidth=0.5, color='orange', label='Intrachromosomal, dna_exons_intra')
    ax.plot(x_ax[rna_introns_intra_y_ax_mask], rna_introns_intra_y_ax_np[rna_introns_intra_y_ax_mask], 
            linewidth=0.5, color='green', label='Intrachromosomal, rna_introns_intra')
    ax.plot(x_ax[rna_exons_intra_y_ax_mask], rna_exons_intra_y_ax_np[rna_exons_intra_y_ax_mask], 
            linewidth=0.5, color='red', label='Intrachromosomal, rna_exons_intra')   
    
    
    print [6] * 10
    for x in x_ax:
        print len([x])
        
        print len ([x] * len(names_statistics["dnaintronsintra"]["values_to_plot_bootstrap"]))
        
        print len(names_statistics["dnaintronsintra"]["values_to_plot_bootstrap"])
        
        ax.scatter([x] * len(names_statistics["dnaintronsintra"]["values_to_plot_bootstrap"]),
                   names_statistics["dnaintronsintra"]["values_to_plot_bootstrap"],
                   alpha=0.1, color='blue', s=1)
        
        ax.scatter([x] * len(names_statistics["dnaexonsintra"]["values_to_plot_bootstrap"]),
                   names_statistics["dnaexonsintra"]["values_to_plot_bootstrap"],
                   alpha=0.1, color='orange', s=1)
        
        ax.scatter([x] * len(names_statistics["rnaintronsintra"]["values_to_plot_bootstrap"]),
                   names_statistics["rnaintronsintra"]["values_to_plot_bootstrap"],
                   alpha=0.1, color='green', s=1)
        
        ax.scatter([x] * len(names_statistics["rnaexonsintra"]["values_to_plot_bootstrap"]),
                   names_statistics["rnaexonsintra"]["values_to_plot_bootstrap"],
                   alpha=0.1, color='red', s=1)        
    

    ax.plot(x_ax, dna_introns_inter_y_ax, 
            '--', linewidth=0.5, color='blue', label='Intrachromosomal, dna_introns_inter')
    ax.plot(x_ax, dna_exons_inter_y_ax, 
            '--', linewidth=0.5, color='orange', label='Intrachromosomal, dna_exons_inter')
    ax.plot(x_ax, rna_introns_inter_y_ax, 
            '--', linewidth=0.5, color='green', label='Intrachromosomal, rna_introns_inter')
    ax.plot(x_ax, rna_exons_inter_y_ax, 
            '--', linewidth=0.5, color='red', label='Intrachromosomal, rna_exons_inter')
    
    ax.legend(loc = 'upper right', prop={'size': 6}, frameon=False)
    
    #plt.savefig(path_pics+"whole_genome_combined_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6))
    ax.set_yscale("log", nonposy='clip')
    #ax.set_ylim(0.000000001, 0.01)
    #plt.savefig(path_pics+"whole_genome_combined_ylog_" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    ax.set_xscale("symlog")
    plt.savefig(path_pics+"whole_genome_combined_ylog_xlog_BS" + str(bin_size)+ ".png", dpi=300, figsize=(6, 6)) 
    plt.close()       
    '''
    
    
    '''        
    dna_introns_inter={}
    dna_exons_inter={}
    rna_introns_inter={}
    rna_exons_inter={}        
    

    
    for chr in chrms:
    #for chr in [""]:
        #### 1-1
        dna_introns_inter[chr]={}
        with open(path + "k562_dna_" + chr + "_introns_raw_values_by_chr_inter_" +str(bin_size)+".tab", 'r+') as introns_inter_data:
            for i, l in enumerate(introns_inter_data):
                ln=l.rstrip().split("\t")
                dna_introns_inter[chr][ln[0]] = [int(ln[1]), float(ln[2]), float(ln[3])]
        #### 1-2
        dna_exons_inter[chr]={}
        with open(path + "k562_dna_" +chr + "_exons_raw_values_by_chr_inter_" +str(bin_size)+".tab", 'r+') as exons_inter_data:
            for i, l in enumerate(exons_inter_data):
                ln=l.rstrip().split("\t")
                dna_exons_inter[chr][ln[0]] = [int(ln[1]), float(ln[2]), float(ln[3])]                    
        #### 1-3               
        rna_introns_inter[chr]={}
        with open(path + "k562_rna_" +chr + "_introns_raw_values_by_chr_inter_" +str(bin_size)+".tab", 'r+') as introns_inter_data:
            for i, l in enumerate(introns_inter_data):
                ln=l.rstrip().split("\t")
                rna_introns_inter[chr][ln[0]] = [int(ln[1]), float(ln[2]), float(ln[3])]
        #### 1-4               
        rna_exons_inter[chr]={}
        with open(path + "k562_rna_" + chr + "_exons_raw_values_by_chr_inter_" +str(bin_size)+".tab", 'r+') as exons_inter_data:
            for i, l in enumerate(exons_inter_data):
                ln=l.rstrip().split("\t")
                rna_exons_inter[chr][ln[0]] = [int(ln[1]), float(ln[2]), float(ln[3])]       
            
    #### 2-2-1
    print "DNA introns"
    tmp_numerator=0
    tmp_denominator=0    
    tmp_denominator2=0
    for chr in chrms:
        chr_length=0        
        for chr_2 in dna_introns_inter[chr]:            
            if chr!=chr_2:
                tmp_numerator+=float(dna_introns_inter[chr][chr_2][1])
            else:                
                chr_length=dna_introns_inter[chr][chr][0]
                #tmp_numerator+=float(dna_introns_inter[chr][chr_2][1])
                tmp_denominator2+=float(dna_introns_inter[chr][chr][2])
        tmp_denominator+=float((num_bins_genome - chr_length) * chr_length) 
         
        print chr ,
        print float(dna_introns_inter[chr][chr_2][1]) ,
        print str(num_bins_genome - chr_length) + " " + str(chr_length)  + " = ",
        print float((num_bins_genome - chr_length) * chr_length)
        print tmp_numerator ,
        print tmp_denominator  ,
        print tmp_numerator/tmp_denominator         
    for x in x_ax:        
        #dna_introns_inter_y_ax.append(tmp_numerator/(1.0 * tmp_denominator * num_bins_genome * num_bins_genome ))     
        dna_introns_inter_y_ax.append(tmp_numerator/(1.0 * tmp_denominator2 * tmp_denominator ))     
        
    #### 2-2-2
    print "DNA exons"
    tmp_numerator=0
    tmp_denominator=0    
    tmp_denominator2=0
    
    for chr in chrms:
        chr_length=0        
        for chr_2 in dna_exons_inter[chr]:
            if chr!=chr_2:
                tmp_numerator+=float(dna_exons_inter[chr][chr_2][1])
            else:
                chr_length=dna_exons_inter[chr][chr][0]
                #tmp_numerator+=float(dna_exons_inter[chr][chr_2][1])
                tmp_denominator2+=float(dna_exons_inter[chr][chr][2])                
        
        tmp_denominator+=float((num_bins_genome - chr_length) * chr_length) 
        
        print chr ,
        print float(dna_exons_inter[chr][chr_2][1]) ,
        print str(num_bins_genome - chr_length) + " " + str(chr_length)  + " = ",
        print float((num_bins_genome - chr_length) * chr_length)
        print tmp_numerator ,
        print tmp_denominator  ,
        print tmp_numerator/tmp_denominator       
    for x in x_ax:        
        #dna_exons_inter_y_ax.append(tmp_numerator/(1.0 * tmp_denominator * num_bins_genome * num_bins_genome )) 
        dna_exons_inter_y_ax.append(tmp_numerator/(1.0 * tmp_denominator2 * tmp_denominator )) 
    
    
    #### 2-2-3
    print "RNA introns"
    tmp_numerator=0
    tmp_denominator=0  
    tmp_denominator2=0
    
    for chr in chrms:
        chr_length=0        
        for chr_2 in rna_introns_inter[chr]:
            if chr!=chr_2:
                tmp_numerator+=float(rna_introns_inter[chr][chr_2][1])
            else:
                chr_length=rna_introns_inter[chr][chr][0]
                #tmp_numerator+=float(rna_introns_inter[chr][chr_2][1])
                tmp_denominator2+=float(rna_introns_inter[chr][chr][2])    
                
        tmp_denominator+=float((num_bins_genome - chr_length) * chr_length)    
        
        print chr ,
        print float(rna_introns_inter[chr][chr_2][1]) ,
        print str(num_bins_genome - chr_length) + " " + str(chr_length)  + " = ",
        print float((num_bins_genome - chr_length) * chr_length)
        print tmp_numerator ,
        print tmp_denominator  ,
        print tmp_numerator/tmp_denominator       
        
    for x in x_ax:        
        #rna_introns_inter_y_ax.append(tmp_numerator/(1.0 * tmp_denominator * num_bins_genome * num_bins_genome ))     
        rna_introns_inter_y_ax.append(tmp_numerator/(1.0 * tmp_denominator2 * tmp_denominator ))     
        
    #### 2-2-4
    print "RNA exons"
    tmp_denominator=0
    tmp_numerator=0
    tmp_denominator2=0
    for chr in chrms:                
        chr_length=0        
        for chr_2 in rna_exons_inter[chr]:
            if chr!=chr_2:
                tmp_numerator+=float(rna_exons_inter[chr][chr_2][1])
                
            else:
                chr_length=rna_exons_inter[chr][chr][0]
                #tmp_numerator+=float(rna_exons_inter[chr][chr_2][1])
                tmp_denominator2+=float(rna_exons_inter[chr][chr][2])  
                
        tmp_denominator+=float((num_bins_genome - chr_length) * chr_length)    
        
        print chr ,
        print float(rna_exons_inter[chr][chr_2][1]) ,
        print str(num_bins_genome - chr_length) + " " + str(chr_length)  + " = ",
        print float((num_bins_genome - chr_length) * chr_length)
        print ""
        print tmp_numerator ,
        print tmp_denominator  ,
        print tmp_numerator/tmp_denominator 
        
    for x in x_ax:        
        #rna_exons_inter_y_ax.append(tmp_numerator/(1.0 * tmp_denominator * num_bins_genome * num_bins_genome ))     
        rna_exons_inter_y_ax.append(tmp_numerator/(1.0 * tmp_denominator2 * tmp_denominator ))     
    ####           
    '''
    '''      
    for (a,b,c) in list(itertools.product(ax1, ax2, ["intra"])):
        print (a,b,c)
        print a+b+c
    #print names['dnaintronsintra']["values_to_plot"]
    '''    
#
#
def Workflow_chromosomes():
    #print "split_mrnas_by_chr_of_rna_part"
    #split_mrnas_by_chr_of_rna_part(k562_rna, "k562_rna")
    #split_dna_contacts_by_chr_of_first_dna_part(k562_dna, "k562_dna")
    print "\n\nload chrms lengths"
    load_length(k562_stat)
    print "\n\nchrms_ordered_by_length"
    print chrms_ordered_by_length
    print "\n\nload_coordinates introns: chrm - number of introns"
    load_coordinates(introns_file, introns, "", introns_by_chr_by_bin)
    #print introns.keys()
    print "\n\nload_coordinates exons: chrm - number of exons"    
    load_coordinates(exons_file, exons, "exons", exons_by_chr_by_bin)
    #print exons.keys()
    '''
    print "REading chr1"
    
    lines=[]
    with open("k562_dna_chr1.tab", 'r') as inlines:
        for i, l in enumerate(inlines):
            lines.append(l)
    
    print len(lines)
    '''
    print "\n\nparse_rna&dna_chrms_data_vs_annotated_coordinates"
    #for chr in chrms_flat:
    tmp_chrms=[]
    tmp_chrms.extend(chrms_flat[:])
    #tmp_chrms.remove('chr1')
    
    number_of_cpu=8
    #parse_chrms_data_vs_annotated_coordinates("", introns_by_chr_by_bin, "introns", "DNA", number_of_cpu)
    #parse_chrms_data_vs_annotated_coordinates("", exons_by_chr_by_bin, "exons", "DNA", number_of_cpu)
    
    #parse_chrms_data_vs_annotated_coordinates("", introns_by_chr_by_bin, "introns", "RNA", number_of_cpu)
    #parse_chrms_data_vs_annotated_coordinates("", exons_by_chr_by_bin, "exons", "RNA", number_of_cpu)    
    
    '''
    for chr in tmp_chrms:
        print chr 
        #parse_dna_chrms_data_vs_annotated_coordinates("k562_dna_" + chr + ".tab", "k562_dna_" + chr + "_introns.tab", introns, chr)
        #parse_dna_chrms_data_vs_annotated_coordinates("k562_dna_" + chr + ".tab", "k562_dna_" + chr + "_exons.tab", exons, chr)
        #parse_rna_chrms_data_vs_annotated_coordinates("k562_rna_" + chr + ".tab", "k562_rna_" + chr + "_introns.tab", introns, chr)
        #parse_rna_chrms_data_vs_annotated_coordinates("k562_rna_" + chr + ".tab", "k562_rna_" + chr + "_exons.tab", exons, chr)
    '''
    
    
    '''
    print "\nmake_intrachromosomal_calculation"
    
    path_values="values/"
    tmp_chrms=[]
    tmp_chrms.extend(chrms_flat[:])
    tmp_chrms.remove('chrY')    
    
    chmrs_needed=[]
    chmrs_needed.extend(chrms_flat[:])
    chmrs_needed.remove('chrY')
    
    statistics=open(path_values + "statistics_"  + str(bin_size)+ ".tab", "w+")    

    for chr in tmp_chrms:
        print chr 
        
        fl1=path_values + "k562_dna_" + chr + "_introns_values_intra_" + str(bin_size)+ ".tab"
        fl2=path_values + "k562_dna_" + chr + "_exons_values_intra_" + str(bin_size)+ ".tab"
        fl3=path_values + "k562_rna_" + chr + "_introns_values_intra_" + str(bin_size)+ ".tab"
        fl4=path_values + "k562_rna_" + chr + "_exons_values_intra_" + str(bin_size)+ ".tab"   
        
        fl1_raw_intra=path_values + "k562_dna_" + chr + "_introns_raw_values_intra_" + str(bin_size)+ ".tab"
        fl2_raw_intra=path_values + "k562_dna_" + chr + "_exons_raw_values_intra_" + str(bin_size)+ ".tab"
        fl3_raw_intra=path_values + "k562_rna_" + chr + "_introns_raw_values_intra_" + str(bin_size)+ ".tab"
        fl4_raw_intra=path_values + "k562_rna_" + chr + "_exons_raw_values_intra_" + str(bin_size)+ ".tab"          
        
        fl1_inter=path_values + "k562_dna_" + chr + "_introns_values_inter_" + str(bin_size)+ ".tab"
        fl2_inter=path_values + "k562_dna_" + chr + "_exons_values_inter_" + str(bin_size)+ ".tab"
        fl3_inter=path_values + "k562_rna_" + chr + "_introns_values_inter_" + str(bin_size)+ ".tab"
        fl4_inter=path_values + "k562_rna_" + chr + "_exons_values_inter_" + str(bin_size)+ ".tab" 
        
        fl1_raw_inter=path_values + "k562_dna_" + chr + "_introns_raw_values_by_chr_inter_" + str(bin_size)+ ".tab"
        fl2_raw_inter=path_values + "k562_dna_" + chr + "_exons_raw_values_by_chr_inter_" + str(bin_size)+ ".tab"
        fl3_raw_inter=path_values + "k562_rna_" + chr + "_introns_raw_values_by_chr_inter_" + str(bin_size)+ ".tab"
        fl4_raw_inter=path_values + "k562_rna_" + chr + "_exons_raw_values_by_chr_inter_" + str(bin_size)+ ".tab" 
        

        make_intrachromosomal_and_interchromosomal_calculation ("RNA", chr, "k562_rna_dna", "exons", "k562_rna_" + chr + "_exons.tab", fl4, fl4_raw_intra, fl4_inter, fl4_raw_inter, statistics, chmrs_needed)
        make_intrachromosomal_and_interchromosomal_calculation ("RNA", chr, "k562_rna_dna", "introns", "k562_rna_" + chr + "_introns.tab", fl3, fl3_raw_intra,  fl3_inter, fl3_raw_inter, statistics, chmrs_needed)
        
        make_intrachromosomal_and_interchromosomal_calculation ("DNA", chr, "k562_dna_dna", "introns", "k562_dna_" + chr + "_introns.tab", fl1, fl1_raw_intra, fl1_inter, fl1_raw_inter, statistics, chmrs_needed)
        make_intrachromosomal_and_interchromosomal_calculation ("DNA", chr, "k562_dna_dna", "exons", "k562_dna_" + chr + "_exons.tab", fl2, fl2_raw_intra, fl2_inter, fl2_raw_inter, statistics, chmrs_needed)        

        make_combined_graph(fl1, fl2, fl3, fl4, fl1_inter, fl2_inter, fl3_inter, fl4_inter, "combined_intra_inter" + chr, chr)
              
    statistics.close()
    '''
    tmp_chrms=[]
    tmp_chrms.extend(chrms_flat[:])
    tmp_chrms.remove('chrY')
       
    make_combined_graph_for_whole_genome("values/", "pics/", tmp_chrms)
#
#
def Workflow_genes():

    load_length(k562_stat)
    print "chrms ordered by length"
    print chrms_ordered_by_length
    print "load_coordinates"
    print "Number of introns per chromosome"   
    load_coordinates(introns_file, introns, "", introns_by_chr_by_bin)
    print "Number of exons per chromosome"   
    load_coordinates(exons_file, exons, "exons", exons_by_chr_by_bin)
    print "Loading gene data"
    load_gene_data(genes_file, genes)    

    print "\nparse_rna&dna_chrms_data_vs_annotated_coordinates\n"
    tmp_chrms=[]
    tmp_chrms.extend(chrms_flat[:])
    
    number_of_cpu=4
    flank=100001
    
    #annotate_DNA_contacts_left_part("", "introns", number_of_cpu)    
    '''
    annotate_DNA_contacts_left_part("", "exons", number_of_cpu)    
    '''
    
    #parse_gene_data_vs_annotated_coordinates ("", "introns", "DNA_exon_intron", number_of_cpu, flank)
    
    
    #parse_gene_data_vs_annotated_coordinates ("", "exons", "DNA_exon_intron", number_of_cpu, flank)
    
    
    #parse_gene_data_vs_annotated_coordinates ("", "introns", "RNA_exon_intron", number_of_cpu, flank)
    
    
    #parse_gene_data_vs_annotated_coordinates ("", "exons", "RNA_exon_intron", number_of_cpu, flank)
    
    print "\nmake_intragene_calculation\n"
    
    path_values="values/"
    tmp_chrms=[]
    tmp_chrms.extend(chrms_flat[:])
    tmp_chrms.remove('chrY')
    
    
    chmrs_needed=[]
    chmrs_needed.extend(chrms_flat[:])
    chmrs_needed.remove('chrY')
    
    
    '''
    fl1=path_values + "k562_dna_introns_values_intragene_" + str(bin_size)+ ".tab"
    fl2=path_values + "k562_dna_exons_values_intragene_" + str(bin_size)+ ".tab"
    fl3=path_values + "k562_rna_introns_values_intragene_" + str(bin_size)+ ".tab"
    fl4=path_values + "k562_rna_exons_values_intragene_" + str(bin_size)+ ".tab"   
    
    fl1_by_gene=path_values + "k562_dna_introns_values_intragene_by_gene_" + str(bin_size)+ ".tab"
    fl2_by_gene=path_values + "k562_dna_exons_values_intragene_by_gene_" + str(bin_size)+ ".tab"
    fl3_by_gene=path_values + "k562_rna_introns_values_intragene_by_gene_" + str(bin_size)+ ".tab"
    fl4_by_gene=path_values + "k562_rna_exons_values_intragene_by_gene_" + str(bin_size)+ ".tab"   
                    
    fl1_raw_intra=path_values + "k562_dna_introns_raw_values_intragene_" + str(bin_size)+ ".tab"
    fl2_raw_intra=path_values + "k562_dna_exons_raw_values_intragene_" + str(bin_size)+ ".tab"
    fl3_raw_intra=path_values + "k562_rna_introns_raw_values_intragene_" + str(bin_size)+ ".tab"
    fl4_raw_intra=path_values + "k562_rna_exons_raw_values_intragene_" + str(bin_size)+ ".tab"          
    
    
    fl1_table=path_values + "k562_dna_introns_values_table_" + str(bin_size)+ ".tab"
    fl2_table=path_values + "k562_dna_exons_values_table_" + str(bin_size)+ ".tab"    
    fl3_table=path_values + "k562_rna_introns_values_table_" + str(bin_size)+ ".tab"
    fl4_table=path_values + "k562_rna_exons_values_table_" + str(bin_size)+ ".tab"
    
    fl1_ratio=path_values + "k562_dna_introns_values_ratio_right_to_left_" + str(bin_size)+ ".tab"
    fl2_ratio=path_values + "k562_dna_exons_values_ratio_right_to_left_" + str(bin_size)+ ".tab"    
    fl3_ratio=path_values + "k562_rna_introns_values_ratio_right_to_left_" + str(bin_size)+ ".tab"
    fl4_ratio=path_values + "k562_rna_exons_values_ratio_right_to_left_" + str(bin_size)+ ".tab"    
    
    '''
    '''
    make_intragene_calculation (mode,  title, df_name, 
                                fl_in, 
                                fl_out_intra, fl_out_intra_by_gene, fl_out_raw_intra, statistics , chrms_needed)
   
    k562_dna_all_genes_introns.tab
    k562_dna_all_genes_exons.tab
    
    k562_rna_all_genes_introns.tab
    k562_rna_all_genes_exons.tab
    '''
    
    
    statistics=open(path_values + "statistics_genes_"  + str(bin_size)+ ".tab", "w+") 
    '''
    make_intragene_calculation_relative_bin ("RNA",  "k562_rna_dna_intragene", "introns", 
                                "k562_rna_all_genes_introns_" + str(flank) + ".tab", 
                                fl3, fl3_by_gene, fl3_raw_intra, statistics, chmrs_needed)
    
    make_intragene_calculation_relative_bin ("RNA",  "k562_rna_dna_intragene", "exons", 
                                "k562_rna_all_genes_exons_" + str(flank) + ".tab", 
                                fl4, fl4_by_gene, fl4_raw_intra, statistics, chmrs_needed)    
    '''
    
    '''
    make_intragene_calculation ("DNA",  "k562_dna_dna_intragene", "introns",
                                "k562_dna_all_genes_introns"  + ".tab", 
                                fl1, fl1_by_gene, fl1_raw_intra, statistics, chmrs_needed, fl1_table)
   
       
    make_intragene_calculation ("DNA",  "k562_dna_dna_intragene", "exons", 
                                "k562_dna_all_genes_exons"  + ".tab", 
                                fl2, fl2_by_gene, fl2_raw_intra, statistics, chmrs_needed, fl2_table)   

    make_intragene_calculation ("RNA",  "k562_rna_dna_intragene", "exons", 
                            "k562_rna_all_genes_exons_" + str(flank) + ".tab", 
                            fl4, fl4_by_gene, fl4_raw_intra, statistics, chmrs_needed, fl4_table) 
    make_intragene_calculation ("RNA",  "k562_rna_dna_intragene", "introns", 
                                "k562_rna_all_genes_introns_" + str(flank) + ".tab", 
                                fl3, fl3_by_gene, fl3_raw_intra, statistics, chmrs_needed, fl3_table)
    
     

        


    
    statistics.close()
    '''
    #make_combined_graph_genes(fl1, fl2, fl3, fl4, "_")

    
 
   
    fl2=path_values + "k562_dna_exons_values_intragene_wo_last_exon_" + str(bin_size)+ ".tab"
    fl4=path_values + "k562_rna_exons_values_intragene_wo_last_exon_" + str(bin_size)+ ".tab"   
    
    fl2_by_gene=path_values + "k562_dna_exons_values_intragene_by_gene_wo_last_exon_" + str(bin_size)+ ".tab"
    fl4_by_gene=path_values + "k562_rna_exons_values_intragene_by_gene_wo_last_exon_" + str(bin_size)+ ".tab"   
                    
    fl2_raw_intra=path_values + "k562_dna_exons_raw_values_intragene_wo_last_exon_" + str(bin_size)+ ".tab"
    fl4_raw_intra=path_values + "k562_rna_exons_raw_values_intragene_wo_last_exon_" + str(bin_size)+ ".tab"   
    
    fl2_table_wo_exon=path_values + "k562_dna_exons_values_table_wo_exon_" + str(bin_size)+ ".tab"    
    fl4_table_wo_exon=path_values + "k562_rna_exons_values_table_wo_exon_" + str(bin_size)+ ".tab"
    
    '''
    filter_out_last_exon("")
    
    '''
    
    '''
    statistics=open(path_values + "statistics_genes_wo_last_exon_"  + str(bin_size)+ ".tab", "w+")        
    
    make_intragene_calculation ("RNA",  "k562_rna_dna_intragene_wo_last_exon", "exons", 
                                "k562_rna_all_genes_exons_wo_last_exon.tab", 
                                fl4, fl4_by_gene, fl4_raw_intra, statistics, chmrs_needed, fl4_table_wo_exon)
    
    make_intragene_calculation ("DNA",  "k562_dna_dna_intragene_wo_last_exon", "exons", 
                                "k562_dna_all_genes_exons_wo_last_exon.tab", 
                                fl2, fl2_by_gene, fl2_raw_intra, statistics, chmrs_needed, fl2_table_wo_exon)   
        
    statistics.close()
    '''
    
    make_combined_graph_genes_wo_last_exon( fl2, fl4, "combined_intragene_wo_last_exon")
    
 
#
#
def test():
    
    
    chr='chr1'
    string_count=0
    fl_in="k562_dna_" + chr + ".tab"
    with gzip.open(fl_in, 'r') as fl_raw:
        for i, l in enumerate(fl_raw):
            ln=l.rstrip().split("\t")
            string_count+=1
            if ln[1]!=chr:
                print "Warning: parse_rna_chrms_data_vs_annotated_coordinates bad input file "
                return           
    print string_count
    
    fl_in="k562_rna_" + chr + ".tab"
    string_count=0
    with open(fl_in, 'r') as fl_raw:
        for i, l in enumerate(fl_raw):
            ln=l.rstrip().split("\t")
            string_count+=1
            if ln[0]!=chr:
                print "Warning: parse_rna_chrms_data_vs_annotated_coordinates bad input file "
                return           
    print string_count
    
    
    #
    '''
    #def core_parse_chrms_data_vs_annotated_coordinates(lines, fl_out, mode, df_starts, df_ends, df_strands):    
    def core_parse_chrms_data_vs_annotated_coordinates( fl_out, mode, df_starts, df_ends, df_strands):    
        import gzip
        import numpy as np
        number_of_entries = 0
        if mode=='RNA':
            with open (fl_out, 'w+') as fl_proc:                         
                for i, ln in enumerate(inLines): 
                    number_of_entries+=1
                    #ln=l.rstrip().split("\t")
                    coord1=int(ln[1])
                    chr=ln[0]            
                    rna_strand=ln[8]   
                    if np.any((df_starts <= coord1) & (df_ends >= coord1) & (df_strands==rna_strand)):
                        fl_proc.write(l)
            return fl_in + "\t"+str(number_of_entries)
        elif mode=='DNA':
            with open (fl_out, 'wt') as fl_proc:             
                for i, ln in enumerate(inLines): 
                    number_of_entries+=1
                    #ln=l.rstrip().split("\t") 
                    coord1=int(ln[2])
                    chr=ln[1]  
                    if np.any((df_starts <= coord1) & (df_ends >= coord1)):
                        fl_proc.write(l)
            return fl_in + "\t"+str(number_of_entries)
            
        else:
            return "Warning: core_parse_chrms_data_vs_annotated_coordinates unknown mode"   
    '''
    
    
    
    '''
    
    def parse_chrms_data_vs_annotated_coordinates(folder, df, df_name, number_of_cpu):
         
     
        #client = ipyparallel.Client(context=pyz.Context(),  location="10.0.1.1") 
        #client = ipyparallel.Client(context=pyz.Context()) 
        #client = ipyparallel.Client() 
        client = ipyparallel.Client(profile='natalia') 
        dview = client[:]
        print "clients ids ",
        print client.ids
        dview.block = False
        #print chrms_flat
        
        filenames_dna_by_sublist=[]
        filenames_rna_by_sublist=[]
        tmp_rna=[]
        tmp_dna=[]
        
        for i in range(len(chrms_ordered_by_length)):
            if chrms_ordered_by_length[i]=="chrY":
                continue
            if len(tmp_rna) == number_of_cpu-1:            
                tmp_rna.append("k562_rna_" + chrms_ordered_by_length[i]+".tab")
                tmp_dna.append("k562_dna_" + chrms_ordered_by_length[i]+".tab")
                filenames_rna_by_sublist.append(tmp_rna)
                filenames_dna_by_sublist.append(tmp_dna)
                tmp_rna=[]
                tmp_dna=[]
            else:
                tmp_rna.append("k562_rna_" + chrms_ordered_by_length[i]+".tab")
                tmp_dna.append("k562_dna_" + chrms_ordered_by_length[i]+".tab")
        if len(tmp_rna)!=0:             
            filenames_rna_by_sublist.append(tmp_rna)
            filenames_dna_by_sublist.append(tmp_dna)    
        
           
        filenames_rna=[] 
        filenames_dna=[]
        for i in range(len(chrms_ordered_by_length)):
            if chrms_ordered_by_length[i]=="chrY":
                continue
            filenames_rna.append("k562_rna_" + chrms_ordered_by_length[i]+".tab")
            filenames_dna.append("k562_dna_" + chrms_ordered_by_length[i]+".tab")
            
    
        results=[None]*number_of_cpu    
        for fl in filenames_dna:
            print fl
            client = ipyparallel.Client(profile='natalia') 
            dview = client[:]
            print "clients ids ",
            print client.ids
            dview.block = False
            with dview.sync_imports():
                import numpy as np        
        
        for i in range(len(filenames_dna_by_sublist)): 
            print filenames_dna_by_sublist[i]
            lines=[None]*number_of_cpu
            client = ipyparallel.Client(profile='natalia') 
            dview = client[:]
            print "clients ids ",
            print client.ids
            dview.block = False
            
            with dview.sync_imports():
                import numpy as np
                
                
            
            with open (os.path.join(path , fl), 'r') as inLines,\
                 open (os.path.join(path , fl.split(".tab")[0] + "_" + df_name + ".tab"), 'w+') as out_file:
                
                for k, l in enumerate (inLines):
                    for engine in range(number_of_cpu):
                        results[j] = dview.apply(core_parse_chrms_data_vs_annotated_coordinates, 
                                                 l,
                                                 "DNA",
                                                 df_starts, 
                                                 df_ends, 
                                                 df_strands)
                                
                    
                    print '\nWaiting dnas...\n'
                    dview.wait(results)                            
                    for j in range(number_of_cpu):
                        if results[j][0]!="":
                            out_file.write(results[j][0])
            client.close()
            gc.collect()                
                
                
            for j in range(len(filenames_dna_by_sublist[i])):
                chr=str(filenames_dna_by_sublist[i][j].split("_")[2].split(".tab")[0])
                if chr not in chrms_flat:
                    print "Warning: parse_chrms_data_vs_annotated_coordinates unknown chromosome ",                
                    print chr
                    return
                
                print chr ,
                print filenames_dna_by_sublist[i][j]
                df_starts=np.array(df[chr][:,0],  dtype='int')
                df_ends=np.array(df[chr][:,1], dtype='int')
                df_strands=[] 
                dview.targets = [j] 
                
                
                
                lines[j]=[]
                with open (os.path.join(path , filenames_dna_by_sublist[i][j]), 'r') as inLines:
                    for k, l in enumerate (inLines):
                        lines[j].append(l.rstrip().split("\t"))  
                print "Reading is done"
                                    
                print "Sending to core is done"
                results[j] = dview.apply(core_parse_chrms_data_vs_annotated_coordinates, 
                                         lines[j],
                                         os.path.join(path , filenames_dna_by_sublist[i][j].split(".tab")[0] + "_" + df_name + ".tab"),
                                         "DNA",
                                         df_starts, df_ends, df_strands)
                
            print '\nWaiting dnas...\n'
            dview.wait(results)
            for j in range(len(filenames_dna_by_sublist[i])):
                print results[j][0]
            client.close()
            gc.collect()
            
            #
            # -----------------------------------------RNA--------------------------------------------------
            #
            
            client = ipyparallel.Client(profile='natalia') 
            dview = client[:]
            print "clients ids ",
            print client.ids
            dview.block = False        
            print filenames_rna_by_sublist[i] 
            lines=[None]*number_of_cpu
            for j in range(len(filenames_rna_by_sublist[i])):
                chr=filenames_rna_by_sublist[i][j].split("_")[2].split(".tab")[0]
                if chr not in chrms_flat:
                    print "Warning: parse_chrms_data_vs_annotated_coordinates unknow chromosome ",
                    print chr
                    return         
                print chr ,
                df_starts=np.array(df[chr][:,0],  dtype='int')
                df_ends=np.array(df[chr][:,1], dtype='int')
                df_strands=np.array(df[chr][:,2], dtype='str')                 
                dview.targets = [j]   
                lines[j]=[]
                with open (os.path.join(path , filenames_rna_by_sublist[i][j]), 'r') as inLines:
                    for k, l in enumerate (inLines):
                        lines[j].append(l.rstrip().split("\t"))
                dview[j]['lines']=lines
                
                results[j] = dview.apply(core_parse_chrms_data_vs_annotated_coordinates, 
                                         #lines[j],
                                         os.path.join(filenames_rna_by_sublist[i][j].split(".tab")[0] + "_" + df_name + ".tab"),
                                         "RNA",
                                         df_starts, df_ends, df_strands)   
                gc.collect()
            print '\nWaiting rnas...\n'
            dview.wait(results)  
            for j in range(len(filenames_rna_by_sublist[i])):
                print results[j][0]
            client.close()
    
        
    #
    
    
    '''
    
    '''
    pattern_dna = re.compile("^k562\_dna\_chr\d*\.tab$")
    pattern_rna = re.compile("^k562\_rna\_chr\d*\.tab$")
    filenames_dna = [j for j in os.listdir(folder) if pattern_dna.match(j)]
    filenames_rna = [j for j in os.listdir(folder) if pattern_rna.match(j)]
    '''      
    
    '''
    def parse_rna_chrms_data_vs_annotated_coordinates(fl_in, fl_out, df, ch):
        print fl_in
        print fl_out
        string_count=0
        with open(fl_in, 'r') as fl_raw:
            for i, l in enumerate(fl_raw):
                ln=l.rstrip().split("\t")
                string_count+=1
                if ln[0]!=ch:
                    print "Warning: parse_rna_chrms_data_vs_annotated_coordinates bad input file "
                    return           
        print string_count
                
        df_starts=np.array(df[ch][:,0],  dtype='int')
        df_ends=np.array(df[ch][:,1], dtype='int')
        df_strands=np.array(df[ch][:,2], dtype='str')    
        
        print df_starts
        print df_ends
        print df_strands
        print df[ch]
        with open(fl_in, 'r') as fl_raw,\
             open (fl_out, 'w+') as fl_proc:
            for i, l in enumerate(fl_raw): 
                if i % 50000==0: 
                    print "\r" + str(i),
                ln=l.rstrip().split("\t")
                coord1=int(ln[1])
                chr=ln[0]            
                rna_strand=ln[8]   
                if np.any((df_starts <= coord1) and (df_ends >= coord1) and (df_strands==rna_strand)):
                    fl_proc.write(l)
                    
    
    
    
    #
    #
    def parse_dna_chrms_data_vs_annotated_coordinates(fl_in, fl_out, df, ch):
        print fl_in , 
        print fl_out
        string_count=0
        with gzip.open(fl_in, 'r') as fl_raw:
            for i, l in enumerate(fl_raw):
                ln=l.rstrip().split("\t")
                string_count+=1
                if ln[1]!=ch:
                    print "Warning: parse_rna_chrms_data_vs_annotated_coordinates bad input file "
                    return           
        print string_count
                
        df_starts=np.array(df[ch][:,0],  dtype='int')
        df_ends=np.array(df[ch][:,1], dtype='int')
        
        print df_starts
        print df_ends
        print df[ch]
        
        with gzip.open(fl_in, 'r') as fl_raw,\
             gzip.open (fl_out, 'wt') as fl_proc:
            for i, l in enumerate(fl_raw): 
                ln=l.rstrip().split("\t") 
                if i % 50000==0: 
                    print "\r" + str(i),
                ln=l.rstrip().split("\t")
                coord1=int(ln[2])
                chr=ln[1]           
    
                if np.any((df_starts <= coord1) and (df_ends >= coord1)):
                    fl_proc.write(l)
                    
                    
                    
                    
                    
        for x in x_ax:
        y_ax_intra_tmp=[]
        for gene in dist_intra.keys():
            
            #dist_by_bins_intra[id][d][binRNA]=+1
            
            length=abs(binId(genes_length[gene][0]) - binId(genes_length[gene][1]))
            if x<0:
                norm_value_num_of_pos=abs(x)
            else:
                norm_value_num_of_pos=length + 1 - x
                
            if x in dist_intra[gene]:     
                y_ax_intra_tmp.append((dist_intra[gene][x]*1.0) / (norm_value[gene] * (norm_value_num_of_pos)))  
                
                values_intra_by_gene.write(str(x) + "\t" + str(gene) + "\t" 
                                   + str((dist_intra[gene][x]*1.0) / (norm_value[gene] * (norm_value_num_of_pos)) ) + "\n")                
                
                values_raw_intra.write(str(x) + "\t" + str(gene) + "\t" + str(binId(dist_intra[gene][x])) 
                                       + "\t" + str(norm_value[gene]) + "\t" + str(norm_value_num_of_pos) + "\n")
            else:
                y_ax_intra_tmp.append(None)  
                
                values_intra_by_gene.write(str(x) + "\t" + str(gene) + "\t" + str(0) + "\n")
                
                values_raw_intra.write(str(x) + "\t" + str(gene) + "\t" + str(0) + "\t" + str(norm_value[gene]) 
                                       + "\t" + str(norm_value_num_of_pos) + "\n")
        tmp_list=[]
        
'''     
    
    
    

  







#
#

'''
bin_size=1000000
Workflow_chromosomes()
'''
#bin_size=100000
#Workflow_chromosomes()


bin_size=1000
Workflow_genes()

'''
bin_size=100
Workflow_genes()
'''

'''
bin_size=10000
Workflow_genes()

'''



