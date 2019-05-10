# -*- coding: utf-8 -*- 

#import ipyparallel
import collections
import os
#import zmq as pyz
import gc
import numpy as np
#import matplotlib.pyplot as plt
#plt.switch_backend('agg')
import matplotlib.ticker as ticker
from collections import Counter

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
mpl.rcParams['agg.path.chunksize'] = 30000
plt.switch_backend('agg')

#def init_names():
    
#os.system("ipcluster stop --profile=natalia &")
#os.system("ipcluster start --profile=natalia -n 4 &")

path = ""
path2 = ""
path3 = ""


path3_febr = ""
paths_bg_febr = ""
paths_rnas_febr = ""

lines = []
paths_rnas = {}
rnas_init = []

paths_bg = {}
bg_list = []

binId = lambda : None
bin_size=-1
bin_sizes = {}
chrms_stat = {}
chrms = []
chrms_flat = []
chrms_length = {}
chrms_length_file_name = ""
fibr_bg_stat = ""
k562_bg_stat = ""
    
fibr = "fibr.annot.full.tab"
fibr_mRNA = "fibr.annot.mRNA.tab"
raw_firb = "fibr.dna.back.bed"
fibr_overlap_mRNA_raw = "fibr.annot.mRNA_raw_overlap.tab"

raw_fibr_binned = ""
fibr_norm_genome_mean = ""
fibr_norm_genome_mean_smoothed = ""

k562 = "K562.annot.full.tab"
k562_mRNA = "K562.annot.mRNA.tab"

raw_k562 = "K562.dna.back.bed"
k562_overlap_mRNA_raw = "K562.annot.mRNA_raw_overlap.tab"

raw_k562_binned = ""
raw_k562_binned_grid = ""

k562_norm_genome_mean = ""
k562_norm_genome_mean_smoothed = ""



k562_norm_grid = path3 + path2 + "K562.dna.back_binned_norm_grid.tab"
k562_bin_test_grid = path3 + path2 + "K562.dna.back_bin_test_grid.tab"
k562_number_of_distinct_rnas = path3 + path2 + "K562.dna.back_number_of_distinct_rnas_grid.tab"
k562_chrm_sum_per_bin = path3+ path2 + "K562.dna.back_chrm_sum_per_bin_grid.tab"
k562_mean_per_bin = path3 + path2 + "K562.dna.back_mean_per_bin_grid.tab"
k562_ratio_value_sum = path3 + path2 + "K562.dna.back_ratio_value_sum_grid.tab"
k562_norm_smoothed_grid = path3 + path2 + "K562.dna.back_binned_norm_smoothed_grid.tab"

k562_norm_chrm_mean = path + path3 + "K562.dna.back_binned_norm_chrm_mean.tab"
k562_norm_no = path + path3 + "K562.dna.back_binned_norm_no.tab"
k562_norm_chrm_mean_smoothed = path3 + "K562.dna.back_binned_norm_chrm_mean_smoothed.tab"
k562_norm_no_smoothed = path3 + "K562.dna.back_binned_norm_no_smoothed.tab"


raw_fibr_binned_grid = ""
fibr_norm_grid = ""
fibr_norm_smoothed_grid = ""

fibr_norm_grid_preview=path3+ path2 + "fibr.dna.back_binned_norm_grid_preview.tab"
raw_fibr_binned_grid_preview=path3+ path2+"fibr.dna.back_binned_grid_preview.tab"
fibr_bin_test_grid=path3+ path2+ "fibr.dna.back_bin_test_grid.tab"
fibr_number_of_distinct_rnas=path3+ path2+ "fibr.dna.back_number_of_distinct_rnas_grid.tab"
fibr_chrm_sum_per_bin=path3+ path2+ "fibr.dna.back_chrm_sum_per_bin_grid.tab"
fibr_mean_per_bin=path3+ path2+ "fibr.dna.back_mean_per_bin_grid.tab"
fibr_ratio_value_sum=path3+ path2+ "fibr.dna.back_ratio_value_sum_grid.tab"

raw_firb_preview = path3 + "fibr.dna.back_preview.bed"
fibr_overlap_mRNA_raw_preview = path3 + "fibr.annot.mRNA_raw_overlap_preview.tab"
raw_fibr_binned_preview = path + path3 + "fibr.dna.back_binned_preview.tab"
fibr_norm_chrm_mean = path + path3 + "fibr.dna.back_binned_norm_chrm_mean.tab"
fibr_norm_no = path + path3 + "fibr.dna.back_binned_norm_no.tab"
fibr_norm_genome_mean_preview = path + path3 + "fibr.dna.back_binned_norm_genome_mean_preview.tab"
fibr_norm_chrm_mean_preview = path+path3+ "fibr.dna.back_binned_norm_chrm_mean_preview.tab"
fibr_norm_no_preview = path + path3 + "fibr.dna.back_binned_norm_no_preview.tab"
fibr_norm_chrm_mean_smoothed = path3 + "fibr.dna.back_binned_norm_chrm_mean_smoothed.tab"
fibr_norm_no_smoothed = path3 + "fibr.dna.back_binned_norm_no_smoothed.tab"
fibr_norm_genome_mean_smoothed_preview = path3 + "fibr.dna.back_binned_norm_genome_mean_smoothed_preview.tab"
fibr_norm_chrm_mean_smoothed_preview = path3 + "fibr.dna.back_binned_norm_chrm_mean_smoothed_preview.tab"
fibr_norm_no_smoothed_preview = path3 + "fibr.dna.back_binned_norm_no_smoothed_preview.tab"



k562_cut_Alexey="K562.annot.full_cut_Alexey.tab"
k562_cut_Alexey_preview="K562.annot.full_cut_Alexey_preview.tab" 





def init(bin_size_local):
    global path, path2, path3
    global bins, bin_size, binId, bin_sizes 
    global fibr_bg_stat, k562_bg_stat, chrms_length, chrms_length_file_name
    global lines, paths_rnas, rnas_init
    global paths_bg, bg_list
    global chrms_stat, chrms, chrms_flat
    
    global raw_fibr_binned, fibr_norm_genome_mean, fibr_norm_genome_mean_smoothed
    global raw_fibr_binned_grid, fibr_norm_grid, fibr_norm_smoothed_grid 
    global raw_k562_binned, raw_k562_binned_grid 
    global k562_norm_genome_mean, k562_norm_genome_mean_smoothed
    global path3_febr, paths_bg_febr, paths_rnas_febr
        
    #path="/home/n/Documents/work/RNA-DNA/"
    path = ""
    path2 = ""
    #path2="grid_genuine/"
    #path2="grid_low_bins_filtered_out/"
    #path2="grid_cut_bin_number/"    
    #path3="100kb/"
    path3="Alexey2/"
    
    bins = collections.defaultdict(list)
    bin_size = bin_size_local
    binId = lambda x: int(x/bin_size)    
    
    fibr_bg_stat = path + path3 + "stat_fibr_bg_" + str(bin_size)
    k562_bg_stat = path + path3 + "stat_K562_bg_" + str(bin_size)    
    
    chrms_length_file_name = path + "hg19_chrms_length"
    chrms_length = {}
    with open(chrms_length_file_name, 'r') as chrms_length_file:
        for i, l in enumerate(chrms_length_file):
            ln = l.rstrip().split("\t")
            chrms_length[ln[0]] = int(ln[1])
    #print chrms_length
    lines = ["fibr", "K562"]    
    paths_rnas = {}
    rnas_init = ["GAPDH", "XIST", "FIRRE", "RUNX1", "AGAP1", 
                 "KCNQ1OT1", "MALAT1" , "AC016205.1"]
    for line in lines:
        for rna in rnas_init:
            paths_rnas[line + "_" + rna] = path + line + "_" + rna \
                + ".full.tab"
            paths_rnas[line + "_" + rna + "_binned"] = path + path3 \
                + line + "_" + rna + ".full_binned.tab"
            paths_rnas[line + "_" + rna + "_binned_normalized"] = path \
                + path3 + line + "_" + rna + ".full_binned_normalized.tab"
            paths_rnas[line + "_" + rna + "_binned_normalized_fold"] = path \
                + path3 + line + "_" + rna + ".full_binned_normalized_fold.tab"
            paths_rnas[line + "_" + rna \
                       + "_binned_normalized_fold_filtered"] = path \
                + path3 + line + "_" \
                + rna + ".full_binned_normalized_fold_filtered.tab"
            paths_rnas[line + "_" + \
                       rna + "_binned_normalized_fold_filtered_smoothed"] = \
                path + path3 + line + "_" + rna \
                + ".full_binned_normalized_fold_filtered_smoothed.tab"            
            paths_rnas[line + "_" + rna + "_binned_normalized_grid"] = \
                path2 + path3 + line + "_" + rna \
                + ".full_binned_normalized_grid.tab"
            paths_rnas[line + "_" + rna + "_binned_normalized_fold_grid"] = \
                path3 + path2 + line + "_" + rna \
                + ".full_binned_normalized_fold_grid.tab"
            paths_rnas[line + "_" + rna \
                       + "_binned_normalized_fold_filtered_grid"] = \
                path3 + path2 + line + "_" + rna \
                + ".full_binned_normalized_fold_filtered_grid.tab"
            paths_rnas[line + "_" + rna \ 
                       + "_binned_normalized_fold_filtered_smoothed_grid"] = \
                path3 + path2 + line + "_" + rna \ 
                + ".full_binned_normalized_fold_filtered_smoothed_grid.tab"
     
    paths_bg = {}    
    bg_list = ["fibr", "K562"]
    
    for bg in bg_list:            
        paths_bg["raw_" + bg + "_binned"] = path + path3 + bg \
            + ".dna.back_binned" + str(bin_size) + ".tab"
        paths_bg[bg + "_norm_genome_mean"] = path + path3 + bg \
            + ".dna.back_binned" + str(bin_size) + "_norm_genome_mean.tab"
        paths_bg[bg + "_norm_genome_mean_smoothed"] = path + path3 + bg \
            + ".dna.back_binned" + str(bin_size) + "_norm_genome_mean_smoothed.tab"
    
        paths_bg["raw_" + bg + "_binned_grid"] = path + path3 + bg \
            + ".dna.back_binned" + str(bin_size) + "_grid.tab"
        paths_bg[bg + "_norm_grid"] = path + path3 + path2 + bg \
            + ".dna.back_binned" + str(bin_size) + "_norm_grid.tab"
        paths_bg[bg + "_norm_smoothed_grid"] = path + path3 + path2 + bg \
            + ".dna.back_binned" + str(bin_size) + "_norm_smoothed_grid.tab"
    
        
    paths_bg_febr = {}
    bg_list_febr = ["Dawns", "adneur", "neurSC", "Hela_DRB", 
                    "Hela_G1", "Hela_M", "K562", "fibr"]
    
    bin_sizes = {"Dawns" : [50000], "adneur" : [50000], "neurSC" : [50000], 
                 "Hela_DRB" : [250000], "Hela_G1" : [250000], 
                 "Hela_M" : [250000], "K562" : [100000, 50000, 20000] , 
                 "fibr" : [100000, 50000, 20000]}
    
    path3_febr = ""    
    path4_febr = "calculations_bg/"
    
    path5_febr = "calculations_rnas/"
    
    cell_lines_febr = ["Dawns", "adneur", "neurSC", 
                       "Hela_DRB", "Hela_G1", 
                       "Hela_M", "K562", "fibr"]   
    paths_rnas_febr = {}
    
    rnas_init_febr = ["MIR3648", "MIR3687", 
                      "Xrna_K562_20585", "MIR3687Xrna", 
                      "SSU-rRNA_Hsa_RepM_10963"]
    
    for bg in bg_list_febr:
        for bn in bin_sizes[bg]:
            paths_bg_febr[bg + "_annot_table"] = path + path3_febr + bg \
                + ".annot.full.tab"        
            paths_bg_febr["raw_" + bg] = path + path3_febr + bg \
                + ".back.bed"
            paths_bg_febr["raw_" + bg + "_binned"] = path + path3_febr \
                + path4_febr + bg \
                + ".dna.back_binned_" + str(bn) + ".tab"
            paths_bg_febr[bg + "_binned_norm_genome_mean"] = path + \
                path3_febr + path4_febr + bg + ".dna.back_binned_" \
                + str(bn) + "_norm_genome_mean.tab"
            paths_bg_febr[bg + "_binned_norm_genome_mean_smoothed"] = path \
                + path3_febr + path4_febr + bg + ".dna.back_binned_" \
                + str(bn) + "_norm_genome_mean_smoothed.tab"
            paths_bg_febr[bg + "_bg_stat"] = path + path3_febr + \
                "stat_" + bg + "_bg_" + str(bn)
    
            paths_bg_febr[bg + "_binned_norm_chrm_mean"] = path + \
                path3_febr + path4_febr + bg + ".dna.back_binned_" \
                + str(bn) + "_norm_chrm_mean.tab"
            paths_bg_febr[bg + "_binned_norm_no"] = path + path3_febr \
                + path4_febr + bg + ".dna.back_binned_" \
                + str(bn) + "_norm_no.tab"
            
    

    for line in cell_lines_febr:
        for rna in rnas_init_febr:
            for bn in bin_sizes[line]:
                paths_rnas_febr[line + "_" + rna] = path + path3_febr + \
                    path5_febr + line + "_" + rna + ".full_" + str(bn) + ".tab"
                paths_rnas_febr[line + "_" + rna + "_raw_DNA_parts"] = \
                    path + path3_febr + path5_febr + line + "_" \
                    + rna + ".full_raw_DNA_parts_" + str(bn) + ".bed"
                
                paths_rnas_febr[line + "_" + rna + "_binned"] = path \
                    + path3_febr + path5_febr + line + "_" + rna \
                    + ".full_binned_" + str(bn) + ".tab"
                paths_rnas_febr[line + "_" + rna + "_binned_last_bin_cor"] = \
                    path + path3_febr + path5_febr + line + "_" + rna \
                    + ".full_binned_last_bin_cor_" + str(bn) + ".bedGraph"
                
                paths_rnas_febr[line + "_" + rna + "_binned_to_ckeck"] = \
                    path + path3_febr + path5_febr + line + "_" + rna \
                    + ".full_binned_to_ckeck_" + str(bn) + ".tab"            
                paths_rnas_febr[line + "_" + rna + "_binned_normalized"] = \
                    path + path3_febr + path5_febr + line + "_" + rna \
                    + ".full_binned_normalized_" + str(bn) + ".tab"
                paths_rnas_febr[line + "_" + rna \
                                + "_binned_normalized_fold"] = path \
                    + path3_febr + path5_febr + line +"_" + rna \
                    + ".full_binned_normalized_fold_" + str(bn) + ".tab"
                paths_rnas_febr[line + "_" + rna \
                                + "_binned_normalized_fold_filtered"] = \
                    path + path3_febr + path5_febr + line + "_" \
                    + rna + ".full_binned_normalized_fold_filtered_" \
                    + str(bn) + ".tab"
                paths_rnas_febr[line + "_" + rna \
                                + "_binned_normalized_fold_filtered_smoothed"] = \
                    path + path3_febr + path5_febr + line + "_" + rna \
                    + ".full_binned_normalized_fold_filtered_smoothed_" \
                    + str(bn) + ".tab"
                paths_rnas_febr[line + "_" + rna \
                                + "_binned_normalized_fold_filtered_smoothed_last_bin_cor"] = \
                    path + path3_febr + path5_febr + line + "_" + rna \
                    + ".full_binned_normalized_fold_filtered_smoothed_last_bin_cor_" \
                    + str(bn) + ".bedGraph"

    raw_fibr_binned = path + path3 + "fibr.dna.back_binned" \
        + str(bin_size) + ".tab"
    fibr_norm_genome_mean = path + path3 + "fibr.dna.back_binned" \
        + str(bin_size) + "_norm_genome_mean.tab"
    fibr_norm_genome_mean_smoothed = path + path3 + "fibr.dna.back_binned" \
        + str(bin_size) + "_norm_genome_mean_smoothed.tab"
    
    raw_fibr_binned_grid=path + path3 + "fibr.dna.back_binned" \
        + str(bin_size) + "_grid.tab"
    fibr_norm_grid=path2 + path3 + "fibr.dna.back_binned" \
        + str(bin_size) + "_norm_grid.tab"
    fibr_norm_smoothed_grid=path2 + path3 + "fibr.dna.back_binned" \
        + str(bin_size) + "_norm_smoothed_grid.tab"
    
    
    raw_k562_binned = path + path3 + "K562.dna.back_binned" \
        + str(bin_size) + ".tab"
    raw_k562_binned_grid = path + path3 + "K562.dna.back_binned" \
        + str(bin_size) + "_grid.tab"
    
    k562_norm_genome_mean = path + path3 + "K562.dna.back_binned" \
        + str(bin_size) + "_norm_genome_mean.tab"
    k562_norm_genome_mean_smoothed = path + path3 +  "K562.dna.back_binned" \
        + str(bin_size) + "_norm_genome_mean_smoothed.tab"
    
    
    

    
    
    chrms_stat={}
    chrms=[]
    chrms_flat=[]
    
    tmp=[]
    for j in range (0, 5):    
        tmp=[]
        for i in range (1, 5):
            tmp.append('chr' + str(j*4 + i))
        chrms.append(tmp)
    
    tmp=[]  
    
    tmp.append('chr21')
    tmp.append('chr22')
    tmp.append('chrX')
    tmp.append('chrY')
    chrms.append(tmp)
    
    chrms_flat=[chr for tmp in chrms for chr in tmp]
    
    #print chrms_flat
    #chrms_flat=["chr"+str(i) for i in range(1,23)]
    #chrms_flat.append("chrX")
    #chrms_flat.append("chrY")
    
    #print chrms    


#
#
def create_binnded_background_for_bed_files(track_name, bed_f, bg_f):
    print bed_f
    print bg_f
    print chrms
    bed = open(bed_f, 'r')
    bg = open(bg_f, 'w+')
    fname = track_name
    #print fname
    
    bg.write("track type=bedGraph name=\"rna_" \
             + fname + "\"description=\"rna track " 
             + fname + " binned\"\n")
     
    chrs = {}
 
    #chrms_flat = [item for sublist in chrms for item in sublist]
    #print chrms_flat
    for chr in chrms_flat:
        chrs[chr] = collections.defaultdict(list)
    for line in bed:
        ln = line.rstrip().split("\t")
        if len (ln)<3:
            print "Warning: create_binnded_background_for_bed_files \
            length of the input line is less than 4 elements"
            print ln
            continue
        chr = ln[0]
        
        if int(ln[1])==27738001 or int(ln[1])==27738479:
            continue
        bin1 = binId(int(ln[1]))
        
        bn = bin1
        if bn in chrs[chr]:
            chrs[chr][bn]+=1
        else:
            chrs[chr][bn]=1
            
        
    for chr in chrms_flat:
        bins = sorted (chrs[chr].keys())
        print chr + " ",
        print len(bins)
        for bn in bins:            
            bg.write(chr + "\t" + str(bn * bin_size) \
                     + "\t" + str(bn * bin_size+bin_size) \
                     + "\t" + str(chrs[chr][bn]) + "\n")
    
    bg.close()
    bed.close()


#
#
def create_binnded_background_from_annotation_table(track_name, 
                                                    annotation_table_file, 
                                                    background_file):
    print annotation_table_file
    print background_file
    print chrms
    chrs = {} 
    reads_count = 0
    annot = open(annotation_table_file, 'r')
    bg = open(background_file, 'w+')
    
    bg.write("track type=bedGraph name=\"bg_" \
             + track_name + "\"description=\"bg track " \
             + track_name + " binned\"\n")         

    for chr in chrms_flat:
        chrs[chr] = collections.defaultdict(list)
        
    for line in annot:
        ln = line.rstrip().split("\t")
        
        assert len (ln)==24, \
               "Warning: create_binnded_background_from_annot_tables \
               length of the input line is less than 24 elements %r" % ln
        
        if ln[11]!="protein_coding":
            continue
        
        if ln[7]==ln[13]:
            continue
        
        chr = ln[13]      
        bin1 = binId(int(ln[14]))        
        reads_count+=1        
        if bin1 in chrs[chr]:
            chrs[chr][bin1]+=1
        else:
            chrs[chr][bin1] = 1            
    
    print reads_count
    
    for chr in chrms_flat:
        print chr ,
        print len(chrs[chr].keys())
        for bn in sorted(chrs[chr].keys()):            
            bg.write(chr + "\t" + str(bn * bin_size) \
                     + "\t" + str(bn * bin_size+bin_size) 
                     + "\t" + str(chrs[chr][bn]) + "\n")
    
    bg.close()
    annot.close()

def create_binnded_background_from_annotation_table_old_format(track_name,                                                                 annotation_table_file,                                                                background_file):
    print annotation_table_file
    print background_file
    print chrms
    chrs = {} 
    reads_count = 0
    annot = open(annotation_table_file, 'r')
    bg = open(background_file, 'w+')
    
    bg.write("track type=bedGraph name=\"bg_" \
             + track_name + "\"description=\"bg track " \
             + track_name + " binned\"\n")         

    for chr in chrms_flat:
        chrs[chr] = collections.defaultdict(list)
        
    for line in annot:
        ln = line.rstrip().split("\t")
        
        assert len (ln)==24, "Warning: create_binnded_background_from_annot_tables \
        length of the input line is less than 24 elements %r" % ln
        
        if ln[11]!="protein_coding":
            continue
        
        if ln[7]==ln[13]:
            continue
        
        chr = ln[13]      
        bin1 = binId(int(ln[14]))        
        reads_count+=1        
        if bin1 in chrs[chr]:
            chrs[chr][bin1]+=1
        else:
            chrs[chr][bin1] = 1            
    
    print reads_count
    
    for chr in chrms_flat:
        print chr ,
        print len(chrs[chr].keys())
        for bn in sorted(chrs[chr].keys()):            
            bg.write(chr + "\t" + str(bn * bin_size)+ "\t" \
                     + str(bn * bin_size+bin_size) + "\t" \
                     + str(chrs[chr][bn]) + "\n")
    
    bg.close()
    annot.close()
    
    
#
#
def simple_normalization_for_bg(fl_in, fl_norm_genome_mean_out, 
                                fl_norm_chrm_mean_out, fl_norm_no_out, 
                                fl_bg_stat):
    bg=open(fl_in, "r")
    bg_norm_genome_mean=open(fl_norm_genome_mean_out,'w+')
    bg_norm_chrm_mean=open(fl_norm_chrm_mean_out,'w+')
    bg_norm_no=open(fl_norm_no_out,'w+')
    bg_stat=open(fl_bg_stat, 'w+')
    
    chr=""
    chr_size=0
    total=0.0
    total_truncated=0.0
    total_genome=0.0
    lines=[]
    norm_bg=[]
    number_of_bins=0
    
    chrms_max={}
    chrms_filled={}
    chrms_sum={}
    genome_mean_total=0.0
    genome_sum=0.0
    bins_total=0.0
    bins_filled=0.0
    
    for l in bg:
        ln=l.rstrip().split("\t")  
        if "bedGraph" in l:
            continue
                    
        if ln[0] in chrms_sum:
            #if int(ln[2]) > chrms_max[ln[0]]:
                #chrms_max[ln[0]]=int(ln[2]) 
            chrms_sum[ln[0]]+=int(ln[3]) 
            chrms_filled[ln[0]]+=1
        else:
            #chrms_max[ln[0]]=int(ln[2]) 
            chrms_sum[ln[0]]=int(ln[3])    
            chrms_filled[ln[0]]=1
        
        
        if ln[0] in chrms_length and ln[0] not in chrms_max:
            chrms_max[ln[0]]=chrms_length[ln[0]]
        elif ln[0] not in chrms_length:
            print "Warning: simple_normalization_for_bg chr not in chrm_length dict " + ln[0]
            
    #print chrms_length
    #print chrms_max
    
    print 'chr chrm_sum_signal chrms_filled_bins chrm_length chrm_mean' 
    for chr in sorted(chrms_sum.keys()):
        print chr + "\t" + str(chrms_sum[chr]) + "\t" + str (chrms_filled[chr]) \
              + "\t" + str (chrms_max[chr]) + "\t" \
              + str(((chrms_sum[chr] * 1.0)/(1.0 * binId(chrms_max[chr]))))
        bg_stat.write(chr + "\t" + str(chrms_sum[chr]) + "\t" \
                      + str(chrms_filled[chr]) + "\t" + str(chrms_max[chr]) \
                      + "\t" +  str((chrms_sum[chr] * 1.0)/(1.0 * binId(chrms_max[chr]))) + "\n")
        genome_sum+=chrms_sum[chr]
        bins_total+=binId(chrms_max[chr])
        bins_filled+=chrms_filled[chr]
    
    print 'total_bins_total ',
    print bins_total    
    bg_stat.write('bins_total\t' + str(bins_total)  + "\n")
    
    print 'total_bins_filled ',
    print bins_filled
    bg_stat.write('bins_filled\t' + str(bins_filled) + "\n" )
    
    print 'genome_mean bins total ',
    genome_mean_total = (genome_sum * 1.0)/(1.0 * bins_total)
    print genome_mean_total
    bg_stat.write('genome_mean bins total\t' + str(genome_mean_total) + "\n")
    
    print 'genome_mean bins filled ',
    genome_mean_filled = (genome_sum * 1.0)/(1.0 * bins_filled)
    print genome_mean_filled
    bg_stat.write('genome_mean bins filled\t' + str(genome_mean_filled) + "\n") 
    
    bg.close()
    bg = open(fl_in, "r")        
    prev_chrm = ""
    prev_end = 0
    for l in bg:         
        if "bedGraph" in l:
            continue
        ln = l.rstrip().split("\t") 
        bin1 = binId(int(ln[1]))
        bin2 = binId(int(ln[2]))
        
        if chr=="" or prev_chrm!=ln[0]:
            chr = ln[0]
            for i in range (1, bin1 + 1):
                bg_norm_genome_mean.write(ln[0] + "\t" + str( bin_size*(i-1)) + "\t" + str(bin_size*(i)) + "\t" + "0\n")
                bg_norm_chrm_mean.write(ln[0] + "\t" + str( bin_size*(i-1)) + "\t" + str(bin_size*(i)) + "\t" + "0\n")
                bg_norm_no.write(ln[0] + "\t" + str( bin_size*(i - 1)) + "\t" + str(bin_size * (i)) + "\t" + "0\n")
                
        if  ln[0]==chr :
            for i in range (1, 1 + bin1 - prev_end):
                bg_norm_genome_mean.write(ln[0] + "\t"+ str(prev_end * bin_size + bin_size*(i - 1)) \
                                          + "\t" + str(prev_end*bin_size+ bin_size*(i)) + "\t" + "0\n")
                bg_norm_chrm_mean.write(ln[0] + "\t" + str(prev_end * bin_size + bin_size * (i - 1)) \
                                        + "\t" + str(prev_end * bin_size + bin_size * (i)) + "\t" + "0\n")
                bg_norm_no.write(ln[0] + "\t" + str(prev_end * bin_size + bin_size * (i - 1)) \
                                 + "\t" + str(prev_end*bin_size+bin_size * (i)) + "\t" + "0\n")
            for i in range(1, (1 + bin2 - bin1)):
                bg_norm_genome_mean.write(ln[0] + "\t" + str(bin1 * bin_size + bin_size * (i - 1)) \
                                          + "\t" + str(bin1 * bin_size + bin_size * (i)) \
                                          + "\t" + str(float(ln[3])/genome_mean_total) + "\n")           
                bg_norm_chrm_mean.write(ln[0] + "\t" + str(bin1 * bin_size + bin_size * (i - 1)) \
                                        + "\t"+  str(bin1 * bin_size + bin_size * (i)) \
                                        + "\t" + str(float(ln[3])/((chrms_sum[chr] * 1.0)/(1.0 * binId(chrms_max[chr])))) + "\n") 
                bg_norm_no.write(ln[0] + "\t" + str(bin1 * bin_size + bin_size * (i - 1)) \
                                 + "\t" + str(bin1 * bin_size + bin_size * (i)) + "\t" + str(ln[3]) + "\n") 
        prev_end = bin2
        prev_chrm = ln[0]         

        
    bg.close()
    bg_norm_genome_mean.close()
    bg_norm_chrm_mean.close()
    bg_norm_no.close() 
    bg_stat.close()


#
#filter full tab files for only protein coding RNA parts
#
def filter_mRNA_from_full_annot_tab(fl_bg_tab, fl_bg_mRNA_tab):
    tab_mRNA_fl = open(fl_bg_mRNA_tab, 'w+')
    with open(fl_bg_tab, 'r') as tab_fl:
            for i, l in enumerate(tab_fl):
                if 'protein_coding' in l:
                    ln = l.rstrip().split("\t")
                    #tab[ln[4]]=[ln[9],ln[10],ln[-1]]
                    tab_mRNA_fl.write(ln[4] + "\t" + ln[9] + "\t" + ln[10] + "\t" + ln[-1] + "\n")   
                if i % 200000==0:
                    print '{0}'.format(str(i)), 
    print i  
    tab_mRNA_fl.close()
#
#
def create_bed_file_for_dna_part_for_rna(fl_in, fl_out):
    bed = open(fl_out, 'w+')
    with open(fl_in, 'r') as rna:
        for i, l in enumerate(rna):
            if i==0:
                bed.write("track type=bed name=\"" 
                          + fl_out.split(".bed")[0] 
                          + "\"description=\"" 
                          + fl_out.split(".bed")[0] + "\"\n")                                            
            ln = l.rstrip().split("\t")        
            bed.write(ln[11] + "\t" + ln[12]+ "\t" + ln[13] + "\n")
    bed.close()
#
#
def filter_rna_new_format(file_in, file_out, mask, name):
    fl_in = open (file_in, 'r')
    fl_out = open (file_out,'w+')
    if 'region' in mask:
        mk = mask.split("_")        
        chr = mk[1]
        start = int(mk[2])
        end = int(mk[3])
    else:
        mk = []
        start = -1
        end = -1
    print "filter_rna : " ,    
    print "mask = " ,
    print mk
    for ln in fl_in:
        lna = ln.rstrip().split("\t")
        
        if lna[11]=='chr19' and int(ln[12])==27738479:            
            continue
        
        if mask=='name' and lna[12]==name :
            fl_out.write(ln)
        if mask == 'pc' and 'protein_coding' in ln:
            fl_out.write(ln)
        if mask == 'pc,trans' and 'protein_coding' in ln and lna[7]!=lna[13]:
            fl_out.write(ln)  
        if mask=='all,trans' and lna[7]!=lna[13]:
            fl_out.write(ln)
        if 'region' in mask:
            if lna[0]==chr and int(lna[1])>=start and int(lna[2])<=end:
                fl_out.write(ln)
            
    fl_in.close()
    fl_out.close()
#
#
#filter annotation tabel files with mask
#
def filter_rna(file_in, file_out, mask, name):
    fl_in = open (file_in, 'r')
    fl_out = open (file_out,'w+')
    if 'region' in mask:
        mk = mask.split("_")
        
        chr = mk[1]
        start = int(mk[2])
        end = int(mk[3])
    else:
        mk = []
        start=-1
        end=-1
    print "filter_rna : " ,    
    print "mask = " ,
    print mk
    for ln in fl_in:
        lna = ln.rstrip().split("\t")
        print lna ,
        print lna[12]
        
        if lna[11]=='chr19' and int(ln[12])==27738479:            
            continue
        
        if mask=='name' and lna[10]==name :
            fl_out.write(ln)
        if mask == 'pc' and 'protein_coding' in ln:
            fl_out.write(ln)
        if mask == 'pc,trans' and 'protein_coding' in ln and lna[5]!=lna[11]:
            fl_out.write(ln)  
        if mask=='all,trans' and lna[5]!=lna[11]:
            fl_out.write(ln)
        if 'region' in mask:
            if lna[0]==chr and int(lna[1])>=start and int(lna[2])<=end:
                fl_out.write(ln)
            
    fl_in.close()
    fl_out.close()
#
#check if tab files with protein coding parts contain all reads in background files created by Nastya in July
#
def check_overlap_for_reads(fl_bg_bed,fl_bg_mRNA_tab,fl_overlap):    
    overlap = open(fl_overlap,'w+')
    tab = {}
    lines = []
    with open(fl_bg_mRNA_tab, 'r') as tab_mRNA_fl:
        print 'Hi'
        for i, l in enumerate(tab_mRNA_fl):
            tab[l[:l.index("\t")]] = l.rstrip()
            lines.append(l.rstrip())
 
            if i % 50000 ==0:
                print '{0}'.format(str(i)), 
                del lines[:]
                del lines                
                gc.collect()
                lines = []
    print i
            
    found = 0
    not_found = 0
    
    with open(fl_bg_bed, 'r') as bed_fl:
        print 'Hi2'
        for i, l in enumerate(bed_fl):
            ln = l.rstrip().split("\t")    
            if ln[3] in tab :
                overlap.write (l.rstrip() + "\t" + tab[ln[3]] + "\n")
                found+=1
            else:
                not_found+=1
            if i % 200000 ==0:
                print '{0}'.format(str(i)),                
            

    print found
    print not_found
    overlap.close()
    




#
def core_smoth_with_window(lines_list, chr, smoothed):
    lines = np.asarray(lines_list,dtype=object)  
    #print chr
    #print lines[0:11,3]
    proximity = [ float(x) for x in lines[0:11,3] ] 
    
    for j in range(len(lines)):
        if j<=5:
            mn = np.mean(proximity)
        elif j>len(lines) - 5 - 1:
            mn = np.mean(proximity)
        else:    
            del proximity[0]
            proximity.append(float(lines[j+5][3]))    
            mn = np.mean(proximity)                            
            if chr=='chr10' and lines[j][1]=='179500':
                print proximity
                print np.mean(proximity)
        smoothed.write(lines[j][0] + "\t" + lines[j][1] + "\t" + lines[j][2] + "\t" +str(mn) + "\n")   
    del lines
    
#smooth files with N-bins window. Files contain all expanded beans
#
def smooth_with_window(fl_bg, fl_bg_smoothed):
    smoothed = open(fl_bg_smoothed, 'w+')
    lines_list = []
    chr = ""
    with open(fl_bg, 'r') as bg:
        for i, l in enumerate(bg):
            ln = l.rstrip().split("\t")
            if chr=="":
                chr = ln[0]
                lines_list.append(ln)
            elif ln[0]==chr:
                lines_list.append(ln)
            else:
                core_smoth_with_window(lines_list,chr, smoothed)                 
                chr = ln[0]
                del lines_list                    
                gc.collect()                    
                lines_list = []
                lines_list.append(ln)
        core_smoth_with_window(lines_list,chr, smoothed) 
           
    #print i  
    smoothed.close()    






#
#
def core_smoth_with_window_rna(lines_list, chr, smoothed):
    lines = np.asarray(lines_list, dtype=object)  
    checked_bins = []
    
    for i, ln in enumerate(lines):
        bin = binId(int(ln[1]))
        proximity1 = [0] * 11
        proximity2 = [0] * 11
        proximity3 = [0] * 11
        #if binId(int(ln[1]))>=5 and binId(int(ln[1])) < chrms_stat[ln[0]][1] -5:    
        for j in range(-5,6):
            if i+j>=0 and i+j<len(lines) and abs(binId(int(lines[i+j][1])) - binId(int(ln[1])))<=5:
                proximity1[5 + binId(int(lines[i + j][1])) - binId(int(ln[1]))] = (float(lines[i + j][3])) 
                proximity2[5 + binId(int(lines[i + j][1])) - binId(int(ln[1]))] = (float(lines[i + j][3])) 
                proximity3[5 + binId(int(lines[i + j][1])) - binId(int(ln[1]))] = (float(lines[i + j][3])) 
        mn2 = np.mean(proximity2)              
                                       
        tmp = {}                  
        for k in range(-1,-6, -1): 
            del proximity1[-1]
            proximity1 = [0.0] + proximity1
            if bin==36182:
                print proximity1
            mn1 = np.mean(proximity1)
            tmp[k] = mn1
            
        for k in range(-5,0, +1):            
            if bin + k not in checked_bins and bin+k>=0:
                smoothed.write(ln[0] + "\t" + str((bin + k) * bin_size) \
                               + "\t" + str((bin + k + 1) * bin_size) + "\t" + str(tmp[k]) + "\n")
                checked_bins.append(bin + k)          
        #
        #  
        if bin not in checked_bins:
            smoothed.write(ln[0] + "\t" + ln[1] + "" + "\t" + ln[2] + "\t" + str(mn2) + "\n")     
        checked_bins.append(bin)            
        #
        #          
        for k in range(1,6):
            del proximity3[0]
            found = 0
            for m in range(0, 11):
                #if i+m<len(lines) and bin+5<=binId(int(lines[i+m][1])) and bin+k+5>=binId(int(lines[i+m][1])):
                if i + m<len(lines) and bin + 5 + k==binId(int(lines[i + m][1])):
                    proximity3.append(float(lines[i + m][3])) 
                    found = 1   
                    break
            if found==0:
                proximity3.append(0.0)  
                
            if bin + k not in checked_bins and bin + k <= binId(chrms_length[ln[0]]):
                mn3 = np.mean(proximity3) 
                tmp_count = 0
                if ln[0]=='chr10' and ((bin+k)*bin_size==1700000 
                                       or (bin+k)*bin_size==700000 
                                       or (bin+k)*bin_size==5800000):
                    print (bin+k)*bin_size
                    print proximity3
                    tmp_count+=1
                    if tmp_count==3:
                        return
                smoothed.write(ln[0] + "\t" + str((bin+k)*bin_size) + "\t" + str((bin+k+1)*bin_size) + "\t" +str(mn3) + "\n")                
                checked_bins.append(bin+k)              
            
            if i+1 <len(lines) and bin+k+1 == binId(int(lines[i+1][1])) :
                break            
    del lines
    
   
#
#
def smooth_with_window_rna(fl_rna, fl_rna_smoothed):
    smoothed = open(fl_rna_smoothed, 'w+')    
    if ".tab" in fl_rna_smoothed:        
        track_name = fl_rna_smoothed.split("/")[-1].split(".tab")[0]
    elif ".bedGraph" in fl_rna_smoothed:
        track_name = fl_rna_smoothed.split("/")[-1].split(".bedGraph")[0]
    else:
        track_name = ""
        print "Warning: normalize_single_rna problem with track name"  
    smoothed.write("track type=bedGraph name=\"" + track_name + "\" description=\"" + track_name + "\"\n")                                             
       
    lines_list = []
    chr = ""
    checked_bins = []
    with open(fl_rna, 'r') as rna:
        for i, l in enumerate(rna):
            ln = l.rstrip().split("\t")
            if len(ln)<4:
                print "Warning: smooth_with_window_rna short line " + l
                #smoothed.write(l)
                continue
            if chr=="":
                chr = ln[0]
                lines_list.append(ln)
                
            elif ln[0]==chr:
                lines_list.append(ln)
            else:
                core_smoth_with_window_rna(lines_list,chr, smoothed)                 
                chr = ln[0]
                del lines_list                    
                gc.collect()                    
                lines_list = []
                lines_list.append(ln)
                checked_bins = []
        if len(lines_list)!=0:
            core_smoth_with_window_rna(lines_list,chr, smoothed)        
    smoothed.close()        
#
#
def load_stat(fl_in):
    global chrms_stat
    chrms_stat = {}
    
    with open(fl_in, 'r') as tab_fl:
        for i, l in enumerate(tab_fl):
            ln = l.rstrip().split("\t")
            if len(ln)>2:
                print ln
                #chrms_stat[ln[0]]  =  [int(ln[-2]), binId(int(ln[-2])), int(ln[2])]
                chrms_stat[ln[0]] = [chrms_length[ln[0]], binId(chrms_length[ln[0]]), int(ln[2])]
#
#
def calculate_single_RNA_coverage_new_format(fl_annot_tab, fl_out):
    total_coverage = 0
    rna = {}
    with open(fl_annot_tab, 'r') as tab_fl:
        for i, l in enumerate(tab_fl):
            ln = l.rstrip().split("\t")            
            bin1 = binId(int(ln[14]))
            bin2 = binId(int(ln[15]))            
            if bin2-bin1>1: 
                print "Warning: calculate_single_RNA_coverage very long RNA " + l            
            if ln[13] in rna:
                rna[ln[13]][bin1] = rna[ln[13]][bin1]+1 if bin1 in rna[ln[13]] else 1
                total_coverage+=1                
            else:
                rna[ln[13]] = {}
                rna[ln[13]][bin1] = 1
                total_coverage+=1
    print  fl_out ,  
    end = 0
    print ' total RNA coverage ',
    print total_coverage
    with open(fl_out, 'w+') as rna_bin:
        if ".tab" in fl_out:
            track_name = fl_out.split(".tab")[0]
        elif ".bedGraph" in fl_out:
            track_name = fl_out.split(".bedGraph")[0]
        else:
            track_name = ""
            print "Warning: normalize_single_rna problem with track name"  
        
        rna_bin.write("track type=bedGraph name=\"" + track_name + "\"description=\"" + track_name + "\"\n")                                            
        #rna_bin.write("track type=bedGraph name=\" " +fl_out.split(".bedGraph")[0] + " \"description=\" track " + fl_out.split(".bedGraph")[0] + " \"\n")                                            
      
        for chr in sorted(rna.keys()):
            for bin in sorted(rna[chr].keys()):
                if (bin+1)*bin_size > chrms_stat[chr][0] : 
                    print chr ,
                    print chrms_stat[chr][0] ,
                    print (bin+1)*bin_size
                    #end = chrms_stat[chr][0]
                    end = (bin+1)*bin_size
                else:
                    end = (bin+1)*bin_size 
                #rna_bin.write(chr+"\t"+str(bin*bin_size) +"\t"+ str((bin+1)*bin_size ) +"\t"+ str(rna[chr][bin]) + "\n")
                rna_bin.write(chr+"\t"+str(bin*bin_size) +"\t"+ str(end) +"\t"+ str(rna[chr][bin]) + "\n")
                
#
#
def calculate_single_RNA_coverage(fl_annot_tab, fl_out):
    total_coverage = 0
    rna = {}
    with open(fl_annot_tab, 'r') as tab_fl:
        for i, l in enumerate(tab_fl):
            ln = l.rstrip().split("\t")            
            bin1 = binId(int(ln[12]))
            bin2 = binId(int(ln[13]))            
            if bin2-bin1>1: 
                print "Warning: calculate_single_RNA_coverage very long RNA " + l            
            if ln[11] in rna:
                rna[ln[11]][bin1] = rna[ln[11]][bin1]+1 if bin1 in rna[ln[11]] else 1
                total_coverage+=1                
            else:
                rna[ln[11]] = {}
                rna[ln[11]][bin1] = 1
                total_coverage+=1
    print  fl_out ,  
    end = 0
    print ' total RNA coverage ',
    print total_coverage
    with open(fl_out, 'w+') as rna_bin:
        if ".tab" in fl_out:
            track_name = fl_out.split(".tab")[0]
        elif ".bedGraph" in fl_out:
            track_name = fl_out.split(".bedGraph")[0]
        else:
            track_name = ""
            print "Warning: normalize_single_rna problem with track name"  
        
        rna_bin.write("track type=bedGraph name=\"" + track_name + "\"description=\"" + track_name + "\"\n")                                            
        #rna_bin.write("track type=bedGraph name=\" " +fl_out.split(".bedGraph")[0] + " \"description=\" track " + fl_out.split(".bedGraph")[0] + " \"\n")                                            
      
        for chr in sorted(rna.keys()):
            for bin in sorted(rna[chr].keys()):
                if (bin+1)*bin_size > chrms_stat[chr][0] : 
                    print chr ,
                    print chrms_stat[chr][0] ,
                    print (bin+1)*bin_size
                    #end = chrms_stat[chr][0]
                    end = (bin+1)*bin_size
                else:
                    end = (bin+1)*bin_size 
                #rna_bin.write(chr+"\t"+str(bin*bin_size) +"\t"+ str((bin+1)*bin_size ) +"\t"+ str(rna[chr][bin]) + "\n")
                rna_bin.write(chr+"\t"+str(bin*bin_size) +"\t"+ str(end) +"\t"+ str(rna[chr][bin]) + "\n")
    

#
#
def normalize_single_rna(fl_rna_binned, fl_rna_norm):
    rna = {}
    total_coverage = []
    with open(fl_rna_binned, 'r') as rna_bin:
        for i, l in enumerate(rna_bin):
            ln = l.rstrip().split("\t")    
            if len(ln) <4:
                print "Warning: normalize_single_rna short line " + l
                continue
            if ln[0] not in rna: rna[ln[0]] = []
            rna[ln[0]].append(ln)
            total_coverage.append(float(ln[-1]))
    
    print 'total coverage ',
    print sum(total_coverage)    
    
    
    number_of_bins = sum(chrms_stat[chr][1] for chr in chrms_stat)
    
    #print chrms_stat
    mean_total = sum(total_coverage)/(number_of_bins*1.0)
    
    print 'mean by all bins in the genome',
    print mean_total
    
    
    print 'mean by only filled bins ',
    mean_filled = np.mean(total_coverage)    
    print mean_filled
    
    with open(fl_rna_norm, 'w+') as rna_norm:
        if ".tab" in fl_rna_norm:
            track_name = fl_rna_norm.split(".tab")[0]
        elif ".bedGraph" in fl_rna_norm:
            track_name = fl_rna_norm.split(".bedGraph")[0]
        else:
            track_name = ""
            print "Warning: normalize_single_rna problem with track name"  
        rna_norm.write("track type=bedGraph name=\"" + track_name + "\"description=\"" + track_name + "\"\n")                                            
        
        for chr in sorted(rna.keys()):
            for ln in rna[chr]:                
                rna_norm.write(chr+"\t"+str(int(ln[1])) +"\t"+ str(int(ln[2])) +"\t"+ str((float(ln[3])*1.0)/mean_total) + "\n")
    
    
    




#
def calculate_fold_enrichment_for_single_RNA(fl_bg, fl_rna,fl_out):
    fold_enrichment =  open(fl_out, 'w+')
    if ".tab" in fl_out:
        track_name = fl_out.split(".tab")[0]
    elif ".bedGraph" in fl_out:
        track_name = fl_out.split(".bedGraph")[0]
    else:
        track_name = ""
        print "Warning: normalize_single_rna problem with track name"  
    
    fold_enrichment.write("track type=bedGraph name=\"" +track_name + "\"description=\"" + track_name + "\"\n")                                             

    '''
    fold_enrichment.write("track type=bedGraph name=" )
    track_name = fl_out.split("_")[0]+"_"+fl_out.split("_")[1]
    fold_enrichment.write(track_name)
    fold_enrichment.write(" description = \"decription\" autoScale=on viewLimits=0.0:25.0\n")
    '''
    bg = {}
    eps = 0.0001
    total_coverage = []
    with open(fl_bg, 'r') as f_bg:
        for i, l in enumerate(f_bg):

            ln = l.rstrip().split("\t")         
            if ln[0] not in bg: bg[ln[0]] = {}
            bg[ln[0]][int(ln[1])] = float(ln[-1])
    rna = {}
    count_zeros = 0
    with open(fl_rna, 'r') as f_rna:
        for i, l in enumerate(f_rna):
            ln = l.rstrip().split("\t")    
            if len(ln)<4:
                print "Warning: calculate_fold_enrichment_for_single_RNA short line " + l
                continue            
            #bin1 =  int(ln[1]) if int(ln[1]) % bin_size ==0 else int(ln[1])-1
            bin1 = int(ln[1])
            if ln[0] not in rna: rna[ln[0]] = {}
            if bin1 not in rna[ln[0]] : rna[ln[0]][bin1] = []
            #rna[ln[0]][bin1].append(float(ln[-1]))
            if bin1 not in bg[ln[0]]: 
                print "Warning: calculate_fold_enrichment_for_single_RNA bin not in bg"
                print ln[0] + "\t" + str(bin1) 
                continue
            fold = (float(ln[-1]))/(bg[ln[0]][bin1]) if bg[ln[0]][bin1]!=0  else (-1) * float(ln[-1])
            if bg[ln[0]][bin1]==0:
                count_zeros+=1
            #rna[ln[0]][bin1].append(fold)
            fold_enrichment.write(ln[0] + "\t" + ln[1] + "\t" +ln[2] + "\t" +str(fold) + "\n")
    print 'zeros ',        
    print count_zeros
    fold_enrichment.close()
    

#
#
def filtering_enrichment_signal_for_RNA(fl_RNA,fl_out):
    filtered = open(fl_out,"w+")
   
    if ".tab" in fl_out:
        track_name = fl_out.split(".tab")[0]
    elif ".bedGraph" in fl_out:
        track_name = fl_out.split(".bedGraph")[0]
    else:
        track_name = ""
        print "Warning: normalize_single_rna problem with track name"  
    filtered.write("track type=bedGraph name=\"" +track_name + "\"description=\"" + track_name + "\"\n")                                             
    
    list_of_coordinates = {}
    with open(fl_RNA, 'r') as fl_RNA:
        for i, l in enumerate(fl_RNA): 
            ln = l.rstrip().split("\t")
            if len(ln)<4:                 
                print "Warning: filtering_enrichment_signal_for_RNA short line " + l
                continue            
            if ln[0] not in list_of_coordinates: list_of_coordinates[ln[0]] = []
            if abs(float(ln[3]))<2:
                continue
            list_of_coordinates[ln[0]].append(ln)
    for chr in sorted(list_of_coordinates.keys()):
        for i,ln in enumerate(list_of_coordinates[chr]):
            bin1 = binId(int(ln[1]))
            start = i-5 if i>4 else 0
            end = i+6 if i<len(list_of_coordinates[chr])-5 else len(list_of_coordinates[chr])
            count = 0
            for j in range(start,end):
                if abs(binId(int(list_of_coordinates[chr][j][1]))-bin1)<=5:
                    count+=1
            if count >=3:
                filtered.write(ln[0]+"\t"+ ln[1]+"\t"+ ln[2]+"\t"+ ln[3]+"\n")
    filtered.close()
       
         




#
def plot_RNA(fl_RNA,fl_bg):
    fig = plt.figure(figsize = (8, 4))
    fig, ax = plt.subplots(figsize = (8, 4))
    fname = fl_RNA.split(".tab")[0]
    x = []
    x_ax = []
    y_ax = [] 
    x_ax2 = []
    y_ax2 = [] 
    
    x_ticks = []
    x_minor_ticks = []
    y = []
    my_xticks = []
    chr_indent = {}
    
    chrms_data = {}
    prev_bin = 0
    prev_chr = ''
    bg_data = {}
    with open(fl_bg, 'r') as fl_bg:
        for i, l in enumerate(fl_bg): 
            ln = l.rstrip().split("\t")            
            if ln[0] not in bg_data: bg_data[ln[0]] = {}  
            bg_data[ln[0]][binId(int(ln[1]))] = float(ln[3])                    
         
        
    with open(fl_RNA, 'r') as fl_RNA:
        for i, l in enumerate(fl_RNA): 
            if i <1:
                continue            
            ln = l.rstrip().split("\t")
            if ln[0] not in chrms_data: chrms_data[ln[0]] = {}   
            chrms_data[ln[0]][binId(int(ln[1]))] = float(ln[3])                    
    
    for chr in sorted(chrms_data.keys()):
        bins = sorted(chrms_data[chr].keys())
        if bins[0]>1:
            chrms_data[chr][bins[0]-1] = 0
        for i,bin in enumerate(bins):
            if i>0 and i<len(bins)-1:
                if bins[i]-bins[i-1] >1: chrms_data[chr][bins[i]-1] = 0
                if bins[i+1]-bins[i] >1: chrms_data[chr][bins[i]+1] = 0    
          
       
    prev_indent = 0
    for i, chr in enumerate(chrms_flat):
        #print i, chr
        if (i+1) % 2 == 0:
            my_xticks.append("\n" +chr.split("chr")[1])
        else:
            my_xticks.append(chr.split("chr")[1])
        #print chrms_stat[chr]
        x_ticks.append(prev_indent+chrms_stat[chr][1])  
        x_minor_ticks.append(prev_indent+chrms_stat[chr][1]/2)
        
        if chr in chrms_data: 
            for bin in sorted(chrms_data[chr]):
                x_ax.append(prev_indent+bin)  
                y_ax.append(chrms_data[chr][bin])
        
        if chr in bg_data: 
            for bin in sorted(bg_data[chr]):
                x_ax2.append(prev_indent+bin)  
                y_ax2.append(bg_data[chr][bin])        
        
        if i==0:
            chr_indent[chr] = 0
        else:
            chr_indent[chr] = prev_indent
        prev_indent+=chrms_stat[chr][1]
        
    x = np.array(x_ax)  
    y = np.array(y_ax)
    
    
    x2 = np.array(x_ax2) 
    y2 = np.array(y_ax2)
    
    x_ticks_np = np.array(x_ticks)
    
    
    #plt.plot(x, y,linewidth = 0.5)
    #ax.fill_between(x, 0, y)
    
    ax.plot(x, y,linewidth = 0.1)
    #ax.plot(x, y,'r', markersize = 0.1)
    for i in range(len(x_ticks_np)):
        if (i+1)% 2 ==0:
            ax.axvspan(x_ticks_np[i-1], x_ticks_np[i], alpha = 0.2, color = 'grey')
    
    #ax.xaxis.set_major_locator(x_ticks)
    #ax.xaxis.set_minor_locator(x_minor_ticks)
    
    ax.xaxis.set_major_locator(ticker.FixedLocator(x_ticks))
    ax.xaxis.set_minor_locator(ticker.FixedLocator(x_minor_ticks))
    
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(my_xticks))    
    
    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)
        tick.label1.set_horizontalalignment('center')
        tick.label1.set_fontsize(7) 
    #plt.xticks(fontsize = 8, rotation = 0)
    #plt.xticks(x_ticks_np, my_xticks)
    #plt.ylim(ymin = 0)
    plt.savefig(path+fname+".png", dpi = 300, figsize = (8, 4)) 
    
    ax.plot(x2, y2, alpha = 0.3, color = 'red')
    plt.savefig(path+fname+"_with_bg.png", dpi = 300, figsize = (8, 4))
    plt.close()
       

#
#
def plot_RNA_signal(fl_RNA, fl_bg):
    fig = plt.figure(figsize = (8, 4))
    fig, ax = plt.subplots(figsize = (8, 4))
    fname = fl_RNA.split(".tab")[0] 
    
    wdth = 10
    x = []
    x_ax = []
    y_ax = [] 
    x_ax2 = []
    y_ax2 = [] 
    
    x_ticks = []
    x_minor_ticks = []
    y = []
    my_xticks = []
    chr_indent = {}
    
    chrms_data = {}
    prev_bin = 0
    prev_chr = ''
    bg_data = {}
    
    with open(fl_bg, 'r') as fl_bg:
        for i, l in enumerate(fl_bg): 
            ln = l.rstrip().split("\t")            
            if ln[0] not in bg_data: bg_data[ln[0]] = {}  
            bg_data[ln[0]][binId(int(ln[1]))] = float(ln[3])                    
         
        
    with open(fl_RNA, 'r') as fl_RNA:
        for i, l in enumerate(fl_RNA): 
            if i<1:
                continue            
            ln = l.rstrip().split("\t")
            if ln[0] not in chrms_data: chrms_data[ln[0]] = {}   
            chrms_data[ln[0]][binId(int(ln[1]))] = float(ln[3])                    
    
    '''
    for chr in sorted(chrms_data.keys()):
        bins = sorted(chrms_data[chr].keys())
        if bins[0]>1:
            chrms_data[chr][bins[0]-1] = 0
        for i, bin in enumerate(bins):
            if i>0 and i<len(bins)-1:
                if bins[i]-bins[i-1] >1: chrms_data[chr][bins[i]-1] = 0
                if bins[i+1]-bins[i] >1: chrms_data[chr][bins[i]+1] = 0    
    '''
    
    prev_indent = 0
    
    
    for i, chr in enumerate(chrms_flat):
        if (i+1) % 2 == 0:
            #my_xticks.append("\n" + chr.split("chr")[1])
            my_xticks.append(chr.split("chr")[1])
        else:
            my_xticks.append(chr.split("chr")[1])
            
        x_ticks.append(prev_indent + chrms_stat[chr][1] * wdth)  
        x_minor_ticks.append(prev_indent + chrms_stat[chr][1] * wdth/2)
        
        #for bin in sorted(chrms_data[chr]):
        #print chrms_stat[chr]        
        for bin in range(binId(chrms_stat[chr][0]) + 1):
            x_ax.append(prev_indent + bin * wdth)  
            if chr in chrms_data and bin in chrms_data[chr]:
                y_ax.append(chrms_data[chr][bin])
            else:
                y_ax.append(0)
        
        
        #for bin in sorted(bg_data[chr]):
        for bin in range(binId(chrms_stat[chr][0]) + 1):
            x_ax2.append(prev_indent + bin * wdth)  
            if chr in bg_data and bin in bg_data[chr]:
                y_ax2.append(bg_data[chr][bin])       
            else:
                y_ax2.append(0)
                
        
        if i==0:
            chr_indent[chr] = 0
        else:
            chr_indent[chr] = prev_indent
        prev_indent+=chrms_stat[chr][1] * wdth
        
    x = np.array(x_ax)  
    y = np.array(y_ax)
    
    
    x2 = np.array(x_ax2) 
    y2 = np.array(y_ax2)
    
    x_ticks_np = np.array(x_ticks)
    
    
    #plt.plot(x, y,linewidth = 0.5)
    #ax.fill_between(x, 0, y)
    
    #ax.plot(x, y, linewidth = 0.1)
    ax.bar(x, y, width = wdth*10, color  =  'black')
    ax.tick_params(axis = 'both', which = 'major', labelsize = 7)
    ax.tick_params(axis = 'both', which = 'minor', labelsize = 7)
    print x
    print y
    #ax.plot(x, y,'r', markersize=0.1)
    '''
    for i in range(len(x_ticks_np)):
        if (i+1) % 2==0:
            ax.axvspan(x_ticks_np[i-1], x_ticks_np[i], alpha = 0.2, color = 'white')
    '''
    #ax.xaxis.set_major_locator(x_ticks)
    #ax.xaxis.set_minor_locator(x_minor_ticks)
    
    ax.xaxis.set_major_locator(ticker.FixedLocator(x_ticks))
    ax.xaxis.set_minor_locator(ticker.FixedLocator(x_minor_ticks))
    
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(my_xticks))    
    
    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)
        tick.label1.set_horizontalalignment('center')
        tick.label1.set_fontsize(4) 
    #plt.xticks(fontsize = 8, rotation = 0)
    #plt.xticks(x_ticks_np, my_xticks)
    #plt.ylim(ymin = 0)
    
    
    
    plt.savefig(path + fname + ".png", dpi = 300, figsize = (8, 4)) 
    
    ax.plot(x2, y2, alpha = 0.3, color = 'red')
    plt.savefig(path + fname + "_with_bg.png", dpi = 300, figsize = (8, 4))
    plt.close()
#
#
def plot_bg_by_chromosome(fl_bg, chr_name):
    fig = plt.figure(figsize = (8, 4))
    fig, ax = plt.subplots(figsize = (8, 4))

    fname = fl_bg.split(".tab")[0] +"_"+chr_name
    x = []
    x_ax = []
    y_ax = [] 
    x_ax2 = []
    y_ax2 = [] 

    x_ticks = []
    x_minor_ticks = []
    y  =  []
    my_xticks=[]
    chr_indent = {}

    chrms_data = {}
    prev_bin = 0
    prev_chr = ''
    bg_data = {}
    with open(fl_bg, 'r') as fl_bg:
        for i, l in enumerate(fl_bg): 
            ln = l.rstrip().split("\t") 
            if ln[0]!=chr_name:
                continue
            if float (ln[3])>240000:
                print l
            if len (ln)<3:
                print "Warning: plot_bg length of the input line is less than 4 elements ",
                print fl_bg
                print l + " line number " + str(i)
                continue            
            if ln[0] not in bg_data: bg_data[ln[0]] = {}  

            bg_data[ln[0]][binId(int(ln[1]))] = float(ln[3])                    


    prev_indent = 0
    for i, chr in enumerate(chrms_flat):
        if chr!=chr_name:
            prev_indent = 0
            continue
        #print i, chr
        if (i+1) % 2 == 0:
            my_xticks.append("\n" +chr.split("chr")[1])
        else:
            my_xticks.append(chr.split("chr")[1])
        #print chrms_stat[chr]
        x_ticks.append(prev_indent+chrms_stat[chr][1])  
        x_minor_ticks.append(prev_indent+chrms_stat[chr][1]/2)



        if chr in bg_data: 
            for bin in sorted(bg_data[chr]):
                x_ax2.append(prev_indent+bin)  
                y_ax2.append(bg_data[chr][bin])        

        if i==0:
            chr_indent[chr] = 0
        else:
            chr_indent[chr] = prev_indent
        prev_indent+=chrms_stat[chr][1]

    x = np.array(x_ax)  
    y = np.array(y_ax)


    x2 = np.array(x_ax2) 
    y2 = np.array(y_ax2)

    x_ticks_np = np.array(x_ticks)


    #ax.fill_between(x2, 0, y2)
    ax.plot(x2, y2,'ro', markersize = 1)
    #ax.plot(x2, y2,linewidth = 0.5)
    for i in range(len(x_ticks_np)):
        if (i+1)% 2 ==0:
            ax.axvspan(x_ticks_np[i-1], x_ticks_np[i], alpha = 0.2, color = 'grey')


    ax.xaxis.set_major_locator(ticker.FixedLocator(x_ticks))
    ax.xaxis.set_minor_locator(ticker.FixedLocator(x_minor_ticks))

    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(my_xticks))    

    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)
        tick.label1.set_horizontalalignment('center')
        tick.label1.set_fontsize(7) 

    plt.savefig(path+fname+".png", dpi = 300, figsize = (8, 4)) 
    plt.close()

    
#
def plot_bg(fl_bg):
    fig = plt.figure(figsize = (8, 4))
    fig, ax = plt.subplots(figsize = (8, 4))
   
    fname = fl_bg.split(".tab")[0]
    x = []
    x_ax = []
    y_ax = [] 
    x_ax2 = []
    y_ax2 = [] 
    
    x_ticks = []
    x_minor_ticks = []
    y = []
    my_xticks = []
    chr_indent = {}
    
    chrms_data = {}
    prev_bin = 0
    prev_chr = ''
    bg_data = {}
    with open(fl_bg, 'r') as fl_bg:
        for i, l in enumerate(fl_bg): 
            ln = l.rstrip().split("\t")     
            if len (ln)<3:
                print "Warning: plot_bg length of the input line is less than 4 elements ",
                print fl_bg
                print l + " line number " + str(i)
                continue            
            if ln[0] not in bg_data: bg_data[ln[0]] = {}  
            
            bg_data[ln[0]][binId(int(ln[1]))] = float(ln[3])                    
         

    prev_indent = 0
    for i, chr in enumerate(chrms_flat):
        #print i, chr
        if (i+1) % 2 == 0:
            my_xticks.append("\n" +chr.split("chr")[1])
        else:
            my_xticks.append(chr.split("chr")[1])
        #print chrms_stat[chr]
        x_ticks.append(prev_indent+chrms_stat[chr][1])  
        x_minor_ticks.append(prev_indent+chrms_stat[chr][1]/2)
        

        
        if chr in bg_data: 
            for bin in sorted(bg_data[chr]):
                x_ax2.append(prev_indent+bin)  
                y_ax2.append(bg_data[chr][bin])        
        
        if i==0:
            chr_indent[chr] = 0
        else:
            chr_indent[chr] = prev_indent
        prev_indent+=chrms_stat[chr][1]
        
    x = np.array(x_ax)  
    y = np.array(y_ax)
    
    
    x2 = np.array(x_ax2) 
    y2 = np.array(y_ax2)
    
    x_ticks_np = np.array(x_ticks)
    
    
    #ax.fill_between(x2, 0, y2)
    ax.plot(x2, y2,'ro', markersize = 0.1)
    #ax.plot(x2, y2,linewidth = 0.5)
    for i in range(len(x_ticks_np)):
        if (i+1)% 2 ==0:
            ax.axvspan(x_ticks_np[i-1], x_ticks_np[i], alpha = 0.2, color = 'grey')

    
    ax.xaxis.set_major_locator(ticker.FixedLocator(x_ticks))
    ax.xaxis.set_minor_locator(ticker.FixedLocator(x_minor_ticks))
    
    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(my_xticks))    
    
    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)
        tick.label1.set_horizontalalignment('center')
        tick.label1.set_fontsize(7) 

    plt.savefig(path+fname+".png", dpi = 300, figsize = (8, 4)) 
    plt.close()


           
#
#
def create_binnded_background_for_bed_files_keeping_reads_ids(bed_f, fl_overlap, bg_f):
    print bed_f
    print bg_f
    #print chrms
    bed = open(bed_f, 'r')
    bg = open(bg_f, 'w+')
    chrs = {}
    chrs_w_ids = {}
    #print chrms_flat
    table = {}
    print 'overlap'
    with open(fl_overlap, 'r') as overlap_table:
        for i, l in enumerate(overlap_table):
            ln = l.rstrip().split("\t")
            if ln[3] not in table : table[ln[3]] = []
            if ln[6] not in table[ln[3]] : table[ln[3]].append(ln[6])
    
    gc.collect()            
    for chr in chrms_flat:
        chrs[chr] = collections.defaultdict(list)
        chrs_w_ids[chr] = collections.defaultdict(list)
    count_not_in_annot_table = 0
    count_total_lines_in_bed = 0
    print 'bed'    
    for line in bed:
        ln = line.rstrip().split("\t")
        chr = ln[0]
        count_total_lines_in_bed+=1
        bin1 = binId(int(ln[1]))
        #bin2 = binId(int(ln[2]))
        
        chrs[chr][bin1] = chrs[chr][bin1]+1 if bin1 in chrs[chr]  else 1
        if bin1 not in chrs_w_ids[chr]: chrs_w_ids[chr][bin1] = [[],[]]
        if ln[3] not in chrs_w_ids[chr][bin1][0]: 
            chrs_w_ids[chr][bin1][0].append(ln[3]) 
            #chrs_w_ids[chr][bin1][1].append([name for name in table[ln[3]] if name not in chrs_w_ids[chr][bin1][1]]) 
            if ln[3] not in table:
                count_not_in_annot_table+=1
                continue
            chrs_w_ids[chr][bin1][1].extend(table[ln[3]]) 
        '''
        if bin1!=bin2:             
            chrs[chr][bin2] = chrs[chr][bin2]+1 if bin2 in chrs[chr]  else 1
            if bin2 not in chrs_w_ids[chr]: chrs_w_ids[chr][bin2] = [[],[]]
            if ln[3] not in chrs_w_ids[chr][bin2]: 
                chrs_w_ids[chr][bin2][0].append(ln[3])
                #chrs_w_ids[chr][bin2][1].extend([name for name in table[ln[3]] if name not in chrs_w_ids[chr][bin2][1]])  
                chrs_w_ids[chr][bin2][1].extend(table[ln[3]])   
        '''        
    print "count_not_in_annot_table ", 
    print count_not_in_annot_table
    print "count_total_lines_in_bed ",
    print count_total_lines_in_bed
    print 'writing'            
    for chr in chrms_flat:
        bins = sorted (chrs[chr].keys())
        print chr +" ",
        print len(bins)
        for bn in bins:
            bg.write(chr + "\t" + str(bn*bin_size)+ "\t" + str(bn*bin_size+bin_size) + "\t" + str(chrs[chr][bn]) + "\t")
            for id in chrs_w_ids[chr][bn][0]:                
                bg.write(id+";")
            bg.write("\t")    
            for id in chrs_w_ids[chr][bn][1]:                
                bg.write(id+";")            
            bg.write("\n")
    bg.close()
    bed.close()


#
def normalize_binnded_background_grid(fl_bg, fl_norm_out, fl_bin_test_grid, fl_number_of_distinct_rnas, fl_chrm_sum_per_bin, fl_mean_per_bin, fl_ratio_value_sum):
    norm = open(fl_norm_out, 'w+')
    
    binned = open(fl_bin_test_grid, 'w+')
    distinct_rnas = open(fl_number_of_distinct_rnas, 'w+')
    chrm_sum_per_bin = open(fl_chrm_sum_per_bin,'w+')
    mean_per_bin =  open(fl_mean_per_bin, 'w+')
    ratio_value_sum =  open(fl_ratio_value_sum, 'w+')
    
    count_to_test = 100
    count_bins_single_bin_on_chr = 0
    count_bins_infrequent_bin_on_chr = 0
    data_bin = {}
    data_name = {}
    print 'agregating'
    with open(fl_bg, 'r') as fl_bg:
        for i, l in enumerate(fl_bg):
            ln = l.rstrip().split("\t")   
            chr = ln[0]
            bin1 = binId(int(ln[1]))
            if chr not in data_bin:
                data_bin[chr] = {}
                data_name[chr] = {}
            
            count = Counter(ln[5].split(";"))
            unique_elem = []
            for name in count.keys():
                if name!="":
                    unique_elem.append(name)
                    if name not in data_name[chr]: data_name[chr][name] = {}
                    data_name[chr][name][bin1] = count[name]
            data_bin[chr][bin1] = [int(ln[3]),unique_elem]
            
    print 'calculating and writing'        
    
    number_of_occupied_bins = []
    prev_chrm = ""
    prev_end = 0    
    
    for chr in sorted(data_bin.keys()):
    #for chr in ['chr1']:
        count_bins_single_bin_on_chr = 0
        count_bins_infrequent_bin_on_chr = 0        
        for bin in sorted(data_bin[chr]):
            m_grid = 0
            bin_value = data_bin[chr][bin][0]
            chrm_sum = 0
            number_of_occupied_bins = []
            value = 0
            for name in data_bin[chr][bin][1]:
                chrm_sum_true = 0
                
                #if len(data_name[chr][name])<0:
                if len(data_name[chr][name])<50:
                    bin_value-=data_name[chr][name][bin]
                    continue
                
                number_of_occupied_bins.append(len(data_name[chr][name]))
                m_grid+=1
                for bn in data_name[chr][name]:
                    chrm_sum_true+=data_name[chr][name][bn]
                    chrm_sum+=data_name[chr][name][bn]
                value+=(data_name[chr][name][bin]*1.0)/(1.0 * chrm_sum_true)    
            
            bin1 = bin        
            if prev_chrm=="" or prev_chrm!=chr:
                prev_chrm = chr
                prev_end = bin1
                for i in range (1, bin1 + 1):
                    norm.write(chr + "\t" + str( bin_size * (i-1)) + "\t" + str(bin_size*(i)) + "\t" + "0\n")
                    binned.write(chr + "\t" + str( bin_size * (i-1)) +"\t" +  str(bin_size*(i)) + "\t" + "0\n")
                    distinct_rnas.write(chr + "\t"+ str( bin_size * (i-1)) + "\t" +  str(bin_size*(i)) + "\t" + "0\n")
                    chrm_sum_per_bin.write(chr + "\t"+ str( bin_size * (i-1)) +"\t" +  str(bin_size*(i)) + "\t" + "0\n")
                    mean_per_bin.write(chr + "\t" + str( bin_size * (i-1)) + "\t" + str(bin_size*(i)) + "\t" + "0\n")
                    ratio_value_sum.write(chr + "\t"+ str( bin_size * (i-1)) + "\t" + str(bin_size*(i)) + "\t" + "0\n")
            for i in range (1, 1 + bin1 - prev_end):
                norm.write(chr + "\t"+ str(prev_end * bin_size + bin_size * (i - 1)) 
                           + "\t" + str(prev_end * bin_size + bin_size*(i)) + "\t" + "0\n")
                binned.write(chr + "\t" + str(prev_end * bin_size + bin_size * (i - 1)) 
                             + "\t" + str(prev_end * bin_size + bin_size*(i)) + "\t" + "0\n")
                distinct_rnas.write(chr + "\t" + str(prev_end * bin_size + bin_size * (i - 1)) 
                                    + "\t" + str(prev_end * bin_size+ bin_size*(i)) + "\t" + "0\n")
                chrm_sum_per_bin.write(chr + "\t" + str(prev_end * bin_size + bin_size * (i - 1)) 
                                       + "\t" + str(prev_end * bin_size + bin_size * (i)) + "\t" + "0\n")
                mean_per_bin.write(chr + "\t" + str(prev_end * bin_size + bin_size * (i - 1)) 
                                   + "\t" + str(prev_end * bin_size + bin_size * (i)) + "\t" + "0\n")
                ratio_value_sum.write(chr + "\t" + str(prev_end * bin_size + bin_size * (i - 1)) 
                                      + "\t" + str(prev_end * bin_size + bin_size * (i)) + "\t" + "0\n")
            #m_grid=len(data_bin[chr][bin][1])  bin_value
            #final_value= (1.0/(1.0 *m_grid)) * ((1.0 * data_bin[chr][bin][0] * chrms_stat[chr][1]) / (1.0 *chrm_sum))
            
            if chrm_sum==0 or m_grid==0:
                final_value = 0
                ratio_value_sum.write (chr + "\t" + str(bin * bin_size) + "\t" 
                                       +  str(bin * bin_size + bin_size) + "\t" + "0" + "\n")         
                
            else:
                #final_value = (1.0/(1.0 *m_grid)) * ((1.0 * bin_value * chrms_stat[chr][1]) / (1.0 * chrm_sum))
                final_value = (1.0/(1.0 * m_grid)) * ((1.0 * value * chrms_stat[chr][1]))
                ratio_value_sum.write (chr + "\t" + str(bin * bin_size) + "\t" 
                                       + str(bin * bin_size + bin_size) + "\t" 
                                       + str((1.0 * data_bin[chr][bin][0])/(1.0 *chrm_sum)) + "\n")         
                
            
            if data_bin[chr][bin][0]==chrm_sum :
                final_value = 0
                count_bins_single_bin_on_chr+=1
            
            #if (1.0*data_bin[chr][bin][0])/ chrm_sum >0.1:
                #final_value = 0
            tmp = 0
            for name in number_of_occupied_bins:
                if name<10:
                    tmp+=1
            
            if tmp==len(number_of_occupied_bins) and tmp!=0:
                count_bins_infrequent_bin_on_chr+=1
                final_value = 0
                
            if final_value==249225.0 and count_to_test>0:
                count_to_test-=100
                print chr + " bin = " + str(bin * bin_size)
                print "m_grid = " + str(m_grid)
                print "data_bin[chr][bin][0] = " + str(data_bin[chr][bin][0])
                print "chrms_stat[chr][1] = " + str(chrms_stat[chr][1])
                print "chrm_sum = " + str(chrm_sum)                
                print "value = " + str(value)
                
            norm.write (chr + "\t" + str(bin * bin_size) 
                        + "\t" +  str(bin*bin_size + bin_size) 
                        + "\t" + str(final_value) + "\n")       
            binned.write(chr + "\t" + str(bin * bin_size) 
                         + "\t" +  str(bin * bin_size + bin_size) 
                         + "\t" + str(data_bin[chr][bin][0]) + "\n") 
            distinct_rnas.write(chr + "\t" + str(bin*bin_size) 
                                 + "\t" +  str(bin * bin_size + bin_size) 
                                 + "\t" + str(len(data_bin[chr][bin][1])) + "\n")                     
            chrm_sum_per_bin.write(chr + "\t" + str(bin * bin_size) 
                                    + "\t" +  str(bin * bin_size + bin_size) 
                                    + "\t" + str(chrm_sum) + "\n")         
            mean_per_bin.write(chr + "\t" + str(bin * bin_size) 
                               + "\t" +  str(bin * bin_size + bin_size) 
                               + "\t" + str((1.0 * chrm_sum)/(1.0 * chrms_stat[chr][1])) + "\n")         
            
            prev_end = bin1 + 1
                        
            if chr=='chr1' and (bin * bin_size)<=1000:
                print chr + " bin = " + str(bin * bin_size)
                print "len(data_bin[chr][bin][1]) = " + str(len(data_bin[chr][bin][1]))
                print "data_bin[chr][bin][0] = " + str(data_bin[chr][bin][0])
                print "chrms_stat[chr][1] = " + str(chrms_stat[chr][1])
                print "chrm_sum = " + str(chrm_sum)                
                print "value = " + str(final_value)
     
        print chr,        
        print " bins number where chr sum resides in single bin ",
        print count_bins_single_bin_on_chr        
        print " bins number where chr sum resides in small number bin ",
        print count_bins_infrequent_bin_on_chr
        print " total number of occupied bins ",
        print chrms_stat[chr][2]        
        print " total number of bins ",
        print chrms_stat[chr][1]
    norm.close()
    binned.close()
    distinct_rnas.close()
    chrm_sum_per_bin.close()
    ratio_value_sum.close()
#
#
def normalize_rna_grid(fl_in, fl_out):
    total_coverage = 0
    total_coverage_by_chr = {}
    rna = {}
    with open(fl_in, 'r') as fl_rna:
        for i, l in enumerate(fl_rna):
            ln = l.rstrip().split("\t")   
                      
            bin1 = binId(int(ln[12]))
            bin2 = binId(int(ln[13]))
            chr = ln[11]
            if bin2 - bin1>1: 
                print "Warning: normalize_rna_grid very long RNA " + l            
            if chr in rna:
                rna[chr][bin1] = rna[chr][bin1] + 1 if bin1 in rna[chr] else 1
                total_coverage+=1  
                total_coverage_by_chr[chr]+=1
            else:
                rna[chr] = {}
                rna[chr][bin1] = 1
                total_coverage_by_chr[chr] = 1
                total_coverage+=1
    print  fl_out ,         
    print ' total RNA coverage ',
    print total_coverage
    
    with open(fl_out, 'w+') as norm:
        for chr in sorted(rna.keys()):
            for bin in sorted(rna[chr].keys()):
                norm.write(chr + "\t" + str(bin * bin_size) 
                           + "\t" + str((bin + 1) * bin_size) 
                           + "\t" + str((1.0 * rna[chr][bin] * chrms_stat[chr][1])/total_coverage_by_chr[chr] ) + "\n")
    
            
    
#
#
#
#
def preview_bg(fl_in, fl_out):
    fl = open(fl_in, 'r')
    fl_out = open(fl_out, 'w+')
    count=200000
    for l in fl:        
        if count>0 :
            if 'chr19' in l:
                ln = l.rstrip().split("\t")
                if float(ln[3])>170:
                    count-=1
                    fl_out.write(l)
        else:
            break
    fl.close()    
    fl_out.close()
    

#
def test_binned_bg(fl_in, fl_out):
    fl = open(fl_in,'r')
    fl_out = open(fl_out,'w+')
    count = 200000
    for l in fl:        
        if count>0 :
            if 'chr19' in l:
                ln=l.rstrip().split("\t")
                if int(ln[1]) >= 27738000 and int(ln[2]) <=27739000:
                    #count-=1
                    fl_out.write(l)
        else:
            break
    fl.close()    
    fl_out.close()

#
def test (view,data):
    num=4    
    ar=[None]*num    
    t = time.time()
    #pr_list = [view.apply_async(calculate_counts_for_TU_plots_wide_intergenic_regions,"random_"+str(i) + "_peaks_genes", rand_peaks[i],genes,wndw) for i in range(num)]
    #view.wait(pr_list)    
    t3 = time.time()-t    
    #print('%f secs (multicore)' % t3)    
    #print pr_list
    t = time.time()
    for i in range (0, num):
        view.targets = [i]
        view.block = False
        ar[i] = view.apply(calculate_counts_for_TU_plots_wide_intergenic_regions,"random_"+str(i) + "_peaks_genes", rand_peaks[i],genes,wndw)	
        #if ar.ready() :
            #ar.get()	
    view.wait(ar)
    for i in range(20):
        raw_counts_g[i]=ar[i].get()[0]      
    t3 = time.time()-t    
    #print('%f secs (multicore)' % t3)  

#
def convert_space_to_tab(fl_in,fl_out):
    out = open(fl_out, "w+")
    with open(fl_in, 'r') as fl_in:
        for i, l in enumerate(fl_in):   
            ln = l.rstrip().split(" ")
            out.write(ln[0])
            for t in ln[1:]:
                if t!="":
                    out.write("\t" + t)
            out.write("\n")
    out.close()
#
def test_annot_table(f_in, f_out):
    with open(f_in, 'r') as fl_in,\
         open (f_out,'w+') as out:
        for i, l in enumerate(fl_in): 
            ln = l.rstrip().split("\t")
            if ln[10]!='RUNX1':
                continue
            if i<5:
                out.write("F\t" + l)
            if len(ln)<14:
                out.write("L\t" +l)
            if int (ln[2]) - int(ln[1]) > 1000:
                out.write(str(int (ln[2]) - int(ln[1])) + "\t" + l)
                
#
#
def split_RNA_parts_annot_table_by_chrom(f_in, fl_0, fl_1, fl_2, fl_3, fl_4):
    out1 = open(fl_1,"w+")
    out2 = open(fl_2,"w+")
    out3 = open(fl_3,"w+")
    out4 = open(fl_4,"w+") 
    out0 = open(fl_0, "w+")
    with open(f_in, 'r') as fl_in:
        for i, l in enumerate(fl_in):
            ln = l.rstrip().split("\t")
            chr = ln[0].split("chr")[1]
            out0.write(ln[0] + "\t" + ln[1] + "\t" + ln[2] + "\n")
            if chr in ["1", "2", "3", "4", "5"]:
                out1.write(ln[0] + "\t" + ln[1] + "\t" + ln[2] + "\n")
            if chr in ["6", "7", "8", "9", "10"]:
                out2.write(ln[0] + "\t" + ln[1] + "\t" + ln[2] + "\n")
            if chr in ["11", "12", "13", "14", "15", "16"]:
                out3.write(ln[0] + "\t" + ln[1] + "\t" + ln[2] + "\n")   
            if chr in ["17", "18", "19", "20", "21", "22", "X", "Y"]:
                out4.write(ln[0] + "\t" + ln[1] + "\t" + ln[2] + "\n")             
    out1.close()
    out2.close()
    out3.close()
    out4.close()    
    out0.close()
#            
#
def split_DNA_parts_annot_table_by_chrom(f_in, fl_0, fl_1, fl_2, fl_3, fl_4):
    out1 = open(fl_1,"w+")
    out2 = open(fl_2,"w+")
    out3 = open(fl_3,"w+")
    out4 = open(fl_4,"w+")     
    out0 = open(fl_0, "w+")
    with open(f_in, 'r') as fl_in:
        for i, l in enumerate(fl_in):
            ln = l.rstrip().split("\t")
            chr = ln[11].split("chr")[1]
            out0.write(ln[11] + "\t" + ln[12] + "\t" + ln[13] + "\n")
            if chr in ["1","2","3", "4"]:
                out1.write(ln[11] + "\t" + ln[12] + "\t" +ln[13] + "\n")
            if chr in ["5","6", "7","8","9" ]:
                out2.write(ln[11]+"\t" + ln[12] + "\t" +ln[13] + "\n")
            if chr in ["10", "11","12", "13","14","15"]:
                out3.write(ln[11]+"\t" + ln[12] + "\t" +ln[13] + "\n")   
            if chr in ["16", "17","18", "19","20","21", "22","X","Y"]:
                out4.write(ln[11]+"\t" + ln[12] + "\t" +ln[13] + "\n")             
    out1.close()
    out2.close()
    out3.close()
    out4.close()                
    out0.close()


#
#
def cut_annotation_table_RNA_parts(fl_in_name, fl_out_name):
    if ".tab" in fl_out_name:
        track_name=fl_out_name.split(".tab")[0]
    elif ".bedGraph" in fl_out_name:
        track_name=fl_out_name.split(".bedGraph")[0]
    else:
        track_name=""
        print "Warning: normalize_single_rna problem with track name"      
    
    with open(fl_in_name, 'r') as fl_in,\
         open(fl_out_name, 'w+') as fl_out:
        fl_out.write("track type=bedGraph name=\" " +track_name + " \"description=\" track " + track_name + " \"\n") 
        for i, l in enumerate(fl_in):
            ln=l.rstrip().split("\t")            
            fl_out.write(ln[0]+"\t" + ln[1] + "\t" +ln[2] + "\n")        
        
#
#for Alexey split RNA parts into 4 groups depending on condition on RNA and DNA part
#
def cut_annotion_table(gene_name, chr, start, end, f_in, fl_1, fl_2, fl_3, fl_4 ):
    out1=open(fl_1,"w+")
    out2=open(fl_2,"w+")
    out3=open(fl_3,"w+")
    out4=open(fl_4,"w+")   
    
    out1.write("track type=bedGraph name=\"rna_"+gene_name+"_1\"description=\"rna track "+gene_name+"_1 not binned\"\n")
    out2.write("track type=bedGraph name=\"rna_"+gene_name+"_2\"description=\"rna track "+gene_name+"_2 not binned\"\n")
    out3.write("track type=bedGraph name=\"rna_"+gene_name+"_3\"description=\"rna track "+gene_name+"_3 not binned\"\n")
    out4.write("track type=bedGraph name=\"rna_"+gene_name+"_4\"description=\"rna track "+gene_name+"_4 not binned\"\n")
    
    count1=0
    count2=0
    count3=0
    count4=0
    
    
    with open(f_in, 'r') as fl_in:
        for i, l in enumerate(fl_in):
            ln=l.rstrip().split("\t")
            if ln[10]!=gene_name:
                continue
            
            out1.write(ln[0]+"\t" + ln[1] + "\t" +ln[2] + "\n")
            count1+=1
            if ln[11]==chr and start<=int(ln[12]) and end>=int(ln[13]):
                out2.write(ln[0]+"\t" + ln[1] + "\t" +ln[2] + "\n")
                count2+=1
            if ln[11]==chr and (start>int(ln[12]) or end<int(ln[13])) :
                out3.write(ln[0]+"\t" + ln[1] + "\t" +ln[2] + "\n")        
                count3+=1
            if ln[11]!=chr:
                out4.write(ln[0]+"\t" + ln[1] + "\t" +ln[2] + "\n")                 
                count4+=1
    
    print count1
    print count2
    print count3
    print count4
    out1.close()
    out2.close()
    out3.close()
    out4.close()

#
#
def create_binnded_background_for_bed_files_only_left_end(bed_f, bg_f):
    print bed_f
    print bg_f
    print chrms
    bed=open(bed_f, 'r')
    bg=open(bg_f, 'w+')
    fname=bed_f.split("k562_")[1]
    print fname
    
    bg.write("track type=bedGraph name=\"rna_" +fname+ "\"description=\"rna track " + fname + " binned\"\n")
     
    chrs={}
    #chrms_flat=[item for sublist in chrms for item in sublist]
    print chrms_flat
    for chr in chrms_flat:
        chrs[chr]=collections.defaultdict(list)
    for line in bed:
        ln=line.rstrip().split("\t")
        if len (ln)<3:
            print "Warning: create_binnded_background_for_bed_files length of the input line is less than 4 elements"
            print ln
            continue
        chr=ln[0]
        bin1=binId(int(ln[1]))      

        if bin1 in chrs[chr]:
            chrs[chr][bin1]+=1
        else:
            chrs[chr][bin1]=1
                
    for chr in chrms_flat:
        bins=sorted (chrs[chr].keys())
        print chr +" ",
        print len(bins)
        for bn in bins:
            bg.write(chr + "\t" + str(bn*bin_size)+ "\t" + str(bn*bin_size+bin_size) + "\t" + str(chrs[chr][bn]) + "\n")
    bg.close()
    bed.close()
#
#
def correct_last_bin_chmrs_end(fl_in_name, fl_out_name):
    with open(fl_in_name, 'r') as fl_in, \
         open(fl_out_name, 'w+') as fl_out:
        for i, l in enumerate(fl_in):
            ln=l.rstrip().split("\t")
            if "track type=" in l:
                fl_out.write(l)
            if len (ln)<4:
                print "Warning: correct_last_bin_chmrs_end short line " + l
                continue  
            
            if int(ln[2]) > chrms_length[ln[0]]:
                fl_out.write(ln[0]+ "\t" + ln[1] + "\t" + str(chrms_length[ln[0]]) + "\t" + ln[3] + "\n")
            else:
                fl_out.write(l)
    
#
#
def fibr_simple():    
    
    create_binnded_background_for_bed_files("fibr.", raw_firb, raw_fibr_binned)
    simple_normalization_for_bg(raw_fibr_binned, fibr_norm_genome_mean, fibr_norm_chrm_mean, fibr_norm_no, fibr_bg_stat)
    
    load_stat(fibr_bg_stat)  
    
    smooth_with_window(fibr_norm_genome_mean, fibr_norm_genome_mean_smoothed)
    
    
    #rnas=["KCNQ1OT1","MALAT1"]
    #rnas=["XIST","GAPDH","AGAP1","RUNX1","FIRRE","KCNQ1OT1","MALAT1"]
    #rnas=["XIST","GAPDH","AGAP1","RUNX1","FIRRE","KCNQ1OT1","MALAT1"]
    #rnas=["AGAP1"]
    #filter_rna(fibr, fibr_FIRRE,"name","FIRRE")
    #filter_rna(fibr, fibr_RUNX1,"name","RUNX1")
    #filter_rna(fibr, fibr_AGAP1,"name","AGAP1")
    
    
    
    for rna in rnas:    
        
        calculate_single_RNA_coverage(paths_rnas["fibr_" +rna], paths_rnas["fibr_"+rna+"_binned"])
        normalize_single_rna(paths_rnas["fibr_"+rna+"_binned"], paths_rnas["fibr_"+rna+"_binned_normalized"])
        calculate_fold_enrichment_for_single_RNA(fibr_norm_genome_mean_smoothed, paths_rnas["fibr_"+rna+"_binned_normalized"], paths_rnas["fibr_"+rna+"_binned_normalized_fold"])
        filtering_enrichment_signal_for_RNA(paths_rnas["fibr_"+rna+"_binned_normalized_fold"],paths_rnas["fibr_"+rna+"_binned_normalized_fold_filtered"])
        smooth_with_window_rna(paths_rnas["fibr" +"_"+rna+"_binned_normalized_fold_filtered"], paths_rnas["fibr" +"_"+rna+"_binned_normalized_fold_filtered_smoothed"])    
        
        plot_RNA(paths_rnas["fibr_"+rna+"_binned"], fibr_norm_genome_mean_smoothed)
        plot_RNA(paths_rnas["fibr_"+rna+"_binned_normalized"], fibr_norm_genome_mean_smoothed)
        plot_RNA(paths_rnas["fibr_"+rna+"_binned_normalized_fold"], fibr_norm_genome_mean_smoothed)  
        plot_RNA(paths_rnas["fibr_"+rna+"_binned_normalized_fold_filtered_smoothed"], fibr_norm_genome_mean_smoothed)  
        
        
    '''
    plot_bg(raw_fibr_binned)
    gc.collect() 
    plot_bg(fibr_norm_genome_mean)
    gc.collect() 
    plot_bg(fibr_norm_genome_mean_smoothed)
    '''
    
    '''
    gc.collect() 
    plot_bg(fibr_norm_chrm_mean)
    gc.collect() 
    plot_bg(fibr_norm_chrm_mean_smoothed)
    gc.collect() 
    plot_bg(fibr_norm_no)
    gc.collect() 
    plot_bg(fibr_norm_no_smoothed)
    '''

#
#
def fibr_grid():
    #filter_mRNA_from_full_annot_tab(fibr,fibr_mRNA)
    #check_overlap_for_reads(raw_firb, fibr_mRNA,fibr_overlap_mRNA_raw)
    
    #rnas=["XIST","GAPDH","AGAP1","RUNX1","FIRRE","KCNQ1OT1","MALAT1"]    
    #rnas=["AGAP1"]
    
    create_binnded_background_for_bed_files_keeping_reads_ids(raw_firb, fibr_overlap_mRNA_raw, raw_fibr_binned_grid)
    normalize_binnded_background_grid(raw_fibr_binned_grid, fibr_norm_grid, fibr_bin_test_grid, fibr_number_of_distinct_rnas, fibr_chrm_sum_per_bin,fibr_mean_per_bin,fibr_ratio_value_sum)
    load_stat(fibr_bg_stat)
    smooth_with_window(fibr_norm_grid, fibr_norm_smoothed_grid)
    
    #plot_bg_by_chromosome(fibr_norm_grid, "chr1")
    #plot_bg(fibr_norm_grid)
    #plot_bg(fibr_ratio_value_sum)
    #plot_bg(fibr_mean_per_bin)    
    #plot_bg(fibr_bin_test_grid)
    #plot_bg(fibr_number_of_distinct_rnas)
    #plot_bg(fibr_chrm_sum_per_bin)    
    #plot_bg(fibr_norm_smoothed_grid)
    
    
    for rna in rnas:  
        normalize_rna_grid(paths_rnas["fibr_" +rna], paths_rnas["fibr_"+rna+"_binned_normalized_grid"])        
        calculate_fold_enrichment_for_single_RNA(fibr_norm_smoothed_grid, paths_rnas["fibr_"+rna+"_binned_normalized_grid"], paths_rnas["fibr_"+rna+"_binned_normalized_fold_grid"])
        filtering_enrichment_signal_for_RNA(paths_rnas["fibr_"+rna+"_binned_normalized_fold_grid"],paths_rnas["fibr_"+rna+"_binned_normalized_fold_filtered_grid"])        
        smooth_with_window_rna(paths_rnas["fibr_"+rna+"_binned_normalized_fold_filtered_grid"], paths_rnas["fibr_"+rna+"_binned_normalized_fold_filtered_smoothed_grid"])
        
        plot_RNA(paths_rnas["fibr_"+rna+"_binned_normalized_grid"], fibr_norm_smoothed_grid)
        plot_RNA(paths_rnas["fibr_"+rna+"_binned_normalized_fold_grid"], fibr_norm_smoothed_grid)
        plot_RNA(paths_rnas["fibr_"+rna+"_binned_normalized_fold_filtered_grid"], fibr_norm_smoothed_grid)    
        plot_RNA(paths_rnas["fibr_"+rna+"_binned_normalized_fold_filtered_smoothed_grid"], fibr_norm_smoothed_grid)    
     




#
#
def k562_simple():
    
    create_binnded_background_for_bed_files("k562.",raw_k562, raw_k562_binned)
    simple_normalization_for_bg(raw_k562_binned, k562_norm_genome_mean,k562_norm_chrm_mean, k562_norm_no, k562_bg_stat)    
    load_stat(k562_bg_stat)    
    smooth_with_window(k562_norm_genome_mean, k562_norm_genome_mean_smoothed)
    
    
    #rnas=["KCNQ1OT1","MALAT1"]
    #rnas=["XIST","GAPDH","AGAP1","RUNX1","FIRRE","KCNQ1OT1","MALAT1"]
    #rnas=["AGAP1"]
    #filter_rna(k562, paths_rnas["K562_FIRRE"],"name","FIRRE")
    #filter_rna(k562, paths_rnas["K562_RUNX1"],"name","RUNX1")
    #filter_rna(k562, paths_rnas["K562_AGAP1"],"name","AGAP1")
    
    for rna in rnas:    
        calculate_single_RNA_coverage(paths_rnas["K562_" +rna], paths_rnas ["K562_"+rna+"_binned"])
        normalize_single_rna(paths_rnas["K562_"+rna+"_binned"], paths_rnas["K562_"+rna+"_binned_normalized"])
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed, paths_rnas["K562_"+rna+"_binned_normalized"], paths_rnas["K562_"+rna+"_binned_normalized_fold"])
        filtering_enrichment_signal_for_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold"],paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered"])
        smooth_with_window_rna(paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered"], paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered_smoothed"])            
        
        plot_RNA(paths_rnas["K562_"+rna+"_binned"], k562_norm_genome_mean_smoothed)
        plot_RNA(paths_rnas["K562_"+rna+"_binned_normalized"], k562_norm_genome_mean_smoothed)
        plot_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold"], k562_norm_genome_mean_smoothed)
        plot_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered"], k562_norm_genome_mean_smoothed)    
    
    plot_bg(raw_k562_binned)
    gc.collect() 
    plot_bg(k562_norm_genome_mean)
    gc.collect() 
    plot_bg(k562_norm_genome_mean_smoothed)
    '''
    gc.collect()     
    plot_bg(k562_norm_chrm_mean)
    gc.collect() 
    plot_bg(k562_norm_chrm_mean_smoothed)
    gc.collect() 
    plot_bg(k562_norm_no)
    gc.collect() 
    plot_bg(k562_norm_no_smoothed)
    '''

#
#
def k562_grid():
    #filter_mRNA_from_full_annot_tab(k562,k562_mRNA)
    #check_overlap_for_reads(raw_k562, k562_mRNA,k562_overlap_mRNA_raw)
    
    rnas=["XIST","GAPDH","AGAP1","RUNX1","FIRRE","KCNQ1OT1","MALAT1"]
    #rnas=["AGAP1"]
    
    create_binnded_background_for_bed_files_keeping_reads_ids(raw_k562, k562_overlap_mRNA_raw, raw_k562_binned_grid)
    normalize_binnded_background_grid(raw_k562_binned_grid, k562_norm_grid, k562_bin_test_grid, k562_number_of_distinct_rnas, k562_chrm_sum_per_bin,k562_mean_per_bin,k562_ratio_value_sum)
    load_stat(k562_bg_stat)    
    smooth_with_window(k562_norm_grid, k562_norm_smoothed_grid)
    
    #plot_bg_by_chromosome(k562_norm_grid, "chr1")
    #plot_bg(k562_norm_grid)
    #plot_bg(k562_ratio_value_sum)
    #plot_bg(k562_mean_per_bin)
    
    #plot_bg(k562_bin_test_grid)
    #plot_bg(k562_number_of_distinct_rnas)
    #plot_bg(k562_chrm_sum_per_bin)
    
    #plot_bg(k562_norm_smoothed_grid)
    
    
    for rna in rnas:  
        normalize_rna_grid(paths_rnas["K562_" +rna], paths_rnas["K562_"+rna+"_binned_normalized_grid"])        
        calculate_fold_enrichment_for_single_RNA(k562_norm_smoothed_grid, paths_rnas["K562_"+rna+"_binned_normalized_grid"], paths_rnas["K562_"+rna+"_binned_normalized_fold_grid"])
        filtering_enrichment_signal_for_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold_grid"],paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered_grid"])
        smooth_with_window_rna(paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered_grid"], paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered_smoothed_grid"])
            
        plot_RNA(paths_rnas["K562_"+rna+"_binned_normalized_grid"], k562_norm_smoothed_grid)
        plot_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold_grid"], k562_norm_smoothed_grid)
        plot_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered_grid"], k562_norm_smoothed_grid)    
        plot_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered_smoothed_grid"], k562_norm_smoothed_grid)    
        
    


#
#
def convert_all_to_bedgraph():
    
    '''
    rnas_init=["MALAT1"]
    for bg in bg_list:
        for rna in rnas_init:
            smooth_with_window_rna(paths_rnas[bg +"_"+rna+"_binned_normalized_fold_filtered"], paths_rnas[bg +"_"+rna+"_binned_normalized_fold_filtered_smoothed"])            
            smooth_with_window_rna(paths_rnas[bg +"_"+rna+"_binned_normalized_fold_filtered_grid"], paths_rnas[bg +"_"+rna+"_binned_normalized_fold_filtered_smoothed_grid"])
    '''
    
    for path in paths_rnas.values():
        if "XIST" not in path and "AC016205" not in path:
            continue
        if '_binned' in path and '_grid' not in path:
            new_file=path.split(".tab")[0]+".bedGraph"
            with open(path, 'r') as fl, \
                 open (new_file, 'w+') as out:  
                count=0
                for i, l in enumerate(fl):   
                    if i==0 and "bedGraph" not in l:
                        out.write("track type=bedGraph name=\"" +path.split(".tab")[0] + "\"description=\"" + path.split(".tab")[0] + "\"\n")                                            
                    if i==0 and "bedGraph" in l:
                        out.write("track type=bedGraph name=\"" +path.split(".tab")[0] + "\"description=\"" + path.split(".tab")[0] + "\"\n")                                            
                        #out.write(l)
                        continue
                    if "_smoothed" in path:
                        ln=l.rstrip().split("\t")
                        #print ln
                        if "*" in ln[1]:
                            pass
                            #print ln[1].split("*")[0]
                        out.write(ln[0] + "\t" + ln[1].split("*")[0]  + "\t" +  ln[2]  + "\t" +  ln[3] + "\n")                        
                    else:
                        out.write(l)
                    count+=1
            if count==1 or count==0:
                print path
                os.system("rm -f " + new_file)
                       
    for path in paths_bg.values():
        if '_binned_norm' in path and '_grid' in path:
            new_file=path.split(".tab")[0]+".bedGraph"
            #print path.split(".tab")[0]
            with open(path, 'r') as fl, \
                 open (new_file, 'w+') as out:                    
                out.write("track type=bedGraph name=\" " + path.split(".tab")[0] + " \"description=\" track " + path.split(".tab")[0] + " \"\n")
                for i, l in enumerate(fl):   
                    out.write(l)
                
    



#
#
def Workflow_September1():
    line="fibr"
    load_stat(fibr_bg_stat)
    
    print "RUNX1"
    name="RUNX1"
    chr="chr21"
    start=36160098
    end=37376965
    
    #the RNA parts!
    cut_annotion_table(name, chr, start, end, fibr,"Alexey/"+line+"_"+name+"_1.tab","Alexey/"+line+"_"+name+"_2.tab","Alexey/"+line+"_"+name+"_3.tab","Alexey/"+line+"_"+name+"_4.tab")
    #"RNA" parts coverage
    create_binnded_background_for_bed_files(line+"_"+name+"_1","Alexey/"+line+"_"+name+"_1.tab", "Alexey/"+line+"_"+name+"_1_bin_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+name+"_2","Alexey/"+line+"_"+name+"_2.tab", "Alexey/"+line+"_"+name+"_2_bin_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+name+"_3","Alexey/"+line+"_"+name+"_3.tab", "Alexey/"+line+"_"+name+"_3_bin_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+name+"_4", "Alexey/"+line+"_"+name+"_4.tab", "Alexey/"+line+"_"+name+"_4_bin_" +str(bin_size)+".tab")
        
    print "MALAT1"
    name="MALAT1"
    chr="chr11"
    start=65265233
    end=65273987
    
    cut_annotion_table(name, chr, start, end, fibr,"Alexey/"+line+"_"+name+"_1.tab","Alexey/"+line+"_"+name+"_2.tab","Alexey/"+line+"_"+name+"_3.tab","Alexey/"+line+"_"+name+"_4.tab")
    
    create_binnded_background_for_bed_files(line+"_"+name+"_1","Alexey/"+line+"_"+name+"_1.tab", "Alexey/"+line+"_"+name+"_1_bin_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+name+"_2","Alexey/"+line+"_"+name+"_2.tab", "Alexey/"+line+"_"+name+"_2_bin_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+name+"_3","Alexey/"+line+"_"+name+"_3.tab", "Alexey/"+line+"_"+name+"_3_bin_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+name+"_4","Alexey/"+line+"_"+name+"_4.tab", "Alexey/"+line+"_"+name+"_4_bin_" +str(bin_size)+".tab")
    
    print "NEAT1"
    name="NEAT1"
    chr="chr11"
    start=65190245
    end=65213011
    
    cut_annotion_table(name, chr, start, end, fibr,"Alexey/"+line+"_"+name+"_1.tab","Alexey/"+line+"_"+name+"_2.tab","Alexey/"+line+"_"+name+"_3.tab","Alexey/"+line+"_"+name+"_4.tab")
    
    create_binnded_background_for_bed_files(line+"_"+name+"_1","Alexey/"+line+"_"+name+"_1.tab", "Alexey/"+line+"_"+name+"_1_bin_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+name+"_2","Alexey/"+line+"_"+name+"_2.tab", "Alexey/"+line+"_"+name+"_2_bin_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+name+"_3","Alexey/"+line+"_"+name+"_3.tab", "Alexey/"+line+"_"+name+"_3_bin_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+name+"_4","Alexey/"+line+"_"+name+"_4.tab", "Alexey/"+line+"_"+name+"_4_bin_" +str(bin_size)+".tab")
        
    #test_annot_table(fibr,"Alexey/"+line+"_test_RUNX1.tab")
    
    split_RNA_parts_annot_table_by_chrom(fibr, "Alexey/"+line+"_RNA_chr_all.tab", "Alexey/"+line+"_RNA_chr_1_5.tab", "Alexey/"+line+"_RNA_chr_6_10.tab", "Alexey/"+line+"_RNA_chr_11_16.tab", "Alexey/"+line+"_RNA_chr_17_Y.tab")
    create_binnded_background_for_bed_files(line+"_RNA_chr_all", "Alexey/"+line+"_RNA_chr_all.tab", "Alexey/"+line+"_RNA_chr_all_bin.tab")
    
    split_DNA_parts_annot_table_by_chrom(fibr, "Alexey/"+line+"_DNA_chr_all.tab", "Alexey/"+line+"_DNA_chr_1_4.tab", "Alexey/"+line+"_DNA_chr_5_9.tab", "Alexey/"+line+"_DNA_chr_10_15.tab", "Alexey/"+line+"_DNA_chr_16_Y.tab")
    create_binnded_background_for_bed_files(line+"_DNA_chr_all_","Alexey/"+line+"_DNA_chr_all.tab", "Alexey/"+line+"_DNA_chr_all_bin.tab")    
    
    plot_bg("Alexey/"+line+"_DNA_chr_all_bin.tab")




#
#
def Workflow_September2():
    line="K562"
    
    rna1="MIR3648"   
    rna1_out="Alexey2/K562_"+ rna1+".full.tab"
    #filter_rna(k562, rna1_out,"name",rna1)
    
    rna2="MIR3687"
    rna2_out="Alexey2/K562_"+ rna2+".full.tab"
    #filter_rna(k562, rna2_out,"name",rna2)    
    
    rna3="Xrna_31259"
    print rna3
    rna3_out="Alexey2/K562_"+ rna3+".full.tab"
    #filter_rna(k562, rna3_out,"name",rna3)
    
    
    rna_region="chr21_9000000_10000000"
    print rna_region
    rna_region_out="Alexey2/K562_region_"+ rna_region+".full.tab"
    #filter_rna(k562, rna_region_out,"region_"+rna_region,"")
    
    '''
    #100Kb
    bin_size=100000    
    create_binnded_background_for_bed_files(line+"_"+rna1, "Alexey2/K562_"+ rna1 +".full.tab", "Alexey2/K562_"+ rna1 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna2, "Alexey2/K562_"+ rna2 +".full.tab", "Alexey2/K562_"+ rna2 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna3, "Alexey2/K562_"+ rna3 +".full.tab", "Alexey2/K562_"+ rna3 +".full_binned_" +str(bin_size)+".tab")
    #1000Kb    
    bin_size=1000000
    create_binnded_background_for_bed_files(line+"_"+rna1, "Alexey2/K562_"+ rna1 +".full.tab", "Alexey2/K562_"+ rna1 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna2, "Alexey2/K562_"+ rna2 +".full.tab", "Alexey2/K562_"+ rna2 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna3, "Alexey2/K562_"+ rna3 +".full.tab", "Alexey2/K562_"+ rna3 +".full_binned_" +str(bin_size)+".tab")
    #10bp
    '''
    
    rnas=[rna1,rna2,rna3]
    
    #100Kb
    print #100Kb
    init(100000)
    #bin_size=100000    
    
    create_binnded_background_for_bed_files(raw_k562_binned.split(".tab")[0], raw_k562, raw_k562_binned)
    simple_normalization_for_bg(raw_k562_binned, k562_norm_genome_mean, k562_norm_chrm_mean, k562_norm_no, k562_bg_stat)    
    load_stat(k562_bg_stat)    
    smooth_with_window(k562_norm_genome_mean, k562_norm_genome_mean_smoothed)
    
    for rna in rnas:    
        calculate_single_RNA_coverage("Alexey2/K562_"+ rna +".full.tab", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph")
        normalize_single_rna("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph" )
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed,  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph",  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph")
        filtering_enrichment_signal_for_RNA("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph")
        smooth_with_window_rna("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph")            
        correct_last_bin_chmrs_end ("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed_last_bin_cor.bedGraph")
        correct_last_bin_chmrs_end ("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_last_bin_cor.bedGraph")

        
    rnas_test= ["XIST","AC016205.1"]
    for rna in rnas_test:    
        calculate_single_RNA_coverage(paths_rnas["K562_" +rna], paths_rnas ["K562_"+rna+"_binned"])
        normalize_single_rna(paths_rnas["K562_"+rna+"_binned"], paths_rnas["K562_"+rna+"_binned_normalized"])
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed, paths_rnas["K562_"+rna+"_binned_normalized"], paths_rnas["K562_"+rna+"_binned_normalized_fold"])
        filtering_enrichment_signal_for_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold"],paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered"])
        smooth_with_window_rna(paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered"], paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered_smoothed"])            
        #correct_last_bin_chmrs_end()
        
    #1000Kb
    print #1000Kb
    #bin_size=1000000     
    init(1000000)
    
    create_binnded_background_for_bed_files(raw_k562_binned.split(".tab")[0], raw_k562, raw_k562_binned)
    simple_normalization_for_bg(raw_k562_binned, k562_norm_genome_mean, k562_norm_chrm_mean, k562_norm_no, k562_bg_stat)    
    load_stat(k562_bg_stat)    
    smooth_with_window(k562_norm_genome_mean, k562_norm_genome_mean_smoothed)
    
    for rna in rnas:    
        calculate_single_RNA_coverage("Alexey2/K562_"+ rna +".full.tab", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph")
        normalize_single_rna("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph" )
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed,  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph",  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph")
        filtering_enrichment_signal_for_RNA("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph")
        smooth_with_window_rna("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph")            
        correct_last_bin_chmrs_end ("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed_last_bin_cor.bedGraph")
        correct_last_bin_chmrs_end ("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_last_bin_cor.bedGraph")
        
    #bin_size=10
    init(10)
    rnas=[rna_region]
    for rna in rnas:    
        cut_annotation_table_RNA_parts("Alexey2/K562_region_"+ rna +".full.tab", "Alexey2/K562_region_"+ rna +".full_rna_parts_not_binned.bedGraph")
        create_binnded_background_for_bed_files(line + "_region_" + rna_region + "binned" + str(bin_size), "Alexey2/K562_region_"+ rna +".full_rna_parts_not_binned.bedGraph", "Alexey2/K562_region_"+ rna +".full_rna_parts_binned" +str(bin_size)+".bedGraph")    
        correct_last_bin_chmrs_end("Alexey2/K562_region_"+ rna +".full_rna_parts_binned" +str(bin_size)+".bedGraph", "Alexey2/K562_region_"+ rna +".full_rna_parts_binned" +str(bin_size)+"_last_bin_cor.bedGraph")
    
#
#
def Workflow_April19():
    global bin_size
    
    bin_size={"GAPDH" : [1000000, 100000], "MALAT1" : [1000000, 1000],  "XIST" : [1000000, 100000]}
    init(-1)
    #cell_lines=["Dawns", "adneur", "neurSC", "Hela_DRB", "Hela_G1", "Hela_M", "K562", "fibr"]    
    
    cell_lines=[ "K562", "fibr"]    
    rnas_init=["GAPDH", "MALAT", "XIST"]
    
    # k562 сырой фона, бинированный 1Mb; 
    # k562 сырой GAPDH, бинированный 1Mb; нормированный GAPDH бина 100 Kb
    # k562 сырой каунт MALAT1 бин 1Mb; грид-нормированный сигнал MALAT для бина 1 Kb
    # k562 cырой каунт XIST бин 1Mb; грид-нормированный сигнал XIST для бина 100 Kb
    
    # fibr cырой каунт XIST бин 1Mb; грид-нормированный сигнал XIST для бина 100 Kb
    
    
    for line in cell_lines:
        print bin_sizes[line]
        
        '''
        #bin_size=bin_sizes[line]
        bin_sizes_range=[100000, 50000, 20000]
        
        if bin_sizes[line] not in bin_sizes_range:
            bin_sizes_range.append(bin_sizes[line])
        '''
        
        for bsize in bin_sizes[line]:
            bin_size=bsize
            print line ,
            print "bin size = " ,
            print bin_size
            
            create_binnded_background_for_bed_files(line, paths_bg_febr["raw_" + line], \
                                                    paths_bg_febr["raw_" + line + "_binned"])
            
            simple_normalization_for_bg(paths_bg_febr["raw_" + line + "_binned"], \
                                        paths_bg_febr[line + "_binned_norm_genome_mean"], \
                                        paths_bg_febr[line + "_binned_norm_chrm_mean"], \
                                        paths_bg_febr[line + "_binned_norm_no"], \
                                        paths_bg_febr[line + "_bg_stat"])
            
            load_stat(paths_bg_febr[ line + "_bg_stat"]) 
            
            smooth_with_window(paths_bg_febr[line + "_binned_norm_genome_mean"], \
                               paths_bg_febr[line + "_binned_norm_genome_mean_smoothed"])
            
            
            #load_stat("stat_K562_bg_100000") 
            
            for rna in rnas_init_febr:    
                print "====================================================" ,
                print line ,
                print rna ,
                print "===================================================="
                
                
                filter_rna(paths_bg_febr[line + "_annot_table"], \
                           paths_rnas_febr[line + "_" + rna], \
                           "name", \
                           rna)
                
                if rna=="Xrna_K562_20585":
                    os.system("cat " + paths_rnas_febr[line + "_" + "Xrna_K562_20585"] + " " + paths_rnas_febr[line + "_" + "MIR3687"] + " > " + paths_rnas_febr[line + "_" + "MIR3687Xrna"])
                    rna="MIR3687Xrna"
                    
                create_bed_file_for_dna_part_for_rna(paths_rnas_febr[line + "_" + rna], \
                                                     paths_rnas_febr[line + "_" + rna + "_raw_DNA_parts"])            
                
                calculate_single_RNA_coverage(paths_rnas_febr[line + "_" + rna], \
                                              paths_rnas_febr[line + "_" + rna + "_binned"])
                
                correct_last_bin_chmrs_end(paths_rnas_febr[line + "_"  + rna + "_binned"], \
                                           paths_rnas_febr[line + "_"  + rna + "_binned_last_bin_cor"])
                
                normalize_single_rna(paths_rnas_febr[line + "_" + rna + "_binned"], \
                                     paths_rnas_febr[line + "_" + rna + "_binned_normalized"])
                
                calculate_fold_enrichment_for_single_RNA(paths_bg_febr[line + "_binned_norm_genome_mean_smoothed"], \
                                                         paths_rnas_febr[line + "_"  + rna + "_binned_normalized"], \
                                                         paths_rnas_febr[line + "_" + rna + "_binned_normalized_fold"])
                
                filtering_enrichment_signal_for_RNA(paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold"], \
                                                    paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered"])
                
                smooth_with_window_rna(paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered"], \
                                       paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered_smoothed"])            
    
                correct_last_bin_chmrs_end(paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered_smoothed"], \
                                           paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered_smoothed_last_bin_cor"])
                                 
        
       
    
    '''            
    binned left part
    create_binnded_background_for_bed_files(line + "_" + rna + "function:create_binned_background", \
                                            paths_rnas_febr[line + "_" + rna], \
                                            paths_rnas_febr[line + "_" + rna + "_binned_to_ckeck"])
           
    rna1="MIR3648"   
    rna1_out="Alexey2/K562_"+ rna1+".full.tab"
    #filter_rna(k562, rna1_out,"name",rna1)
    
    rna2="MIR3687"
    rna2_out="Alexey2/K562_"+ rna2+".full.tab"
    #filter_rna(k562, rna2_out,"name",rna2)    
    
    rna3="Xrna_31259"
    print rna3
    rna3_out="Alexey2/K562_"+ rna3+".full.tab"
    #filter_rna(k562, rna3_out,"name",rna3)
    
    
    rna_region="chr21_9000000_10000000"
    print rna_region
    rna_region_out="Alexey2/K562_region_"+ rna_region+".full.tab"
    #filter_rna(k562, rna_region_out, "region_" + rna_region, "")
    '''
    '''
    #100Kb
    bin_size=100000    
    create_binnded_background_for_bed_files(line+"_"+rna1, "Alexey2/K562_"+ rna1 +".full.tab", "Alexey2/K562_"+ rna1 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna2, "Alexey2/K562_"+ rna2 +".full.tab", "Alexey2/K562_"+ rna2 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna3, "Alexey2/K562_"+ rna3 +".full.tab", "Alexey2/K562_"+ rna3 +".full_binned_" +str(bin_size)+".tab")
    #1000Kb    
    bin_size=1000000
    create_binnded_background_for_bed_files(line+"_"+rna1, "Alexey2/K562_"+ rna1 +".full.tab", "Alexey2/K562_"+ rna1 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna2, "Alexey2/K562_"+ rna2 +".full.tab", "Alexey2/K562_"+ rna2 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna3, "Alexey2/K562_"+ rna3 +".full.tab", "Alexey2/K562_"+ rna3 +".full_binned_" +str(bin_size)+".tab")
    #10bp
    '''
    '''
    rnas=[rna1,rna2,rna3]
    
    #100Kb
    print #100Kb
    init(100000)
    #bin_size=100000    
    
    create_binnded_background_for_bed_files(raw_k562_binned.split(".tab")[0], raw_k562, raw_k562_binned)
    simple_normalization_for_bg(raw_k562_binned, k562_norm_genome_mean, k562_norm_chrm_mean, k562_norm_no, k562_bg_stat)    
    load_stat(k562_bg_stat)    
    smooth_with_window(k562_norm_genome_mean, k562_norm_genome_mean_smoothed)
    
    for rna in rnas:    
        calculate_single_RNA_coverage("Alexey2/K562_" + rna + ".full.tab", "Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + ".bedGraph")
        normalize_single_rna("Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + ".bedGraph", "Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + "_normalized.bedGraph")
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed,  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph",  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph")
        filtering_enrichment_signal_for_RNA("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph")
        smooth_with_window_rna("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph")            
        correct_last_bin_chmrs_end ("Alexey2/K562_" + rna + ".full_binned" + str(bin_size)+ "_normalized_fold_filtered_smoothed.bedGraph", "Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + "_normalized_fold_filtered_smoothed_last_bin_cor.bedGraph")
        correct_last_bin_chmrs_end ("Alexey2/K562_" + rna + ".full_binned" + str(bin_size)+ ".bedGraph", "Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + "_last_bin_cor.bedGraph")

        
    rnas_test= ["XIST", "AC016205.1"]
    for rna in rnas_test:    
        calculate_single_RNA_coverage(paths_rnas["K562_" + rna], paths_rnas["K562_" + rna + "_binned"])
        normalize_single_rna(paths_rnas["K562_" + rna + "_binned"], paths_rnas["K562_" + rna + "_binned_normalized"])
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed, paths_rnas["K562_"+rna+"_binned_normalized"], paths_rnas["K562_"+rna+"_binned_normalized_fold"])
        filtering_enrichment_signal_for_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold"],paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered"])
        smooth_with_window_rna(paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered"], paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered_smoothed"])            
        #correct_last_bin_chmrs_end()
        
    #1000Kb
    print #1000Kb
    #bin_size=1000000     
    init(1000000)
    
    create_binnded_background_for_bed_files(raw_k562_binned.split(".tab")[0], raw_k562, raw_k562_binned)
    simple_normalization_for_bg(raw_k562_binned, k562_norm_genome_mean, k562_norm_chrm_mean, k562_norm_no, k562_bg_stat)    
    load_stat(k562_bg_stat)    
    smooth_with_window(k562_norm_genome_mean, k562_norm_genome_mean_smoothed)
    
    for rna in rnas:    
        calculate_single_RNA_coverage("Alexey2/K562_"+ rna +".full.tab", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph")
        normalize_single_rna("Alexey2/K562_" + rna + ".full_binned" + str(bin_size)+".bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph" )
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed,  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph",  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph")
        filtering_enrichment_signal_for_RNA("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph")
        smooth_with_window_rna("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph")            
        correct_last_bin_chmrs_end ("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed_last_bin_cor.bedGraph")
        correct_last_bin_chmrs_end ("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_last_bin_cor.bedGraph")
        
    #bin_size=10
    init(10)
    rnas=[rna_region]
    for rna in rnas:    
        cut_annotation_table_RNA_parts("Alexey2/K562_region_" + rna +".full.tab", "Alexey2/K562_region_"+ rna +".full_rna_parts_not_binned.bedGraph")
        create_binnded_background_for_bed_files(line + "_region_" + rna_region + "binned" + str(bin_size), "Alexey2/K562_region_"+ rna +".full_rna_parts_not_binned.bedGraph", "Alexey2/K562_region_"+ rna +".full_rna_parts_binned" +str(bin_size)+".bedGraph")    
        correct_last_bin_chmrs_end("Alexey2/K562_region_"+ rna +".full_rna_parts_binned" +str(bin_size)+".bedGraph", "Alexey2/K562_region_"+ rna +".full_rna_parts_binned" +str(bin_size)+"_last_bin_cor.bedGraph")
    '''
#
#
def Workflow_February19():
    global bin_size
    init(-1)
    #cell_lines=["Dawns", "adneur", "neurSC", "Hela_DRB", "Hela_G1", "Hela_M", "K562", "fibr"]    
    cell_lines=["Dawns", "adneur", "neurSC", "Hela_DRB", "Hela_G1", "Hela_M", "K562", "fibr"]    
    #cell_lines=["Dawns"]    
    rnas_init_febr=["MIR3648", "MIR3687", "Xrna_K562_20585", "SSU-rRNA_Hsa_RepM_10963"]
    
    
    for line in cell_lines:
        print bin_sizes[line]
        '''
        #bin_size=bin_sizes[line]
        bin_sizes_range=[100000, 50000, 20000]
        
        if bin_sizes[line] not in bin_sizes_range:
            bin_sizes_range.append(bin_sizes[line])
        '''
        for bsize in bin_sizes[line]:
            bin_size=bsize
            print line ,
            print "bin size = " ,
            print bin_size
            
            create_binnded_background_for_bed_files(line, paths_bg_febr["raw_" + line], \
                                                    paths_bg_febr["raw_" + line + "_binned"])
            
            simple_normalization_for_bg(paths_bg_febr["raw_" + line + "_binned"], \
                                        paths_bg_febr[line + "_binned_norm_genome_mean"], \
                                        paths_bg_febr[line + "_binned_norm_chrm_mean"], \
                                        paths_bg_febr[line + "_binned_norm_no"], \
                                        paths_bg_febr[line + "_bg_stat"])
            
            load_stat(paths_bg_febr[ line + "_bg_stat"]) 
            
            smooth_with_window(paths_bg_febr[line + "_binned_norm_genome_mean"], \
                               paths_bg_febr[line + "_binned_norm_genome_mean_smoothed"])
            
            
            #load_stat("stat_K562_bg_100000") 
            
            for rna in rnas_init_febr:    
                print "====================================================" ,
                print line ,
                print rna ,
                print "===================================================="
                
                
                filter_rna(paths_bg_febr[line + "_annot_table"], \
                           paths_rnas_febr[line + "_" + rna], \
                           "name", \
                           rna)
                
                if rna=="Xrna_K562_20585":
                    os.system("cat " + paths_rnas_febr[line + "_" + "Xrna_K562_20585"] + " " + paths_rnas_febr[line + "_" + "MIR3687"] + " > " + paths_rnas_febr[line + "_" + "MIR3687Xrna"])
                    rna="MIR3687Xrna"
                    
                create_bed_file_for_dna_part_for_rna(paths_rnas_febr[line + "_" + rna], \
                                                     paths_rnas_febr[line + "_" + rna + "_raw_DNA_parts"])            
                
                calculate_single_RNA_coverage(paths_rnas_febr[line + "_" + rna], \
                                              paths_rnas_febr[line + "_" + rna + "_binned"])
                
                correct_last_bin_chmrs_end(paths_rnas_febr[line + "_"  + rna + "_binned"], \
                                           paths_rnas_febr[line + "_"  + rna + "_binned_last_bin_cor"])
                
                normalize_single_rna(paths_rnas_febr[line + "_" + rna + "_binned"], \
                                     paths_rnas_febr[line + "_" + rna + "_binned_normalized"])
                
                calculate_fold_enrichment_for_single_RNA(paths_bg_febr[line + "_binned_norm_genome_mean_smoothed"], \
                                                         paths_rnas_febr[line + "_"  + rna + "_binned_normalized"], \
                                                         paths_rnas_febr[line + "_" + rna + "_binned_normalized_fold"])
                
                filtering_enrichment_signal_for_RNA(paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold"], \
                                                    paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered"])
                
                smooth_with_window_rna(paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered"], \
                                       paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered_smoothed"])            
    
                correct_last_bin_chmrs_end(paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered_smoothed"], \
                                           paths_rnas_febr[line + "_"  + rna + "_binned_normalized_fold_filtered_smoothed_last_bin_cor"])
                                 
        
       
    
    '''            
    binned left part
    create_binnded_background_for_bed_files(line + "_" + rna + "function:create_binned_background", \
                                            paths_rnas_febr[line + "_" + rna], \
                                            paths_rnas_febr[line + "_" + rna + "_binned_to_ckeck"])
           
    rna1="MIR3648"   
    rna1_out="Alexey2/K562_"+ rna1+".full.tab"
    #filter_rna(k562, rna1_out,"name",rna1)
    
    rna2="MIR3687"
    rna2_out="Alexey2/K562_"+ rna2+".full.tab"
    #filter_rna(k562, rna2_out,"name",rna2)    
    
    rna3="Xrna_31259"
    print rna3
    rna3_out="Alexey2/K562_"+ rna3+".full.tab"
    #filter_rna(k562, rna3_out,"name",rna3)
    
    
    rna_region="chr21_9000000_10000000"
    print rna_region
    rna_region_out="Alexey2/K562_region_"+ rna_region+".full.tab"
    #filter_rna(k562, rna_region_out, "region_" + rna_region, "")
    '''
    '''
    #100Kb
    bin_size=100000    
    create_binnded_background_for_bed_files(line+"_"+rna1, "Alexey2/K562_"+ rna1 +".full.tab", "Alexey2/K562_"+ rna1 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna2, "Alexey2/K562_"+ rna2 +".full.tab", "Alexey2/K562_"+ rna2 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna3, "Alexey2/K562_"+ rna3 +".full.tab", "Alexey2/K562_"+ rna3 +".full_binned_" +str(bin_size)+".tab")
    #1000Kb    
    bin_size=1000000
    create_binnded_background_for_bed_files(line+"_"+rna1, "Alexey2/K562_"+ rna1 +".full.tab", "Alexey2/K562_"+ rna1 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna2, "Alexey2/K562_"+ rna2 +".full.tab", "Alexey2/K562_"+ rna2 +".full_binned_" +str(bin_size)+".tab")
    create_binnded_background_for_bed_files(line+"_"+rna3, "Alexey2/K562_"+ rna3 +".full.tab", "Alexey2/K562_"+ rna3 +".full_binned_" +str(bin_size)+".tab")
    #10bp
    '''
    '''
    rnas=[rna1,rna2,rna3]
    
    #100Kb
    print #100Kb
    init(100000)
    #bin_size=100000    
    
    create_binnded_background_for_bed_files(raw_k562_binned.split(".tab")[0], raw_k562, raw_k562_binned)
    simple_normalization_for_bg(raw_k562_binned, k562_norm_genome_mean, k562_norm_chrm_mean, k562_norm_no, k562_bg_stat)    
    load_stat(k562_bg_stat)    
    smooth_with_window(k562_norm_genome_mean, k562_norm_genome_mean_smoothed)
    
    for rna in rnas:    
        calculate_single_RNA_coverage("Alexey2/K562_" + rna + ".full.tab", "Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + ".bedGraph")
        normalize_single_rna("Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + ".bedGraph", "Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + "_normalized.bedGraph")
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed,  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph",  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph")
        filtering_enrichment_signal_for_RNA("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph")
        smooth_with_window_rna("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph")            
        correct_last_bin_chmrs_end ("Alexey2/K562_" + rna + ".full_binned" + str(bin_size)+ "_normalized_fold_filtered_smoothed.bedGraph", "Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + "_normalized_fold_filtered_smoothed_last_bin_cor.bedGraph")
        correct_last_bin_chmrs_end ("Alexey2/K562_" + rna + ".full_binned" + str(bin_size)+ ".bedGraph", "Alexey2/K562_" + rna + ".full_binned" + str(bin_size) + "_last_bin_cor.bedGraph")

        
    rnas_test= ["XIST", "AC016205.1"]
    for rna in rnas_test:    
        calculate_single_RNA_coverage(paths_rnas["K562_" + rna], paths_rnas["K562_" + rna + "_binned"])
        normalize_single_rna(paths_rnas["K562_" + rna + "_binned"], paths_rnas["K562_" + rna + "_binned_normalized"])
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed, paths_rnas["K562_"+rna+"_binned_normalized"], paths_rnas["K562_"+rna+"_binned_normalized_fold"])
        filtering_enrichment_signal_for_RNA(paths_rnas["K562_"+rna+"_binned_normalized_fold"],paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered"])
        smooth_with_window_rna(paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered"], paths_rnas["K562_"+rna+"_binned_normalized_fold_filtered_smoothed"])            
        #correct_last_bin_chmrs_end()
        
    #1000Kb
    print #1000Kb
    #bin_size=1000000     
    init(1000000)
    
    create_binnded_background_for_bed_files(raw_k562_binned.split(".tab")[0], raw_k562, raw_k562_binned)
    simple_normalization_for_bg(raw_k562_binned, k562_norm_genome_mean, k562_norm_chrm_mean, k562_norm_no, k562_bg_stat)    
    load_stat(k562_bg_stat)    
    smooth_with_window(k562_norm_genome_mean, k562_norm_genome_mean_smoothed)
    
    for rna in rnas:    
        calculate_single_RNA_coverage("Alexey2/K562_"+ rna +".full.tab", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph")
        normalize_single_rna("Alexey2/K562_" + rna + ".full_binned" + str(bin_size)+".bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph" )
        calculate_fold_enrichment_for_single_RNA(k562_norm_genome_mean_smoothed,  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized.bedGraph",  "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph")
        filtering_enrichment_signal_for_RNA("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph")
        smooth_with_window_rna("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph")            
        correct_last_bin_chmrs_end ("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed.bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_normalized_fold_filtered_smoothed_last_bin_cor.bedGraph")
        correct_last_bin_chmrs_end ("Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+".bedGraph", "Alexey2/K562_"+ rna +".full_binned"  +str(bin_size)+"_last_bin_cor.bedGraph")
        
    #bin_size=10
    init(10)
    rnas=[rna_region]
    for rna in rnas:    
        cut_annotation_table_RNA_parts("Alexey2/K562_region_" + rna +".full.tab", "Alexey2/K562_region_"+ rna +".full_rna_parts_not_binned.bedGraph")
        create_binnded_background_for_bed_files(line + "_region_" + rna_region + "binned" + str(bin_size), "Alexey2/K562_region_"+ rna +".full_rna_parts_not_binned.bedGraph", "Alexey2/K562_region_"+ rna +".full_rna_parts_binned" +str(bin_size)+".bedGraph")    
        correct_last_bin_chmrs_end("Alexey2/K562_region_"+ rna +".full_rna_parts_binned" +str(bin_size)+".bedGraph", "Alexey2/K562_region_"+ rna +".full_rna_parts_binned" +str(bin_size)+"_last_bin_cor.bedGraph")
    '''
#
#
def Worflow_April19_test():

    bin_size_local = 100000
    
    init(bin_size_local)
    
    '''
    create_binnded_background_from_annotation_table("fibr_background", 
                                                    "fibr.full.tab", 
                                                    "fibr.background.binned" + str(bin_size_local) + ".tab")
    
    
    simple_normalization_for_bg("fibr.background.binned" + str(bin_size_local) + ".tab", 
                                "fibr.background.binned" + str(bin_size_local) + ".norm_genome_mean.tab", 
                                "fibr.background.binned" + str(bin_size_local) + ".norm_chrm_mean.tab", 
                                "fibr.background.binned" + str(bin_size_local) + ".norm_no.tab", 
                                "fibr.background.binned" + str(bin_size_local) + ".stat")
    
     
    
    smooth_with_window("fibr.background.binned" + str(bin_size_local) + ".norm_genome_mean.tab", 
                       "fibr.background.binned" + str(bin_size_local) + ".norm_genome_mean.smoothed.tab")
    
    '''
    load_stat("fibr.background.binned" + str(bin_size_local) + ".stat")
    
    rnas = ["XIST", "MALAT1"]
    line = "fibr"
    for rna in rnas:    
        print "====================================================" ,
        print line ,
        print rna ,
        print "===================================================="
        
        '''
        filter_rna_new_format("fibr.full.tab", 
                   "fibr." + rna + ".tab", 
                   "name", 
                   rna)
          
        
        calculate_single_RNA_coverage_new_format("fibr." + rna + ".tab", 
                                      "fibr." + rna + ".binned" + str(bin_size_local) + ".tab")
        
       
        normalize_single_rna("fibr." + rna + ".binned" + str(bin_size_local) + ".tab", 
                             "fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.tab")
        
        calculate_fold_enrichment_for_single_RNA("fibr.background.binned" + str(bin_size_local) + ".norm_genome_mean.smoothed.tab", 
                                                 "fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.tab", 
                                                 "fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.tab")
        
        filtering_enrichment_signal_for_RNA("fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.tab",
                                            "fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.tab")
        
        smooth_with_window_rna("fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.tab", 
                               "fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.smoothed.tab")            
    
        correct_last_bin_chmrs_end("fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.smoothed.tab", 
                                   "fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.smoothed.corrected.tab")
        
        
        
        '''
        plot_RNA_signal("fibr." + rna + ".binned" + str(bin_size_local) + ".tab", 
                        "fibr.background.binned" + str(bin_size_local) + ".norm_genome_mean.smoothed.tab")
        
        plot_RNA_signal("fibr." + rna + ".binned" + str(bin_size_local) + ".normalized.FE.filtered.smoothed.corrected.tab", 
                        "fibr.background.binned" + str(bin_size_local) + ".norm_genome_mean.smoothed.tab")
    
    
    '''
    temp_path = "/home/ubuntu/RNAdata/Normalization/September/"
    
    
    create_binnded_background_for_bed_files("fibr_background", temp_path + "fibr.dna.back.bed", temp_path + "fibr.background.binned.tab")
    
    simple_normalization_for_bg(temp_path + "fibr.background.binned.tab", 
                                temp_path + "fibr.background.binned.norm_genome_mean.tab", 
                                temp_path + "fibr.background.binned.norm_chrm_mean.tab", 
                                temp_path + "fibr.background.binned.norm_no.tab", 
                                temp_path + "fibr.background.binned.stat")
    
    load_stat(temp_path + "fibr.background.binned.stat") 
    
    smooth_with_window(temp_path + "fibr.background.binned.norm_genome_mean.tab", 
                       temp_path + "fibr.background.binned.norm_genome_mean.smoothed.tab")
    
    '''
    


#
#


#
#
#
#load_stat(k562_bg_stat)

#filter_rna(fibr, paths_rnas["fibr_AC016205.1"],"name","AC016205.1")
#filter_rna(k562, paths_rnas["K562_AC016205.1"],"name","AC016205.1")

#rnas=["XIST","GAPDH","AGAP1","RUNX1","FIRRE","KCNQ1OT1", "AC016205.1"]
#rnas=["XIST", "AC016205.1"]
#fibr_simple()
#fibr_grid()

#filter_rna(fibr, paths_rnas["fibr_AC016205.1"],"name","AC016205.1")
#filter_rna(k562, paths_rnas["K562_AC016205.1"],"name","AC016205.1")

#k562_simple()
#k562_grid()

#preview_bg(raw_fibr_binned_grid, raw_fibr_binned_grid_preview)
#test_grid_bg(fibr)
#test_binned_bg(raw_firb,raw_firb_preview)
#covert_space_to_tab("fibr.dna.back.wo_ids.sorted.tab.bed","fibr.dna.back.wo_ids.sorted.tab.bed")

#convert_all_to_bedgraph()

#Workflow_September1()

#WWorkflow_September2()

#Workflow_February19()

Worflow_April19_test()

