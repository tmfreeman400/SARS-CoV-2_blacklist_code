# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 20:33:15 2021

This script contains all of the python scripts in the correct sequence to
generate all of the inputs required for the associated R script afterwards.

All pre-requisites are listed at the top of the script in the modules section,
and all paths for files and directories should be changed to match those of
your own system before running this script. In addition, the VT program should
be installed and the path to the installed VT program should be specified

@author: tim
"""

#Import necessary modules, including PyVCF module and matlotlib modules
import subprocess
import vcf

import numpy as np
import random
import statistics
import os.path
from progress.bar import IncrementalBar

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.patches as mpatches

#Specify location of main project directory
main_project_directory_path = "/media/tim/DATA/MovedDownloads/coronavirus_experiments/"

#Specify folder containing the Incremental database files
#These should be of the form: meanAF_(basecaller)_primertrimmed.txt and
# sdAF_(basecaller)_primertrimmed.txt
IncDB_folder = '/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/'

#Specify location of vt program
vt_path = "/media/tim/DATA/MovedDownloads/coronavirus_experiments/vt/vt"

#Specify locations of medaka VCF files for each basecaller
#Use 4 separate directories containing the VCFs for those basecallers
#with medaka only and no other VCFs
medaka_flipflop_files = "/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/medaka/flipflop/SHEF*.vcf"
medaka_rle_files = "/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/medaka/rle/SHEF*.vcf"
medaka_hac4_files = "/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/medaka/hac/SHEF*.vcf"
medaka_hac3_files = "/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/medaka/hac3/SHEF*.vcf"

#Create list of filepaths for all VCF directories for all combinations of
#variant callers and basecallers
vc_bc_dirs_list = [
"/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/medaka/flipflop/",
"/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/medaka/rle/",
"/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/medaka/hac/",
"/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/medaka/hac3/",
"/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/nanopolish/flipflop/",
"/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/nanopolish/rle/",
"/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/nanopolish/hac/",
"/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/nanopolish/hac3/"
]

#Specify file path for the list of all patient identifiers
patient_identifiers_file = "/media/tim/DATA/MovedDownloads/coronavirus_experiments/reduced_884_patient_list.txt"

#Read in reference FASTA file for coronavirus
fasta_path = '/media/tim/DATA/MovedDownloads/coronavirus_experiments/nCoV-2019.reference.fasta'

#Specify the paths to the directories containing the VCF files and depth files
#respecitvely. Both of these directories should have these files split across
#a nested directory structure indicating the basecaller and variant caller, e.g.
# 'medaka/flipflop/VCFfile' or 'nanopolish/hac/depthfile' within them.
VCFs_directory = '/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/vcfs_884/'
depths_directory = '/media/tim/DATA/MovedDownloads/coronavirus_experiments/IncDB_code-master/nanomed_884/depths_884/'

#List names of all basecallers and variant callers used for the above directory
#structure sub-folders:
basecaller_folders = ["hac","hac3","flipflop","rle"]
variantcaller_folders = ["nanopolish", "medaka"]

#Specify file endings used only for nanopolish and medaka VCF files of interest,
#from which figures are plotted in R later
nanopolish_VCF_file_ending = ".merged.vcf"
medaka_VCF_file_ending = ".merged_test.dc.vcf"


#Note: May need to edit "primertrimmed" ending of meanAF and sdAF IncDBs in
#accordance with output from IncDB generation section

#%%
##No more edits should be carried out below this line

##

#Manually set parameter values used for the MC model in our pipeline

#If not running as script, should manually assign parameters in this code chunk

#Set number of simulated patients
npat = 884 

#Number of Monte Carlo simulations for calculating SD threshold
nsims = 20000

#Assumed error rate for model, this is ~5% in ONT SARS-CoV-2 sequencing data
error_rate = 0.05 

#Confidence interval for standard deviation (must be between 0-1).
#Standard deviation values that fall below the lower SD bound for their
#observed allelic fraction will be annotated as "suspect" at this confidence
#interval
conf_interval = 0.999

#Save plot values in output folder path?
save_plot_values = True

#Representative sample of mean coverage across sequenced loci 
#single numeric value
depth = 30
repr_depths = depth

#Minimum coverage depth filter for calling genomic variants
depth_filt = 20 

#Set ploidy for MC model (2 for diploid genomes, 1 for haploid etc.)
ploidy = 1

#Output file folder path
output_folder_path = main_project_directory_path

#%%
##Set other values that should not be changed

#Aggregate allelic fraction values to simulate
allelic_f_vec = [x/100 for x in range(101)] 

#Number of different allelic fraction values to simulate
niter = len(allelic_f_vec)

#Create list of zeros to store mean, standard deviation values later
mean_vec = np.zeros(niter)
std_vec = np.zeros(niter)
#Create matrix of zeros to store standard deviation values later
std_mat = np.zeros((nsims, niter))

#Get representative coverage depths for sequenced loci that are above the
#minimum coverage depth threshold
if type(repr_depths) == "list":
      repr_depths = [value for value in repr_depths if value > depth_filt]


#%%

#Run MC simulation

#Start progress bar
bar = IncrementalBar('MC simulation progress', max = niter )

#For each aggregate allelic fraction value, calculate standard deviation range
# of simulated values
for k in range(niter):
    
    #Update progress bar. Finishes at k=101.
    bar.next()
    
    #Nested for loop below calculates standard deviation at aggregate allelic
    #fraction across multiple samples, repeated for n simulations
    for i in range(nsims):
        
        #Coverage depth c sampled for simulated genomic locus from vector of
        #mean locus coverage depths representative of the sequencing data for
        #this chromosome, rounded up to the nearest whole number. 
        if type(repr_depths) == "list":
            c = random.sample(repr_depths, 1)[0]
        else:
            c = repr_depths
        
        #Initialise list of all simulated patient allelic fractions
        all_simulated_observed_afs = np.zeros(npat)
        
        #For each simulated patient, simulate observed allelic fraction
        for j in range(npat):
            
            #Simulate genotype for patient by drawing from binomial
            #distribution where p is the simulated population allelic fraction
            single_pat_gt = np.random.binomial(n=ploidy, p=allelic_f_vec[k])
        
            #From the simulated genotype, simulate number of reads observed 
            #providing evidence for the variant given coverage depth c and
            #given error rate.
            single_pat_reads = np.random.binomial(
            n=c, p=single_pat_gt*(1-error_rate)/ploidy +
            error_rate*(1-single_pat_gt)/ploidy)
        
            #Divide simulated number of reads by total depth to get simulated
            #observed allelic fraction for each patient 
            all_simulated_observed_afs[j] = single_pat_reads/c
            
        std_mat[i,k] = statistics.stdev(all_simulated_observed_afs)
       
bar.finish()

#%%    
#Calculate lower and upper bounds for standard deviation at each allelic 
#fraction from MC simulation results
min_sd_bound = np.zeros(len(allelic_f_vec))
for i in range(1, (len(allelic_f_vec)-1)):
    min_sd_bound[i] = np.quantile(std_mat[:,i], (1 - conf_interval))

max_sd_bound = np.zeros(len(allelic_f_vec))
for i in range(1, (len(allelic_f_vec)-1)):
    max_sd_bound[i] = np.quantile(std_mat[:,i], conf_interval)
        
#Round values to 4 decimal places
min_sd_bound = np.round(min_sd_bound, 4)
max_sd_bound = np.round(max_sd_bound, 4)

#Write function to generate non-overwriting filepath names for ".txt" files
#This function does not work for proposed filenames containing "(" or ")".
def nonoverwrite_output_path(proposed_output_filepath, rep_num=0):
    #Check path name doesn't contain "(" or ")"
    if (len(proposed_output_filepath.split(")")) + len(
    proposed_output_filepath.split("(")) > 2) and (rep_num == 0) :
        print("Function not suitable for file names containing '(' or ')'")        
        return None
    #If file already exists, alter stem and check again recursively   
    if os.path.isfile(proposed_output_filepath):
        path_stem = proposed_output_filepath[:-4]
        if path_stem[-1] == ")":
            rep_num = int(path_stem.split("(")[1].split(")")[0])
            path_stem = path_stem.split("(")[0]
        return(nonoverwrite_output_path(
        (path_stem+"("+str(rep_num+1)+").txt"), (rep_num+1)))
    #If file path doesn't already exist, return file path
    else:
        return(proposed_output_filepath)

#Write minimum and maximum thresholds for standard deviation to file
output_filepath = nonoverwrite_output_path(
    output_folder_path+"min_sd_"+str(npat)+"_samples.txt")
with open(output_filepath, "w") as f:
    for line in min_sd_bound.tolist():
        f.write(str(line)+"\n")
    print("Lower standard deviation threshold saved in "+output_filepath)
        
output_filepath2 = nonoverwrite_output_path(
    output_folder_path+"max_sd_"+str(npat)+"_samples.txt")
with open(output_filepath2, "w") as f:
    for line in max_sd_bound.tolist():
        f.write(str(line)+"\n")
    print("Upper standard deviation threshold saved in "+output_filepath2)

#%%
#Save plot values for density plot?
if save_plot_values:
    plotx = allelic_f_vec*nsims
    ploty = np.round(np.reshape(std_mat, np.size(std_mat)), 4)
    output_filepath3 = (
    output_folder_path+"MC_plotvalues_"+str(npat)+"_samples.txt")
    output_filepath3 = nonoverwrite_output_path(output_filepath3)
    with open(output_filepath3, "w") as f:
        for line in range(len(plotx)):
            f.write(str(plotx[line])+","+str(ploty[line])+"\n")
        f.write(str(npat))
    print("Values for plotting MC model saved in "+output_filepath3)


#%%This is the script for generating the plot using the above arguments

#Set parameter values that can be adjusted

#Path to MC_plot_values, lower and upper standard deviation boundaries files
#These values were saved in the previous 2 code sections
input_file = output_folder_path+"MC_plotvalues_"+str(npat)+"_samples.txt"
input_file_min = output_folder_path+"min_sd_"+str(npat)+"_samples.txt"
input_file_max = output_folder_path+"max_sd_"+str(npat)+"_samples.txt"

#Load plot values generated by get_mc_model_for_sd.py
total_lines = len(open(input_file, "r").readlines())
with open(input_file, "r") as f:
    plotx = [0] * (total_lines - 1)
    ploty = [0] * (total_lines - 1)
    for line_no in range((total_lines - 1)):      
        value = f.readline().split(",")            
        plotx[line_no] = float(value[0])
        ploty[line_no] = float(value[1])
    Nsamp = f.readline()

#Create lighter version of binary color scale for clearer visualisation
#Define cmap_map function from matplotlib for this
def cmap_map(function, cmap):
    """
    Applies function to alter vectors of shape 3 [r, g, b] on colormap cmap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ("red", "green", "blue"):
        step_dict[key] = list(map(lambda x: x[0], cdict[key]))
    step_list = sum(step_dict.values(), [])
    step_list = np.array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : np.array(cmap(step)[0:3])
    old_LUT = np.array(list(map(reduced_cmap, step_list)))
    new_LUT = np.array(list(map(function, old_LUT)))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i, key in enumerate(["red","green","blue"]):
        this_cdict = {}
        for j, step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j, i]
            elif new_LUT[j,i] != old_LUT[j, i]:
                this_cdict[step] = new_LUT[j, i]
        colorvector = list(map(lambda x: x + (x[1], ), this_cdict.items()))
        colorvector.sort()
        cdict[key] = colorvector

    return colors.LinearSegmentedColormap("colormap",cdict,1024)
    
light_binary = cmap_map(lambda x: x/2 + 0.5, cm.binary)

#Plot density plot            
n_temp = plt.hist2d(plotx, ploty, bins = 40)
plt.close()   
n = plt.hist2d(
    plotx, ploty, cmap = light_binary, norm=colors.LogNorm(
    vmin=1, vmax=n_temp[0].max()), bins = 40)
    
#Set axes limits and add axes labels and title
plt.xlim((0, 1.0))
plt.ylim((0, 0.6))
plt.xlabel("Aggregate allelic fraction")
plt.ylabel("Standard deviation")
plt.title("HW model of expected aggregate allelic fraction vs."+"\n"+
           "standard deviation for "+str(Nsamp)+" simulated patients")

#%%
#Read in lower and upper boundaries for the confidence interval used to 
#annotate suspect loci
with open(input_file_min, "r") as f:
    min_sd_bound_list = f.readlines()

min_sd_bound = [float(value) for value in min_sd_bound_list]

with open(input_file_max, "r") as f:
    max_sd_bound_list = f.readlines()

max_sd_bound = [float(value) for value in max_sd_bound_list]

#Plot lower and upper bounds for standard deviation confidence interval
line1 = plt.plot(
    [x/100 for x in range(101)], min_sd_bound, linewidth=2.0,
     label="Model lower threshold")
line2 = plt.plot(
    [x/100 for x in range(101)], max_sd_bound, linewidth=2.0,
     label="Model upper threshold")

plt.legend(loc = "upper right")

#%%
#Check if filename already exists, and add (n) to filename if so to avoid
#overwriting previous file. Use recursive function for this.
exists = False
file_num = 0

def safe_savefig(output_folder, Nsamples, already_exists, file_number):    
    if already_exists:
        plot_filepath = (
        output_folder_path+"MC_plot_"+str(Nsamples)
        +"("+str(file_number)+").png")
    elif already_exists == False:
        plot_filepath = output_folder_path+"MC_plot_"+str(Nsamples)+".png"
        
    if os.path.isfile(plot_filepath):
        #If file already exists, we must update the variable plot_filepath
        #and try again
        already_exists = True
        file_number += 1
        safe_savefig(
        output_folder, Nsamples, already_exists, file_number)
    else:
        #If file does not exist, we are safe to save file without overwriting
        plt.savefig(plot_filepath, bbox_inches="tight")
        print("Figure saved as "+plot_filepath)

#Save plot in output folder
safe_savefig(
    output_folder_path, Nsamp, exists, file_num)
plt.close()

#%%Now use the MC model to calculate suspect loci from the IncDBs

#Read in reference FASTA file for coronavirus
ref_sequence = []
with open(fasta_path, "r") as f:
    f.readline()
    for line_no in range(500):
        ref_sequence = ref_sequence+[letter for letter in f.readline().split('\n')[0]]

#Read in IncDB for coronavirus
number_of_patients = npat

input_file_min = output_folder_path+"min_sd_"+str(npat)+"_samples.txt"
input_file_max = output_folder_path+"max_sd_"+str(npat)+"_samples.txt"

IncDB_names = basecaller_folders 

for seq in range(4):
    
    meanDB_path = IncDB_folder+'meanAF_'+IncDB_names[seq]+'_primertrimmed.txt'
    sdDB_path = IncDB_folder+'sdAF_'+IncDB_names[seq]+'_primertrimmed.txt'     
    
    total_lines = len(open(meanDB_path, "r").readlines())
    with open(meanDB_path, "r") as f:
        means = [0] * (total_lines)    
        for line_no in range(total_lines):      
            value = f.readline().split()
            means[line_no] = value

    #Rows that do not sum to >0.95 are not universally covered across
    # all sequences and should not be included. Filter these out by replacing
    # mean and sd values with [0, 0, 0, 0] so not detected as suspect
    rowsums = [
    sum([float(line) for line in means[i]]) for i in range(len(means)) ]
    not_covered_coordinates = np.where(
    [ (value < 0.95) for value in rowsums])[0].tolist()
    for i in not_covered_coordinates:
        means[i] = ['0'] * 6
 

    meanA = [float(line[0]) for line in means]
    meanC = [float(line[1]) for line in means]
    meanG = [float(line[2]) for line in means]
    meanT = [float(line[3]) for line in means]
    meanDEL = [float(line[4]) for line in means]
    meanINS = [float(line[5]) for line in means]
    
    #Filter out mean values where it matches reference if this setting is on
    filter_reference = True
    if filter_reference:
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'A']:
            meanA[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'C']:
            meanC[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'G']:
            meanG[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'T']:
            meanT[index] = 0

    mean_all = meanA+meanC+meanG+meanT+meanINS+meanDEL
    
    total_lines = len(open(sdDB_path, "r").readlines())
    with open(sdDB_path, "r") as f:
        sds = [0] * (total_lines)    
        for line_no in range(total_lines):      
            value = f.readline().split()
            sds[line_no] = value

    for i in not_covered_coordinates:
        sds[i] = ['0'] * 6
    
    sdA = [float(line[0]) for line in sds]
    sdC = [float(line[1]) for line in sds]
    sdG = [float(line[2]) for line in sds]
    sdT = [float(line[3]) for line in sds]
    sdDEL = [float(line[4]) for line in sds]
    sdINS = [float(line[5]) for line in sds]
    
        #Filter out sd values where it matches reference if this setting is on
    filter_reference = True
    if filter_reference:
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'A']:
            sdA[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'C']:
            sdC[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'G']:
            sdG[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'T']:
            sdT[index] = 0
            
    sd_all = sdA+sdC+sdG+sdT+sdINS+sdDEL
    
    n_temp = plt.hist2d(mean_all, sd_all, bins = 40)
    plt.close()   
    n = plt.hist2d(
        mean_all, sd_all, cmap = 'binary', norm=colors.LogNorm(
        vmin=0.1, vmax=n_temp[0].max()), bins = 40)
        
    #Set axes limits and add axes labels and title
    plt.xlim((0, 1.0))
    plt.ylim((0, 0.6))
    plt.xlabel("Aggregate allelic fraction")
    plt.ylabel("Standard deviation")
    plt.title("Observed aggregate allelic fraction vs."+"\n"+
               "standard deviation for "+str(number_of_patients)+" simulated patients")
               
    #%%
    #Read in lower and upper boundaries for the confidence interval used to 
    #annotate suspect loci
    with open(input_file_min, "r") as f:
        min_sd_bound_list = f.readlines()
    
    min_sd_bound = [float(value) for value in min_sd_bound_list]
    
    with open(input_file_max, "r") as f:
        max_sd_bound_list = f.readlines()
    
    max_sd_bound = [float(value) for value in max_sd_bound_list]
    
    #Plot lower and upper bounds for standard deviation confidence interval
    line1 = plt.plot(
        [x/100 for x in range(101)], min_sd_bound, linewidth=2.0,
         label="Model lower threshold")
    line2 = plt.plot(
        [x/100 for x in range(101)], max_sd_bound, linewidth=2.0,
         label="Model upper threshold")
    
    plt.legend(loc = "upper right")
    
    #%%
    #Check if filename already exists, and add (n) to filename if so to avoid
    #overwriting previous file. Use recursive function for this.
    exists = False
    file_num = 0
    
    def safe_savefig(output_folder, sequencer, already_exists, file_number):    
        if already_exists:
            plot_filepath = (
            output_folder+"suspect_plot_"+sequencer
            +"("+str(file_number)+").png")
        elif already_exists == False:
            plot_filepath = output_folder+"suspect_plot_"+sequencer+".png"
            
        if os.path.isfile(plot_filepath):
            #If file already exists, we must update the variable plot_filepath
            #and try again
            already_exists = True
            file_number += 1
            safe_savefig(
            output_folder, sequencer, already_exists, file_number)
        else:
            #If file does not exist, we are safe to save file without overwriting
            plt.savefig(plot_filepath, bbox_inches="tight")
            print("Figure saved as "+plot_filepath)
    
    #Save plot in output folder
    safe_savefig(
        output_folder_path, IncDB_names[seq], exists, file_num)
    plt.close()

    #%%
    #For each sequencing method, write a file containing the coordinates of
    #all suspect loci


    allelic_f_vec = [x/100 for x in range(101)] 

    for allelenum in range(6):
    
        xvec = [float(value[allelenum]) for value in means]
        yvec = [float(value[allelenum]) for value in sds]
        alleles = ['A', 'C', 'G', 'T', 'DEL', 'INS']
        allele = alleles[allelenum]
        #Filter out reference alleles
        finitexindex = np.where([value < 0.7 for value in xvec])
        xindexissuspect = [0] * len(xvec)
        xindexissuspect2 = [0] * len(xvec)
        suspectindex = []
        for i in finitexindex[0].tolist():
            x1 = xvec[i]
            y1 = yvec[i]
            x1 = np.round(x1, 2)
            y2 = min_sd_bound[int(np.where(allelic_f_vec == x1)[0].tolist()[0])]
            y2h = max_sd_bound[int(np.where(allelic_f_vec == x1)[0].tolist()[0])]
            if y1 < y2:
                xindexissuspect[i] = 1
            if y1 > y2h:
                xindexissuspect2[i] = 1    

        suspectindex = np.where(xindexissuspect)[0].tolist()
        suspectindex2 = np.where(xindexissuspect2)[0].tolist()
        
        #For all coordinates in suspectindex, get mean AF, SD, allele and
        #Write file in format: "chr1", "0 coordinate", "actual coordinate",
        #mean, sd, m-sd, -2sd, +sd, +2sd, saving under file name
        # allele/chr1_with_z.bed
        suspect_filename = allele+'/'+allele+'chr1_with_z_'+IncDB_names[seq]+'.bed'
        actual = [(value + 1) for value in suspectindex]
        mean = [line[allelenum] for line in [means[value] for value in suspectindex]]
        sd = [line[allelenum] for line in [sds[value] for value in suspectindex]]
        min1sd = [ (float(mean[i]) - float(sd[i])) for i in range(len(suspectindex)) ]
        min2sd = [ (float(mean[i]) - 2*float(sd[i])) for i in range(len(suspectindex)) ]
        min1sd = [ 0 if (value < 0) else value for value in min1sd ]
        min2sd = [ 0 if (value < 0) else value for value in min2sd ]
        plus1sd = [ (float(mean[i]) + float(sd[i])) for i in range(len(suspectindex)) ]
        plus2sd = [ (float(mean[i]) + 2*float(sd[i])) for i in range(len(suspectindex)) ]
        plus1sd = [ 1 if (value > 1) else value for value in plus1sd ]
        plus2sd = [ 1 if (value > 1) else value for value in plus2sd ]
        
        if os.path.isdir(str(output_folder_path+allele)) == False:
            os.mkdir(str(output_folder_path+allele))
        with open(str(output_folder_path+suspect_filename), 'w') as f:
            for line_no in range(len(suspectindex)):
                f.write('chr1\t'+str(suspectindex[line_no])+'\t'
                +str(actual[line_no])+'\t'+str(mean[line_no])+'\t'
                +str(sd[line_no])+'\t'+str(min1sd[line_no])+'\t'
                +str(min2sd[line_no])+'\t'+str(plus1sd[line_no])
                +'\t'+str(plus2sd[line_no])+'\n')
        
        print(len(suspectindex))


#%% Code from get_sd_ratio_884.py inserted here to plot scatterplots of
#suspect loci

IncDB_names = basecaller_folders 

for seq in range(4):
    
    meanDB_path = IncDB_folder+'meanAF_'+IncDB_names[seq]+'_primertrimmed.txt'
    sdDB_path = IncDB_folder+'sdAF_'+IncDB_names[seq]+'_primertrimmed.txt'
    
    total_lines = len(open(meanDB_path, "r").readlines())
    with open(meanDB_path, "r") as f:
        means = [0] * (total_lines)    
        for line_no in range(total_lines):      
            value = f.readline().split()
            means[line_no] = value

    #Rows that do not sum to >0.95 are not universally covered across
    # all sequences and should not be included. Filter these out by replacing
    # mean and sd values with [0, 0, 0, 0] so not detected as suspect
    rowsums = [
    sum([float(line) for line in means[i]]) for i in range(len(means)) ]
    not_covered_coordinates = np.where(
    [ (value < 0.95) for value in rowsums])[0].tolist()
    for i in not_covered_coordinates:
        means[i] = ['0'] * 6 

    meanA = [float(line[0]) for line in means]
    meanC = [float(line[1]) for line in means]
    meanG = [float(line[2]) for line in means]
    meanT = [float(line[3]) for line in means]
    meanDEL = [float(line[4]) for line in means]
    meanINS = [float(line[5]) for line in means]
    
    #Filter out mean values where it matches reference if this setting is on
    filter_reference = True
    if filter_reference:
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'A']:
            meanA[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'C']:
            meanC[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'G']:
            meanG[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'T']:
            meanT[index] = 0

    mean_all = meanA+meanC+meanG+meanT+meanINS+meanDEL
    
    total_lines = len(open(sdDB_path, "r").readlines())
    with open(sdDB_path, "r") as f:
        sds = [0] * (total_lines)    
        for line_no in range(total_lines):      
            value = f.readline().split()
            sds[line_no] = value

    for i in not_covered_coordinates:
        sds[i] = ['0'] * 6
    
    sdA = [float(line[0]) for line in sds]
    sdC = [float(line[1]) for line in sds]
    sdG = [float(line[2]) for line in sds]
    sdT = [float(line[3]) for line in sds]
    sdDEL = [float(line[4]) for line in sds]
    sdINS = [float(line[5]) for line in sds]
    
        #Filter out sd values where it matches reference if this setting is on
    filter_reference = True
    if filter_reference:
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'A']:
            sdA[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'C']:
            sdC[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'G']:
            sdG[index] = 0
        for index in [i for i in range(len(ref_sequence)) if ref_sequence[i] == 'T']:
            sdT[index] = 0
            
    sd_all = sdA+sdC+sdG+sdT+sdINS+sdDEL
    
    #Define function to filter out values with meanAF <5% for both mean and sd
    def filt5(meanvec, sdvec):
        gt5_index = [i for i in range(len(meanvec)) if meanvec[i] > 0.05]
        mean5vec = [meanvec[i] for i in gt5_index]
        sd5vec = [sdvec[i] for i in gt5_index]
        return(mean5vec,sd5vec)
    
    #Plot scatterplot (# out 2D density plot) of loci with min 5% AF   
    (mean_all_5,sd_all_5)=filt5(mean_all,sd_all)
    (meanA_5,sdA_5)=filt5(meanA,sdA)
    (meanC_5,sdC_5)=filt5(meanC,sdC)
    (meanG_5,sdG_5)=filt5(meanG,sdG)
    (meanT_5,sdT_5)=filt5(meanT,sdT)
    (meanDEL_5,sdDEL_5)=filt5(meanDEL,sdDEL)
    (meanINS_5,sdINS_5)=filt5(meanINS,sdINS)
    lens = [len(meanA_5), len(meanC_5), len(meanG_5), len(meanT_5), len(meanDEL_5), len(meanINS_5)]
    [value*100/sum(lens) for value in lens]

        
    colorlist = [
      [70, 168, 108] for i in range(len(meanA_5))]+[
      [46, 169, 176] for i in range(len(meanC_5))]+[
      [242, 147, 48] for i in range(len(meanG_5))]+[
      [234, 93, 78] for i in range(len(meanT_5))]+[
      [102, 103, 173] for i in range(len(meanINS_5))]

    mean_all_5 = meanA_5+meanC_5+meanG_5+meanT_5+meanINS_5
    sd_all_5 = sdA_5+sdC_5+sdG_5+sdT_5+sdINS_5
    
    #Scramble values into non-sequential plotting order
    #so that no color is favoured
    scrambled_index = []
    for i in range(0,101):
        temp_scrambled = np.arange(start=i, stop=len(mean_all_5), step=101)
        scrambled_index = [
           value for value in scrambled_index]+[
           value for value in temp_scrambled]

    mean_all_5_scrambled = [mean_all_5[i] for i in scrambled_index]
    sd_all_5_scrambled = [sd_all_5[i] for i in scrambled_index]
    colorlist_scrambled = [colorlist[i] for i in scrambled_index]
    
    C = np.array(colorlist_scrambled)
    for i in range(len(mean_all_5_scrambled)):
        plt.scatter(
        mean_all_5_scrambled[(i*1):(i*1 + 1)],
        sd_all_5_scrambled[(i*1):(i*1 + 1)],
        c = C[(i*1):(i*1 + 1)]/255.0, marker='.', s=20)
        
    #Set axes limits and add axes labels and title
    plt.xlim((0, 1.0))
    plt.ylim((0, 0.6))
    plt.xlabel("Aggregate allelic fraction")
    plt.ylabel("Standard deviation")
    plt.title("Observed aggregate allelic fraction vs."+"\n"+
               "standard deviation for "+str(number_of_patients)+
               " patients (-DEL)")
    
    hA = mpatches.Patch(color=('#%02x%02x%02x' % (70, 168, 108)), label='A')
    hC = mpatches.Patch(color=('#%02x%02x%02x' % (46, 169, 176)), label='C')
    hG = mpatches.Patch(color=('#%02x%02x%02x' % (242, 147, 48)), label='G')
    hT = mpatches.Patch(color=('#%02x%02x%02x' % (234, 93, 78)), label='T')
    hDEL = mpatches.Patch(color=('#%02x%02x%02x' % (0, 0, 0)), label='DEL')
    hINS = mpatches.Patch(color=('#%02x%02x%02x' % (102, 103, 173)), label='INS')   
    plt.legend(handles=[hA,hC,hG,hT,hINS])    
    
    #%%
    #Read in lower and upper boundaries for the confidence interval used to 
    #annotate suspect loci
    with open(input_file_min, "r") as f:
        min_sd_bound_list = f.readlines()
    
    min_sd_bound = [float(value) for value in min_sd_bound_list]
    
    with open(input_file_max, "r") as f:
        max_sd_bound_list = f.readlines()
    
    max_sd_bound = [float(value) for value in max_sd_bound_list]
    
    #Plot lower and upper bounds for standard deviation confidence interval
    line1 = plt.plot(
        [x/100 for x in range(101)], min_sd_bound, '#193E72',
         linewidth=2.0, label="Model lower threshold")
    line2 = plt.plot(
        [x/100 for x in range(101)], max_sd_bound, '#FED47A',
         linewidth=2.0, label="Model upper threshold")
    
    #Save plot using names of the format:
    #suspect_plot_scatter_basecaller_primertrimmed_nonreferenceonly_noDEL
    plot_filepath = (output_folder_path+'suspect_plot_scatter_'+
    IncDB_names[seq]+'_nonREF_noDEL')
    plt.savefig(plot_filepath+'.png', bbox_inches="tight")
    plt.savefig(plot_filepath+'.pdf', bbox_inches="tight")
    plt.close()

    #%%Repeat the above plotting code, but with deletions included for contrast
    
    colorlist = [
      [70, 168, 108] for i in range(len(meanA_5))]+[
      [46, 169, 176] for i in range(len(meanC_5))]+[
      [242, 147, 48] for i in range(len(meanG_5))]+[
      [234, 93, 78] for i in range(len(meanT_5))]+[
      [0, 0, 0] for i in range(len(meanDEL_5))]+[
      [102, 103, 173] for i in range(len(meanINS_5))]

    mean_all_5 = meanA_5+meanC_5+meanG_5+meanT_5+meanDEL_5+meanINS_5
    sd_all_5 = sdA_5+sdC_5+sdG_5+sdT_5+sdDEL_5+sdINS_5
    
    #Scramble values into non-sequential plotting order
    #so that no color is favoured
    scrambled_index = []
    for i in range(0,101):
        temp_scrambled = np.arange(start=i, stop=len(mean_all_5), step=101)
        scrambled_index = [
           value for value in scrambled_index]+[
           value for value in temp_scrambled]

    mean_all_5_scrambled = [mean_all_5[i] for i in scrambled_index]
    sd_all_5_scrambled = [sd_all_5[i] for i in scrambled_index]
    colorlist_scrambled = [colorlist[i] for i in scrambled_index]
    
    C = np.array(colorlist_scrambled)
    for i in range(len(mean_all_5_scrambled)):
        plt.scatter(
        mean_all_5_scrambled[(i*1):(i*1 + 1)],
        sd_all_5_scrambled[(i*1):(i*1 + 1)],
        c = C[(i*1):(i*1 + 1)]/255.0, marker='.', s=20)
        
    #Set axes limits and add axes labels and title
    plt.xlim((0, 1.0))
    plt.ylim((0, 0.6))
    plt.xlabel("Aggregate allelic fraction")
    plt.ylabel("Standard deviation")
    plt.title("Observed aggregate allelic fraction vs."+"\n"+
               "standard deviation for "+str(number_of_patients)+
               " patients (+DEL)")
    
    hA = mpatches.Patch(color=('#%02x%02x%02x' % (70, 168, 108)), label='A')
    hC = mpatches.Patch(color=('#%02x%02x%02x' % (46, 169, 176)), label='C')
    hG = mpatches.Patch(color=('#%02x%02x%02x' % (242, 147, 48)), label='G')
    hT = mpatches.Patch(color=('#%02x%02x%02x' % (234, 93, 78)), label='T')
    hDEL = mpatches.Patch(color=('#%02x%02x%02x' % (0, 0, 0)), label='DEL')
    hINS = mpatches.Patch(color=('#%02x%02x%02x' % (102, 103, 173)), label='INS')   
    plt.legend(handles=[hA,hC,hG,hT,hDEL,hINS])    
    
    #%%
    #Read in lower and upper boundaries for the confidence interval used to 
    #annotate suspect loci
    with open(input_file_min, "r") as f:
        min_sd_bound_list = f.readlines()
    
    min_sd_bound = [float(value) for value in min_sd_bound_list]
    
    with open(input_file_max, "r") as f:
        max_sd_bound_list = f.readlines()
    
    max_sd_bound = [float(value) for value in max_sd_bound_list]
    
    #Plot lower and upper bounds for standard deviation confidence interval
    line1 = plt.plot(
        [x/100 for x in range(101)], min_sd_bound, '#193E72',
         linewidth=2.0, label="Model lower threshold")
    line2 = plt.plot(
        [x/100 for x in range(101)], max_sd_bound, '#FED47A',
         linewidth=2.0, label="Model upper threshold")

    #Save plot using names of the format:
    #suspect_plot_scatter_basecaller_primertrimmed_nonreferenceonly_withDEL
    plot_filepath = (output_folder_path+'suspect_plot_scatter_'+
    IncDB_names[seq]+'_nonREF_withDEL')
    plt.savefig(plot_filepath+'.png', bbox_inches="tight")
    plt.savefig(plot_filepath+'.pdf', bbox_inches="tight")
    plt.close()

    #%%
    #For each sequencing method, write a file containing the coordinates of
    #all suspect loci


    allelic_f_vec = [x/100 for x in range(101)] 

    for allelenum in range(6):
    
        xvec = [float(value[allelenum]) for value in means]
        yvec = [float(value[allelenum]) for value in sds]
        alleles = ['A', 'C', 'G', 'T', 'DEL', 'INS']
        allele = alleles[allelenum]
        #Filter out reference alleles
        finitexindex = np.where([value < 1.1 for value in xvec])
        sdratio = [0] * len(xvec)
        sdratio2 = [0] * len(xvec)
        for i in finitexindex[0].tolist():
            x1 = xvec[i]
            y1 = yvec[i]
            x1 = np.round(x1, 2)
            y2 = min_sd_bound[int(np.where(allelic_f_vec == x1)[0].tolist()[0])]
            y2h = max_sd_bound[int(np.where(allelic_f_vec == x1)[0].tolist()[0])]
            if y2 > 0:
                sdratio[i] = (y1/y2)
                sdratio2[i] = (y1/y2h)
            else:
                sdratio[i] = 0
                sdratio2[i] = 0

        #Write SD ratio vectors to file so these ratios can be plotted
        filename = allele+'/'+allele+'_minSDratios_'+IncDB_names[seq]+'.bed'
        actual = [(value + 1) for value in range(len(sdratio))]
        mean = [line[allelenum] for line in [means[value] for value in range(len(actual))]]
        sd = [line[allelenum] for line in [sds[value] for value in range(len(actual))]]
        
        if os.path.isdir(str(output_folder_path+allele)) == False:
            os.mkdir(str(output_folder_path+allele))
        with open(str(output_folder_path+filename), 'w') as f:
            for line_no in range(len(sdratio)):
                f.write(str(actual[line_no])+'\t'+str(mean[line_no])+'\t'
                +str(sd[line_no])+'\t'+str(sdratio[line_no])+'\t'
                +str(sdratio2[line_no])+'\n')
        

##

#%%

#Write BASH code to get a list of all medaka VCF filepaths as the input for
#the next step
cmd1 = '''
#Create list of all medaka input files that do not already contain *.dc.*
ls '''+medaka_flipflop_files+''' | grep -v '.dc.' > '''+main_project_directory_path+'''/medaka_vcf_paths.txt;
ls '''+medaka_rle_files+''' | grep -v '.dc.' >> '''+main_project_directory_path+'''/medaka_vcf_paths.txt;
ls '''+medaka_hac4_files+''' | grep -v '.dc.' >> '''+main_project_directory_path+'''/medaka_vcf_paths.txt;
ls '''+medaka_hac3_files+''' | grep -v '.dc.' >> '''+main_project_directory_path+'''/medaka_vcf_paths.txt;
'''

#Run first BASH command 
subprocess.run(cmd1,
    shell=True, check=True,
    executable='/bin/bash')


#For first to last file on the above list, read in the input filename,
# generate a corresponding output filename, and then run the following code
# to decompose multiallele rows:

cmd2 = '''
last_file_num=`wc '''+main_project_directory_path+'''/medaka_vcf_paths.txt | awk '{print $1}'`
for filenum in $(seq 1 $last_file_num)
do
INPUT=`head -n $filenum '''+main_project_directory_path+'''/medaka_vcf_paths.txt | tail -n 1`
OUTPUT=`echo $INPUT | sed 's/\.vcf/\.dc&/'`
'''+vt_path+''' decompose_blocksub -a $INPUT -o $OUTPUT
done
'''

#Run second BASH command to use VT program to split multiallelic variants in
#all of the above specified medaka VCFs, renaming the output files produced
#with .dc.vcf ending to indicate that this has occurred
subprocess.run(cmd2,
    shell=True, check=True,
    executable='/bin/bash')
    


#Read in the list of all patient identifiers to iterate over
with open(patient_identifiers_file, 'r') as l:
    patient_list = l.readlines()

patient_list = [p.rstrip() for p in patient_list]

ref_sequence = []
with open(fasta_path, "r") as f:
    f.readline()
    for line_no in range(500):
        ref_sequence = ref_sequence+[letter for letter in f.readline().split('\n')[0]]
        
#Read in VCF file and corresponding coverage depth file for that patient

for patient in patient_list:
  for aligner in basecaller_folders: # Or other basecaller
    for vcaller in variantcaller_folders: # "medaka" Or "nanopolish"
      ending = nanopolish_VCF_file_ending
      if vcaller == "medaka":
          ending = medaka_VCF_file_ending
      
      vcf_path = (VCFs_directory+'/'+
      vcaller+'/'+aligner+'/'+patient+ending)
      depths_path = (depths_directory+'/'+
      vcaller+'/'+aligner+'/'+patient+".depth")
      
      total_lines = len(open(depths_path, "r").readlines())
      with open(depths_path, 'r') as f:
          depths = [0] * (total_lines)    
          for line_no in range(total_lines):      
              value = f.readline().split()
              depths[line_no] = value
      
      #Keep same first 7 columns, which already includes QUAL and FILTER
      
      vcf_reader = vcf.Reader(open(vcf_path, 'r'))
      col1 = []
      col2 = []
      col3 = []
      col4 = []
      col5 = []
      col6 = []
      col7 = []
      for record in vcf_reader:
        col1.append(record.CHROM)
        col2.append(record.POS)
        col3.append(record.ID)
        col4.append(record.REF)
        col5.append(record.ALT)
        col6.append(record.QUAL)
        col7.append(record.FILTER)
      
      #If there are multiple rows showing a variant call for the same variant,
      #just retain the higher quality variant
      
      catcol = [""]*len(col1)
      checkCOL = ["KEEP"]*len(col1)
      #For each row, check if it has duplicates in columns 2,4&5
      for rownum in range(len(checkCOL)):
          catcol[rownum] = str(col2[rownum])+str(col4[rownum])+str(col5[rownum])
      for rownum in range(len(checkCOL)):
          if (catcol.count(catcol[rownum]) > 1) & (checkCOL[rownum] == "KEEP"):
              #Get indices of all matching variant rows
              match_index = [i for i, x in enumerate(catcol) if x == catcol[rownum]]
              #Get index of duplicate row with highest quality
              match_quals = [col6[i] for i in match_index]
              bestqual_index = match_quals.index(max(match_quals))
              #Set checkCOL value to "DISCARD" for other indices
              match_index.pop(bestqual_index)
              for i in match_index:
                  checkCOL[i] = "DISCARD"
  
      #Discard all rows where checkCOL == "DISCARD"
      col1 = [col1[i] for i in range(len(col1)) if checkCOL[i]!="DISCARD"]
      col2 = [col2[i] for i in range(len(col2)) if checkCOL[i]!="DISCARD"]
      col3 = [col3[i] for i in range(len(col3)) if checkCOL[i]!="DISCARD"]
      col4 = [col4[i] for i in range(len(col4)) if checkCOL[i]!="DISCARD"]
      col5 = [col5[i] for i in range(len(col5)) if checkCOL[i]!="DISCARD"]
      col6 = [col6[i] for i in range(len(col6)) if checkCOL[i]!="DISCARD"]
      col7 = [col7[i] for i in range(len(col7)) if checkCOL[i]!="DISCARD"]
      
      #Remove [ ] around col5 values to reflect standard appearance in VCF
      col5 = [col5[i][0] for i in range(len(col5))]    
      
      #Fix REF, POS and ALT values at any positions where an INDEL occurs
      # combined with starting alleles that do not match
      fix_index = [index for index in range(len(col5)) if (
      ((len(col4[index]) > 1) or (len(col5[index]) > 1)) and (
      str(col4[index])[0]!=str(col5[index])[0]))]
            
      fixed_pos = [int(col2[i])-1 for i in fix_index]
      fix_prebase = [ref_sequence[num-1] for num in fixed_pos]
      val0 = 0
      for j in fix_index:
          col2[j] = fixed_pos[val0]
          col4[j] = str(fix_prebase[val0])+str(col4[j])
          col5[j] = str(fix_prebase[val0])+str(col5[j])
          val0 += 1
          
      #Swap in DEL or INS for deletions and insertions respectively and note
      #original sequence to swap back in at end
      DELindex = [index for index in range(len(col5)) if [
      len(col5[i]) for i in range(len(col5))][index]<[
      len(col4[i]) for i in range(len(col4))][index]
      ]
      DELlist = [col5[i] for i in DELindex]
      INSindex = [index for index in range(len(col5)) if [
      len(col5[i]) for i in range(len(col5))][index]>[
      len(col4[i]) for i in range(len(col4))][index]
      ]
      INSlist = [col5[i] for i in INSindex]
      
    
      #Fix col3 and col7 values if these are missing 
      #(vcf module removes "." and "PASS")
      for i in range(len(col3)):
        if col3[i]==None:
              col3[i]='.'
      
      for i in range(len(col7)):
          if col7[i]==[]:
              col7[i]='PASS'
          else:
              col7[i]=col7[i][0]
        
      #8th column-A,C,G,T,DEL,INS coverage (need to select the corresponding line)
         
      col8 = [','.join(depths[i-1]) for i in col2]
      #The above needs to be shifted along 1 bp if the reference or variant is multiple nucleotides
      j=0    
      for i in col2:
          if ([len(value) for value in col4][j]>1 or [len(value) for value in col5][j]>1):
              col8[j]=','.join(depths[i])
          j+=1
    
      for i in DELindex:
          col5[i] = "DEL"
      for i in INSindex:
          col5[i] = "INS"
      
      #9th column-total depth (should be sum of 8th column)
      col9=[sum([int(v) for v in col8[i].split(',')]) for i in range(len(col8))]
      
      #10th column-variant allelic fraction
      #col10 = [(int(depths[col2[i]-1][['A','C','G','T','DEL','INS'].index(col5[i])]))/col9[i] for i in range(len(col2))]  
      #The above needs to be shifted along 1 bp if the variant is an inserted
      #base to the right, in order to get the INS fraction, or a DEL at a ref
      #of more than 1 nucleotide
      col10=[None]*len(col9)      
      for i in range(len(col2)):
          temp=(int(depths[col2[i]-1][['A','C','G','T','DEL','INS'].index(col5[i])]))          
          if temp==0:
              col10[i]=0
          else:
              col10[i]=temp/col9[i]
          if ([len(value) for value in col4][i]>1 or [len(value) for value in col5][i]>1):
              col10[i]=(int(depths[col2[i]][['A','C','G','T','DEL','INS'].index(col5[i])]))/col9[i]
    
      #5th column: We can now put the exact indel names back in
      # instead of DEL and INS
      j=0    
      for i in DELindex:
          col5[i] = DELlist[j]
          j+=1
  
      j=0    
      for i in INSindex:
          col5[i] = INSlist[j]
          j+=1
  
      #Get VCF header to re-add into output VCF file by reading in all lines
      #that start with the "#" symbol
      with open(vcf_path, 'r') as h:
          header = h.readlines()
      
      header = [line for line in header if line.startswith('##')]
      header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tACGTDELINSCOUNTS\tTOTALDEPTH\tVARFRAC\n")
      
      #Write VCF file with all of these columns
      
      with open(vcf_path+'.summary', 'w') as f2:
                  for line in header:
                      f2.write(line)
                  for line_no in range(len(col1)):
                      f2.write(col1[line_no]+'\t'+str(col2[line_no])+'\t'
                      +str(col3[line_no])+'\t'+str(col4[line_no])+'\t'
                      +str(col5[line_no])+'\t'+str(col6[line_no])+'\t'
                      +str(col7[line_no])+'\t'+str(col8[line_no])+'\t'
                      +str(col9[line_no])+'\t'+str(col10[line_no])+'\n')
                      
  
      #Also write VCF file in the format required for PyClone
  
      #1st column: mutation_id: A unique identifier for the mutation.
      #Concatenate 4/5/2nd columns together to get unique variant ID
      newcol1 = [str(col4[i])+str(col2[i])+str(col5[i]) for i in range(len(col2))]
      
      #Swap in DEL or INS for deletions and insertions respectively and note
      #original sequence to swap back in at end
      DELindex2 = [index for index in range(len(col4)) if [len(col4[i]) for i in range(len(col4))][index]<1]
      DELlist2 = [col4[i] for i in DELindex2]
      for i in DELindex2:
          col4[i] = "INS"
      INSindex2 = [index for index in range(len(col4)) if [len(col4[i]) for i in range(len(col4))][index]>1]
      INSlist2 = [col4[i] for i in INSindex2]
      for i in INSindex2:
          #Use second allele of multiple allele REF to get REF depth, since it
          #is this position where DEL or INS occurs
          tempword = col4[i]
          col4[i] = str(tempword)[1]
          
      DELindex3 = [index for index in range(len(col5)) if [
      len(col5[i]) for i in range(len(col5))][index]<[
      len(col4[i]) for i in range(len(col4))][index]
      ]
      DELlist3 = [col5[i] for i in DELindex3]
      for i in DELindex3:
          col5[i] = "DEL"
      INSindex3 = [index for index in range(len(col5)) if [
      len(col5[i]) for i in range(len(col5))][index]>[
      len(col4[i]) for i in range(len(col4))][index]
      ]
      INSlist3 = [col5[i] for i in INSindex3]
      for i in INSindex3:
          col5[i] = "INS"
          
          
      #2nd column: ref_counts: The number of reads supporting the reference.
      newcol2=[None]*len(newcol1)      
      for i in range(len(col2)):        
          newcol2[i]=(int(depths[col2[i]-1][['A','C','G','T','DEL','INS'].index(col4[i])]))
          if col4[i]=='INS':
              newcol2[i]=(int(depths[col2[i]][['A','C','G','T','DEL','INS'].index(col4[i])]))
  
      #3rd column: var_counts: The number of reads supporting the variant.
      newcol3=[None]*len(newcol1)      
      for i in range(len(col2)):        
          newcol3[i]=(int(depths[col2[i]-1][['A','C','G','T','DEL','INS'].index(col5[i])]))
          if col5[i]=='INS':
              newcol3[i]=(int(depths[col2[i]][['A','C','G','T','DEL','INS'].index(col5[i])]))
                            
      #4th column: normal_cn: The copy number of the locus in SARS-CoV-2.
      newcol4=[1]*len(newcol1)
                            
      #5th column: minor_cn: This equals 0 in SARS-CoV-2.
      newcol5=[0]*len(newcol1)                      
                            
      #6th column: major_cn: This equals 1 in SARS-CoV-2.
      newcol6=[1]*len(newcol1)
      
      #New header only conserves names of 6 columns
      newheader = "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\n"
      
      #Write PyClone TSV file using these new columns
      with open(vcf_path+'.pyclone.tsv', 'w') as f3:
                  for line in newheader:
                      f3.write(line)
                  for line_no in range(len(newcol1)):
                      f3.write(newcol1[line_no]+'\t'+str(newcol2[line_no])+'\t'
                      +str(newcol3[line_no])+'\t'+str(newcol4[line_no])+'\t'
                      +str(newcol5[line_no])+'\t'+str(newcol6[line_no])+'\n')
                      


#Read in the list of all patient identifiers to iterate over
with open(patient_identifiers_file, 'r') as l:
    patient_list = l.readlines()

patient_list = [p.rstrip() for p in patient_list]

#Read in VCF summary file and corresponding suspect loci files for that patient

for patient in patient_list: #E.g. "SHEF-C09F2" or "SHEF-C0989"
  for aligner in basecaller_folders: # Or other basecaller
    for vcaller in variantcaller_folders: # "medaka" Or "nanopolish"
      ending = nanopolish_VCF_file_ending+".summary" 
      if vcaller == "medaka":
          ending = medaka_VCF_file_ending+".summary" 
      
      vcf_path = (VCFs_directory+'/'+
      vcaller+'/'+aligner+'/'+patient+ending)
      
      #Read in all 10 columns and separate out header
      with open(vcf_path, 'r') as f:
          vcfrows = f.readlines()
      
      header = [v for v in vcfrows if v.startswith('##')]
      corerows = [v for v in vcfrows if not v.startswith('#')]

      #Split corerows into respective columns
      col1 = []
      col2 = []
      col3 = []
      col4 = []
      col5 = []
      col6 = []
      col7 = []
      col8 = []
      col9 = []
      col10= []
      for i in range(len(corerows)):
        splitrows = corerows[i].split('\t')
        col1.append(splitrows[0])
        col2.append(splitrows[1])
        col3.append(splitrows[2])
        col4.append(splitrows[3])
        col5.append(splitrows[4])
        col6.append(splitrows[5])
        col7.append(splitrows[6])
        col8.append(splitrows[7])
        col9.append(splitrows[8])
        col10.append(splitrows[9].rstrip('\n'))
      
      #Swap in DEL or INS for deletions and insertions respectively and note
      #original sequence to swap back in at end
      DELindex = [index for index in range(len(col5)) if [len(col5[i]) for i in range(len(col5))][index]<1]
      DELlist = [col5[i] for i in DELindex]
      for i in DELindex:
          col5[i] = "DEL"
      INSindex = [index for index in range(len(col5)) if [len(col5[i]) for i in range(len(col5))][index]>1]
      INSlist = [col5[i] for i in INSindex]
      for i in INSindex:
          col5[i] = "INS"
      
      #Load in suspect loci for each allele respectively
      suspect_loci_file_paths = [
          output_folder_path+'/'+
          allele+'/'+allele+'chr1_with_z_'+aligner+'.bed' for allele in ['A','C','G','T','DEL','INS']]
      for allelenum in range(6):
          with open(suspect_loci_file_paths[allelenum], 'r') as fALT:
              temprows=fALT.readlines()
          if allelenum == 0:
              Arows=temprows
          if allelenum == 1:
              Crows=temprows
          if allelenum == 2:
              Grows=temprows
          if allelenum == 3:
              Trows=temprows
          if allelenum == 4:
              DELrows=temprows
          if allelenum == 5:
              INSrows=temprows
          
      
      #Do calculations for each row of summary file
      zs = []
      for rownum in range(len(corerows)):
          ALTcol1 = []
          ALTcol2 = []
          ALTcol3 = []
          ALTcol4 = []
          ALTcol5 = []
          ALTcol6 = []
          ALTcol7 = []
          ALTcol8 = []
          ALTcol9 = []
          ALT = col5[rownum]
          if ALT=='A':
              ALTrows=Arows
          if ALT=='C':
              ALTrows=Crows
          if ALT=='G':
              ALTrows=Grows
          if ALT=='T':
              ALTrows=Trows
          if ALT=='DEL':
              ALTrows=DELrows
          if ALT=='INS':
              ALTrows=INSrows
          for i in range(len(ALTrows)):
            splitrows = ALTrows[i].split('\t')
            ALTcol1.append(splitrows[0])
            ALTcol2.append(splitrows[1])
            ALTcol3.append(splitrows[2])
            ALTcol4.append(splitrows[3])
            ALTcol5.append(splitrows[4])
            ALTcol6.append(splitrows[5])
            ALTcol7.append(splitrows[6])
            ALTcol8.append(splitrows[7])
            ALTcol9.append(splitrows[8].rstrip('\n'))
          #If variant is in suspect loci file, do calculation for z values
          try:
              ALTindex = ALTcol3.index(col2[rownum])
              varAF=float(col10[rownum])
              sysAF=float(ALTcol4[ALTindex])
              sysSD=float(ALTcol5[ALTindex])
              zMETRIC=(varAF-sysAF)/sysSD
          #Otherwise, mark corresponding position as "NS" (not suspect)
          except ValueError:
              zMETRIC="NS"
          zs.append(zMETRIC)
      
      #5th column: We can now put the exact indel names back in
      # instead of DEL and INS, prior to writing the output file
      j=0    
      for i in DELindex:
          col5[i] = DELlist[j]
          j+=1
  
      j=0    
      for i in INSindex:
          col5[i] = INSlist[j]
          j+=1
  
      #Append additional line to header:
      header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tACGTDELINSCOUNTS\tTOTALDEPTH\tVARFRAC\tZMETRIC\n")
      
      #Write VCF file with all of these columns
      
      with open(vcf_path+'2', 'w') as f2:
                  for line in header:
                      f2.write(line)
                  for line_no in range(len(col1)):
                      f2.write(col1[line_no]+'\t'+str(col2[line_no])+'\t'
                      +str(col3[line_no])+'\t'+str(col4[line_no])+'\t'
                      +str(col5[line_no])+'\t'+str(col6[line_no])+'\t'
                      +str(col7[line_no])+'\t'+str(col8[line_no])+'\t'
                      +str(col9[line_no])+'\t'+str(col10[line_no])+'\t'
                      +str(zs[line_no])+'\t'+patient+'\n')
                      


#Run for loop for each of the directories in the VCF directory list
for i in range(len(vc_bc_dirs_list)):
    
    #BASH code starts here:
    cmd_i = '''
    #Concatenate all results for each pipeline together, preserving patient IDs
    cat '''+vc_bc_dirs_list[i]+'''/*2 | egrep -v "^#" | awk '($9 >= 20)' > '''+vc_bc_dirs_list[i]+'''/ALL_summary3
    
    #Use awk to filter for variants with z values below e.g. 1,
    #to consider as suspect variants (use awk to remove lines with NS in 
    #column 11, then awk on column 11)
    awk '($11 != "NS")' '''+vc_bc_dirs_list[i]+'''/ALL_summary3  | awk '($11 < 1)' > '''+vc_bc_dirs_list[i]+'''/ALL_suspect_variants.vcf
    
    #Generate table of suspect variants, listing how many patients they 
    #occurred in at the z cutoff (Use awk to only take columns 1-5, then sort 
    #and uniq to get number of copies of each unique row, then sort again in 
    #descending order of most to fewest patients per suspect variant)
    awk -vOFS="\t" '{print $1,$2,$3,$4,$5,".","."}' '''+vc_bc_dirs_list[i]+'''/ALL_suspect_variants.vcf | sort | uniq -c | sort -r -n | awk -vOFS="\t" '{print $2,$3,$4,$5,$6,$7,$8,$1}' > '''+vc_bc_dirs_list[i]+'''/ALL_suspect_variants_rankedtable.vcf
    '''
    
    #Run BASH command to do all of the above
    subprocess.run(cmd_i,
        shell=True, check=True,
        executable='/bin/bash')


