#!/usr/bin/python

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 12:45:18 2021

sarscov2vcfblacklister.py as described in Freeman et al. 2021
This python script takes a VCF input on the command line with -i
Input the path to the Freeman_full_blacklist.xlsx file with -b
Outputs an annotated version of the VCF stating if variants are blacklisted
Output file has altered ending ".bl.vcf" instead of ".vcf" to avoid overwrite

The following python modules are required: vcf, csv, sys, getopt

@author: tim
"""

#Import modules
import vcf
import csv
import sys, getopt

#%%This section indicates which arguments to take for the script
##The section can be omitted and values input manually if not being run from
##the command line
def main(argv):
   vcffile = ""
   blfile = ""
   try:
      opts, args = getopt.getopt(
      argv,"hi:b:",
      ["INPUT_vcf=","BLACKLIST_file="])
   except getopt.GetoptError:
      print(
      "sarscov2vcfblacklister.py -i <vcffile> -b <blfile>")
      sys.exit(2)
   for opt, arg in opts:
       if opt == "-h":
           print(
           "sarscov2vcfblacklister.py -i <vcffile> -b <blfile>")
           sys.exit()
       elif opt in ("-i", "--INPUT_vcf"):
           vcffile = arg
       elif opt in ("-b", "--BLACKLIST_file"):
           blfile = arg
   print("Input VCF path is "+vcffile)
   print("Blacklist file path is "+blfile)
   return([vcffile, blfile])

#Set parameter values for script from input arguments
if __name__ == "__main__":
   [vcffile, blfile] = main(sys.argv[1:])


#%%Define paths
vcf_path = vcffile #E.g "something.vcf"
blacklist_path = blfile # E.g. "path/to/Freeman_full_blacklist.xlsx"
output_vcf_path = vcf_path[:-4]+'.bl.vcf'

#Read in VCF file
vcf_reader = vcf.Reader(open(vcf_path, 'r'))
col1 = []
col2 = []
col3 = []
col4 = []
col5 = []
col6 = []
col7 = []
col8 = []
col9 = []
for record in vcf_reader:
  col1.append(record.CHROM)
  col2.append(record.POS)
  col3.append(record.ID)
  col4.append(record.REF)
  col5.append(record.ALT)
  col6.append(record.QUAL)
  col7.append(record.FILTER)
  col8.append(record.INFO)
  col9.append(record.FORMAT)
        
#Read in blacklist file
bcol1 = []
bcol2 = []
bcol3 = []
bcol4 = []
bcol5 = []
bcol6 = []
bcol7 = []
bcol8 = []
bcol9 = []
bcol10 = []
bcol11 = []
bcol12 = []
bcol13 = []
with open(blacklist_path, "r") as bl:
  reader = csv.reader(bl, delimiter = "\t")
  for row in reader: 
    bcol1.append(row[0])    
    bcol2.append(row[1])
    bcol3.append(row[2])
    bcol4.append(row[3])
    bcol5.append(row[4])
    bcol6.append(row[5])
    bcol7.append(row[6])
    bcol8.append(row[7])
    bcol9.append(row[8])
    bcol10.append(row[9])
    bcol11.append(row[10])
    bcol12.append(row[11])
    bcol13.append(row[12])

#For each variant in the VCF file:
#Annotate if variant is "caution", "mask" or unlisted in blacklist
col10 = []
col11 = []
for i in range(len(col1)):
  pos = col2[i]
  match_in_b = [j for j in range(1,len(bcol1)) if (
    col2[i]==int(bcol1[j]) and
    col4[i]==bcol2[j] and
    str(col5[i])=='['+bcol3[j]+']')]
#If variant is blacklisted, add information on why in another column
  if len(match_in_b)==0 :
    col10.append("not_blacklisted")
    col11.append("NA")
    continue
  if len(match_in_b)==1 :
    match_in_b_int = match_in_b[0]
    col10.append(bcol13[match_in_b_int])
    col11.append("|")
    if bcol4[match_in_b_int]!='NA':
      col11[i]+=(bcol4[0]+':'+bcol4[match_in_b_int]+"|")
    if bcol5[match_in_b_int]!='FALSE':
      col11[i]+=(bcol5[0]+':'+bcol5[match_in_b_int]+"|")
    if bcol6[match_in_b_int]!='TRUE':
      col11[i]+=(bcol6[0]+':'+bcol6[match_in_b_int]+"|")
    if bcol7[match_in_b_int]!='FALSE':
      col11[i]+=(bcol7[0]+':'+bcol7[match_in_b_int]+"|")
    if bcol8[match_in_b_int]=='FALSE':
      col11[i]+=(bcol8[0]+':'+bcol8[match_in_b_int]+"|")
    if bcol9[match_in_b_int]=='FALSE':
      col11[i]+=(bcol9[0]+':'+bcol9[match_in_b_int]+"|")
    if bcol10[match_in_b_int]!='FALSE':
      col11[i]+=(bcol10[0]+':'+bcol10[match_in_b_int]+"|")
    if bcol11[match_in_b_int]!='FALSE':
      col11[i]+=(bcol11[0]+':'+bcol11[match_in_b_int]+"|")
    if bcol12[match_in_b_int]!='FALSE':
      col11[i]+=(bcol12[0]+':'+bcol12[match_in_b_int]+"|")
    continue
  print("Error: Number of matches is not 0 or 1 at i=="+str(i))
  break
      
#Get VCF header to save
with open(vcf_path, 'r') as h:
  header = h.readlines()
  header = [line for line in header if line.startswith('##')]
  header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBLACKLIST?\tREASON\n")

#Save annotated VCF with altered file ending
with open(str(output_vcf_path), 'w') as f:
  for line in header:
    f.write(line)
  for line_no in range(len(col1)):
    f.write(str(col1[line_no])+'\t'+str(col2[line_no])+'\t'
    +str(col3[line_no])+'\t'+str(col4[line_no])+'\t'
    +str(col5[line_no])+'\t'+str(col6[line_no])+'\t'
    +str(col7[line_no])+'\t'+str(col8[line_no])+'\t'
    +str(col9[line_no])+'\t'+str(col10[line_no])+'\t'
    +str(col11[line_no])+'\n')




