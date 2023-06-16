#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
This script summarize the output log files from NGS clean
Authors: Liangjiao Xue liangjiao.xue@gmail.com
         Jake Reeves jake.reeves2013@gmail.com 
Copy Right : Tsai Lab, the University of Georgia, cjtsai@uga.edu
             Liangjiao Xue, liangjiao.xue@gmail.com
"""


from sys import argv
import os
import re

#################
### FUNCTIONS ###
#################


# Get list of log files
def get_log_files(data_dir):
    all_files = os.listdir(data_dir)
    files = []
    for file in all_files:
        if file.endswith("Final.log.txt"):
            files.append(file)
    return files


# Pull proper percent of reads from each file in files
def analyze_file(filenames, data_dir, file_out):

    out_stream = open(file_out, 'w')

    # ADD ANY NEW MARKERS HERE
    headers = ["File",
               "Read length",
               "Read qual foramt",
               "Total Read Number",
               "After trimming",
               "rRNA reads",
               "HPT reads",
               "NPT reads",
               "IRP reads"]
    out_stream.write("\t".join(headers) + "\n")

    for file in filenames:
        with open(data_dir + '/' + file, 'r') as infile:
            headline = []
            dataline = []
            dataline.append(file)
            #headline.append("file")
            total = 0        
                          
            for aline in infile:
                aline = str(aline)
                if aline.startswith('Read length'):
                    m = re.search('\d+', aline).group(0)
                    dataline.append(str(m))
                    #headline.append('Read length')
                elif aline.startswith('Read Qual foramt'):
                    m = re.search('\d+', aline).group(0)
                    dataline.append(str(m))
                    #headline.append('Read Qual foramt')
                elif aline.startswith('Trimming'):
                    m = re.search('\d+', aline).group(0)
                    dataline.append(str(m))
                    #headline.append('Total Read Number')
                elif re.search('Number of input reads', aline):
                    m = re.search('\d+', aline).group(0)
                    dataline.append(str(m))
                    #headline.append('After trimming')
                elif re.search('Uniquely mapped reads number', aline):
                    m = re.search('\d+', aline).group(0)
                    total += int(m)
                elif re.search('Number of reads mapped to multiple loci', aline):
                    m = re.search('\d+', aline).group(0)
                    total += int(m)
                    dataline.append(str(total))
                    #headline.append('rRNA reads')
                elif aline.startswith('HPT'):
                    m = re.search('\d+', aline).group(0)
                    dataline.append(str(m))
                    #headline.append('HPT reads')
                elif aline.startswith('NPT'):
                    m = re.search('\d+', aline).group(0)
                    dataline.append(str(m))
                    #headline.append('NPT reads')
                elif aline.startswith('IRP'):
                    m = re.search('\d+', aline).group(0)
                    dataline.append(str(m))
                    #headline.append('Irp reads')

        #out_stream.write("\t".join(headers) + "\n") #uncomment and delete headers if this blows up
        out_stream.write("\t".join(dataline) + "\n")
    out_stream.close()


############
### MAIN ###
############

if __name__ == '__main__':
    script, data_dir = argv
    out_file = data_dir + '/trim_summary_all.txt'
    files = get_log_files(data_dir)
    analyze_file(files, data_dir,out_file)