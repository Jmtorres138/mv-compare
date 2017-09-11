#!/usr/bin/python -O
# Jason Matthew Torres
'''
Run enrichment jobs for each metabolite:ontology group
Usage: python 02.1_run-jobs
'''
# libraries

import sys,os,gzip

work_dir = "/well/got2d/jason/projects/mv-compare/"
in_dir = work_dir + "input_files/"
out_dir = work_dir + "output_files/"
job_dir = work_dir + "jobs/"
log_dir = work_dir + "logs/"
sig_file = in_dir + "metaxcan_full_results_0.01_only.txt"

def get_metab_list():
    metab_list = []
    fin = open(sig_file,'r')
    head_list = fin.readline().strip().split()
    ind = head_list.index("metabolite")
    for line in fin:
        l = line.strip().split()
        metab = l[ind]
        metab_list.append(metab)
    fin.close()
    metab_list = list(set(metab_list))
    return(metab_list)

def run_job(metab,ont):
    job_file = job_dir+"job_"+metab+"_"+ont+".sh"
    fout=open(job_file,'w')
    command_list  = ["Rscript","--vanilla",work_dir+"02.0_go-enrich-metab.R",metab,ont]
    command = " ".join(command_list)
    script='''
#$ -N %s
#$ -pe shmem 1
#$ -P mccarthy.prjc
#$ -q short.qc
#$ -e %s.error
#$ -o %s.out
echo "start time" `date`

module load R/3.3.1

%s
echo "end time" `date`
    ''' % (metab+"_"+ont, log_dir+metab+"_"+ont,log_dir+metab+"_"+ont, command)
    fout.write(script)
    fout.close()
    call = ["qsub", job_file]
    sp.check_call(call))

def main():
    metab_list = get_metab_list()
    ont_list = ["BP","MF","CC"]
    for metab in metab_list:
        for ont in ont_list:
            run_job(metab,ont)

if (__name__=="__main__"): main()
