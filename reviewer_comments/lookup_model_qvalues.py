import os

lookup_dir = "/well/got2d/apayne/GTEx_v7/metabolite_paper/lasso_cluster_files/"
work_dir = "/well/mccarthy/users/jason/projects/mv-compare/reviewer_comments/"
tiss_list =  os.listdir(lookup_dir)

fout = open(work_dir + "lasso_models.txt",'w')
metab_file = lookup_dir + "Liver" + "/final_tables/" + "M00054.txt"
fin = open(metab_file,'r')
head_line =  fin.readline().strip().split(",")
head_list = head_line[0:2] + head_line[6:12] + ["Tissue"]
fout.write("\t".join(head_list)+"\n")
fin.close()

print tiss_list
for tiss in tiss_list:
    try:
        print tiss
        metab_file = lookup_dir + tiss + "/final_tables/" + "M00054.txt"
        #print metab_file
        fin = open(metab_file,'r')
        fin.readline() # header
        for line in fin:
            l = line.strip().split(",")
            write_list = l[0:2] + l[6:12] + [tiss]
            fout.write("\t".join(write_list)+"\n")
        fin.close()
    except:
        pass
fout.close()
