Edit Run_intersection.py file to change the fst file that you are using.

Change the following line:

'''
fst_file_name = f"{cwd}/two_pool_5%_0.05_5kb_merged.fst"
'''

The fst bed file has the following format:

chrom,start,end,mean_fst
NC_062234.1,12001,13000,0.072
NC_062234.1,126001,128000,0.082
NC_062234.1,282001,283000,0.057
NC_062234.1,334001,335000,0.086
NC_062234.1,349001,350000,0.066


The chromosome names are based on GCF_905475465.1_ilPieNapi1.2 annotation. 


