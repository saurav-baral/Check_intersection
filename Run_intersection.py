import os

cwd = os.getcwd()

fst_file_name = f"{cwd}/two_pool_5%_0.05_5kb_merged.fst"

gff_file_name = f"{cwd}/GCF_905475465.1_ilPieNapi1.2_genomic.gff"

test_gene_set = f"{cwd}/extended_positive_rho_genes_networked.txt"



### Setting Othogroups

import pandas as pd
import os

orthogroup_dictionary  = {}

with open(f"{cwd}/OrthologousGroups.tsv", 'r') as orthogroups_file:
	for lines in orthogroups_file.readlines()[1:]:
		lines_split = lines.strip().split("\t")
		orthogroup_dictionary[lines_split[1][:-3]] = lines_split[0]


#### Read GFF file

import pandas as pd

with open(gff_file_name,'r') as gff_file:
	gff_lines = gff_file.readlines()

output = "Gene_id\tChromosome\tGene_Start\tGene_End"
ortho_id_list = []

for row_details in gff_lines:
	if row_details.startswith("#"):
		continue

	row_split = row_details.strip().split("\t")

	if "gene" in row_split[2]:
		gene_id = row_split[8][3:]
		# print(gene_id)
		if gene_id not in orthogroup_dictionary:
			continue
		ortho_id = orthogroup_dictionary[gene_id]
		ortho_id_list.append(ortho_id)
		gene_start = min(int(row_split[3]),int(row_split[4]))
		gene_end = max(int(row_split[3]),int(row_split[4]))
		scaff = row_split[0]
		output += f"\n{ortho_id}\t{scaff}\t{int(max(0,gene_start-10000))}\t{int(gene_end+10000)}"

with open(f"{gff_file_name[:-4]}_gene_padded_10k.tsv", 'w') as out_file:
	out_file.write(output)



with open(f"{cwd}/rerconverge_tested_gene_list.txt",'r') as rer_file:
	busco_gene_list_full_list = rer_file.readlines()

busco_gene_list_full = [file_name.strip() for file_name in busco_gene_list_full_list]

busco_gene_list = []
for gene_id in busco_gene_list_full:
	if gene_id in ortho_id_list:
		busco_gene_list.append(gene_id)


### Read Fst file

fst_dictionary = {}
with open(fst_file_name,'r') as fst_file:
	for fst_lines in fst_file.readlines()[1:]:
		line_split = fst_lines.strip().split(",")
		
		
		fst_dictionary.setdefault(line_split[0],[])
		fst_dictionary[line_split[0]].append([int(line_split[1]), int(line_split[2])])

### Read padded Genes

gene_loc_dic = {}
with open(f"{gff_file_name[:-4]}_gene_padded_10k.tsv", 'r') as gene_loc_file:
	for gene_line in gene_loc_file.readlines()[1:]:
		line_split = gene_line.strip().split("\t")
		gene_loc_dic[line_split[0]] = [line_split[1],line_split[2],line_split[3]]


### Intersecting

intersecting_genes = []
missing_count = 0

with open(test_gene_set, 'r') as rho_set_file:
	for rho_set_count,lines in enumerate(rho_set_file.readlines()):
		rho_gene_name = lines.strip().split("\t")[0]
		
		try:
			chromosome,start,end = gene_loc_dic[rho_gene_name]
		except KeyError:
			# print(rho_gene_name,"Missing")
			missing_count += 1
			continue
		
		if chromosome not in fst_dictionary:
			continue
		for sections in (fst_dictionary[chromosome]):
			# print(sections)
			s1, e1 = int(start), int(end)
			s2, e2 = sections
			if s1 <= e2 and s2 <= e1:
				# print(lines, end = "")
				intersecting_genes.append(rho_gene_name)
			
print("\n",len(intersecting_genes), rho_set_count+1, len(set((intersecting_genes))), missing_count)
final_gene_number = len(set((intersecting_genes)))

#### Run randomization


rho_set_size = rho_set_count+1 -  missing_count
import random
intersecting_genes_length_list = []
for i in range(1000):
	if i % 10 == 0:
		print("run ",i)
	test_gene_list = random.sample(busco_gene_list, rho_set_size)
	
	intersecting_genes = []
	for gene_name in test_gene_list:
		try:
			chromosome,start,end = gene_loc_dic[gene_name]
	
		except KeyError:
			# print(gene_name, "missing")
			pass
		# print(chromosome,start,end)
		# try:
		if chromosome not in fst_dictionary:
			continue
		for sections in (fst_dictionary[chromosome]):
			# print(sections)
			s1, e1 = int(start), int(end)
			s2, e2 = sections
			if s1 <= e2 and s2 <= e1:
				intersecting_genes.append(gene_name)
	
	intersecting_genes_length_list.append(len(set(intersecting_genes)))
print((intersecting_genes_length_list))

import scipy.stats as stats

# Updated list of observations
observations = intersecting_genes_length_list


value_of_interest = final_gene_number
# print(value_of_interest)
count_extreme = sum(x >= value_of_interest for x in observations)
p_value = count_extreme / len(observations)

# p_value


import seaborn as sns
import matplotlib.pyplot as plt
plt.figure(figsize=(6, 12))
sns.kdeplot(y=intersecting_genes_length_list, fill=True, linewidth=1, color="blue")

# Vertical line at defined place (example: x=10)
defined_value = value_of_interest
plt.axhline(defined_value, color='red', linestyle='--', linewidth=5, label=f"x = {defined_value}\nP = {p_value}")
sns.despine(top=True, right=True)
plt.ylabel("Number of Genes")
plt.xlabel("Density")
# plt.title("Smooth Density Plot with Vertical Line")
plt.legend(fontsize=40)
# plt.grid(axis='y', linestyle='--', alpha=0.7)


plt.savefig(f"{cwd}/density_plot.pdf", format="pdf", bbox_inches="tight")
plt.show()