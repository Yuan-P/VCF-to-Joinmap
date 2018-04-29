#!/usr/bin/python

import sys
import os

if __name__ == "__main__":
	
	if (len(sys.argv) <= 1) or ('help' in sys.argv):
		print 'Usage: python [vcf]multiVCF_to_joinmap.py <vcf_file> <name_of_joinmap_input> <population_type> <heterozygote_ratio> <threshold_for_missing_data_ratio> <minimum_DP_threshold_for_genotype>'
	else:
		vcf_list	= open(sys.argv[1],'r').readlines()
		out_file	= open(sys.argv[1]+'.locus','w')
		out_h_file	= open(sys.argv[1]+'.head','w')
#		location_list	= []
		name	= sys.argv[2]
		popt	= sys.argv[3]
		hetero_ratio	= sys.argv[4]
		miss_ratio	= sys.argv[5]
		minimum_dp_threshold	= sys.argv[6]
		header	= [x for x in vcf_list if '#CHROM' in x][0].strip().split()
		nind	= str(len(header)-9)
		nloc	= 0
	
		for line in vcf_list:
			string	= line.strip()
			stack	= ''
			gt_list	= []
		
			if '#' in string:
				continue
			else:
				cell	= string.split('\t')
				value_order	= cell[8]
				genotype_list	= range(9,len(cell))
#				print len(genotype_list)
				stack	= cell[0].split('|')[0] +':'+ cell[1]
		
				for position,item in enumerate(value_order.split(':')):
					if item	== 'GT':
						gt_val_position	= position
					elif item	== 'DP':
						dp_val_position	= position
		
				for genotype_num in genotype_list:
					genotype	= cell[genotype_num]
		
					GT	= genotype.split(':')[gt_val_position]

					if GT	== './.':
						joinmapGT	= '\t'+'-'
						gt_list.append('-')

					else:
						DP	= int(genotype.split(':')[dp_val_position])
						if DP	< int(minimum_dp_threshold): #minimum depth threshold for genotype
							joinmapGT	= '\t'+'-'
							gt_list.append('-')
						else:
							if GT	== '0/0':
								joinmapGT	= '\t'+'a'
								gt_list.append('a')
							elif GT	== '0/1':# or '1/0':
								joinmapGT	= '\t'+'h'
								gt_list.append('h')
							elif GT	== '1/1':
								joinmapGT	= '\t'+'b'
								gt_list.append('b')
							else:
								joinmapGT	= '\t'+'-'
								gt_list.append('-')
	
					stack	+= joinmapGT
	
#				if ('h' not in gt_list) and 
				if (gt_list.count('h')	< len(genotype_list)*float(hetero_ratio)) and (len(set(gt_list))	> 1) and (gt_list.count('-')  < len(genotype_list)*float(miss_ratio)):
					out_file.write(stack+'\n')
					nloc	+= 1
		out_h_file.write('name = '+ name +'\n'+ 'popt = '+ popt +'\n'+ 'nloc = '+ str(nloc) +'\n'+ 'nind = '+ nind +'\n')
##		os.system('cat '+ sys.argv[1] +'.head '+ sys.argv[1] +'.locus > '+ sys.argv[1] +'.loc') # ; rm '+ sys.argv[1] +'.head')
	os.systems('cat 'sys.argv[1]+'.head '+sys.argv[1]+'.locus > '+ sys.argv[1]+'.loc' )
