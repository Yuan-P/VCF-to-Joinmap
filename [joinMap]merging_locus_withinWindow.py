#!/usr/bin/python

import sys
args	= sys.argv

if len(args)	== 1:
	print 'This script is used at merging variation genotype within size window from *.loc file'
	print 'Usage: python '+ args[0] +' [loc_file]'

else:
	loc_list	= open(args[1],'r').xreadlines()
	out_file	= open(args[1]+'.merged.loc','w')
	tempGT_list	= []
	info_list	= []
	stack_list	= []
	loci_string	= ''
	temp_loc_string	= ''
	stackGT_string	= ''

	for loc in loc_list:
		loc_string	= loc.strip()

		if '=' not in loc_string:
			loc_cell	= loc_string.split('\t')
			loc_loci	= loc_cell[0]
			loc_linkage	= loc_loci.split(':')[0]
			loc_pos	= int(loc_loci.split(':')[1])
			loc_gt	= loc_cell[1:]

			if temp_loc_string	== '':
				temp_loc_string	= loc_string
				temp_loc_cell	= temp_loc_string.split('\t')
				temp_loc_loci   = temp_loc_cell[0]
				temp_loc_linkage        = temp_loc_loci.split(':')[0]
				temp_loc_pos    = int(temp_loc_loci.split(':')[1])
				temp_loc_gt     = temp_loc_cell[1:]
				loci_string	= temp_loc_loci
				for gt in temp_loc_gt:
					tempGT_list.append(list(gt))
				continue
			
			temp_loc_cell	= temp_loc_string.split('\t')
			temp_loc_loci	= temp_loc_cell[0]
			temp_loc_linkage	= temp_loc_loci.split(':')[0]
			temp_loc_pos	= int(temp_loc_loci.split(':')[1])
			temp_loc_gt	= temp_loc_cell[1:]
#			elif (temp_loc_string != '') and (loc_linkage == temp_loc_linkage) and (loc_pos <= temp_loc_pos+10000):

			if temp_loc_string	!= '':
				if (loc_linkage == temp_loc_linkage) and (loc_pos <= temp_loc_pos+10000):
					for pos,gt in enumerate(loc_gt):
						tempGT_list[pos].append(gt)
					continue
				else:
					stackGT_string += loci_string
					for gt_list in tempGT_list:
						gt_num	= len(gt_list)
						count_a	= gt_list.count('a')
						count_b	= gt_list.count('b')
						count_h	= gt_list.count('h')
						if (float(count_a)/gt_num) > 0.9:
							stackGT_string	+= '\t'+'a'
						if (float(count_b)/gt_num) > 0.9:
							stackGT_string	+= '\t'+'b'
						if (float(count_h)/gt_num) > 0.9:
							stackGT_string	+= '\t'+'h'
						else:
							stackGT_string	+= '\t'+'-'
					stack_list.append(stackGT_string)
					tempGT_list     = []
					stackGT_string  = ''
					temp_loc_string = loc_string
					temp_loc_cell   = temp_loc_string.split('\t')
					temp_loc_loci   = temp_loc_cell[0]
					temp_loc_linkage        = temp_loc_loci.split(':')[0]
					temp_loc_pos    = int(temp_loc_loci.split(':')[1])
					temp_loc_gt     = temp_loc_cell[1:]
					loci_string     = temp_loc_loci
					for gt in temp_loc_gt:
						tempGT_list.append(list(gt))
					continue
		else:
			info_list.append(loc_string)

	for info in info_list:
		if 'name' in info:
#			out_file.write(info +'\n')
			print info
		if 'popt' in info:
#			out_file.write(info +'\n')
			print info
		if 'nloc' in info:
#			out_file.write('nloc = '+str(len(stack_list)) +'\n')
			print 'nloc = '+str(len(stack_list))
		if 'nind' in info:
#			out_file.write(info +'\n')
			print info

	for item in stack_list:
#		out_file.write(item +'\n')
		print item
