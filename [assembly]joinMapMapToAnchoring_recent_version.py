import sys
import numpy as np
from scipy import stats
from shim import fasta_dict

delim = sys.argv[1]

def contigGrouping(JM_list,reformed_JM_file):
	contig_count_dict	= {}
	contig_group_dict	= {}
	split_contig	= []
	for i in JM_list:
		string  = i.strip()
		if string       != '':
			if 'group' in string:
				key     = string.replace(' ','_')
			else:
				out_string	= key +'\t'+ '\t'.join(string.split()[0].split('_')[0:2]) +'\t'+ string.split()[1]
				open(reformed_JM_file,'a').write(out_string+'\n')
				contig	= string.split()[0].split('_')[0]
				if contig not in contig_count_dict.keys():
					contig_count_dict[contig]	= {key:1}
				else:
					if key in contig_count_dict[contig].keys():
						contig_count_dict[contig][key]	+= 1
					else:
						contig_count_dict[contig][key]	= 1
	open(reformed_JM_file).close()
	for j in contig_count_dict.keys():
		if len(contig_count_dict[j].keys())	== 1:
			contig_group_dict[j]	= contig_count_dict[j].keys()[0]
		else:
			split_contig.append(j)
			group	= contig_count_dict[j].keys()
			count	= [ float(contig_count_dict[j][x]) for x in contig_count_dict[j].keys()]
			tot	= np.array(count).sum()
			ratio	= np.array(count)/tot
			for k in range(len(ratio)):
				if ratio[k] >= 0.9:
					#print j +'\t'+ group[k]
					contig_group_dict[j] = group[k]
					split_contig.remove(j)
				else:
					continue
	
	group_contig_dict	= {}

	for i in contig_group_dict.keys():
		if contig_group_dict[i] not in group_contig_dict.keys():
			group_contig_dict[contig_group_dict[i]]	= [i]
		else:
			group_contig_dict[contig_group_dict[i]].append(i)

 	return group_contig_dict,split_contig

def getMapInfo(reformed_JM_file,group_contig_dict):
	JM_out_list	= open(reformed_JM_file,'r').readlines()
        dic_contig    = {}

        for jm in JM_out_list:
		cell	= jm.strip().split()
		group	= cell[0]
                contig	= cell[1]
                p_position      = int(cell[2])
                g_position      = float(cell[3])

		if contig in group_contig_dict[group]:
	                if contig in dic_contig.keys():
        	                dic_contig[contig][0].append(p_position)
				dic_contig[contig][1].append(g_position)
	                else:
				dic_contig[contig]  = [[p_position],[g_position]]
        return dic_contig

def getLinearRegression(dic_mapinfo):
        new_dic = dic_mapinfo
        info_dic       = {}
        for scaf in dic_mapinfo.keys():
                if len(dic_mapinfo[scaf][0])    <= 1:
                        p_pos_med       = dic_mapinfo[scaf][0][0]
                        g_pos_med       = dic_mapinfo[scaf][1][0]
                        slope   = 0
                else:
                        p_pos   = dic_mapinfo[scaf][0]
                        g_pos   = dic_mapinfo[scaf][1]
                        slope, intercept, r_value, p_value,std_err      = stats.linregress(p_pos,g_pos)
                        p_pos_med       = np.median(np.array(p_pos))
                        g_pos_med       = slope * p_pos_med + intercept

                #new_dic[scaf].append([slope, p_pos_med, g_pos_med])
                info_dic[scaf]	= [slope, p_pos_med, g_pos_med]
        return info_dic

def getSuperscaffolds(group,contig_of_group,contig_fa,info_dic_fromLinear):
        contig_dic        = fasta_dict(contig_fa)
        map_info_dic	= {x:info_dic_fromLinear[x] for x in contig_of_group}
	order_info_dic	= {info_dic_fromLinear[x][2]:[x,info_dic_fromLinear[x][0]] for x in contig_of_group}

        #print order_info_dic
        order_list      = order_info_dic.keys()
        order_list.sort()
        Linkage_group   = '>'+ group +'\n'
        number  = 0
        anchored_list   = []
        for order in order_list:
                number  += 1
                anchored_list.append(order_info_dic[order][0])
                #print order_info_dic[order][0] 
                if number       == len(order_list):
                        if order_info_dic[order][1] >= 0:
                                Linkage_group += contig_dic[order_info_dic[order][0]].upper()
                        else:
                                Linkage_group += revcom(complement(contig_dic[order_info_dic[order][0]].upper()))
                else:
                        if order_info_dic[order][1] >= 0:
                                Linkage_group += contig_dic[order_info_dic[order][0]].upper() + 'n'*500
                        else:
                                Linkage_group += revcom(complement(contig_dic[order_info_dic[order][0]].upper())) + 'n'*500
        return Linkage_group, order_info_dic, anchored_list

def complement(s):
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        letters = list(s)
        letters = [basecomplement[base] for base in letters]
        return ''.join(letters)

def revcom(s):
        result  = complement(s[::-1])
        return result

if __name__=="__main__":
	''' python script _(delimeter) JM_out JM_out_reformed(automatically_made) contig.fa file_out '''
	JM_list = open(sys.argv[2],'r').readlines()
	reformed_JM_out_file	= sys.argv[3]
	contig_fa	= open(sys.argv[4],'r').read()
	file_out	= open(sys.argv[5],'a')
	file_unanchored_out	= open(sys.argv[5]+'.unanchored.fa','a')
	log_out	= open(sys.argv[5]+'.log','a')
	group_contig_dict,split_contig	= contigGrouping(JM_list,reformed_JM_out_file)
	info_dic	= getLinearRegression(getMapInfo(reformed_JM_out_file,group_contig_dict))
	group	= group_contig_dict.keys()
	contig_dic	= fasta_dict(contig_fa)
	contig_key_list	= contig_dic.keys()

	for LG in group:
		contig_of_group	= group_contig_dict[LG]
		superscaffold, log_info, anchored_list	= getSuperscaffolds(LG,contig_of_group,contig_fa,info_dic)
		file_out.write(superscaffold +'\n')
		log_out.write(LG +'\t'+ str(log_info) +'\n')
		
		for anchored in anchored_list:
			if anchored.strip() in contig_key_list:
				contig_key_list.remove(anchored.strip())

	number_unachored	= 0
	for contig in contig_key_list:
		number_unachored	+= 1
		file_unanchored_out.write('>'+ contig +'\n'+ contig_dic[contig] +'\n')
	log_out.write('Unanchored scaffolds :'+'\t'+ str(number_unachored) +'\n'+ 'Duplicated :'+'\t'+ str(split_contig))
