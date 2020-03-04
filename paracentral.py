# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 13:28:04 2020

paracentral - if the foveal, or central sight, describes a 5 degree arc the 
paracentral cone is the first periphery of the 8 degree arc which surronds the 
foveal cone.

_v0.2 Internal Beta (I fully approve of who I am, even as i get better.) 
    - output up, down flanking seqeunce of target as well as whole sequence

@author: pspealaman
"""

import argparse
import subprocess
import json

parser = argparse.ArgumentParser()

parser.add_argument('-i',"--input_file")
parser.add_argument('-o',"--output_file")

#python paracentral.py -g2f -i DGY2310_assembly.gfa -o DGY2310_assembly.fa
parser.add_argument('-g2f','--gfa_to_fa', action='store_true')

#python paracentral.py -b -i DGY2310_assembly.gfa -s LacZ.fsa -o test.file
parser.add_argument('-b','--blastn', action='store_true')
#parser.add_argument('-q','--query_fa')
parser.add_argument('-s','--subject_fa')
parser.add_argument('-f','--flanking_size')

args = parser.parse_args()

def output_handler(output):
    if len(output.strip()) > 1:
        print(output)

def gfa_to_fa(infile_name, outfile_name):
    infile = open(infile_name)
    outfile = open(outfile_name,'w')  
    
    for line in infile:
        if line[0]=='S':
            name = line.split('\t')[1].strip()
            sequence = line.split('\t')[2].strip()
    
            outline = ('>{}\n{}\n').format(name,sequence)
            outfile.write(outline)
    
    infile.close()
    outfile.close()

def gff_from_json(json_file_name, max_eval=5e-05):
    
    '''
    contig, contig_base = these are 'chromosome' hits in the subject fasta
    '''
    active = False
    query_deets_dict = {}
    contig_dict = {}
    lookup_dict = {}
    try:
        data = json.load(open(json_file_name))
        outfile = open(json_file_name.split('.json')[1]+('.gff'), 'w')
        
    except:
        print('json file ' + json_file_name + ' not found')
        return()


    for report_index in range(len(data["BlastOutput2"])):
        data_dict = (data["BlastOutput2"][report_index])
        for each_report in data_dict.items():
            for a_key, a_value in enumerate(each_report):
                if type(a_value)==dict:                   
                    for b_key, b_value in a_value.items():
                        if type(b_value)==dict:
                            for c_key, c_value in b_value.items():
                                if ('bl2seq' in c_key) and (type(c_value)==list):
                                    hit_dict = c_value[0]
                                    
                                    for d_key, d_value in hit_dict.items():                                        
                                        q_title = str(hit_dict['query_title'])
                                        is_str = str(hit_dict)
                                        if ('hits' in d_key) and (type(d_value)==list) and (len(d_value)>0) and 'No hits found' not in is_str:
                                            for each_hits in d_value:
                                                print('hits', d_value, str(each_hits['num']))
                                                for e_key, e_value in each_hits.items():
                                                
                                                    base = q_title + '.'+str(each_hits['num']) 
                                                    contig = each_hits['description']
                                                    contig = contig[0]
                                                    contig_dict[base] = str(contig['id'])
                                                        
                                                    if (e_key == 'hsps') and (type(e_value)==list):
                                                        for e_index in range(len(e_value)):
                                                            each_hsps = e_value[e_index]
    
                                                            numb = str(base)+'.'+str(each_hsps['num'])
                                                            
                                                            if len(numb)>1:
                                                                active = True
                                                                
                                                                hit_from = int(each_hsps["hit_from"])                                                                
                                                                hit_to = int(each_hsps["hit_to"])                                                                
                                                                query_from = str(each_hsps["query_from"])                                                                
                                                                query_to = str(each_hsps["query_to"])                                                            
                                                                bit_score = float(each_hsps["bit_score"])
                                                                evalue_score = float(each_hsps["evalue"])
                                                                query_strand = str(each_hsps["query_strand"])                                                                
                                                                hit_strand = str(each_hsps["hit_strand"])                                                                
                                                                qseq = str(each_hsps["qseq"])                                                            
                                                                hseq = str(each_hsps["hseq"])
                                                            
                                                            if evalue_score > max_eval:
                                                                active = False
                                                            
                                                            if active:
                                                                active = False
                                                                query_deets_dict[numb] = ['q_id','hit_from','hit_to','query_from','query_to','bit_score','query_strand','hit_strand','qseq', 'hseq','q_title']
                                                                query_deets_dict[numb][0] = base
                                                                query_deets_dict[numb][1] = hit_from
                                                                query_deets_dict[numb][2] = hit_to
                                                                query_deets_dict[numb][3] = query_from
                                                                query_deets_dict[numb][4] = query_to
                                                                query_deets_dict[numb][5] = bit_score
                                                                query_deets_dict[numb][6] = query_strand
                                                                query_deets_dict[numb][7] = hit_strand
                                                                query_deets_dict[numb][8] = qseq
                                                                query_deets_dict[numb][9] = hseq
                                                                query_deets_dict[numb][10] = q_title
                                                                print('is hit', numb, query_deets_dict[numb])
                                                                numb = 0
                                                                

    ct = 0    
    for numb, deets in query_deets_dict.items():
        contig = query_deets_dict[numb][10]
                
        #base = query_deets_dict[numb][0]
                
        #contig = contig_dict[base]
        
        hit_from = int(query_deets_dict[numb][1])
        hit_to = int(query_deets_dict[numb][2])
        
        q_from = int(query_deets_dict[numb][3])
        q_to = int(query_deets_dict[numb][4])
        
        if hit_from < hit_to:
            start = hit_from
            stop = hit_to
        else:
            start = hit_to
            stop = hit_from
                    
        if query_deets_dict[numb][7] == 'Plus':
            sign = '+'
        else:
            sign = '-'
            
        bit_score = float(query_deets_dict[numb][5])
        
        if query_deets_dict[numb][6] != query_deets_dict[numb][7]:
            orient = 'reverse'
        else:
            orient = 'forward'
            
        mod_seq = ('{}_{}').format(q_from, q_to)
        
        node_loci = ('{},{},{},{}').format(contig, query_deets_dict[numb][3], query_deets_dict[numb][4], float(query_deets_dict[numb][5]))
        gff_line = ('{}\terisapfel\tblastn_aligned\t{}\t{}\t.\t{}\t{}\tnode_uid={}; orient={}; from_to={}\n').format(contig, start, stop, sign, int(round(bit_score)), node_loci, orient, mod_seq)
        
        outfile.write(gff_line)
        
        lookup_dict[ct]={'contig': contig, 'hit_from':start, 'hit_to':stop, 'sign':sign, 'len':0}
        
        ct+=1

    return(lookup_dict)
    
def pull_sequences(infile_name):
    contig_seq = {}
    
    infile = open(infile_name)
    
    for line in infile:
        #print(line)
        if line[0] == 'S':
            contig_name = line.split('\t')[1].strip()
            seq = line.split('\t')[2].strip()
            contig_seq[contig_name]=seq

    return(contig_seq)
 
if args.gfa_to_fa:
    infile_name = args.input_file
    outfile_name = args.output_file
    gfa_to_fa(infile_name, outfile_name)
    
if args.blastn:
    infile_name = args.input_file
    fafile_name = args.input_file+'_temp.fa'
    
    if '/' in args.subject_fa:
        subject_name = args.subject_fa.rsplit('/',1)[1].split('.fa')[0]
    else:
        subject_name = args.subject_fa.split('.fa')[0]
    
    if args.flanking_size:
        flanking_size = int(args.flanking_size)
    else:
        flanking_size = 1000
        
    #Future Versions: A better structure would be to call this once instead of each time
    gfa_to_fa(infile_name, fafile_name)
    
    bashCommand = ('blastn -query {query_fa} -subject {subject_fa} -outfmt 15 -out paracentral_temp.json').format(query_fa=fafile_name, subject_fa=args.subject_fa)
    output_handler(subprocess.check_output([bashCommand],stderr=subprocess.STDOUT,shell=True))
    #print(bashCommand)
            
    lookup_dict = gff_from_json('paracentral_temp.json', max_eval=5e-05)
    
    if len(lookup_dict) > 0:
        outfile_name = args.output_file
        outfile = open(outfile_name,'w')
        
        contig_seq = pull_sequences(infile_name)
        
        for ct, deets in lookup_dict.items():
            print(ct, deets)
            contig = deets['contig']
            hit_from = int(deets['hit_from'])
            hit_to = int(deets['hit_to'])
            sign = deets['sign']
            
            if deets['len'] == 0:
                deets['len'] = len(contig_seq[contig])
            
            contig_length = deets['len']
            
            if contig in contig_seq:
                align_name = ('{}_{}_{}').format(subject_name, contig, ct)
                core_seq = contig_seq[contig][hit_from:hit_to+1]
                
                if sign == '+':
                    upstream_seq = contig_seq[contig][max(hit_from-flanking_size,1):min(contig_length, hit_to-1)]
                    downstream_seq = contig_seq[contig][max(hit_from+1, 1):min(contig_length, flanking_size+hit_to+1)]
                    
                else:
                    upstream_seq = contig_seq[contig][max(hit_from+1, 1):min(flanking_size+hit_to+1,contig_length)]
                    downstream_seq = contig_seq[contig][max(hit_from-flanking_size, 1):min(hit_to-1, contig_length)]  
                    
                
                whole_seq = contig_seq[contig][max(hit_from-flanking_size,1):min(flanking_size+hit_to+1,contig_length)]
                    
                outline = ('>{align_name}_upstream\n{up_seq}\n>{align_name}_downstream\n{dn_seq}\{align_name}_whole_fragment\n{wh_seq}\n').format(align_name=align_name, up_seq=upstream_seq, dn_seq=downstream_seq, wh_seq=whole_seq)
                outfile.write(outline)
        
        outfile.close()
        
    
