#!/usr/bin/env python

import operator
import re
from Bio.Blast import NCBIXML
import argparse


## currently optimized for python2.7. Have not tested this with python3 yet


#This script will parse an xml blast file and output a tsv file with the columns denoting query name, length, number of alignments, alignment name, HSP evalue, and HSP length
#recall that there can be multiple HSPs per alignment, and alignment is a hit on a sibject sequence
#HSPs are high scoring pairs within an alignment
#The output will list the query name for each HSP, but will not repeat the other information
#different query sequences are separated by a divider
#If the seqs flag is specified, the program will output the HSP query and reference sequences in the 7th and 8th columns, respectively

#command line arguments
parser=argparse.ArgumentParser(description='This script will parse a blast output file')
parser.add_argument('-infile', type=str, required=True, help='full path blast results file in xml format')
parser.add_argument('-outfile', type=str, required=True, help='full path name of parsed output file')
parser.add_argument('-seqs', action='store_true', help='add option if you want the program to output the HSP sequence') #Note that we now specify a new keyword, action, and give it the value "store_true". This means that, if the option is specified, assign the value True to args.verbose. Not specifying it implies False.
args=parser.parse_args()
(infile, outfile, seqs)=(args.infile, args.outfile, args.seqs)


#make file handle of blast results
result_handle=open(infile)

#make iterator
blast_records=NCBIXML.parse(result_handle)

outfile_handle=open(outfile, "w")

if seqs:
	outfile_handle.write("Query\tQuery_len\tNum_alignments\tAlignment_hit\tHSP_evalue\tHSP_length\tHSP_query_seq\tHSP_subject_seq\n")
else:
	outfile_handle.write("Query\tQuery_len\tNum_alignments\tAlignment_hit\tHSP_evalue\tHSP_length\n")

for blast_record in blast_records:
	query=blast_record.query
	qlen=blast_record.query_length
	#change spaces to underscores in query
	query=query.replace(" ", "_")
	#alignmetns is a list
	num_alignments=len(blast_record.alignments)
	if num_alignments==0:
		if seqs:
			outfile_handle.write(query+"\t"+str(qlen)+"\t"+str(num_alignments)+"\t"+"NA\tNA\tNA\tNA\tNA\n")
		else:
			outfile_handle.write(query+"\t"+str(qlen)+"\t"+str(num_alignments)+"\t"+"NA\tNA\tNA\n")
		divider= "~"*120
		outfile_handle.write(divider+"\n") 
		continue

	first_alignment_iteration=True
	for alignment in blast_record.alignments:
		#the alignment title format is odd, I need to parse it to get the correct name
		match=re.match(r'(gnl\|BL_ORD_ID\|\d+\s)(.+)', alignment.title)
		parsed_alignment=match.group(2)
		parsed_alignment=parsed_alignment.replace(" ", "_")
		qlen_placeholder=" "*len(str(qlen))
		num_alignments_placeholder=" "*len(str(num_alignments))
		parsed_alignment_placeholder=" "*len(parsed_alignment)		
		first_hsp_iteration=True
		for hsp in alignment.hsps:
			if first_hsp_iteration==True and first_alignment_iteration==True:
				if seqs:
					outfile_handle.write(query+"\t"+str(qlen)+"\t"+str(num_alignments)+"\t"+parsed_alignment+"\t"+str(hsp.expect)+"\t"+str(hsp.align_length)+"\t"+hsp.query+"\t"+hsp.sbjct+"\n")
				else:
					outfile_handle.write(query+"\t"+str(qlen)+"\t"+str(num_alignments)+"\t"+parsed_alignment+"\t"+str(hsp.expect)+"\t"+str(hsp.align_length)+"\n")
				first_hsp_iteration=False
				first_alignment_iteration=False
			elif first_hsp_iteration==True and first_alignment_iteration==False:
				if seqs:
					outfile_handle.write(query+"\t"+qlen_placeholder+"\t"+num_alignments_placeholder+"\t"+parsed_alignment+"\t"+str(hsp.expect)+"\t"+str(hsp.align_length)+"\t"+hsp.query+"\t"+hsp.sbjct+"\n")
				else:
					outfile_handle.write(query+"\t"+qlen_placeholder+"\t"+num_alignments_placeholder+"\t"+parsed_alignment+"\t"+str(hsp.expect)+"\t"+str(hsp.align_length)+"\n")
				first_hsp_iteration=False
			else:
				if seqs:
					outfile_handle.write(query+"\t"+qlen_placeholder+"\t"+num_alignments_placeholder+"\t"+parsed_alignment_placeholder+"\t"+str(hsp.expect)+"\t"+str(hsp.align_length)+"\t"+hsp.query+"\t"+hsp.sbjct+"\n")
				else:	
					outfile_handle.write(query+"\t"+qlen_placeholder+"\t"+num_alignments_placeholder+"\t"+parsed_alignment_placeholder+"\t"+str(hsp.expect)+"\t"+str(hsp.align_length)+"\n")

	divider= "~"*120
	outfile_handle.write(divider+"\n")
outfile_handle.close()
result_handle.close()			
