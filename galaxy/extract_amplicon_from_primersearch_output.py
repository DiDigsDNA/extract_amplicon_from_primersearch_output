#!/usr/bin/env python
'''Reads in full-length fasta sequences and corresponding output from the Emboss primersearch tool
for these sequences. Parses primersearch output for sequence id, description and primer locations and
uses these locations to extract amplicon regions from the full-length fasta sequences. If more than one
primersearch amplimer is present in the primersearch output, extracts sequence from the original amplicon
using only the first amplimer. Prints extracted amplicon sequences to fasta file.
Note: Fasta definition lines should be follow typical format: 
	e.g. >KY925925 A/Santo Antonio da Patrulha/LACENRS-2621/2012 2012/08/10 7 (MP)
Not unusual formats:
	e.g. >A_/_H1N1_|_A/WAKAYAMA/163/2016_|_MP|_|_970779
'''

'''usage: python extract_amplicon_from_primersearch_output.py fastaToParse.fasta results.primersearch 
		  output_filename'''
'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory,Feb 2018'''

import sys,string,os, time, Bio, re
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Data.IUPACData

fastaToParse = sys.argv[1] #fasta file of nucleotide sequences to parse
primersearch_results = sys.argv[2] #emboss primersearch results file
outputFastaHandle = sys.argv[3] + ".fasta" #from user=specified output name
outputFasta = open(outputFastaHandle, 'w') #output fasta containing sequences requiring further analysis

def extractIndex(hit_line_string):
	'''Extract and return hit position from Emboss primersearch forward or reverse strand hit line.'''
	regex = re.compile("strand\\ at\\ [\\[]{0,}[0-9]{1,}[\\]]{0,}\\ with")
	matchArray = regex.findall(hit_line_string)
	position = ''
	if len(matchArray)>0: #check for a match
		match = matchArray[0].replace('[','') #remove '[' and ']'
		match = match.replace(']','')
		digitList = [int(word) for word in match.split() if word.isdigit()]
		position = digitList[0]
		#print("Extracted position ", position)
	return position

def extractAmpliconLength(amplimer_length_line):
	'''Extract and return amplicon length from Emboss primersearch amplimer length line.'''
	amplicon_length = ''
	regex = re.compile("[0-9]{1,}")
	matchArray = regex.findall(amplimer_length_line)
	if len(matchArray)>0: #check for a match
		amplicon_length = int(matchArray[0])
	return amplicon_length

#read full-length fasta sequences into dict of SeqrRecords with key = record.id
with open(fastaToParse, 'r') as sequenceFile:
	full_length_dict = SeqIO.to_dict(SeqIO.parse(fastaToParse, "fasta", alphabet = IUPAC.ambiguous_dna))
	number_to_extract = len(full_length_dict) #store the number of full-length seqs

#parse primersearch results file
with open(primersearch_results, 'r') as inputFile:
	extracted_amplicons = [] #empty list of extracted amplicon sequence
	
	for line in inputFile: #read in each line
		if "Amplimer " in line.strip(): #grab number of amplier
			pattern = '[0-9]+$'
			results = re.search(pattern, line)
			if results:
				num = results.group()
		elif "Sequence: " in line: #grab sequence identifier
			id = line.replace("Sequence: ", "").strip()
			description = inputFile.readline().strip() #grab description from next line
		elif "hits forward strand at " in line:
			forward_hit_line = line.strip()
			reverse_hit_line = inputFile.readline().strip() #next line will be reverse hit
			forwardHitPosition = extractIndex(forward_hit_line) #determine F and R hit positions
			amplimer_length_line = inputFile.readline().strip() #grab next line
			amplicon_length = extractAmpliconLength(amplimer_length_line) #determine length of amplicon
			startIndex = forwardHitPosition-1 #primersearch results weren't zero-indexed
			endIndex = startIndex + amplicon_length
			all_extracted_amplicon_ids = [rec.id for rec in extracted_amplicons]
			#if id exists as a key in the dict and isn't already in extracted amplicons list
			if id in full_length_dict and (id not in all_extracted_amplicon_ids):
				full_length_record = full_length_dict[id] #created SeqRecord from amplicon and add to list
				extracted_sequence = full_length_record.seq[startIndex:endIndex]
				ampliconRec = SeqRecord(extracted_sequence,id = id, name = id, description = description)
				extracted_amplicons.append(ampliconRec)
	
	for rec in extracted_amplicons: #print extracted amplicon record id's and lengths to console
		print("Amplicon extracted for %s: %i bp" % (rec.id,len(rec.seq)))
	print(">>> %i amplicons extracted from %i full-length sequences" % (len(extracted_amplicons), len(full_length_dict)))
	SeqIO.write(extracted_amplicons, outputFastaHandle, "fasta")
	print(">>> Output file: %s" % (outputFastaHandle))
	
inputFile.close()
outputFasta.close()