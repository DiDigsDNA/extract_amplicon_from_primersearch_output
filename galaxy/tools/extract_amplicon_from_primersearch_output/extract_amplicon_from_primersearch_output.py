#!/usr/bin/env python
'''Extracts PCR amplicons from full-length nucleotide fasta sequences, using Emboss primersearch output. 
Reads in full-length fasta sequences and corresponding Emboss primersearch output for the sequences. 
Parses primersearch output for sequence id and description and uses primer locations to extract amplicon
regions from full-length fasta sequences. If more than one amplimer is present in the primersearch output, 
sequence is extracted using the first amplimer. Prints extracted amplicon sequences to fasta file.
Note: Fasta definition lines should follow the typical format: 
	e.g. >KY925925 A/Santo Antonio da Patrulha/LACENRS-2621/2012 2012/08/10 7 (MP)
Instead of:
	e.g. >A_/_H1N1_|_A/WAKAYAMA/163/2016_|_MP|_|_970779
'''

'''usage: python extract_amplicon_from_primersearch_output.py fastaToParse.fasta results.primersearch 
		  output_filename'''
'''Author: Diane Eisler, Molecular Microbiology & Genomics, BCCDC Public Health Laboratory,Feb 2018'''

import sys, string, os, time, Bio, re
from Bio import Seq, SeqIO, SeqUtils, Alphabet, SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Bio.Data.IUPACData

fastaToParse = sys.argv[1]  # fasta file of nucleotide sequences to parse
primersearch_results = sys.argv[2]  # emboss primersearch results file
ampliconOutputHandle = "extractedAmplicons.fasta"  # amplicons extracted from primersearch output
failed_seqs_handle = "sequences_failing_amplicon_extraction.fasta"
extracted_amplicons = open(ampliconOutputHandle, 'w')  # output fasta containing sequences requiring further analysis
failed_sequences = open(failed_seqs_handle, 'w') #fasta of sequences failing amplicon extraction

def extractIndex(hit_line_string):
    '''Extract and return hit position from Emboss primersearch forward or reverse strand hit line.'''
    regex = re.compile("strand\\ at\\ [\\[]{0,}[0-9]{1,}[\\]]{0,}\\ with")
    matchArray = regex.findall(hit_line_string)
    position = ''
    if len(matchArray) > 0:  # check for a match
        match = matchArray[0].replace('[', '')  # remove '[' and ']'
        match = match.replace(']', '')
        digitList = [int(word) for word in match.split() if word.isdigit()]
        position = digitList[0]
    # print("Extracted position ", position)
    return position

def extractAmpliconLength(amplimer_length_line):
    '''Extract and return amplicon length from Emboss primersearch amplimer length line.'''
    amplicon_length = ''
    regex = re.compile("[0-9]{1,}")
    matchArray = regex.findall(amplimer_length_line)
    if len(matchArray) > 0:  # check for a match
        amplicon_length = int(matchArray[0])
    return amplicon_length

def determineFailedExtractions(full_length_dict, extracted_amplicon_list):
    unextracted_sequence_list = []
    print("--> determining failed extractions.....")
    for rec_id in full_length_dict:  # compare id of each full-length sequence
        extracted = False  # assume the sequence's amplicon region wasn't extracted
        for amplicon in extracted_amplicon_list:  # look at each amplicon extracted
            if rec_id == amplicon.id:  # if the rec in full-length list found in extracted_amplicon_list
                extracted = True  # we know its amplicon was extracted successfully
        if extracted == False:  # if the sequence failed amplicon extraction
            unextracted_sequence_list.append(rec_id)  # add it to list of unextracted sequences
    #print("%i full-length sequences failed amplicon extraction" % (len(unextracted_sequence_list)))
    return unextracted_sequence_list  # return the list of unextracted sequences

# read full-length fasta sequences into dict of SeqrRecords with key = record.id
with open(fastaToParse, 'r') as sequenceFile:
    full_length_dict = SeqIO.to_dict(SeqIO.parse(fastaToParse, "fasta", alphabet=IUPAC.ambiguous_dna))
    number_to_extract = len(full_length_dict)  # store the number of full-length seqs

# parse primersearch results file
with open(primersearch_results, 'r') as inputFile:
    extracted_amplicon_list = []  # empty list of extracted amplicon sequence

    for line in inputFile:  # read in each line
        if "Amplimer " in line.strip():  # grab number of amplier
            pattern = '[0-9]+$'
            results = re.search(pattern, line)
            if results:
                num = results.group()
        elif "Sequence: " in line:  # grab sequence identifier
            id = line.replace("Sequence: ", "").strip()
            description = inputFile.readline().strip()  # grab description from next line
        elif "hits forward strand at " in line:
            forward_hit_line = line.strip()
            reverse_hit_line = inputFile.readline().strip()  # next line will be reverse hit
            forwardHitPosition = extractIndex(forward_hit_line)  # determine F and R hit positions
            amplimer_length_line = inputFile.readline().strip()  # grab next line
            amplicon_length = extractAmpliconLength(amplimer_length_line)  # determine length of amplicon
            startIndex = forwardHitPosition - 1  # primersearch results weren't zero-indexed
            endIndex = startIndex + amplicon_length
            all_extracted_amplicon_ids = [rec.id for rec in extracted_amplicon_list]
            # if id exists as a key in the dict and isn't already in extracted amplicons list
            if id in full_length_dict and (id not in all_extracted_amplicon_ids):
                full_length_record = full_length_dict[id]  # created SeqRecord from amplicon and add to list
                extracted_sequence = full_length_record.seq[startIndex:endIndex]
                ampliconRec = SeqRecord(extracted_sequence, id=id, name=id, description=description)
                extracted_amplicon_list.append(ampliconRec)

    for rec in extracted_amplicon_list:  # print extracted amplicon record id's and lengths to console
        print("Amplicon extracted for %s: %i bp" % (rec.id, len(rec.seq)))
    # check that amplicons were extracted from all original full-length sequences
    if len(full_length_dict) == len(extracted_amplicon_list):
        print("Amplicons extracted from all (%i/%i) full-length sequences" % (len(extracted_amplicon_list), len(full_length_dict)))
    else:  # determine which sequences failed amplicon extraction and print to console
        failed_amplicon_extraction_list = determineFailedExtractions(full_length_dict, extracted_amplicon_list)
        print('''%i amplicons extracted from %i sequences (%i failed extraction)''' % (len(extracted_amplicon_list), len(full_length_dict), len(failed_amplicon_extraction_list)))
        SeqIO.write(failed_amplicon_extraction_list, failed_sequences, "fasta")
        print(">>> Failed to extract amplicons from: %s" % failed_seqs_handle)
    SeqIO.write(extracted_amplicon_list, extracted_amplicons, "fasta")
    print(">>> Extracted Amplicons: %s" % ampliconOutputHandle)

inputFile.close()
extracted_amplicons.close()
