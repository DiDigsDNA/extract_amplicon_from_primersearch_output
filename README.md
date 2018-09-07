# extract_amplicon_from_primersearch_output
Extracts PCR amplicons from full-length nucleotide fasta sequences, using Emboss primersearch output. Reads in full-length 
fasta sequences and corresponding Emboss primersearch output for the sequences. Parses primersearch output for sequence id 
and description and uses primer locations to extract amplicon regions from full-length fasta sequences. If more than one
amplimer is present in the primersearch output, sequence is extracted using the first amplimer. Prints extracted 
amplicon sequences to fasta file.
Note: Fasta definition lines should follow the typical format: 
	e.g. >KY925925 A/Santo Antonio da Patrulha/LACENRS-2621/2012 2012/08/10 7 (MP)
As opposed to:
	e.g. >A_/_H1N1_|_A/WAKAYAMA/163/2016_|_MP|_|_970779
