To extract amplicon regions from full-length sequences using primersearch output:

1) Create a primer file in the required format for the tool (see included .tsv file), with forward and reverse primers only.

2) Upload the primer.tsv file, along with a nucleotide fasta file containing sequences to search within, to the Emboss
   Primersearch tool, using either a webtool (http://bioinfo.nhri.org.tw/cgi-bin/emboss/primersearch) or Galaxy installation. 

3) Set your allowed percent mistmatch between the primers and fasta sequences to be searched. Example output was obtained
  with an allowed percent mismatch of 20% and appears as follows:

	Primer name Flu_A_Univ
	Amplimer 1
		Sequence: KY925925  
		A/Santo Antonio da Patrulha/LACENRS-2621/2012 2012/08/10 7 (MP)
		GACC[GA]ATCCTGTCACCTCTGAC hits forward strand at 171 with 1 mismatches
		GGGCATT[TC]TGGACAAA[TG]CGTCTACG hits reverse strand at [753] with 1 mismatches
		Amplimer length: 105 bp
	Amplimer 2
		Sequence: KY925930  
		A/Novo Hamburgo/LACENRS-385/2016 2016/04/04 7 (MP)
		GACC[GA]ATCCTGTCACCTCTGAC hits forward strand at 171 with 1 mismatches
		GGGCATT[TC]TGGACAAA[TG]CGTCTACG hits reverse strand at [753] with 1 mismatches
		Amplimer length: 105 bp
		
4) Obtain primersearch output by running primersearch. This output will be used by extract_amplicon_from_primersearch_output.py
   to extract the theoretical amplicon regions from the full-length fasta sequences.