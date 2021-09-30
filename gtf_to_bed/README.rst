gtf_to_be
=======

Written by Guy Bottu for the GenePattern server of VIB BioinforlmaticsCore,
takes as input a GTF file and writes a BED file in 12 column format
with information about transcripts, for use with RSeqC.

The "thick" information is about the coding region, ideally it goes from
start codon to stop codon, but is information is lacking (e.g. because
of missing sequence or missing annotation), we use the CDS information.
For some transcripts there are multiple start or stop codons. We amways
choose the "thick" so that is has maximum length.

If there is no CDS information (as for ncRNA) the "thick" will have just a
repeat of the transcript start position, as per BED convention.

modified for integration under GenePattern

usage : perl gtf_to_bed.pl <GTF file> <output file>
