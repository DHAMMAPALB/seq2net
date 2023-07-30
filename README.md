# seq2net: a tool for evolutionary protein interactions prediction and analysis of bacterial genomes
This pipeline predicts the protein-protein interactions (PPIs) of bacterial genomes using evolutionary approaches. I present here a simple approach to predict bacterial PPIs using  Gene Neighborhood (GN), Phylogenetic Profile (PP), Gene Fusion (GF), Gene-Coevolution (GC) and Interlog (IN) methods. The consistent and already known PPIs are used to build protein interactions networks. A network is further analyzed to enrich nodes with functional annotations and pathways.

# Requirements
- Ncbi-blast+
- Clustalo
- Perl: List::Util subroutines
- Python modules: sys, requests, pandas, matplotlib, biopython 
- R packages: igraph
 
# Installation

```
sudo chmod +x -R seq2net/Scripts/*
echo 'export PATH="your-dir/seq2net/Scripts/:$PATH"' >> ~/.bashrc
source ~/.bashrc

```
# Usage

```
Usage: seq2net [-h] [-f1 <file1>] [-Ortho] [-GN] [-PP] [-f2 <file2>] [-OrthoPara] [-GF] [-GC] [-IN <PPIs_faa> <PPIs.txt>] [-r <reference_genome>] [-e <value1> [<value2>]] [-c <coverage1> [<coverage2>]] [-t <number>] [-kPPIs] [-enrich] [-o <outdir>]

positional arguments:
	-f1 <file1>			: Input file with a list of fasta (.faa) files for Ortho, GN and PP
	-Ortho				: Predict Orthologs of the reference genome
	-GN				: Predict PPIs using Gene Neighborhood method
	-PP				: Predict PPIs using Phylogenetic Profile method
	-f2 <file2>			: Input file with a list of fasta (.faa) files for OrthoPara, GF and GC
	-OrthoPara			: Prediction Orthologs and Paralogs of the reference genome
	-GF				: Predict PPIs using Gene Fusion method
	-GC				: Predict PPIs using Gene Coevolution method
	-IN <PPIs_faa> <PPIs.txt>	: Predict PPIs using Interlog method; requires input protein fasta sequences (.faa) and a list of known PPIs
	-r <reference_genome>		: Reference genome name

optional arguments:
	-e <value1> [<value2>]		: E-value for Orthologs or/and Paralogs respectively
	-c <coverage1> [<coverage2>]	: Sequence coverage for Orthologs or/and Paralogs respectively
	-t <number>			: Number of methods of PPIs being consistently predicted
	-kPPIs				: List the known PPIs from the predict PPIs using String database
	-enrich				: Enrichment of proteins in Panther Slim Gene Ontology and Pathways
	-o <outdir>			: Output directory
	-h, --help			: Shows this help message and exit

Others: Orthologs.pl <file1> <reference_genome>] <value1> <coverage1> <outdir>
	GN.pl <file1> <reference_genome> <outdir>
	PP.pl <reference_genome> <outdir>
	OrthoParalogs.pl <file1> <reference_genome>] <value1> <coverage1> <outdir>
	GF.pl <reference_genome> <outdir>
	GC.pl <file2> <reference_genome> <outdir>
	IN.pl <reference_genome> <PPIs_faa> <PPIs.txt> <value1> <coverage1> <outdir>
	Consistency.pl <number> <outdir>
	Networks.R <PPIs_output_file> <color>
	kPPIs.py <reference_genome> <PPIs_output_filename> <Known_PPIs_output_filename> <outdir>
	Enrich.py <reference_genome> <PPIs_output_filename> <outdir>
```

# Example

```
cd Example
seq2net -f1 Gseq_names1.txt -r Mycobacterium_tuberculosis_H37Rv -Ortho -GN -o Mtb_PINs
seq2net -f1 Gseq_names1.txt -r Mycobacterium_tuberculosis_H37Rv -Ortho -PP -o Mtb_PINs
```
or
```
seq2net -f1 Gseq_names1.txt -r Mycobacterium_tuberculosis_H37Rv -Ortho -GN -PP -o Mtb_PINs
```
Similarly,
```
seq2net -f2 Gseq_names2.txt -r Mycobacterium_tuberculosis_H37Rv -OrthoPara -GF -GC -o Mtb_PINs
```
```
seq2net -r Mycobacterium_tuberculosis_H37Rv -IN DbPPIs/DIP.faa DbPPIs/DIP_PPIs.txt -o Mtb_PINs
```
A single command for the entire analysis including known PPIs prediction and enrichment analysis would be
```
seq2net -f1 Gseq_names1.txt -Ortho -GN -PP -f2 Gseq_names2.txt -OrthoPara -GF -GC -IN DbPPIs/DIP.faa DbPPIs/DIP_PPIs.txt -r Mycobacterium_tuberculosis_H37Rv.faa -e 0.001 0.01 -c 40 30 -t 2 -kPPIs -enrich -o Mtb_PINs
```
* The reference genome name should be same as NCBI taxonomic name with white spaces given underscores. This will ensure the mapping of reference genome to NCBI Taxonomy databse for the identification of known PPIs from String database and for the gene ontologies and pathways from PANTHER classification system. Protein tables (NCBI format) are required to map the gene identifiers with their gene symbols or synonyms; the protein table names should be same as the names of fasta files but with the extension of csv.

#### The example figures:

[Figure 1.png](https://github.com/DHAMMAPALB/seq2net/issues/1)

Figure 1: Network of PPIs consistently predicted by atleast two methods.

[Figure 2.png](https://github.com/DHAMMAPALB/seq2net/issues/2)

Figure 2: Network of PPIs which were consistently predicted by atleast two methods and were already reported in String database.

[Figure 3.png](https://github.com/DHAMMAPALB/seq2net/issues/3)

Figure 3: Enrichment of network proteins in PANTHER GO-Slim Biological Process.

[Figure 4.png](https://github.com/DHAMMAPALB/seq2net/issues/4)

Figure 4: Enrichment of network proteins in PANTHER GO-Slim Molecular Function.

[Figure 5.png](https://github.com/DHAMMAPALB/seq2net/issues/5)

Figure 5: Enrichment of network proteins in PANTHER GO-Slim Cellular Component.

[Figure 6.png](https://github.com/DHAMMAPALB/seq2net/issues/6)

Figure 6: Enrichment of network proteins in PANTHER Pathways.


# References

- Burne D. Investigation of protein interaction networks in mycobacterium tuberculosis using computational approaches [Doctoral thesis]. The University of Hyderabad. 2019. http://hdl.handle.net/10603/342780, https://www.researchgate.net/publication/369943687_Investigation_of_protein_interaction_networks_in_mycobacterium_tuberculosis_using_computational_approaches
- Bharne D, Vindal V. In silico identification of highly interacting proteins in Mycobacterium tuberculosis H37Rv. Conference Proceedings, J Proteins Proteom. 2014;5(3):111
- Sievers F, Wilm A, Dineen D, et al. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol Syst Biol. 2011;7:539
- Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. https://igraph.org

