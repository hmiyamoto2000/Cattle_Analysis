# Methods: Genomic Analysis and Butyrate Biosynthetic Gene Cluster Identification

## Genomic Analysis and Biosynthetic Gene Cluster Prediction

The complete genome sequence of *Caldifermentibacillus hisashii* N11 (AP028807.1) was used for all genomic analyses. Biosynthetic gene clusters (BGCs) were predicted using antiSMASH (v8.0.4; https://antismash.secondarymetabolites.org). The butyrate biosynthetic gene cluster was identified by manual inspection of the genome sequence using the annotated FASTA format genome sequence (.fasta) and predicted protein sequences (.faa) derived from the registered genome (AP028807.1). 

Candidate genes encoding butyrate biosynthetic enzymes were identified by sequence similarity searches using tblastn against the registered genome sequence, and the results were cross-referenced with the genome annotation in the feature table. The analysis workflow consisted of the following steps:

**(1) Reference protein acquisition.** Characterized butyrate biosynthetic proteins were retrieved from the NCBI and UniProt databases, including enoyl-CoA hydratase (BER99916.1, BES00123.1), 3-hydroxybutyryl-CoA dehydrogenase (BES00124.1), acetyl-CoA C-acetyltransferase (BES00121.1), and CoA transferase subunit B (BES00125.1).

**(2) BLAST database construction.** A local BLAST nucleotide database was constructed from the *C. hisashii* genome sequence (AP028807.1) using the makeblastdb utility (NCBI BLAST+ suite).

**(3) Sequence similarity searches.** tblastn searches (e-value threshold < 1 × 10⁻⁵) were performed against the local genome database using reference butyrate biosynthetic protein sequences as queries to identify genomic sequences with high sequence similarity.

**(4) Genomic annotation cross-reference.** All identified genomic sequences were systematically cross-referenced with the NCBI feature table annotation (AP028807.1_feature_table.txt) to confirm gene structure, exact genomic coordinates (base pair positions), open reading frames (ORFs), and functional assignments based on automated annotation.

**(5) Proteome-level confirmation.** Identified genes were cross-referenced with available proteomic data to distinguish genome-encoded genes from those confirmed at the protein expression level.

**(6) Cluster organization analysis.** The genomic organization and synteny of identified butyrate biosynthetic genes were analyzed to determine chromosomal arrangement, gene order, intergenic spacing, and spatial clustering of functionally related genes.

**(7) Cluster visualization.** The genomic organization of the identified butyrate biosynthetic gene cluster was visualized using custom Python scripts (matplotlib library) to generate publication-quality diagrams showing gene orientation, sequence identity percentages, proteome confirmation status, and cluster boundaries.

The identified butyrate cluster comprised five genes located in a contiguous genomic region spanning approximately 4.8 kilobases, including BES00121.1 (acetyl-CoA C-acetyltransferase, 100% sequence identity), BES00122.1 (acyl-CoA dehydrogenase family protein, 100% identity), BES00123.1 (enoyl-CoA hydratase-related protein, 100% identity), BES00124.1 (3-hydroxybutyryl-CoA dehydrogenase, 100% identity), and BES00125.1 (CoA transferase subunit B, 99.5% identity). Genes BES00121.1 and BES00123.1 were confirmed at the proteome level, while the remaining genes were identified at the genomic level only.

[Procedure]
1st step
grep -B 3 -A 10 "target_name" AP028807_feature_table.txt

2nd step
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=BES00121.1&rettype=fasta&retmode=text" > acetyl_BES00121.faa
cat acetyl_BES00121.faa

curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=BES00122.1&rettype=fasta&retmode=text" > acylcoa_BES00122.faa
cat acylcoa_BES00122.faa

curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=BES00123.1&rettype=fasta&retmode=text" > enoyl_BES00123.faa
cat enoyl_BES00123.faa

curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=BES00124.1&rettype=fasta&retmode=text" > hydroxybutyryl_BES00124.faa
cat hydroxybutyryl_BES00124.faa

curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=BES00125.1&rettype=fasta&retmode=text" > ctfB_BES00125.faa
cat ctfB_BES00125.faa

curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=BER99643.1&rettype=fasta&retmode=text" > bfmBB_BER99643.faa
cat bfmBB_BER99643.faa

"curl ""https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=BER99646.1&rettype=fasta&retmode=text"" > lpdA2_BER99646.faa
cat lpdA2_BER99646.faa"                                

"curl ""https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=BER99297.1&rettype=fasta&retmode=text"" > DUF1796_BER99297.faa
cat DUF1796_BER99297.faa"                                

"curl ""https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=BER97788.1&rettype=fasta&retmode=text"" > daibat_BER97788.faa
cat daibat_BER97788.faa"                                

[Visualization]
    python3 plot_butyrate_cluster_v2.py (the code registered on the github)
