# HMM-Based Modeling of Kunitz Domains

This repository outlines the workflow necessary to generate and evaluate a Hidden Markov Model (HMM) based on multiple structural alignment of Kunitz domain-containing proteins. For this scope, the proteins were selected from the Protein Data Bank (PDB) database based on the availability of:
• high-resolution structural data 
• annotation in the reviewed SwissProt database. 
The aim of this work is to build a structure-informed HMM that captures conserved features of the Kunitz domain.

---

## Project Workflow
The workflow includes:
- Structural Data Acquisition and Preprocessing
- Multiple structural alignment
- HMM Profile Construction 
- Compilation of Validation Datasets
- HMM Scanning and Hit Mapping
- Performance Assessment
---

## Structural Data Acquisition and Preprocessing
#### Structural Data Collection:
A set of proteins (N = 160) characterized by the presence of a Kunitz domain and whose structure has been described and collected in the Protein Data Bank was identified performing an Advanced Search. The applied selection criteria were:  
• Data collection resolution ≤ 3.5 Å;    
• PFAM identifier = PF00014);  
• 45 ≤ Polymer Entity Sequence Length ≤ 80.   
A custom report (```rcsb_pdb_custom_report_20250511052702.csv```) was generated from the output of the Advanced Search, with the intent to store informations about:  
• PDB ID;  
• Total Number of Non-polymer Instances;  
• Refinement Resolu-tion (Å);  
• Data Collection Resolution;  
• Sequence;  
• Polymer Entity Sequence Length;  
• Auth Asym ID (the chain identifier);  
• Annotation Identi-fier (PFAM and GO);  
• Entity ID.  
#### Clustering:
A fasta file (```pdb_kunitz.fasta```) containing the PDB IDs and the sequences of the proteins in the collected dataset (N = 160) was generated from the report file running the following command:
```bash
tr -d '"' < rcsb_pdb_custom_report_20250511052702.csv | awk -F ',' '{if (length($2)>0) {name=$2}; print name,$6,$8,$9}' | grep PF00014 | awk '{print ">"$1"_"$3; print $2}' > pdb_kunitz.fasta
```
Using Cd-hit - a tool for clustering and filtering sequences based on sequence identity - with the standard threshold (sequence identity = 90%) allowed the identification of the set of proteins (N = 25) on which to perform the multiple structural alignment while avoiding redundancy. The clustering was performed running the following command: 
```bash
cd-hit -i pdb_kunitz.fasta -o pdb_kunitz_cluster.txt    
```
Two files were generated:  
• A ```.clstr``` file (e.g. ```pdb_kunitz_cluster.txt.clstr```)  
• A ```.txt``` file (e.g. ```pdb_kunitz_cluster.txt``` — your clustered sequences)  
The ```.clstr``` file contains raw cluster info, but its format may result hard to parse directly.
That’s why the ```clstr2txt``` command was applied to convert it into a tabular format:
```bash
clstr2txt.pl pdb_kunitz_cluster.txt.clstr > clusters_table.txt
```
## Multiple Structural Alingment
The multiple sequence alignment was performed through PDBeFold. To do so, a file listing the PDB identifiers and corresponding chain information for the 25 selected proteins (```pdb_kunitz_ids.txt```) was generated running this command:
```bash
grep '^>' pdb_kunitz_cluster.txt | sed 's/^>//' | sed 's/_/:/'  > pdb_kunitz_ids_25.txt
```
To run the alignment with PDBeFold:  
  
<img width="361" alt="Screenshot 2025-05-19 at 19 51 55" src="https://github.com/user-attachments/assets/2d3e0f2e-fcfe-44ef-a55b-1f2c35b6df44" />  
  
1. At this link [https://www.ebi.ac.uk/msd-srv/ssm/], click on ```Launch PDBeFold```
2. Flag ```Submission Form: multiple```;  
3. Select ```Source: List of PDB codes``` and then ```Choose file``` (e.g. ```pdb_kunitz_ids_25.txt```);
4. Click on ```Submit your query```.  
The output file was saved as ```pdb_kunitz_msa25.fasta``` (with 25 indicating the number of  aligned proteins).
Reviewing the results of the alignment, two proteins were taken out of the dataset:
- 2ODY:F (highest number of aligned residues (Nres = 127); lowest Q-score (0.2173));  
- 5JBT:Y (lowest number of aligned residues (Nres = 38); highest RMSD value (2.9238).

Following the same procedure, a second structural alignment was performed on the set of 23 proteins, using the ```pdb_kunitz_ids_23.txt``` file as input file (it was created manually editing the list given the low number of IDs to remove).  
From the analysis of results of this alignment (saved as ```pdb_kunitz_msa23.fasta```), the decision to remove the protein 4BQD:A was taken because of its higher length (Nres = 78) that led it to produce a large initial gap in the alignment.  
The third and final alignment was performed on a set of 22 proteins, using ```pdb_kunitz_ids_22.txt``` as input file. Its results are saved in the file ```pdb_kunitz_msa22.fasta```.

It is possible to review the files describing the results of the three alignments in this GitHub folder. You can find them in the ```msa_results``` directory.

---
## Build the HMM with HMMER
Use hmmbuild to generate an HMM profile from the alignment:

``` bash
hmmbuild pdb_kunitz.hmm pdb_kunitz_msa22_clean.ali
```
---
## Collecting sequences from UniProtKB to test the model performance 
1. Collect **all human proteins containing a kunitz domain** from the UniProtKB database (N = 18) and download the fasta file (e.g. ```kunitz_human.fasta```). Filters are:
    1. Human (*Taxonomy [OC]* = 9606)
    2. *PFAM id* = PF00014
    3. SwissProt reviewed
2. Collect all **not-human proteins containing a kunitz domain** from the UniProtKB (N = 380) and download the fasta file (e.g ```kunitz_not_human.fasta```). Filters are:
    1. Not human (NOT *Taxonomy [OC]* 9606)
    2. *PFAM id* = PF00014
    3. SwissProt reviewed
3. Collect all **not-kunitz proteins** (~) and download the fasta file (e.g. ```not_kunitz.fasta```). Filters are:
    1. Not PFAM id PF00014)
    2. SwissProt
---
