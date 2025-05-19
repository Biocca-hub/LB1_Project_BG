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
The file downloaded directly from PDBeFold, needs some cleaning, because it has some not-necessary new lines and not all the amino acids are represented as upper case characters:  
```
>PDB:1bun:B STRUCTURE OF BETA2-BUNGAROTOXIN: POTASSIUM CHANNEL
------rkRhpdCD-KPPDT--KICqTVVRAFYYKPSAKRCVQFRYG-GCNgNGNHFKSDHLCRCECley
r
>PDB:1dtx:A CRYSTAL STRUCTURE OF ALPHA-DENDROTOXIN FROM THE GR
------epRrklCI-LHRNP--GRCyDKIPAFYYNQKKKQCERFDWSgCGG-NSNRFKTIEECRRTCig-
-
>PDB:1f5r:I RAT TRYPSINOGEN MUTANT COMPLEXED WITH BOVINE PANCR
--------RpdfCL-EPPYT--GPCkARIIRYFYNAKAGLCQTFVYGgCRA-KRNNFKSAEDCMRTCgg-
-
```
To build the HMM profile, the alignment file needs to edited to look like this:
```
>1BUN:B
------rkRhpdCD-KPPDT--KICqTVVRAFYYKPSAKRCVQFRYG-GCNgNGNHFKSDHLCRCECleyr
>1DTX:A
------epRrklCI-LHRNP--GRCyDKIPAFYYNQKKKQCERFDWSgCGG-NSNRFKTIEECRRTCig--
>1F5R:I
--------RpdfCL-EPPYT--GPCkARIIRYFYNAKAGLCQTFVYGgCRA-KRNNFKSAEDCMRTCgg--
```
# WRITE COMMANDS TO EDIT FILE
The cleaned alignment file (```pdb_kunitz_msa22.ali```) can now be used to build the HMM.

It is possible to review the files describing the results of the three alignments in this GitHub folder. You can find them in the ```msa_results``` directory.

---
## HMM Profile Construction based on Structural Alignment 
Use hmmbuild to generate an HMM profile from the alignment:

``` bash
hmmbuild pdb_kunitz.hmm pdb_kunitz_msa22_clean.ali
```
---
## Compilation of Validation Datasets
1. Collect **all proteins containing a kunitz domain** from the UniProtKB database (N = 18) and download the fasta file (e.g. ```kunitz_all.fasta```). Filters are:
    2. *PFAM id* = PF00014
    3. SwissProt reviewed
2. Collect all **not-human proteins containing a kunitz domain** from the UniProtKB (N = 380) and download the fasta file (e.g ```kunitz_not_human.fasta```). Filters are:
    1. Not human (NOT *Taxonomy [OC]* 9606)
    2. *PFAM id* = PF00014
    3. SwissProt reviewed
3. Merge the two Kunitz dataset files to form a unified collection of positive examples to test the HMM with (```kunitz_all.fasta```):
   ```bash
   cat kunitz_human.fasta kunitz_not_human.fasta > kunitz_all.fasta
   ```
5. Remove all sequences with a sequence identity ≥ 95% and a Nres ≥ 50 mapping the aligned sequences set on the positive dataset using the ```blastp``` command:
   ```bash
   # create a blast database with the kunitz proteins from UniProtKB/SwissProt
   makeblastdb -in kunitz_all.fasta -input_type fasta -dbtype prot -out kunitz_all.fasta
   ```
   ```bash
   # run a blastp search on the aligned sequences:
   blastp -query pdb_kunitz_msa22_clean.ali -db kunitz_all.fasta -out kunitz_pdb22.blast -outfmt 7
   ```
   ```bash
   # filter highly similar hits:
   grep -v "^#" kunitz_pdb22.blast | awk '{if ($3>=95 && $4>=50) print $2}' | sort -u > high_match_22.txt
   ```
   ```bash
   # create a file with the ids to remove:
   cut -d'|' -f2 high_match_22.txt > to_remove.txt
   ```
   ```bash
   # extract the unmatched ids for the proteins that will form the positive database:
   comm -23 <(sort kunitz_all.txt) <(sort to_remove.txt) > kunitz_final.txt 
   ```
6.  Collect all **SwissProt reviewed not-kunitz proteins** (573.230) and download the fasta file (e.g. ```uniprot_sprot.fasta```). Filters are:
    1. Not PFAM id PF00014
    2. SwissProt

---
