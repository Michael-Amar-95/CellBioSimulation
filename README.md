# Restriction_enzyme_analysis

<ins>**Part 1**</ins>

Restriction Enzyme Analysis is a script that analyzes restriction enzyme target sites in bacterial gene sequences. It reads a file containing restriction enzymes and their target sequences, then scans multiple FASTA files containing bacterial gene sequences to count how many genes contain each target sequence.

The input consists of an enzymes file and one or more FASTA files. The enzymes file is a tab-delimited text file where each line contains an enzyme name followed by its target sequence. Example: "EcoRI GAATTC", "HindIII AAGCTT", "BamHI GGATCC". The FASTA files contain bacterial gene sequences in standard FASTA format, where each gene starts with a ">" symbol followed by its sequence.

The output prints the number of genes in each bacterial genome that contain the target sequence of each enzyme. Example output: "EcoRI found in 3 genes in genome1.fasta", "HindIII found in 5 genes in genome1.fasta", "BamHI found in 2 genes in genome2.fasta".

<ins>**Part 2**</ins>  

This part of the code defines a function named transcribe(dna_seq) that performs transcription of a given DNA sequence into its corresponding RNA sequence.

Background:
DNA is composed of four nucleotide bases: Adenine (A), Guanine (G), Thymine (T), and Cytosine (C). These bases pair in a specific way:

A pairs with T

G pairs with C

During transcription, the DNA sequence is converted into an RNA sequence based on base pairing rules, where Thymine (T) is replaced by Uracil (U) in the RNA. The RNA bases are: Adenine (A), Uracil (U), Guanine (G), and Cytosine (C). The pairing rules for transcription are:

* A → U

* T → A

* G → C

* C → G

Function requirements:

The function receives a DNA sequence (dna_seq) as input.

It returns the corresponding transcribed RNA sequence, maintaining the proper 5' to 3' directionality.

The output should be returned in uppercase letters.

You can assume the input is valid.

Example:
If the input DNA sequence is:

ATCAAG
The function will return:

UAGUUC

<ins>**Part 3**</ins>  

This section implements a function named translate(rna_seq) that simulates the biological process of translating an RNA sequence into a chain of codons, which correspond to amino acids in a protein.

Background:
Translation begins when a start codon (AUG, which codes for Methionine) is found in the RNA sequence. From that point onward, the sequence is read in triplets (codons), each corresponding to a specific amino acid.
Translation continues until one of the following conditions is met:

A stop codon (UAA, UAG, or UGA) is encountered.

The end of the RNA sequence is reached and fewer than 3 nucleotides remain.

Function requirements:

The function receives an RNA sequence (rna_seq) as input.

It scans all three possible reading frames and selects the one with the longest valid protein that starts from the first AUG found in that frame.

If multiple frames yield proteins of equal length, the first one encountered (lowest frame index) is selected.

If no valid protein is found (i.e., no AUG or no codons after it), the function returns None.

The returned value is a semicolon-separated string of codons in the selected reading frame, starting from the first AUG up to (but not including) the first stop codon, or the end of the sequence.

Example:
For the input RNA sequence:

- AUCAUGAACAUGCAGAUCAA
  
The function will return:

- AUG;AAC;AUG;CAG;AUC

---

# Cell_Simulation

This code is an extension of the "Restriction_enzyme_analysis.py" script.

It implements a multi-part Python program designed to simulate biological processes related to DNA and RNA. All functions include validation to ensure the input is biologically valid, such as checking that DNA sequences contain only the characters A, G, T, and C.

<ins>**Part 1**</ins>  

Enhances the function find_srr(dna_seq), which identifies all simple sequence repeats (SSRs) in a DNA sequence. An SSR is defined as a pattern of 1 to 6 nucleotides repeated at least 3 times consecutively. The function returns all possible SSRs, their counts, and the longest repeat found, or None if no SSRs exist.

Enhances the function transcribe(dna_seq), which converts a DNA sequence into its complementary RNA strand using base-pairing rules. The conversion respects the correct 5' to 3' direction and returns the result in uppercase. Input is validated to ensure only valid DNA bases are processed.

Enhances the function translate(rna_seq), which simulates translation of an RNA sequence into codons. The function finds the first start codon (AUG) and translates the sequence into codons until a stop codon or the end of the sequence. Among all reading frames, the longest valid protein is selected. If no valid protein is found, the function returns None.

<ins>**Part 2**</ins>  

Implements a class called Polymerase, which represents either a DNA or RNA polymerase enzyme but not both of them. The polymerase is initialized with a type and error rate. RNA polymerase can transcribe DNA into RNA using the transcribe function from part 2.

<ins>**Part 3**</ins>  

Implements a class called Ribosome, initialized with a codon-to-amino acid dictionary and a set of start codons. The ribosome translates RNA into protein using the longest valid open reading frame and the logic from the upgraded translate function.

<ins>**Part 4**</ins>  

Implements a Cell class, which contains attributes such as name, genome, number of copies, genetic code, start codons, and division rate. On initialization, the cell creates its own polymerases and ribosome. The class supports cell division via mitosis (returns identical cells) and meiosis (returns two cells with split or complementary genomes). It also includes a repertoire method that returns information on SSRs, RNA transcripts, and protein sequences derived from the genome.

<ins>**Part 5**</ins>  

Implements specific cell types by subclassing the Cell class, including ProkaryoticCell, StemCell, and NeuronCell. Each subclass has predefined properties.


Main program:
The main program takes command-line arguments including the cell type, number of division cycles, maximum number of cells, and the genomic sequences. It initializes a cell culture and simulates cell division while printing the final number of cells and the repertoire of each. It then attempts meiosis and prints the resulting genomes and their repertoires.

---

# Mutation_Simulator

The program is divided into two main parts.

<ins>**Part 1 (Mutation-enhanced cell simulation)**</ins>  

The Polymerase class has been upgraded to support random mutations during transcription. Given an error rate and a sequence length, a number of random mutations are introduced by selecting positions randomly and changing their base to another valid one, depending on whether it's a DNA or RNA polymerase.

A new class called MutantCell is introduced, inheriting from Cell. Mutant cells:

Track the number of accumulated mutations through generations.

Initially contain no mutations (generation zero).

During mitosis, produce one mutated child and one regular child.

When duplicated using the multiplication operator, return only the mutated offspring.

Another class, CancerCell, is introduced. If a MutantCell accumulates 10 or more mutations, it transforms into a CancerCell with a division rate 10 times higher than the original.

<ins>**Part 2 (Protein domain classification using PROSITE patterns)**</ins>  

This part of the code focuses on classifying protein sequences using domain patterns based on PROSITE syntax. These patterns are converted into Python regular expressions.

A function called prosite_to_python converts a dictionary of domain patterns into compiled Python regular expressions. Invalid patterns are either filtered out or raise a ValueError.

Another function, patterns_to_domains, receives a path to a pattern file and returns a dictionary mapping domain names to compiled regex patterns.

A class named SequenceClassifier is responsible for:

- Initializing with a pattern file

- Classifying a list of protein sequences using the known patterns

- Writing the classification results to a CSV file

- Each protein is annotated with all matching domains, separated by semicolons, or "Not Available" if no match is found

The main program receives a JSON configuration file with the following inputs:

- Path to pattern file

- Path to output file

- Number of division cycles

- Maximum number of cells

- Genomic sequences

The simulation begins with a single MutantCell and performs cell divisions up to the specified limits. At the end, it:

- Prints the original cell and the final number of cells

- Reports the total number of mutations

- Identifies and classifies the proteins found in all cells

- Produces a CSV file listing all proteins and their associated domains

The output varies slightly due to the randomness of mutation. A fixed random seed is used to support reproducible testing.

---

# Biological_Replication_Simulator 
This code is an extension of the other scripts and builds on the mutant and cancer cell simulation developed in them. The program includes full input validation, ensuring that DNA sequences only contain valid characters (A, G, T, and C), and is capable of running simulations using genomic sequences 

<ins>**Part 1 (Enhancing Input and Cell Initialization)**</ins>  

The input format is upgraded to accept a FASTA file containing multiple genomic sequences. The Biopython library is used to parse the sequences. The MutantCell class is extended to receive an additional parameter for error rate during initialization. This allows each cell to simulate transcription with a different mutation rate.

<ins>**Part 2 (Simulation of Mutant Cell Cultures)**</ins>  

The goal is to simulate multiple cell cultures with varying mutation rates and division cycles, and to measure the effects on biological outcomes. The simulation is repeated three times for each combination of parameters to account for randomness.

At the end of each simulation, the following results are collected:

- Number of cancerous cells

- Number of mutant cells that are not cancerous

- Number of distinct proteins found in the culture

The mutation rates tested range from 0.01 to 0.2, in steps of 0.01, including an additional run for 0.05. The number of division cycles ranges from 1 to 5.

All simulation results are saved in a CSV file. Each row in the CSV includes:

- Name of the genomic sequence (from the FASTA file)

- Mutation rate

- Number of division cycles

- Number of cancer cells

- Number of non-cancer mutant cells

- Number of unique proteins

Optional: additional user-defined metrics.

All output files and plots are saved in the same directory as the main script, using relative paths. If the CSV already exists, it is overwritten.

<ins>**Part 3 (Visualization)**</ins>  

Two sets of visualizations are generated.

For each genomic sequence, and for each of the three biological measures (cancer cells, mutant cells, protein count), a line graph is plotted showing:

- X-axis: mutation rate

- Y-axis: the biological measure

- Title: name of the sequence

A separate graph (or multiple graphs) shows how protein count varies across mutation rates and division cycles. The goal is to highlight the relationship between error rate and protein diversity in the culture.
