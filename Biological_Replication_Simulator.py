import csv
import random
import re
import sys
from itertools import repeat
from multiprocessing import Pool
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import os.path

amino_acid_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
prosite_letters = ['x', '[', ']', '{', '}', '-', '(', ')', '<', '>']  # and all numbers


def integer_checker(some_char):
    assert some_char.isdigit(), 'content inside () need to be numeric'


# Python3 code to Check for
# balanced parentheses in an expression
open_list = ["[", "{", "("]
close_list = ["]", "}", ")"]

"""Check if the prosite have illegal chars"""


def isValid(prosite):
    for char in prosite:
        # Check if the prosite have illegal chars
        if char in amino_acid_list or char in prosite_letters:
            continue
        else:
            return False
    return True


# Function to check parentheses
def checkParentheses(myStr):
    stack = []
    for m in myStr:
        if m in open_list:
            stack.append(m)
        elif m in close_list:
            pos = close_list.index(m)
            if ((len(stack) > 0) and
                    (open_list[pos] == stack[len(stack) - 1])):
                stack.pop()
            else:
                return "Unbalanced"
    if len(stack) == 0:
        return "Balanced"
    else:
        return "Unbalanced"


regex_prosite = "<?(x|[a-zA-Z]|\[[a-zA-Z]+\]|\{[a-zA-Z]+\})(\((0|[1-9]{1}[0-9]*)\)|\(([1-9]{1}[0-9]*|0)," \
                "{1}[1-9]{1}[0-9]*\))?(-(x|[a-zA-Z]|\[[a-zA-Z]+\]|" \
                "\{[a-zA-Z]+\})(\((0|[1-9]{1}[0-9]*)\)|\(([1-9]{1}[0-9]*|0),{1}[1-9]{1}[0-9]*\))?)*>?"


class SequenceClassifier:
    """Constructor"""

    def __init__(self, pattern_file):
        self.pattern_file = self.__patterns_to_domains(pattern_file)

    """This function read the data from the CSV file"""

    def __patterns_to_domains(self, pattern_file):
        # open the file
        with open(pattern_file, 'r') as file:
            # convert the csv to dictionary
            csv_file = csv.DictReader(file)
            csv_file_dict = {}
            for row in csv_file:
                csv_file_dict[row.get('Pattern')] = row.get('Domain')
            converted_dict = self.__prosite_to_python(csv_file_dict)
        return converted_dict

    """this function get dictionary in pattern site format and returns a new dictionary with keys in python format"""

    def __prosite_to_python(self, pattern_dict):
        dict_prosite_to_python = {}
        # iterate the prosite
        for key in pattern_dict.keys():
            # check the chars in the prosite. Need to be illegal char ot=r an Amino acid
            if not isValid(key):
                continue
            x = re.fullmatch(regex_prosite, key)
            if not x or x[0] != key:
                continue
            python_pattern = key
            # replace the char with the correct char
            python_pattern = python_pattern.replace('-', '')
            python_pattern = python_pattern.replace('x', '.')
            python_pattern = python_pattern.replace('{', '[^')
            python_pattern = python_pattern.replace('}', ']')
            python_pattern = python_pattern.replace('(', '{')
            python_pattern = python_pattern.replace(')', '}')

            # check the valid of the RE
            try:
                re.compile(python_pattern)
            except re.error as error:
                raise ValueError(error)

            dict_prosite_to_python[python_pattern] = pattern_dict.get(key)
        return dict_prosite_to_python

    """This function create csv file with the list of strings """

    def classify(self, seq_list, csv_file):
        # remove all the duplicates from the sequence list
        seq_list = list(dict.fromkeys(seq_list))
        with open(csv_file, 'w', newline='') as file:
            # create the name of the fields
            fieldnames = ['Sequence', 'Domains', ]
            # convert the csv to dic
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            # add titles to file
            writer.writerow({'Sequence': 'Sequence', 'Domains': 'Domains'})
            patterns_keys_dict = self.pattern_file.keys()
            # iterate on the seq
            for seq in seq_list:
                str_domain = ''
                # iterate the regex
                for regex in patterns_keys_dict:
                    # find the pattern
                    result = re.findall(regex, seq)
                    if result:
                        # check if we already have this domain on this seq
                        if str_domain and convert_list_to_str(result) in str_domain:
                            str_domain += ';'
                        str_domain += self.pattern_file.get(regex)
                if not str_domain:
                    writer.writerow({'Sequence': seq, 'Domains': 'NA'})
                else:
                    writer.writerow({'Sequence': seq, 'Domains': str_domain})


"""This class represent polymerase"""


class Polymerase:
    """This function initialize the polymerase"""

    def __init__(self, type, error_rate=0):
        # check the validation of the type
        assert type.isalpha(), 'type of polymerase need to be DNA or RNA '
        type = type.upper()
        assert type == 'DNA' or type == 'RNA', 'type of polymerase need to be DNA or RNA'
        self.type = type
        # check if the error rate is between 0 to 1
        assert 0 <= error_rate <= 1, 'error rate need to be between 0 to 1'
        self.error_rate = error_rate

    """This function transcribe the sequence to rna in the sequence."""

    def transcribe(self, dna_seq):
        """
        :param dna_seq: The dna that the user want to transcribe.
        :return: An rna from the dna.
        """
        dna_seq = dna_seq.upper()
        # Swap the string.
        temp_seq = dna_seq[::-1]
        # This string will contain the rna
        comp_dna = ''
        # creat a key - values fairs
        # check the type of the polymerase
        if self.error_rate == 0:
            if self.type == 'DNA':
                dic_of_fairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            else:
                dic_of_fairs = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        else:
            if self.type == 'DNA':
                temp_seq = dna_seq
                dic_of_fairs = {'A': 'A', 'T': 'T', 'C': 'C', 'G': 'G'}
            else:
                dic_of_fairs = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
            # iterate on the sequence
        for letter in temp_seq:
            comp_dna += dic_of_fairs[letter]

        if self.error_rate > 0:
            # calculate how much mutation i need to get
            k = np.ceil(len(comp_dna) * self.error_rate)
            # get k random index
            mutation_index = random.sample(range(len(comp_dna)), int(k))
            comp_temp_list = convert_str_to_list(comp_dna)

            # run on the indexes
            for index in mutation_index:
                # check if we have dna or rna
                if self.type == 'DNA':
                    nucleotide = ['A', 'C', 'G', 'T']
                    # delete the current nucleotide we have
                    nucleotide.remove(comp_temp_list[index])
                    comp_temp_list[index] = random.choice(nucleotide)
                else:
                    nucleotide = ['A', 'C', 'G', 'U']
                    # delete the current nucleotide we have
                    nucleotide.remove(comp_temp_list[index])
                    comp_temp_list[index] = random.choice(nucleotide)
            comp_dna = convert_list_to_str(comp_temp_list)
        return comp_dna


"""This function convert string to list"""


def convert_str_to_list(string):
    list1 = []
    list1[:0] = string
    return list1


"""This function convert list to string"""


def convert_list_to_str(list_before_convert):
    str_from_list = ""
    return str_from_list.join(list_before_convert)


"""This class represent ribosome"""


class Ribosome:
    """This function initialize the ribosome"""

    def __init__(self, genetic_code, start_codons):
        self.genetic_code = genetic_code
        # save the starts codon.
        self.start_codons = start_codons

    """This function translate the rna"""

    def translate(self, rna_seq):
        """
        :param rna_seq: The rna we got from the transcribe function.
        :return: Return the longest protein, if exist.
        """
        # check the valid of rna seq
        char_list = ['A', 'U', 'G', 'C']
        matched_list = [characters in char_list for characters in rna_seq]
        assert all(matched_list), 'RNA seq need to contain only A,U,G,C chars'
        # Make all the char to be uppercase.
        rna_seq.upper()
        # Save the length of the rna
        n = len(rna_seq)
        # Initialize an empty strings.
        current_seq = {}
        longest_seq = {}
        # This loop run on the index of the rna.
        for i in range(n):
            # Search for a start codon.
            if rna_seq[i: i + 3] in self.start_codons:
                # Save the protein that we found.
                current_seq = self.find_seq(rna_seq[i:n])
            # Save the longest protein between the previous and the one that we found in this iteration.
            if len(current_seq) > len(longest_seq):
                longest_seq = current_seq
        # Check if we found a protein.
        # If we didnt we send None, o.s, we send the protein.
        if longest_seq:
            list_of_codons = []
            for element in longest_seq:
                list_of_codons.append(longest_seq[element])
            return list_of_codons
        else:
            return None

    """This function help me to find frame of the translate"""

    def find_seq(self, rna_seq):
        """
        :param rna_seq: The sequence that start from the start codon we found in this iteration.
        :return: The protein.
        """
        # initialize and save the start codon
        seq = {0: rna_seq[0: 3]}
        # Save the length of the current sequence.
        n = len(rna_seq)
        # This loop jump on each set of three in the frame and check if we have a stop codon.
        # If its a stop codon its will stop, o.w it will continue until the last three that logic.
        for i in range(1, (n // 3)):
            # look if there is an nucleotide to add
            if self.genetic_code[rna_seq[(i * 3):(i * 3) + 3]]:
                seq[i] = rna_seq[i * 3: (i * 3) + 3]
            # o.w we have a stop codon
            else:
                # you said we dont need to return the stop codon
                return seq
        return seq

    """This function synthesize the protein."""

    def synthesize(self, rna_seq):
        # get the list of the nucleotides
        list_of_codons = self.translate(rna_seq)
        # check if there any amino acid
        if list_of_codons:
            protein = ''
            for codon in list_of_codons:
                # add the amino acid to the string
                protein += self.genetic_code[codon]
            # return the protein
            return protein
        # o.w return none
        else:
            return None


"""This class represent cell"""


class Cell:
    """This function initialize the cell"""

    def __init__(self, name, genome, num_copies, genetic_code, start_codons, error_rate, division_rate):
        self.name = name
        # a list of dna sequences
        self.genome = genome
        # check the valid of copies
        assert num_copies > 0 and isinstance(num_copies, int), 'num fo copies need to be a positive integer'
        self.num_copies = num_copies
        self.genetic_code = genetic_code
        self.start_codons = start_codons
        # check the valid of division rate
        assert division_rate > 1 and isinstance(division_rate, int), 'num fo division rate need to be an ' \
                                                                     'integer bigger than 1'
        self.division_rate = division_rate
        # initialize the polymerases
        self.polDna = Polymerase('DNA', error_rate)
        self.polRna = Polymerase('RNA', 0)
        # initialize the ribosome
        self.ribo = Ribosome(self.genetic_code, self.start_codons)
        self.error_rate = error_rate

    """This function duplicate the cell by mitosis"""

    def mitosis(self):
        # create an array of cells
        cells = [self] * self.division_rate
        # duplicate the cell respectively to the division rat
        # return cells
        return cells

    """This function duplicate the cell by meiosis"""

    def meiosis(self):
        # check if we have an even number of copies
        if self.num_copies % 2 == 0:
            # create an array of cells
            cells = []
            # this is the cell with the origin genome
            cell_origin = Cell(self.name, self.genome, int(self.num_copies / 2)
                               , self.genetic_code, self.start_codons, self.error_rate, self.division_rate)
            # create complementary seq
            complementary_seq = []
            for sequence in self.genome:
                # transcribe the seq and append to the list
                complementary_seq.append(self.polDna.transcribe(sequence))
            # this is the cell with the complementary seq
            cell_complementary = Cell(self.name, complementary_seq, int(self.num_copies / 2)
                                      , self.genetic_code, self.start_codons, self.error_rate, self.division_rate)
            # add the cells to the list
            cells.append(cell_origin)
            cells.append(cell_complementary)
            # return the cells
            return cells
        else:
            return None

    """"This function send the repertoire of the cell"""

    def repertoire(self):
        data = []
        for sequence in self.genome:
            # load the ssr in the sequence
            srr_seq = self.find_srr(sequence)
            # check if there are a srr.
            if srr_seq is None:
                srr_seq = 'No simple repeats in DNA sequence'
            # take the dna and transcribe it to a rna
            rna_seq = self.polRna.transcribe(sequence)
            # create protein
            protein = self.ribo.synthesize(rna_seq)
            # look for the longest protein that exist. If their isn't protein i will inform the user.
            if protein is None:
                protein = 'Non-coding RNA'
            data.append((srr_seq, rna_seq, protein))
        return data

    """This function return details about the cell"""

    def __repr__(self):
        return '<{}, {}, {}>'.format(self.name, str(self.num_copies), str(self.division_rate))

    """This function look for srr in the sequence."""

    def find_srr(self, dna_seq):
        """
        :param dna_seq: the dna i get from the user
        :return: All the srr that exist in the sequence.
        """
        n = len(dna_seq)
        # Initialize an empty string
        srr_dic = {}
        # This loop limit the length of the repeated sequence.
        for i in range(1, 7):
            # This loop look for the srr
            for m in range(n - (3 * i) + 1):
                # i want to check if this seq repeated in the strain.
                seq = dna_seq[m:m + i]
                # This var help me to count the repeated seq.
                flag = 0
                # I initialize the seq to be from the current index until the end of the strain.
                temp_seq = dna_seq[m:n + 1]
                # Count how much time it is repeated.
                while temp_seq.startswith(seq) and flag * i < n:
                    flag = flag + 1
                    temp_seq = dna_seq[m + i * flag: n]
                else:
                    # i check if the seq repeat it self more than 3 times .
                    if flag >= 3:
                        # i check if the srr is already in the dictionary
                        if seq in srr_dic:
                            # if the current value is lower than flag i pass
                            if srr_dic[seq] >= flag:
                                continue
                            # o.w i update the value - number of repeat
                            else:
                                srr_dic[seq] = flag
                        # if the seq is not in the dic i put him
                        else:
                            srr_dic[seq] = flag
        # if there are srr.
        if srr_dic:
            # i sort the dictionary`s keys by alphabetic
            srr_list = sorted(srr_dic)
            output = ''
            #     # put the dic in a string
            for item in srr_list:
                output += item + ',' + str(srr_dic[item]) + ';'
            return output[0: len(output) - 1]
        # o.w return none
        else:
            return None


class ProkaryoticCell(Cell):
    """Constructor"""

    def __init__(self, genome):
        start_codons = ['AUG', 'GUG', 'UUG']
        division_rate = 4
        num_copies = 1
        name = 'ProkaryoticCell'
        genetic_code = {
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
            'UAC': 'Y', 'UAU': 'Y', 'UAA': None, 'UAG': None,
            'UGC': 'C', 'UGU': 'C', 'UGA': 'U', 'UGG': 'W'}
        super().__init__(name, genome, num_copies, genetic_code, start_codons, 0, division_rate)


"""This class represent eukaryotic cell"""


class EukaryoticCell(Cell):
    """initialize the eukaryotic"""

    def __init__(self, name, genome, error_rate, division_rate):
        num_copies = 2
        start_codons = ['AUG']
        genetic_code = {
            'AUA': 'I', 'AUC': 'I', 'AUU': 'I', 'AUG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
            'AAC': 'N', 'AAU': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGU': 'S', 'AGA': 'R', 'AGG': 'R',
            'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
            'CAC': 'H', 'CAU': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
            'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
            'GAC': 'D', 'GAU': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
            'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
            'UUC': 'F', 'UUU': 'F', 'UUA': 'L', 'UUG': 'L',
            'UAC': 'Y', 'UAU': 'Y', 'UAA': None, 'UAG': None,
            'UGC': 'C', 'UGU': 'C', 'UGA': None, 'UGG': 'W'}
        super().__init__(name, genome, num_copies, genetic_code, start_codons, error_rate,
                         division_rate)


"""This class represent neuron"""


class NeuronCell(EukaryoticCell):
    """initialize the neuron"""

    def __init__(self, genome):
        # self.division_rate = 2
        super().__init__('NeuronCell', genome, 0, 2)


"""This class represent stem cell"""


class StemCell(EukaryoticCell):
    """Constructor for mutant and cancer cells"""

    def __init__(self, name, genome, error_rate, division_rate=3):
        super().__init__(name, genome, error_rate, division_rate)


"""This class represent mutant cell"""


class MutantCell(StemCell):
    """initialize the mutant cell"""

    def __init__(self, genome, num_mutations=0, error_rate=0.05):
        if num_mutations > 10:
            super().__init__('CancerCell', genome, error_rate, 10)
        else:
            super().__init__('MutantCell', genome, error_rate, 3)
        self.num_mutations = num_mutations

    def mitosis(self):
        # create the correct cell
        if self.name == 'MutantCell':
            cell = MutantCell(self.genome, self.num_mutations)
        # elif self.name == 'CancerCell':
        else:
            cell = CancerCell(self.genome, self.num_mutations)
        # duplicate the cell respectively to the division rate
        cells = (cell.__mul__(cell.division_rate - 1))
        cells.insert(0, cell)
        return cells

    """This function return all the cells that was created from this specific cell expect from it"""

    def __mul__(self, division_rate):
        # create an array of cells
        cells = []
        # how much mutation i have before the mitosis
        num_mutations = self.num_mutations
        all_genome = []
        # create the mutation genome
        for sequence in self.genome:
            mutation_genome = self.polDna.transcribe(sequence)
            # add the genome
            all_genome.append(mutation_genome)
            # update the number of mutation
            num_mutations += np.ceil(len(sequence) * self.error_rate)
        if num_mutations >= 10:
            new_cell = CancerCell(all_genome, num_mutations)
        else:
            new_cell = MutantCell(all_genome, num_mutations)
        # multiply the same cell as time as division rate
        for i in range(division_rate):
            cells.append(new_cell)
        return cells


class CancerCell(MutantCell):
    """initialize the mutant cell"""

    def __init__(self, genome, num_mutations):
        super().__init__(genome, num_mutations)


"""This function help me to find the max mutation"""


def max_mutation(cell_list):
    max_mut = 0
    # look for the cell that contain the max num of mutation
    for cell in cell_list:
        if cell.num_mutations > max_mut:
            max_mut = cell.num_mutations
    return int(max_mut)


"""This function do the mitosis cycle for any group of cells using to  using the limitation of max cycles  """


def mitosis_cycles(cells_list, max_cycles):
    number_of_cells = 1
    current_division_cycles = 0
    # the while loop do the mitosis process, while monitoring the conditions limits that where given
    # each iteration represent a cycles
    while current_division_cycles < max_cycles:
        temp_cells_list = []
        flag = True
        # each iteration in the loop represent one cell mitosis process
        cell_list_length = len(cells_list)
        # for i in range(len(cells_list)):
        i = 0
        while i < cell_list_length:
            # monitoring the limiting conditions
            temp_cells_list += cells_list[i].mitosis()
            number_of_cells += cells_list[i].division_rate - 1
            i += 1
        current_division_cycles += 1
        # this are all the cells that i creates
        cells_list = temp_cells_list
        if not flag:
            break
    # return the cell list
    return cells_list


"""The function read all the info from the fasta file and save each sequence object in a list and return it"""


def read_from_fas(path):
    seq_list = []
    for seq_record in SeqIO.parse(path, "fasta"):
        seq_list.append(seq_record)
    return seq_list


"""Write the information into a csv file"""


def write_to_csv(genomes):
    with open('./output.csv', 'w', newline='') as csv_file:
        # create the name of the fields
        fieldnames = ['Name of sequence', 'Error rate', 'Number of cycles', 'Number of mutant cells',
                      'Number of cancer cells', 'Number of proteins']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        # add titles to file
        writer.writerow({'Name of sequence': 'Name of sequence', 'Error rate': 'Error rate',
                         'Number of cycles': 'Number of cycles',
                         'Number of mutant cells': 'Number of mutant cells',
                         'Number of cancer cells': 'Number of cancer cells',
                         'Number of proteins': 'Number of proteins'})

        for genome in genomes:
            error_rate_list = [num for num in np.arange(0.05, 0.51, 0.05)]
            error_rate_list.insert(0, 0.01)
            with Pool() as pool:  # using multiproccesing
                info = pool.starmap(multi_from_error_rate, zip(error_rate_list, repeat(genome)))
                """extract the data from and given cycle into the csv file.
                 used because you cant send a csv writer object into a gunction while using multiproccesing"""
            for i in info:
                for j in i:
                    writer.writerow({'Name of sequence': j[0], 'Error rate': j[1],
                                     'Number of cycles': j[2], 'Number of mutant cells': j[3],
                                     'Number of cancer cells': j[4], 'Number of proteins': j[5]})


"""this function continue the write_to_csv function for multiprocessing purposes."""


def multi_from_error_rate(error_rate, genome):
    mitosis_cycle = 1
    info = []
    while mitosis_cycle <= 5:  # the loop keeps track of the wanted cycle and if it wasnt passed
        counter = 0
        mutant_cell = MutantCell([str(genome.seq)], 0, error_rate)
        while counter < 3:
            cells_list = [mutant_cell]
            cells_list = mitosis_cycles(cells_list, mitosis_cycle)
            protein_list = get_protin_list(cells_list)  # returns the protein from a given cell list
            protein_list = list(dict.fromkeys(protein_list))  # delete duplicated protein
            mut_number = get_number_of_cells_by_type(cells_list, 'MutantCell')  # the number of mutant cell in a list
            cancer_number = len(cells_list) - mut_number  # number of cell list
            # save the info into a list
            info.append([genome.name, error_rate, mitosis_cycle, mut_number, cancer_number, len(protein_list)])
            counter += 1
        mitosis_cycle += 1
    return info


"""This function return the longest reading frame from every genome of the genomes of every cell in the cell list"""


def get_protin_list(cells_list):
    data = []
    protein_list = []
    # look for all repertoires
    for cell in cells_list:
        data += cell.repertoire()
    for cell in data:
        # delete Non-coding RNA
        if cell[2] != 'Non-coding RNA':
            protein_list.append(cell[2])
    return protein_list


"""This function returns the number of specific type of cell in a list of cells"""


def get_number_of_cells_by_type(cells_list, cell_type):
    counter = 0
    for cell in cells_list:
        if cell.name == cell_type:
            counter += 1
    return counter


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # check if the input is valid
    assert len(sys.argv) == 2, 'need to be two argument'
    path = sys.argv[1]

    # load the fasta file
    fasta_file = read_from_fas(path)
    write_to_csv(fasta_file)
    panda_file = pd.read_csv('./output.csv')

    for seq in fasta_file:
        # Take only the rows with the same name
        info = panda_file.loc[panda_file['Name of sequence'] == str(seq.name)]
        # Take only the rows in cycles 5
        info = info.loc[info['Number of cycles'] == 2]

        # plot number of mutations
        plt.plot(info['Error rate'].tolist(), info['Number of mutant cells'].tolist())
        plt.title(str(seq.name))
        plt.ylabel('Number of mutant cells')
        plt.xlabel('Error rate')
        plt.savefig('./' + str(seq.name) + '_Error rate_Number of mutant cells output.png')
        plt.show()

        # plot number of cancer cells
        plt.plot(info['Error rate'].tolist(), info['Number of cancer cells'].tolist())
        plt.title(str(seq.name))
        plt.ylabel('Number of cancer cells')
        plt.xlabel('Error rate')
        plt.savefig('./' + str(seq.name) + '_Error rate_Number of cancer cells output.png')
        plt.show()

        # plot number of cancer cells
        plt.plot(info['Error rate'].tolist(), info['Number of proteins'].tolist())
        plt.title(str(seq.name))
        plt.ylabel('Number of proteins')
        plt.xlabel('Error rate')
        plt.savefig('./' + str(seq.name) + '_Error rate_Number of proteins output.png')
        plt.show()

        # plot number of protein as function of error rate and cycles
        # Take only the rows with the same name
        info = panda_file.loc[panda_file['Name of sequence'] == str(seq.name)]
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.plot(info['Error rate'].tolist(), info['Number of cycles'].tolist(),
                info['Number of proteins'].tolist(), label=str(seq.name))
        ax.set_xlabel('Error rate')
        ax.set_ylabel('Number of cycles')
        ax.set_zlabel('Number of proteins')
        plt.savefig('./' + str(seq.name) + '_3 dim output.png')
        plt.show()
