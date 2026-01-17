import sys
from math import floor

"""This class represent polymerase"""


class Polymerase:
    """This function initialize the polymerase"""

    def __init__(self, type, error_rate = 0):
        # check the validation of the type
        assert type.isalpha(), 'type of polymerase need to be DNA or RNA '
        # if type.isalpha():
        #     type = type.upper()
        # else:
        #     raise AssertionError('type of polymerase need to be DNA or RNA ')
        type = type.upper()
        assert type == 'DNA' or type == 'RNA', 'type of polymerase need to be DNA or RNA'
        # if type == 'DNA' or type == 'RNA':
        #     self.type = type
        # else:
        #     raise AssertionError('type of polymerase need to be DNA or RNA')
        self.type = type
        # check if the error rate is between 0 to 1
        assert 0 <= error_rate <= 1, 'error rate need to be between 0 to 1'
        # if 0 <= error_rate <= 1:
        #     self.error_rate = error_rate
        # else:
        #     raise AssertionError('error rate need to be between 0 to 1')
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
        if self.type == 'DNA':
            dicOfFairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        else:
            dicOfFairs = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'}
        # iterate on the sequence
        for letter in temp_seq:
            comp_dna += dicOfFairs[letter]
        return comp_dna


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
        # if all(matched_list):
        #     pass
        # else:
        #     raise AssertionError('RNA seq need to contain only A,U,G,C chars')
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
            listOfCodons = []
            for element in longest_seq:
                listOfCodons.append(longest_seq[element])
            return listOfCodons
        else:
            return None

    """This function help me to find frame of the translate"""

    def find_seq(self, rna_seq):
        """
        :param rna_seq: The sequence that start from the start codon we found in this iteration.
        :return: The protein.
        """
        # initialize a dictionary
        seq = {}
        # save the start codon
        seq = {0: rna_seq[0: 3]}
        # Save the length of the current sequence.
        n = len(rna_seq)
        # This loop jump on each set of three in the frame and check if we have a stop codon.
        # If its a stop codon its will stop, o.w it will continue untill the last three that logic.
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
        listOfCodons = self.translate(rna_seq)
        # check if there any amino acid
        if listOfCodons:
            protein = ''
            for codon in listOfCodons:
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

    def __init__(self, name, genome, num_copies, genetic_code, start_codons, division_rate):
        self.name = name
        # a list of dna sequences
        self.genome = genome
        # check the valid of copies
        assert num_copies > 0 and isinstance(num_copies, int), 'num fo copies need to be a positive integer'
        # if num_copies > 0 and isinstance(num_copies, int):
        #     self.num_copies = num_copies
        # else:
        #     raise AssertionError('num fo copies need to be a positive integer')
        self.num_copies = num_copies
        self.genetic_code = genetic_code
        self.start_codons = start_codons
        # check the valid of division rate
        assert division_rate > 1 and isinstance(division_rate, int), 'num fo division rate need to be an ' \
                                                                     'integer bigger than 1'
        # if division_rate > 1 and isinstance(division_rate, int):
        #     self.division_rate = division_rate
        # else:
        #     raise AssertionError('num fo division rate need to be an integer bigger than 1')
        self.division_rate = division_rate
        # initialize the polymerases
        self.polDna = Polymerase('DNA')
        self.polRna = Polymerase('RNA')
        # initialize the ribosome
        self.ribo = Ribosome(self.genetic_code, self.start_codons)

    """This function duplicate the cell by mitosis"""

    def mitosis(self):
        # create an array of cells
        cells = []
        cell = Cell(self.name, self.genome, self.num_copies, self.genetic_code, self.start_codons, self.division_rate)
        # duplicate the cell respectively to the division rate
        for i in range(self.division_rate):
            cells.append(cell)
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
                               , self.genetic_code, self.start_codons, self.division_rate)
            # create complementary seq
            complementarySeq = []
            for i in self.genome:
                # transcribe the seq and append to the list
                complementarySeq.append(self.polDna.transcribe(i))
            # this is the cell with the complementary seq
            cell_complementary = Cell(self.name, complementarySeq, int(self.num_copies / 2)
                                      , self.genetic_code, self.start_codons, self.division_rate)
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
        for i in self.genome:
            # load the ssr in the sequence
            srr_seq = self.find_srr(i)
            # check if there are a srr.
            if srr_seq is None:
                srr_seq = 'No simple repeats in DNA sequence'
            # take the dna and transcribe it to a rna
            rna_seq = self.polRna.transcribe(i)
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
    def __init__(self, genome):
        self.start_codons = ['AUG', 'GUG', 'UUG']
        self.division_rate = 4
        self.num_copies = num_copies = 1
        self.name = 'ProkaryoticCell'
        self.genetic_code = {
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
        super().__init__(self.name, genome, self.num_copies, self.genetic_code, self.start_codons, self.division_rate)


"""This class represent eukaryotic cell"""


class EukaryoticCell(Cell):
    """initialize the eukaryotic"""

    def __init__(self, name, genome, division_rate):
        self.num_copies = 2
        self.start_codons = ['AUG']
        self.genetic_code = {
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
        super().__init__(name, genome, self.num_copies, self.genetic_code, self.start_codons, self.division_rate)


"""This class represent neuron"""


class NeuronCell(EukaryoticCell):
    """initialize the neuron"""

    def __init__(self, genome):
        self.division_rate = 2
        super().__init__('NeuronCell', genome, self.division_rate)


"""This class represent stem cell"""


class StemCell(EukaryoticCell):
    """initialize the neuron"""

    def __init__(self, genome):
        self.division_rate = 3
        super().__init__('StemCell', genome, self.division_rate)


"""This function use as factory"""


def factory(name, genomes):
    if name == 'StemCell':
        cell = StemCell(genomes)
    elif name == 'NeuronCell':
        cell = NeuronCell(genomes)
    elif name == 'ProkaryoticCell':
        cell = ProkaryoticCell(genomes)
    else:
        raise AssertionError('no much between name and type of cell')
    return cell


"""This function help me to find the final number of cells"""


def finalNumOfCell(division, cell, maxCells):
    i = 0
    # find the max of division i can
    while i < division and cell.division_rate ** i <= maxCells:
        i += 1
    # check if we can use all the division
    if i == division and cell.division_rate ** i <= maxCells:
        return int(cell.division_rate ** i)
    # find the rest of cells
    rest = maxCells - (cell.division_rate ** i)
    # number of cell that can divide
    div = floor(rest / (cell.division_rate - 1))
    # number of cells that cant divide
    same = cell.division_rate ** i - div
    return int(div * cell.division_rate + same)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    n = len(sys.argv)
    # name of cell
    cellName = sys.argv[1]
    # check the num of division
    assert sys.argv[2].isdigit() and int(sys.argv[2]) > 0, 'num fo division need to be a positive integer'
    # if sys.argv[2].isdigit() and int(sys.argv[2]) > 0:
    #     # num of division
    #     numOfDivision = int(sys.argv[2])
    # else:
    #     raise AssertionError('num fo division need to be a positive integer')
    # num of division
    numOfDivision = int(sys.argv[2])
    assert sys.argv[3].isdigit() and int(sys.argv[3]) > 1, 'num of division need to be a positive integer'
    # if sys.argv[3].isdigit() and int(sys.argv[3]) > 1:
    #     # max num of cells
    #     maxCells = int(sys.argv[3])
    # else:
    #     raise AssertionError('num of division need to be a positive integer')
    # max num of cells
    maxCells = int(sys.argv[3])
    # all the genomes sequence
    genomes = []
    for i in range(4, n):
        char_list = ['A', 'T', 'G', 'C']
        matched_list = [characters in char_list for characters in sys.argv[i]]
        #check the valid of seq
        assert all(matched_list), 'DNA seq need to contain only A,T,G,C chars'
        # if all(matched_list):
        #     genomes.append(sys.argv[i])
        # else:
        #     raise AssertionError('DNA seq need to contain only A,T,G,C chars')
        genomes.append(sys.argv[i])
    # create the cell
    cell = factory(cellName, genomes)
    # the final number of cells
    final_Number_Of_Cell = finalNumOfCell(numOfDivision, cell, maxCells)
    # print all what we ask for
    print('Original cell: ' + str(cell))
    print('Final number of cells: ' + str(final_Number_Of_Cell))
    print('Repertoire: ' + str(cell.repertoire()))
    meiosis = cell.meiosis()
    if meiosis:
        print('Undergoing meiosis...')
        cell1 = meiosis[0]
        cell2 = meiosis[1]
        print('First cell genome: ' + str(cell1.genome))
        print('First cell repertoire: ' + str(cell1.repertoire()))
        print('Second cell genome: ' + str(cell2.genome))
        print('Second cell repertoire: ' + str(cell2.repertoire()))
    else:
        print('Cannot undergo meiosis')

