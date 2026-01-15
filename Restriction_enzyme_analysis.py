import sys


"""This function look for srr in the sequence."""
def find_srr(dna_seq):
    """
    :param dna_seq: the dna i get from the user
    :return: All the srr that exist in the sequence.
    """
    n = len(dna_seq)
    #Initialize an empty string
    final_seq = ""
    #This loop limit the length of the repeated sequence.
    for i in range(1, 7):
        #This loop look for the srr
        for m in range(n - (3 * i) + 1):
            #this is the seq i want to check if he repeated in the strain.
            seq = dna_seq[m:m + i]
            #This var help me to count how much the seq repeated it self.
            flag = 0
            # I initialize the seq to be from the current index until the end of the strain.
            temp_seq = dna_seq[m:n + 1]
            # Count how much time it is repeated.
            while temp_seq.startswith(seq) and flag * i < n:
                flag = flag + 1
                temp_seq = dna_seq[m + i * flag: n]
            else:
                #if the seq repeat it self more than 3 times than i save him.
                if flag >= 3:
                    final_seq += seq + "," + str(flag) + ";"
    # if i have no srr i send None.
    if final_seq is None:
        return None
    return final_seq[0: len(final_seq) - 1]

"""This function transcribe the sequence to rna in the sequence."""
def transcribe(dna_seq):
    """
    :param dna_seq: The dna that the user want to transcribe.
    :return: An rna from the dna.
    """
    # Swap the string.
    temp_seq = dna_seq[::-1]
    # This string will contain the rna
    comp_dna = ''
    # This loop run on the nucleotide.
    for i in range(len(dna_seq)):
        # Check witch nucleotide i have.
        if temp_seq[i] == 'a' or temp_seq[i] == 'A':
            comp_dna += 'U'
        if temp_seq[i] == 't' or temp_seq[i] == 'T':
            comp_dna += 'A'
        if temp_seq[i] == 'c' or temp_seq[i] == 'C':
            comp_dna += 'G'
        if temp_seq[i] == 'g' or temp_seq[i] == 'G':
            comp_dna += 'C'
    return comp_dna

"""This function translate the rna"""
def translate(rna_seq):
    """
    :param rna_seq: The rna we got from the transcribe function.
    :return: Return the longest protein, if exist.
    """
    # Make all the char to be upper char.
    rna_seq.upper()
    # Save the length of the rna
    n = len(rna_seq)
    # Initialize an empty strings.
    current_seq = ''
    longest_seq = ''
    # Save the start codon.
    start_codon = 'AUG'
    # This loop run on the index of the rna.
    for i in range(n):
        # Search for a start codon.
        if rna_seq[i:n].startswith(start_codon):
            # Save the protein that we found.
            current_seq = find_seq(rna_seq[i:n])
        # Save the longest protein between the previous and the one that we found in this iteration.
        if len(current_seq) > len(longest_seq):
            longest_seq = current_seq
    # Check if we found a protein.
    # If we didnt we send None, o.s, we send the protein.
    if longest_seq is None:
        return None
    else:
        return longest_seq[0: len(longest_seq) - 1]

"""This function help me to find frame of the translate"""
def find_seq(rna_seq):
    """
    :param rna_seq: The sequence that start from the start codon we found in this iteration.
    :return: The protein.
    """
    # Save the start codon.
    seq = 'AUG;'
    # Save the length of the current sequence.
    n = len(rna_seq)
    # This loop jump on each set of three in the frame and check if we have a stop codon.
    # If its a stop codon its will stop, o.w it will continue untill the last three that logic.
    for i in range(1, (n // 3)):
        # Look if we have a stop codon. If we have we will send the protein that we found.
        if rna_seq[(i * 3):n].startswith('UAA') or rna_seq[(i * 3):n].startswith('UGA') \
                or rna_seq[(i * 3):n].startswith('UAG'):
            # You said we dont need to return the stop codon
            return seq
        # O.w i will add the three nucleotides to the string.
        else:
            seq += rna_seq[i * 3: (i * 3) + 3] + ';'
    return seq


"""This is the main function"""
if __name__ == '__main__':
    # I check i have valid input
    if len(sys.argv) != 2:
        print('invalid input')
    else:
        # sequence is the string
        sequence = sys.argv[1]
        # load the ssr in the sequence
        srr_seq = find_srr(sequence)
        # Check if there is a srr. If exist i will print it, o.w , print their isn't.
        if srr_seq is None:
            print('No simple repeats in DNA sequence')
        else:
            print(srr_seq)
        # Take the dna and transcribe it to a rna
        rna_seq = transcribe(sequence)
        print('RNA sequence:' + rna_seq)
        # Look for the longest protein that exist. If their isn't protein i will inform the user.
        if translate(rna_seq) is None:
            print('Non-coding RNA')
        else:
            print('Translation:' + translate(rna_seq))
