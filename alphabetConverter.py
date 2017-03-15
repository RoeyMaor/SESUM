import sys
import multiprocessing
from SSMotifFinderUtils import *

MAXIMUM_FREE_ENERGY_GAP_TO_MFE = 0 #in kcal/mol
NUMBER_OF_SEQUENCES_THRESHOLD_TO_USE_PARALLEL = 20

alphabet_pairing_dictionary = {
            ('A', '.'): 'A',
            ('C', '.'): 'C',
            ('G', '.'): 'G',
            ('U', '.'): 'U',
            ('A', '('): 'B',
            ('C', '('): 'D',
            ('G', '('): 'H',
            ('U', '('): 'V',
            ('A', ')'): 'B',
            ('C', ')'): 'D',
            ('G', ')'): 'H',
            ('U', ')'): 'V',
        }



'''
def get_MFE_pairing_string_for_sequence(sequence):
    rnafold_output = execute_command_and_get_output("RNAfold","",sequence)
    rows = rnafold_output.split("\n")
    output_pairing_string = rows[1].split(" ")[0]
    return [output_pairing_string]
'''

def get_MFE_and_suboptimal_pairing_strings_for_sequence(sequence,max_energy_distance_from_MFE):
    '''

    :param sequence:
    :param max_energy_distance_from_MFE: in units of kcal/mol
    :return: a list of the suboptimal pairing strings of te sequence in the given enery range
    '''
    #
    rnasubopt_output = execute_command_and_get_output("RNAsubopt","-e "+str(max_energy_distance_from_MFE)+" ",sequence)
    rows = rnasubopt_output.split("\n")[1:-1] #without the sequence itself, the MFE, and trailing ''
    if rows: #suboptimal strctures exist
        return [st.split(" ")[0] for st in rows]
    return rows


def convert_sequence_and_pairing_string_to_new_alphabet(sequence,pairing_string):
    if not len(sequence) == len(pairing_string):
        print("fatal error, sequence and pairing strings are not in the same length")
        sys.exit(-1)
    res = ""
    for i in range(len(sequence)):
        sequence_char = sequence[i]
        pairing_char = pairing_string[i]
        res += alphabet_pairing_dictionary[(sequence_char,pairing_char)]
    return res


def get_list_of_sequences_with_new_alphabet_from_file(sequencesList):
    res = []
    for sequence in sequencesList:
        sequence = sequence.replace("T","U")
        res.append([convert_sequence_and_pairing_string_to_new_alphabet(sequence,structure) for structure in
                    get_MFE_and_suboptimal_pairing_strings_for_sequence(sequence,MAXIMUM_FREE_ENERGY_GAP_TO_MFE)])
    return res

'''
def parallel_work_portion(sequencesList,startIndex,endIndex,res):
    for i in range(startIndex,endIndex+1):
        sequence = sequencesList[i]
        sequence = sequence.replace("T","U")
        res[i] = [convert_sequence_and_pairing_string_to_new_alphabet(sequence,structure) for structure in
                    get_MFE_and_suboptimal_pairing_strings_for_sequence(sequence,MAXIMUM_FREE_ENERGY_GAP_TO_MFE)]
'''

if __name__ == "__main__":
    '''
    Testing alphabet converter:
    '''
    fileName = "DGCR8.fasta"
    tonyFile = open("tony.fasta","w")
    for element in get_list_of_sequences_with_new_alphabet_from_file(fileName):
        for sequence in element:
            tonyFile.write(">\n")
            tonyFile.write(sequence+"\n")
    tonyFile.close()