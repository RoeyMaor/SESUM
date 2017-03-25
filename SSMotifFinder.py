import alphabetConverter
import sys
import multiprocessing
import time

from MemeHelper import *
from convert_to_indexed_element_string import *

numberOfCores = 1 #default value

def check_input_is_ok():
    if len(sys.argv) != 3:
        print("usage: python SSMotifFinder <RNA sequences fasta file> <number of cpu cores to utilize>")
        sys.exit(-1)

    path_to_file = sys.argv[1]
    validateFilePathIsLegal(path_to_file) #exists if not a legal path

    try:
        nubmerOfCores = int(sys.argv[2])
    except ValueError:
        print("number of cores must be an int!")
        sys.exit(-1)

    if multiprocessing.cpu_count() < numberOfCores:
        print("Not enough CPU cores")
        exit(-1)

    if numberOfCores<=0:
        print("number of cores must be greater than 0")
        exit(-1)


if __name__ == "__main__":

    check_input_is_ok()
    path_to_file = sys.argv[1]
    numberOfCores = int(sys.argv[2])

    start = time.clock()
    sequencesList = fasta_File_To_List_Of_Sequences(path_to_file)
    end = time.clock()
    print("Time it took to parse fasta file: "+str(end-start)+" seconds")

    start = time.clock()
    new_alphabet_input = alphabetConverter.get_list_of_sequences_with_new_alphabet_from_file(sequencesList)
    end = time.clock()
    print("Time it took to create dataset with new alphabet: " + str(end - start) + " seconds")

    new_alphabet_sequences = [sequence for sequences in new_alphabet_input for sequence in sequences] #flatten list
    new_alphabet_file_path = path_to_file+"-NEWALPHABET"
    create_fasta_file_from_list_of_sequences(new_alphabet_sequences,new_alphabet_file_path)
    new_alphabet_file = open(new_alphabet_file_path,"r")

    start = time.time()
    meme_result = get_Meme_Result_From_File(new_alphabet_file,new_alphabet_file_path,numberOfCores)
    end = time.time()
    print("Time it took to get meme results: " + str(end - start) + " seconds")

    start = time.clock()
    sequences_substructures = [get_classified_lists_for_all_variants_of_structures_from_sequence(sequence) for sequence in sequencesList]
    flattened_sequences_substructures = [structure for structures in sequences_substructures for structure in structures]
    end = time.clock()
    print("Time it took to create detailed structure information: " + str(end - start) + " seconds")
    assert(len(flattened_sequences_substructures)==len(new_alphabet_sequences))
    new_alphabet_file.close()
    start = time.clock()
    fill_maps(meme_result,flattened_sequences_substructures)
    end = time.clock()
    print("Time it took to extract gapped motif information: "+str(end-start) + " seconds")
    print("\n")
    print("Found "+str(len(meme_result.motifs))+" ungapped sequence-strctural motifs.")
    for (m1,m2) in map_from_motif_index_tupples_to_number_of_joint_appearences:
        if m1==m2:
            continue
        value = map_from_motif_index_tupples_to_number_of_joint_appearences[(m1,m2)]
        if value > THRESHOLD_TO_CONSIDER_TWO_MOTIFS_TO_BE_ON_THE_SAME_SUBSTRCTURE*len(new_alphabet_sequences):
            print("Motifs "+str(m1)+" and "+str(m2)+" are part of the same gapped sequence-strctural motif - on the same loop")

