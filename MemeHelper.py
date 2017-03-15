from SSMotifFinderUtils import *
import xml.etree.ElementTree as ET
import itertools

NAME_OF_ALPHABET_FILE="SSMotifAlphabet"
MAX_E_VALUE=0.001 #maximum e-value to be considered as a motif
MAX_NUMBER_OF_MOTIFS=3
MAXSIZE = 300000
MEME_TEXT_OUTPUT=False
THRESHOLD_TO_CONSIDER_TWO_MOTIFS_TO_BE_ON_THE_SAME_SUBSTRCTURE = 0.8

map_from_motif_index_to_substructures = {} # motif_index -> [(sequenceIndex,(loopKind,loopIndex),start_index)]
map_from_substructure_to_motifs = {} # (sequenceIndex,(loopKind,loopIndex)) -> [(motif_index,start_index)]
map_from_motif_index_tupples_to_number_of_joint_appearences = {}# (motif_index,motif_index) -> count
motifs_that_are_with_themselves_in_substrctures = []

class UnGappedMotifLocation:
    def __init__(self):
        self.sequenceIndex = 0 #the sequence of the motif in the input dataset
        self.startIndex = 0

class UnGappedMotif:
    def __init__(self):
        self.eValue = MAX_E_VALUE
        self.width = 0
        self.locations = []
        return

class MemeResult:
    def __init__(self):
        self.motifs = [] #list of UnGappedMotifs
        return

def get_Meme_Result_From_File(fastaFile, fastaFile_path, numberOfCores):
    '''

    :param fastaFile: the actual fasta file object of the file containing the NEW alphabet sequences
    :return: returns a MemeResult object
    '''
    #-nostatus -p processors -evt max_e_value
    #meme -p 4 -alph SSMotifAlphabet -nostatus -evt 0.01 "CARPIN1SET.fasta-NEWALPHABET-CROPPED"
    numberOfLines = getLineCountOfFile(fastaFile)
    command = "/home/roym/meme/bin/meme"
    commandArgs = "-p "+str(numberOfCores)+" -alph "+NAME_OF_ALPHABET_FILE+" -nostatus"+" -evt "+str(MAX_E_VALUE)
    commandArgs += " -minsites "+str((numberOfLines/10))+" -nmotifs "+str(MAX_NUMBER_OF_MOTIFS) + " -wg 1"
    commandArgs += " -maxsize "+str(MAXSIZE)+(" -text" if MEME_TEXT_OUTPUT else "")+" "+fastaFile_path
    res = execute_command_and_get_output(command, commandArgs, "")

    tree = ET.parse('meme_out/meme.xml')
    root = tree.getroot()
    allMotifs = [motif for motif in root.findall('motifs')[0].findall('motif')]
    res = MemeResult()
    for motif in allMotifs:
        ungappedMotif = UnGappedMotif()
        motif_attributes = motif.attrib
        ungappedMotif.eValue = motif_attributes['e_value']
        ungappedMotif.width = int(motif_attributes['width'])
        for site in motif.findall('contributing_sites')[0].findall('contributing_site'):
            ungapped_motif_location = UnGappedMotifLocation()
            site_attributes = site.attrib
            ungapped_motif_location.sequenceIndex = int(site_attributes['sequence_id'].split('_')[1])
            ungapped_motif_location.startIndex = int(site_attributes['position'])
            ungappedMotif.locations.append(ungapped_motif_location)
        res.motifs.append(ungappedMotif)
    return res

def print_meme_result(meme_result):
    print("motifs:")
    i = 0
    for motif in meme_result.motifs:
        print("-motif "+str(i)+":")
        print(" -e-value: "+str(motif.eValue))
        print(" -width: "+str(motif.width))
        print(" -locations:")
        for location in motif.locations:
            print("     -sequence index: "+str(location.sequenceIndex))
            print("     -start index: " + str(location.startIndex))
        i+=1

def fill_maps(meme_result, flattened_sequences_substructures):
    global map_from_substructure_to_motifs, map_from_motif_index_to_substructures
    motif_index = 0

    for motif in meme_result.motifs:

        width = motif.width

        for location in motif.locations:

            sequence_index = location.sequenceIndex
            relevant_subsequence_structure = flattened_sequences_substructures[sequence_index]
            start_index = location.startIndex
            structure_tupples_already_discovered_in_location = []

            for i in range(start_index,start_index+width):

                subsequence_tuple = relevant_subsequence_structure[i]
                if subsequence_tuple in structure_tupples_already_discovered_in_location:
                    continue
                structure_tupples_already_discovered_in_location.append(subsequence_tuple)

                subsequence_type = subsequence_tuple[0]
                if subsequence_type=='i' or subsequence_type=='m':
                    if not ((sequence_index,subsequence_tuple) in map_from_substructure_to_motifs):
                        map_from_substructure_to_motifs[(sequence_index,subsequence_tuple)] = []
                    map_from_substructure_to_motifs[(sequence_index,subsequence_tuple)].append((motif_index,start_index))
                    if not motif_index in map_from_motif_index_to_substructures:
                        map_from_motif_index_to_substructures[motif_index] = []
                    map_from_motif_index_to_substructures[motif_index].append((sequence_index,subsequence_tuple,start_index))
        motif_index+=1
    # we are interested only in substrctures with more than one motif in it:
    map_from_substructure_to_motifs = {k: v for k, v in map_from_substructure_to_motifs.items() if len(v)>1}
    for i in range(len(meme_result.motifs)):
        for j in range(len(meme_result.motifs)):
            if j>i:
                continue
            map_from_motif_index_tupples_to_number_of_joint_appearences[(i,j)] = 0
    for l in map_from_substructure_to_motifs.items():
        l = l[1]
        motif_indexes_in_list = []
        for tupple in l:
            m_index = tupple[0]
            if not m_index in motif_indexes_in_list:
                motif_indexes_in_list.append(m_index)
        cartesian_product = list(itertools.product(motif_indexes_in_list,motif_indexes_in_list))
        for tupple in cartesian_product:
            if tupple[0]<tupple[1]:
                continue
            map_from_motif_index_tupples_to_number_of_joint_appearences[tupple]+=1
    print(map_from_motif_index_to_substructures[2])
    print(map_from_motif_index_to_substructures)
    print(map_from_substructure_to_motifs)
    print(map_from_motif_index_tupples_to_number_of_joint_appearences)

def avg_distance_between_motifs_on_the_same_substrcture(ungappedMotif1,ungappedMotif2):
    '''

    :param ungappedMotif1: UnGappedMotif
    :param ungappedMotif2: UnGappedMotif
    :return: 0 if the motifs are not on the same substrcture, otherwise returns the mean distance in nucleutides between them
    '''

    if ungappedMotif1 == ungappedMotif2:
        pass
    else:
        pass
    avgDistance = 0
    locations1 = ungappedMotif1.locations
    locations2 = ungappedMotif2.locations



if __name__=="__main__":
    pass
