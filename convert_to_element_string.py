from alphabetConverter import get_MFE_and_suboptimal_pairing_strings_for_sequence, MAXIMUM_FREE_ENERGY_GAP_TO_MFE
import forgi.graph.bulge_graph as cgb


def convert_dataset_to_elementstring():
    with open("DGCR8.fasta",'r') as fasta:
        lines = fasta.readlines()
        for seq in lines:
            if seq.startswith(">"):
                continue
            with open("output.elementstring", "a+") as output:

                get_element_string_from_sequence(output, seq)


def get_element_string_from_sequence(output, seq):
    for s in get_MFE_and_suboptimal_pairing_strings_for_sequence(seq, MAXIMUM_FREE_ENERGY_GAP_TO_MFE):
        g = cgb.BulgeGraph()
        g.from_dotbracket(s)
        output.write(g.to_element_string() + "\n")

if __name__ == "__main__":
    convert_dataset_to_elementstring()
