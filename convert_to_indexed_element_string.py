import forgi.graph.bulge_graph as cgb

from alphabetConverter import get_MFE_and_suboptimal_pairing_strings_for_sequence, MAXIMUM_FREE_ENERGY_GAP_TO_MFE


# l = fus.dotbracket_to_pairtable('((((((((((((..((((((((.((((....)).)).))))))))..))))))))))))')
# g = cgb.BulgeGraph()
# g.from_dotbracket('((((((((((((..((((((((.((((....)).)).))))))))..))))))))))))')
# print(g.to_element_string())
# print(l)
# with open("out.try","w+") as out:
#     print(get_element_string_from_sequence(out,'ACGUGCCACGAUUCAACGUGGCACAG'))
#

def get_element_with_id_list_from_dotbracket(dot_bracket_str):
    graph = cgb.BulgeGraph()
    graph.from_dotbracket(dot_bracket_str)
    element_str = graph.to_element_string(True)
    list_result = element_str.split('\n')
    elements = list_result[0]
    ids = list_result[1]
    return zip([str(elem) for elem in elements], [str(elem_id) for elem_id in ids])


def get_classified_lists_for_all_variants_of_structures_from_sequence(sequence):
    all_variants_classified = []
    for s in get_MFE_and_suboptimal_pairing_strings_for_sequence(sequence, MAXIMUM_FREE_ENERGY_GAP_TO_MFE):
        all_variants_classified.append(get_element_with_id_list_from_dotbracket(s))
    return all_variants_classified


if __name__ == "__main__":
    # get_element_with_id_list_from_dotbracket('((((((((((((..((((((((.((((....)).)).))))))))..))))))))))))')
    all_variants = get_classified_lists_for_all_variants_of_structures_from_sequence(
        'GGAGGAGCAGCAAAGTGAAACTCACCAGAACTGTGTCAGTTTCACCCTGCTGCTCCTCC')
    for x in all_variants:
        print x
