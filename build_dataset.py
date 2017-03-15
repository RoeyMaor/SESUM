import random
import subprocess

from alphabetConverter import get_MFE_and_suboptimal_pairing_strings_for_sequence, MAXIMUM_FREE_ENERGY_GAP_TO_MFE

cmd = ["RNAplot", "-o", "svg"]

nucs = ['A', 'G', 'C', 'U']


def random_pick(some_list, probabilities):
    x = random.uniform(0, 1)
    cumulative_probability = 0.0
    for item, item_probability in zip(some_list, probabilities):
        cumulative_probability += item_probability
        if x < cumulative_probability: break
    return item


def get_ran_parts(ktzavot_len, begin_part_len, middle_part_len, end_part_len):
    begin_part = ""
    end_part = ""
    middle_part = ""
    katze_1 = ""
    katze_2 = ""
    nuc_dic = {"A":"U","U":"A","G":"C","C":"G"}
    for i in range(begin_part_len):
        begin_part += random.choice(nucs)
    for i in range(end_part_len):
        end_part += random.choice(nucs)
    for i in range(middle_part_len):
        middle_part += random.choice(nucs)
    for i in range(ktzavot_len):
        c = random.choice(nucs)
        katze_1 += c
    for c in reversed(katze_1):
        katze_2 += nuc_dic[c]
    return katze_1, begin_part, middle_part, end_part, katze_2


def build_rna(i):
    # important_part1 = "GCCCCUUUUUCCCCCG"
    # important_part2 = "CCCCCCUUUUUCCCC"
    important_part1 = "UUUUU"
    important_part2 = "CCCCC"
    nuc_list1 = []
    nuc_list2 = []
    for nuc in important_part1:
        nuc_list1.append(nuc)
    for nuc in important_part2:
        nuc_list2.append(nuc)
    index1 = random.choice(range(5))
    index2 = random.choice(range(5))
    nuc1 = random_pick(nucs, [0.05, 0.05, 0.05, 0.85])
    nuc2 = random_pick(nucs, [0.05, 0.05, 0.85, 0.05])
    nuc_list1[index1] = nuc1
    nuc_list2[index2] = nuc2
    important_part1 = "".join(nuc_list1)
    important_part2 = "".join(nuc_list2)
    begin_part_len = 7
    end_part_len = 7
    middle_part_len = 34
    begin_part_len += i / 1000
    end_part_len += i / 1000
    middle_part_len -= 2 * (i / 2000)
    the_parts = get_ran_parts(6,begin_part_len, middle_part_len, end_part_len)
    return the_parts[0] + the_parts[1] + important_part1 + the_parts[2] + important_part2 + the_parts[3] + the_parts[4]


for i in range(10000):
    seq = build_rna(i)
    with open("helper.in", "w+") as helper:
        with open("seqs4.out", "a+") as seqs:
            # seq = "AAAGGAAAAAGGGCCCCUUUCCCCCGGGAAAAAAAAAUUUUUUUUUCCCCCCCCUUUCCCCCCUUUUCCUUUU"
            seqs.write(seq + "\n")
            if seq.find("UUUUU") != -1 and seq.find("CCCCC") != -1:
                helper.write(">plots4/good/" + str(i) + "\n")
            else:
                helper.write(">plots4/middle/" + str(i) + "\n")

            helper.write(seq + "\n")
            print(i)
            helper.write(
                get_MFE_and_suboptimal_pairing_strings_for_sequence(seq, MAXIMUM_FREE_ENERGY_GAP_TO_MFE)[0] + "\n")

    f = open("helper.in")
    subprocess.call(cmd, stdin=f)
