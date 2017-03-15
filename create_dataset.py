import random


def get_motifed_rna():
    res = ""
    nucs = ['u', 'g', 'c', 'a']
    res += random.choice(nucs)+random.choice(nucs)
    res += "acu__a_c__acgu"
    for i in range(3):
        res += random.choice(nucs)
    res += "gcguuauuggggu"
    res += random.choice(nucs)
    return res

if __name__ == "__main__":
    print(get_motifed_rna())