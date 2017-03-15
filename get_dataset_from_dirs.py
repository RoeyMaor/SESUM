from os import listdir
from os.path import isfile, join
mypath1 = "plots4/middle"
onlyfiles1 = [f for f in listdir(mypath1) if isfile(join(mypath1, f))]
mapped = [int(f.split("_")[0]) for f in onlyfiles1]
goods = sorted(mapped)[:100]

with open("seqs4.out",'r') as seqs:
    for i,s in enumerate(seqs.readlines()):
        if i in goods:
            with open("new sequences/almost_sequences.out",'a+') as almost_seqs:
                almost_seqs.write(s)

