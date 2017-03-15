inputName = "new sequences/good_sequences.out"
outputName = "new sequences/good_sequences.fasta"
if __name__ == "__main__":
    f = open(inputName,"r")
    o = open(outputName,"w")
    index = 0
    for line in f:
        o.write("> "+str(index)+"\n")
        o.write(line)
        index+=1
    f.close()
    o.close()
