import os
import sys
import subprocess

def validateFilePathIsLegal(path):
    if not os.path.isfile(path):
        print("couldn't find file: "+path)
        sys.exit(-1)

def fetchFileFromRelativePath(path):
    validateFilePathIsLegal(path)
    return open(path)


def execute_command_and_get_output(command,args,stdin_for_command):
    command = command + " " + args
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE, shell=True)
    out, err = p.communicate(stdin_for_command)
    return out

def getLineCountOfFile(file):
    i=-1
    for i, l in enumerate(file):
        pass
    return i + 1

def fasta_File_To_List_Of_Sequences(path):
    file = fetchFileFromRelativePath(path)
    content = [x.strip() for x in file.readlines()]
    file.close()
    # now we need to get rid of fasta comments, and join lines that belong to the same sequence:
    joinedSequencesContent = []
    currentJoinedSequence = ""
    for row in content:
        if row.startswith(">"): # fasta comment
            joinedSequencesContent.append(currentJoinedSequence.upper())
            currentJoinedSequence = ""
        else:
            currentJoinedSequence += row
    if not currentJoinedSequence == "":
        joinedSequencesContent.append(currentJoinedSequence.upper())

    return joinedSequencesContent[1:]

def create_fasta_file_from_list_of_sequences(sequences_list,new_file_name):
    file = open(new_file_name, 'w')
    index = 0
    for sequence in sequences_list:
        file.write(">"+str(index)+"\n"+sequence+"\n")
        index+=1
    file.close()
