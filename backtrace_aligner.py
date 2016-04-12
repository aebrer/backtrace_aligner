#!/usr/bin/env python
import sys
import Bio.SubsMat.MatrixInfo as sim_matrix
from Bio import SeqIO as seq
import itertools
import argparse
import os.path
import re

# a modified dictionary class that enables rapid initialization of tree structures
# great for matricies as well
class NestedDict(dict):
    def __missing__(self, key):
        self[key] = NestedDict()
        return self[key]

#checks if there is a valid file at a specified location
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

# importing arguments from the user
parser=argparse.ArgumentParser(
    description=''' Performes pairwise sequence alignments using an affine gap penalty. ''',
    epilog=""""""
    )
parser.add_argument('-i',
                    dest="input", 
                    required=True,
                    metavar="FILE",
                    type=lambda x: is_valid_file(parser, x), 
                    help='This should be a fasta formatted file ' +\
                          'with multiple protein sequences.'
                    )
parser.add_argument('-g',
                    dest="gap_open", 
                    default = 20,
                    type = int,
                    help='Gap open penalty.')
parser.add_argument('-e',
                    dest="gap_extend", 
                    default = 1,
                    type = int,
                    help='Gap extend penalty.')
args=parser.parse_args()



# get the score from the similarity matrix
def get_score(a,b,sim):
    if (a,b) in sim:
        return sim[(a,b)]
    else:
        try:
            return sim[(b,a)]
        except:
            return None
            
    
# sets boundry conditions for local alignments
def init_local_matrix(x,y,align_dict):

    i = 0
    j = 0

    while i < len(x):
        align_dict[i][0]["score"] = 0
        align_dict[i][0]["trace"] = "up"
        i += 1
    while j < len(y):
        align_dict[0][j]["score"] = 0
        align_dict[0][j]["trace"] = "left"
        j += 1

    align_dict[0][0]["trace"] = "stop"

# sets boundry conditions for global alignments
def init_global_matrix(x,y,align_dict):

    i = 0
    j = 0

    while i < len(x):
        align_dict[i][0]["score"] = -1 * i
        align_dict[i][0]["trace"] = "up"
        i += 1
    while j < len(y):
        align_dict[0][j]["score"] = -1 * j
        align_dict[0][j]["trace"] = "left"
        j += 1

    align_dict[0][0]["trace"] = "stop"

# checks if the last position was a gap or not
def in_gap_check(prev_trace):

    if prev_trace == "up" or prev_trace == "left":
        return gap_extend
    else:
        return gap_open


# writes tsv files for matricies
def table_writer(x,y,align_dict,matrix_name, outname):
    
    output = open(outname, "w")
    
    # printing y
    output.write("\t")   
    for char in y:
        output.write(str(char) + "\t")
    output.write("\n")    

    # print scores and x
    for i in range(0,len(x)):
        char = x[i]
        output.write(str(char) + "\t")
        
        for j in range(0,len(y)):
            score = str(align_dict[i][j][matrix_name])
            output.write(score + "\t")
        output.write("\n")

    output.write("\n")
    output.write("\n")
    output.write("\n")
    
    output.close()


# do a global alignment with affine gap penalty
def global_align(x,y,all_scores,align_dict):    
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            
            current_score = 0
            current_trace = "stop"
             
            diag_score = align_dict[i-1][j-1]["score"]
            up_score = align_dict[i-1][j]["score"]
            left_score = align_dict[i][j-1]["score"]
            
            diag_trace = align_dict[i-1][j-1]["trace"]
            up_trace = align_dict[i-1][j]["trace"]
            left_trace = align_dict[i][j-1]["trace"]
            
            # nogap
            if (diag_score + get_score(x[i] , y[j], sim) > \
                 left_score + in_gap_check(left_trace))\
            and (diag_score + get_score(x[i] , y[j], sim) > \
                up_score + in_gap_check(up_trace)):

                current_score = diag_score + get_score(x[i] , y[j], sim)
                all_scores.append(str(current_score))
                current_trace = "diag"
                    
            
            #gaps
            elif up_score + in_gap_check(up_trace) > \
                 left_score + in_gap_check(left_trace):

                current_score = up_score + in_gap_check(up_trace)
                all_scores.append(str(current_score))
                current_trace = "up"
                    
            else:

                current_score = left_score + in_gap_check(left_trace)
                all_scores.append(str(current_score))
                current_trace = "left"

            align_dict[i][j]["score"] = current_score
            align_dict[i][j]["trace"] = current_trace
    return align_dict, all_scores

# do local alignment
def local_align(x,y,all_scores,align_dict):    
    for i in range(1,len(x)):
        for j in range(1,len(y)):
            
            current_score = 0
            current_trace = "stop"
             
            diag_score = align_dict[i-1][j-1]["score"]
            up_score = align_dict[i-1][j]["score"]
            left_score = align_dict[i][j-1]["score"]
            
            diag_trace = align_dict[i-1][j-1]["trace"]
            up_trace = align_dict[i-1][j]["trace"]
            left_trace = align_dict[i][j-1]["trace"]
            
            # match or mismatch
            if (diag_score + get_score(x[i] , y[j], sim) > \
                 left_score + in_gap_check(left_trace))\
            and (diag_score + get_score(x[i] , y[j], sim) > \
                up_score + in_gap_check(up_trace))\
            and (diag_score + get_score(x[i] , y[j], sim) > 0):

                current_score = diag_score + get_score(x[i] , y[j], sim)
                all_scores.append(str(current_score))
                current_trace = "diag"
                    
            
            #gaps
            elif (up_score + in_gap_check(up_trace) > \
                 left_score + in_gap_check(left_trace))\
            and (up_score + in_gap_check(up_trace) > 0):

                current_score = up_score + in_gap_check(up_trace)
                all_scores.append(str(current_score))
                current_trace = "up"
                    
            elif (left_score + in_gap_check(left_trace) > 0):
                current_score = left_score + in_gap_check(left_trace)
                all_scores.append(str(current_score))
                current_trace = "left"
            

            align_dict[i][j]["score"] = current_score
            align_dict[i][j]["trace"] = current_trace
            
    return align_dict, all_scores

    
# this will get the actual sequence of the alignment, and write it to the
# specified file
def backtracer(x,y,align_dict,x_pos,y_pos,x_name,y_name,outname):
    
    trace = align_dict[x_pos][y_pos]["trace"]
    x_align = []
    y_align = []

    while trace != "stop":
        if trace == "diag":
            x_align.append(x[x_pos])
            y_align.append(y[y_pos])
            x_pos -= 1
            y_pos -= 1
        elif trace == "up":
            x_align.append(x[x_pos])
            y_align.append("-")
            x_pos -= 1
        elif trace == "left":
            x_align.append("-")
            y_align.append(y[y_pos])
            y_pos -= 1
        
        trace = align_dict[x_pos][y_pos]["trace"]   

    x_align.reverse()
    y_align.reverse()
    
    output = open(outname, "w")
    output.write(">" + x_name + "\n")
    counter = 1
    for xchar in x_align:        
        output.write(xchar)
        if counter > 80:
            output.write("\n")
            counter = 0
        counter += 1
    output.write("\n")
    output.write(">" + y_name + "\n")
    
    counter = 1
    for ychar in y_align:        
        output.write(ychar)
        if counter > 80:
            output.write("\n")
            counter = 0
        counter += 1
    output.write("\n")
    output.close

















##################################
# Init code
#################################


# scores need to be negative here
# high cost of opening but low cost of extending is typical in protein alignments
gap_open = -1 * args.gap_open
gap_extend = -1 * args.gap_extend
sim = sim_matrix.blosum62


#generates a list with the sequences and a list with the headers, with matching idicies
fasta_head = list()
fasta_seq = list()
for element in seq.parse(args.input, 'fasta'):
    fasta_head.append(element.id)
    fasta_seq.append(element.seq)

#generates a list of all pairwise comparisions possible based on length of multiple alignment
pairwise = itertools.combinations(range(0,len(fasta_head)), 2)

percent_id = NestedDict()
best_global_scores = NestedDict()
best_local_scores = NestedDict()


for pair in pairwise:
    x = fasta_seq[pair[0]].upper()
    y = fasta_seq[pair[1]].upper()
    
    
    x_name = fasta_head[pair[0]]
    y_name = fasta_head[pair[1]]
    
    outname = x_name + "_" + y_name + "_"
    
    # formatting required for the matrix, the sequences start at index 1 after this
    if x[0] != "-":
        x = "-" + x
    if y[0] != "-":
        y = "-" + y



    #################################
    # global alignment
    #################################



    global_align_dict = NestedDict()
    all_global_scores = []
    init_global_matrix(x,y,global_align_dict)
    global_align_dict, all_global_scores = global_align(x,
                                                        y,
                                                        all_global_scores, 
                                                        global_align_dict)
    best_global_score = max(all_global_scores)


    # write score matrix
    table_writer(x,
                 y,
                 global_align_dict,
                 "score", 
                 outname + "global_score.tsv")
    # write trace matrix
    table_writer(x,
                 y,
                 global_align_dict,
                 "trace",
                 outname + "global_trace.tsv")

    # get the alignment
    backtracer(x,
               y,
               global_align_dict,
               len(x)-1,
               len(y)-1,
               x_name,
               y_name, 
               outname + "global_align.fasta")




    ############################################################
    # local alignment
    ############################################################


    local_align_dict = NestedDict()
    all_local_scores = []
    init_local_matrix(x,y,local_align_dict)
    local_align_dict, all_local_scores = local_align(x,
                                                     y,
                                                     all_local_scores, 
                                                     local_align_dict)


    # write score matrix
    table_writer(x,
                 y,
                 local_align_dict,
                 "score", 
                 outname + "local_scores.tsv")
    # write trace matrix
    table_writer(x,
                 y,
                 local_align_dict,
                 "trace", 
                 outname + "local_trace.tsv")

    # get the best score

    best_local_score = 0
    best_x = 0
    best_y = 0

    for i in local_align_dict:
        for j in local_align_dict[i]:
            if local_align_dict[i][j]["score"] > best_local_score:
                best_local_score = local_align_dict[i][j]["score"]
                best_x = i
                best_y = j

    # get the alignment
    backtracer(x,
               y,
               local_align_dict,
               best_x,
               best_y,
               x_name,
               y_name, 
               outname + "local_align.fasta")




    #######################
    # Analyse alignments
    #######################
    
    aligned_seq = []
    aligned_head = []
    
    for element in seq.parse(outname + "global_align.fasta", 'fasta'):
        aligned_head.append(element.id)
        aligned_seq.append(element.seq)
    
    #calculate percent id
    matches = 0.0
    for pos in range(0,(len(aligned_seq[0]))):
        if aligned_seq[1][pos] == aligned_seq[0][pos]:
            matches += 1.0
    
    
    percent_id[aligned_head[0]][aligned_head[1]] = matches/float(len(aligned_seq[0])) * 100
    best_global_scores[x_name][y_name] = best_global_score
    best_local_scores[x_name][y_name] = best_local_score


output = open("percent_identiy_table.tsv", "w")
output.write("x\ty\tPercent Identity\n")
for x in percent_id:
    for y in percent_id[x]:
        output.write(str(x) + "\t" + str(y) + "\t" + str(percent_id[x][y]) + "\n")
output.close()


output = open("pairwise_best_global_scores.tsv", "w")
output.write("x\ty\tScore\n")
for x in best_global_scores:
    for y in best_global_scores[x]:
        output.write(str(x) + "\t" + str(y) + "\t" + str(best_global_scores[x][y]) + "\n")
output.close()

output = open("pairwise_best_local_scores.tsv", "w")
output.write("x\ty\tScore\n")
for x in best_local_scores:
    for y in best_local_scores[x]:
        output.write(str(x) + "\t" + str(y) + "\t" + str(best_local_scores[x][y]) + "\n")
output.close()








