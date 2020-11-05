"""
Methods to perform global sequence alignment with 2 sequences using
an affine gap penalty with input scoring parameters. Parameters for
the function are:
The X sequence and the Y sequence to be aligned.
p1 - The match value. This is the score added in the scoring matrix
when matching characters are aligned.
p2 - The mismatch penalty. This value is substracted when the different
characters are aligned.
g - The gap open penalty, added when a gap is opened.
s - The gap extension penalty, added when a gap is extended.

To use this function, run the program with either 1 or 6 command
line arguments.

If running with input in a text file, use a single command line argument
indicating the input file path. File shoule be 3 lines: line 1 is the 4
scoring parameters, line 2 is the x sequence and line 3 is the y sequence.
Example:
$ python affine_align.py input_file.txt
File:
2 1 0 2
TACGAGTACGA
ACTGACGACTGAC

If running with inputs passed in directly, use 6 command line arguments with
the arguments as follows:
1 - The X sequence
2 - The Y sequence
3 - p1, the match value
4 - p2, the mismatch penalty
5 - g, the gap open penalty
6 - s, the gap extension penalty.
Example:
$ python affine_align.py TACGAGTACGA ACTGACGACTGAC 2 1 0 2.


Output will be printed to the terminal.


Format of scoring matrix used to perform sequence alignment and traceback:
Vertical sequence is y and indexed by i.
Horizontal sequence is x and indexed by j.
M[i][0] is first column
M[0][j] is first row

M:
   x,j ->
       A   C  G  T  .  .  .
y,i A  n   n  n  
|   C  n   n
V   G  n
    T
.
.
.
"""

import math
import sys

def affine_align(x, y, p1, p2, g, s):
    """Method to perform affine gap penalty global alignment. Note that I've switched
    the labelling of the matrix from what we used in class, so Iy is the matrix
    that represents aligning a character in Y with a gap (a vertical move). Ix
    represents aligning a character in X with a gap (a horizontal move). Therefore
    Iy > Ix for resolving ties.
    Scoring parameters:
    p1 = match
    p2 = mismatch
    g = gap open penalty
    e = gap extension penalty
    """
    #Create M, Ix, and Iy as Y x X matrices of 0's
    M = [[0]*(len(x)+1) for i in range(len(y)+1)]
    Ix = [[0]*(len(x)+1) for i in range(len(y)+1)]
    Iy = [[0]*(len(x)+1) for i in range(len(y)+1)]
    #Set up initial values for Ix and Iy
    #M infs along both axes
    for i in range(1, len(y)+1):
        M[i][0] = -math.inf
    for j in range(1, len(x)+1):
        M[0][j] = -math.inf
    #Ix: Aligning X with gap, horizontal move, infs along top row
    for i in range(0, len(y)+1):
        Ix[i][0] = -math.inf
    #Gap penalties along left column
    for j in range(1, len(x)+1):
        Ix[0][j] = -g if Ix[0][j-1] == -math.inf else Ix[0][j-1] - s
    #Iy: Aligning Y with gap, vertical move, infs along left column
    for j in range(0, len(x)+1):
        Iy[0][j] = -math.inf
    #Gap penalties along top row
    for i in range(1, len(y)+1):
        Iy[i][0] = -g if Iy[i-1][0] == -math.inf else Iy[i-1][0] - s
    #Populate remaining cells
    for i in range(1, len(y)+1):
        for j in range(1, len(x)+1):
            M[i][j] = max(M[i-1][j-1] + delta(x[j-1], y[i-1], p1, p2),
                          Ix[i-1][j-1] + delta(x[j-1], y[i-1], p1, p2),
                          Iy[i-1][j-1] + delta(x[j-1], y[i-1], p1, p2))
            Ix[i][j] = max(M[i][j-1] - g,
                           Ix[i][j-1] - s)
            Iy[i][j] = max(M[i-1][j] - g,
                           Iy[i-1][j] - s)
    #TRACEBACK
    x_ret=""; y_ret=""
    i = len(y); j = len(x)
    #Determine start matrix
    align_scores = (M[i][j], Iy[i][j], Ix[i][j])
    matrix_idx = align_scores.index(max(align_scores))
    #matrix_key will track the current matrix through the traceback
    matrix_key = ["M", "Iy", "Ix"][matrix_idx]
    while i > 0 and j > 0:
        #From M: Check diagonal moves back to all three matrices, align characters
        if matrix_key == "M":
            if M[i][j] == M[i-1][j-1] + p1 or M[i][j] == M[i-1][j-1] - p2:
                x_ret = x[j-1] + x_ret
                y_ret = y[i-1] + y_ret
                i -= 1; j -= 1
                matrix_key = "M"
            elif M[i][j] == Iy[i-1][j-1] + p1 or M[i][j] == Iy[i-1][j-1] - p2:
                x_ret = x[j-1] + x_ret
                y_ret = y[i-1] + y_ret
                i -= 1; j -= 1
                matrix_key = "Iy"
            elif M[i][j] == Ix[i-1][j-1] + p1 or M[i][j] == Ix[i-1][j-1] - p2:
                x_ret = x[j-1] + x_ret
                y_ret = y[i-1] + y_ret
                i -= 1; j -= 1
                matrix_key = "Ix"
        #From Iy: Check vertical move to Iy and M, align y character with x gap
        elif matrix_key == "Iy":
            if Iy[i][j] == M[i-1][j] - g:
                x_ret = "_" + x_ret
                y_ret = y[i-1] + y_ret
                i -= 1
                matrix_key = "M"
            elif Iy[i][j] == Iy[i-1][j] - s:
                x_ret = "_" + x_ret
                y_ret = y[i-1] + y_ret
                i -= 1
                matrix_key = "Iy"
        #From Ix: Check horizontal move to Ix and M, align x character with y gap
        elif matrix_key == "Ix":
            if Ix[i][j] == M[i][j-1] - g:
                x_ret = x[j-1] + x_ret
                y_ret = "_" + y_ret
                j -= 1
                matrix_key = "M"
            elif Ix[i][j] == Ix[i][j-1] - s:
                x_ret = x[j-1] + x_ret
                y_ret = "_" + y_ret
                j -= 1
                matrix_key = "Ix"
    #Finish sequence if edge was reached
    #i>0 means mach remaining characters in y with gaps in x
    if i > 0:
        x_ret = ("_"*i) + x_ret
        y_ret = y[0:i] + y_ret
    #j>0 means mach remaining characters in x with gaps in y
    if j > 0:
        x_ret = x[0:j] + x_ret
        y_ret = ("_"*j) + y_ret
    #Return alinged strings
    return (x_ret, y_ret)
    
def delta(char1, char2, p1, p2):
    """Function for matching characters. Returns
    p1 if characters match, p2 if mismatch.
    """
    if char1 == char2:
        return p1
    else:
        return -p2

#MAIN TAG FOR TESTING
if __name__ == "__main__":
    if len(sys.argv) == 2:
        #Read file, write to screen
        with open(sys.argv[1]) as file:
            inputs = file.readlines()
            print("INPUTS:")
            print(inputs[0])
            print(inputs[1])
            print(inputs[2])
            inputs[0] = [int(x) for x in inputs[0].split()]
            outputs = affine_align(inputs[1], inputs[2],
                                   inputs[0][0], inputs[0][1],
                                   inputs[0][2], inputs[0][3])
            print("OUTPUT:")
            for output in outputs:
                print(output)
    elif len(sys.argv) == 7:
        x = sys.argv[1]; y = sys.argv[2]
        p1 = int(sys.argv[3]); p2 = int(sys.argv[4])
        g = int(sys.argv[5]); s = int(sys.argv[6])
        print("INPUT:")
        print(p1, p2, g, s)
        print(x)
        print(y)
        outputs = affine_align(x, y, p1, p2, g, s)
        print("OUTPUT: ")
        for output in outputs:
            print(output)
    else:
        print("Incorrect input . Must run with 1 or 6 arguments")
