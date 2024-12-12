from Bio.Align import substitution_matrices
import numpy as np

# Match score from BLOSUM92 matrix with Bio.Align.substitution_matrices
def match_score(A, B):
    scoring_matrix = substitution_matrices.load('BLOSUM62')
    return scoring_matrix.get((A, B))

# Smith-Waterman Local Alignment Algorithm
def SW(seqX, seqY, op_pen = -5, ex_pen = -1):
    col = len(seqX) + 1
    row = len(seqY) + 1
    
    dia_mat = np.zeros((row, col)) # matrix for moving diagonally
    hor_mat = np.zeros((row, col))  # matrix for moving horizontally
    ver_mat = np.zeros((row, col))  # matrix for moving vertically
    
    # Populating matrices
    for x in range(1, col):
        for y in range(1, row):
            # print(x, y)
            hor_mat[y][x] = max(hor_mat[y][x-1] + ex_pen, dia_mat[y][x-1] + ex_pen + op_pen, 0)
            ver_mat[y][x] = max(ver_mat[y-1][x] + ex_pen, dia_mat[y-1][x] + ex_pen + op_pen, 0)
            dia_mat[y][x] = max(hor_mat[y][x], ver_mat[y][x], dia_mat[y-1][x-1] + match_score(seqY[y-1], seqX[x-1]), 0)
            
    # Tracing back
    y, x = np.unravel_index(np.argmax(dia_mat), dia_mat.shape)
    s = np.max(dia_mat)
    
    curr_mat = 0 # 0 for diagonal, 1 for horizonal, 2 for vertical
    
    trace_x = []
    trace_y = []

    while x > 0 and y > 0 and s > 0:
        # print(x, y, curr_mat, s)
        # print(trace_x)
        # print(trace_y)
        
        if curr_mat == 0:
            if s == dia_mat[y-1][x-1] + match_score(seqY[y-1], seqX[x-1]):
                trace_x.append(seqX[x-1])
                trace_y.append(seqY[y-1])
                x -= 1
                y -= 1
                s = dia_mat[y][x]
                    
            elif s == hor_mat[y][x]:
                curr_mat = 1
                
            elif s == ver_mat[y][x]:
                curr_mat = 2
                
        elif curr_mat == 1:
            if s == hor_mat[y][x-1] + ex_pen:
                trace_x.append(seqX[x-1])
                trace_y.append('-')
                x -= 1
                s = hor_mat[y][x]
                
            elif s == dia_mat[y][x-1] + ex_pen + op_pen:
                trace_x.append(seqX[x-1])
                trace_y.append('-')
                x -= 1
                curr_mat = 0
                s = dia_mat[y][x]
                
        elif curr_mat == 2:
            if s == ver_mat[y-1][x] + ex_pen:
                trace_x.append('-')
                trace_y.append(seqY[y-1])
                y -= 1
                s = ver_mat[y][x]
                
            elif s == dia_mat[y-1][x] + ex_pen + op_pen:
                trace_x.append('-')
                trace_y.append(seqY[y-1])
                y -= 1
                curr_mat = 0
                s = dia_mat[y][x]

    return trace_x[::-1], trace_y[::-1]
    
# Test case
# [algnx, algny] = SW(["F", "K", "E", "R", "A", "F", "F", "Q", "W"],["F", "K", "E", "N", "N", "R", "A", "Q", "W"])