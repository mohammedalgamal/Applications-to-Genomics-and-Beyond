"""
project 4 in algorithmic thinking
"""

def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    """
    Takes as input a set of characters alphabet
    and three scores diag_score, off_diag_score, and dash_score.
    The function returns a dictionary of dictionaries
    whose entries are indexed by pairs of characters in alphabet plus -.
    """
    
    ans = {}
    alpha = list(alphabet)
    alpha.append('-')
    
    for idx in alpha:
        ans[idx] = {}
        
        for idx2 in alpha:
            
            if idx == idx2:
                ans[idx][idx2] = diag_score
                if idx == '-':
                    ans[idx][idx2] = dash_score
            elif idx == '-' or idx2 == '-':
                ans[idx][idx2] = dash_score
            else:
                ans[idx][idx2] = off_diag_score

    return ans


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    """
    Takes as input two sequences seq_x and seq_y
    whose elements share a common alphabet with the scoring matrix.
    The function computes and returns the alignment matrix for seq_x and seq_y.
    """
    
    if not global_flag:
        mul = 0
    else:
        mul = 1
        
    mlen = len(seq_x)
    nlen = len(seq_y)
    table_s = [[0 for dummy_col in range(nlen + 1)] for dummy_row in range(mlen + 1)]

    for idx in range(1, mlen + 1):
        table_s[idx][0] = (table_s[idx - 1][0] + scoring_matrix[seq_x[idx - 1]]['-'])
        if table_s[idx][0] < 0:
            table_s[idx][0] = mul * table_s[idx][0]

    for jdx in range(1, nlen + 1):
        table_s[0][jdx] = (table_s[0][jdx - 1] + scoring_matrix['-'][seq_y[jdx - 1]])
        if table_s[0][jdx] < 0:
            table_s[0][jdx] = mul * table_s[0][jdx]

    for idx in range(1, mlen + 1):
        for jdx in range(1, nlen + 1):
            table_s[idx][jdx] = (max(table_s[idx - 1][jdx] + scoring_matrix[seq_x[idx - 1]]['-'],
                                 table_s[idx][jdx - 1] + scoring_matrix['-'][seq_y[jdx - 1]],
                                 table_s[idx - 1][jdx - 1] + scoring_matrix[seq_x[idx - 1]][seq_y[jdx - 1]]))
            if table_s[idx][jdx] < 0:
                table_s[idx][jdx] = mul * table_s[idx][jdx]
            
                                    
    return table_s


def score(seq_x, seq_y, scoring_matrix):
    """
    returns the score of 2 alignments
    """
    
    rscore = 0
    idx = 0
    
    if len(seq_x) != len(seq_y):
        return
    else:
        while(idx < len(seq_x)):
            rscore += scoring_matrix[seq_x[idx]][seq_y[idx]]
            idx += 1

    return rscore        

def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    This function computes a global alignment of seq_x and seq_y
    using the global alignment matrix.
    """
    
    ilen = len(seq_x)
    jlen = len(seq_y)
    alignx = ''
    aligny = ''
    while ilen != 0 and jlen != 0:
        if alignment_matrix[ilen][jlen] == alignment_matrix[ilen - 1][jlen - 1] + scoring_matrix[seq_x[ilen - 1]][seq_y[jlen - 1]]:
            alignx = seq_x[ilen - 1] + alignx
            aligny = seq_y[jlen - 1] + aligny
            ilen -= 1
            jlen -= 1
        else:
            if alignment_matrix[ilen][jlen] == alignment_matrix[ilen - 1][jlen] + scoring_matrix[seq_x[ilen - 1]]['-']: 
                alignx = seq_x[ilen - 1] + alignx
                aligny = '-' + aligny
                ilen -= 1
            else:
                alignx = '-' + alignx
                aligny = seq_y[jlen - 1] + aligny
                jlen -= 1

    while ilen != 0:
        alignx = seq_x[ilen - 1] + alignx
        aligny = '-' + aligny
        ilen -= 1
        
    while jlen != 0:
        alignx = '-' + alignx
        aligny = seq_y[jlen - 1] + aligny
        jlen -= 1

        
    return (score(alignx, aligny, scoring_matrix), alignx, aligny)

def get_higest_score(alignment_matrix, seq_x, seq_y):
    """
    returns the highest score in alignment_matrix
    based on the lengths of seq_x, seq_y
    """
    
    ilen = 0
    jlen = 0
    max_val = float('-inf')
    
    for idx in range(len(seq_x) + 1):
        for jdx in range(len(seq_y) + 1):
            if alignment_matrix[idx][jdx] > max_val:
                max_val = alignment_matrix[idx][jdx]
                ilen = idx
                jlen = jdx 
       
    return [ilen, jlen, max_val]            

def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    """
    This function computes a local alignment of seq_x and seq_y
    using the local alignment matrix.
    """
    
    lengths = get_higest_score(alignment_matrix, seq_x, seq_y)
    ilen = lengths[0]
    jlen = lengths[1]     
    alignx = ''
    aligny = ''
    while ilen != 0 and jlen != 0 and alignment_matrix[ilen][jlen] != 0:
        if alignment_matrix[ilen][jlen] == alignment_matrix[ilen - 1][jlen - 1] + scoring_matrix[seq_x[ilen - 1]][seq_y[jlen - 1]]:
            alignx = seq_x[ilen - 1] + alignx
            aligny = seq_y[jlen - 1] + aligny
            ilen -= 1
            jlen -= 1
        else:
            if alignment_matrix[ilen][jlen] == alignment_matrix[ilen - 1][jlen] + scoring_matrix[seq_x[ilen - 1]]['-']: 
                alignx = seq_x[ilen - 1] + alignx
                aligny = '-' + aligny
                ilen -= 1
            else:
                alignx = '-' + alignx
                aligny = seq_y[jlen - 1] + aligny
                jlen -= 1
   
    return (score(alignx, aligny, scoring_matrix), alignx, aligny)








