#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from multiprocessing import Pool
import itertools
import pandas as pd
import numpy as np

#scoring function
def calculate_element_score3(match_award, mismatch_penalty, gap_penalty, seq1, seq2):
    m, n = len(seq1), len(seq2)  # length of two sequences

    # Generate DP table and traceback path pointer matrix
    score = np.zeros((m+1, n+1)).astype(int)      # the DP table
    
    max_score = 0        # initial maximum score in DP table
    # Calculate DP table

    scores_list = []

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            score_diagonal = score[i-1][j-1] + match_score(seq1[i-1], seq2[j-1])
            score_up = score[i][j-1] + gap_penalty
            score_left = score[i-1][j] + gap_penalty
            score[i][j] = max(0,score_left, score_up, score_diagonal)
            if score[i][j] >= max_score:
                max_i = i
                max_j = j
                max_score = score[i][j];

    perfect_score = min(len(seq1), len(seq2)) * match_award

    perfect_matches = np.count_nonzero(score == perfect_score)

    first_tolerance_match = np.count_nonzero(score == perfect_score - 1)
    second_tolerance_match = np.count_nonzero(score == perfect_score - 2)
    third_tolerance_match = np.count_nonzero(score == perfect_score - 3)

    final_score = perfect_matches + first_tolerance_match/(perfect_score * 2) + second_tolerance_match/(perfect_score * 4) + third_tolerance_match/(perfect_score * 8)
    
    return final_score

def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    elif alpha == '-' or beta == '-':
        return gap_penalty
    else:
        return mismatch_penalty
    


    
#set other parameters
match_award      = 2
mismatch_penalty = -1
gap_penalty      = -1 # both for opening and extanding

#s is an individual sequence handed to the function
#score is returned as a numerical value

def f(s):
    
    with open('temp.txt', 'r') as file2:
        ele = file2.readline()
        
        
    score = calculate_element_score3(match_award, mismatch_penalty, gap_penalty, ele, s)

    return score



if __name__ == '__main__':
    
    
    #load data

    df = pd.read_csv('final file hopefully no errors.csv')

    dfp = df.loc[df['Gen'] == 'ssRNA(+)']
    dfn = df.loc[df['Gen'] == 'ssRNA(-)']

    dfp = dfp.reset_index(drop = True)
    dfn = dfn.reset_index(drop = True)
    
    
    #make windows

    winlist6 = []

    for i in range(6,7):
        for i in (itertools.product('ATGC', repeat=i)):
            winlist6.append(''.join(i))
            
            
    winlist = winlist6
    
    
    #working dataframe and result df 
    df = pd.concat([dfp,dfn])
    resultdf = pd.DataFrame(index = list(range(len(df))), columns = winlist)

          
    #working loop
    
    for element in winlist:
        pool = Pool(processes=14)
        
        with open('temp.txt', 'w') as file:
            file.write(element)
            
        it = iter(df['Seq'].to_list())
        
        #result returned as a list the length of iterator it
        result = pool.map(f, it)
                
        resultdf[element] = result
        
        pool.close()
                
    resultdf.to_csv('win6 .csv', index = False)

