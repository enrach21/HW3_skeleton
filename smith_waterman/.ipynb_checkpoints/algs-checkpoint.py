import numpy as np
import pandas as pd
import glob
from Bio import SeqIO
import random

def make_matrix(sizex, sizey):
    """Creates a sizex by sizey matrix filled with zeros.
    """
    test = [[0]*sizey for i in  range(sizex)]
    return test



def getSeq(x):
    """ Will return the sequence from sequenceing folder when given file location
    """
    for record in SeqIO.parse(x, "fasta"):
        seq = (record.seq)
    return seq

def readPosPairs():
    """ returns a data frame of the positive pairs
    """
    BLOSUM = glob.glob('Pospairs.txt')
    data = pd.read_table(BLOSUM[0], comment='#', delim_whitespace=True)
    return data

def readNegPairs():
    """ returns a data frame of the negative pairs
    """
    BLOSUM = glob.glob('Negpairs.txt')
    data = pd.read_table(BLOSUM[0], comment='#', delim_whitespace=True)
    return data

def readMatrix(text):
    """ reads in a matrix and converts it into a dictionary for scoring
    """
    BLOSUM = glob.glob(text)
    data = pd.read_table(BLOSUM[0], comment='#', delim_whitespace=True)
    data.index = data.columns.values
    Dict = {}
    for i in data.columns.values:
        for j in data.index.values:
            Dict[(i,j)]=data[i][j]
    return Dict

def get_AA():
    """ Function that will return a list of the Amino acids used in each matrix
    """
    BLOSUM = glob.glob('BLOSUM50')
    data = pd.read_table(BLOSUM[0], comment='#', delim_whitespace=True)
    x = data.columns.values
    return x

def sw(s1, s2, gap, extention, text):
    """
    x is one seq
    y is the other seq
    gap is the gap penalty
    extention is the extention penalty
    text is the matrix
    output, is the score and location
    """
    Dict = readMatrix(text)
    F = make_matrix(len(s1) + 1, len(s2) + 1) # Keeps track of alignment scores
    X = make_matrix(len(s1) + 1, len(s2) + 1) # Keeps sequence X gaps and extentions
    Y = make_matrix(len(s1) + 1, len(s2) + 1) # Keeps sequence Y gaps and extentions
    T = make_matrix(len(s1) + 1, len(s2) + 1) # Keeps sequence Y gaps and extentions
    best = 0
    optloc = (0,0)
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            X[i][j]= max(0,
                        F[i-1][j] - gap, # new gap
                        X[i-1][j] - extention # contiue to extend
                        ) # Makes a gap on the 'Y' strand
            Y[i][j]= max(0,
                        F[i][j-1] - gap, # New Gap
                        Y[i][j-1] - extention # Continue to extend
                        ) # Makes a gap on the 'X' strand
            F[i][j]= max(0,
                        F[i-1][j-1] + Dict[(s1[i-1],s2[j-1])],
                        X[i-1][j-1] + Dict[(s1[i-1],s2[j-1])],
                        Y[i-1][j-1] + Dict[(s1[i-1],s2[j-1])], 
                        )
            if F[i][j] >= best:
                best = F[i][j]
                optloc = (i,j)
            if(max(X[i][j],Y[i][j],F[i][j]) == 0):
                T[i][j]=0
            elif (F[i][j] == max(X[i][j],Y[i][j],F[i][j])):
                T[i][j]=1
            elif (X[i][j] == max(X[i][j],Y[i][j],F[i][j])):
                 T[i][j]=2
            elif (Y[i][j] == max(X[i][j],Y[i][j],F[i][j])):
                T[i][j]=3
                
            
    print(X)
    print(Y)
    print(F)
    print(T)
    
    # Perform Traceback
    i = optloc[0]
    # print(i)
    j = optloc[1]
    # print(j)
    align1 = np.array([])
    align2 = np.array([])
    
        
    
    while T[i][j] != 0: # This would be the start of the sequence  
        if T[i][j] == 1:
            a1 = s1[i-1]
            a2 = s2[j-1]
            i -= 1
            j -= 1
        elif T[i][j] == 2:
            a1 = s1[i-1]
            a2 = '-'
            i -= 1
        elif T[i][j] == 3:
            a1 = '-'
            a2 = s2[j-1]
            j -= 1
        # print(i)
        # print(j)
        align1 = np.append(a1, align1)
        align2 = np.append(a2, align2)
        print(align1)
        print(align2)
    
    # print(align1)
    # print(align2)
    # return the opt score and the best location
    return best, optloc, align1, align2

def sw_new(s1, s2, gap, extention, d):
    """
    x is one seq
    y is the other seq
    gap is the gap penalty
    extention is the extention penalty
    text is the matrix
    output, is the score and location
    """
    Dict = d
    fromF, fromX, fromY = 1, 2, 3
    F = make_matrix(len(s1) + 1, len(s2) + 1) # Keeps track of alignment scores
    X = make_matrix(len(s1) + 1, len(s2) + 1) # Keeps sequence X gaps and extentions
    Y = make_matrix(len(s1) + 1, len(s2) + 1) # Keeps sequence Y gaps and extentions
    T = make_matrix(len(s1) + 1, len(s2) + 1) # Keeps sequence Y gaps and extentions
    
    # Initialize the out part of the matrix
    for i in range(0, len(s1)+1):
        F[i][0] = (0, fromX)
        X[i][0] = (0, fromX)
        Y[i][0] = (0, fromX)
    
    for i in range(0, len(s2)+1):
        F[0][i] = (0, fromY)
        X[0][i] = (0, fromY)
        Y[0][i] = (0, fromY)
    
    # Keep track of highest score and it's location
    best = 0
    optloc = (0,0)
    
    # Go throught the matrix and update with optimal score at each postion
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            X[i][j]= max((0,fromX),
                        ((F[i-1][j])[0]-gap, fromF), # new gap
                        ((X[i-1][j])[0]-extention, fromX) # contiue to extend
                        ) # Makes a gap on the 'Y' strand
            Y[i][j]= max((0,fromY),
                        ((F[i][j-1])[0] - gap,fromF), # New Gap
                        ((Y[i][j-1])[0] - extention,fromY) # Continue to extend
                        ) # Makes a gap on the 'X' strand
            F[i][j]= max((0,fromF),
                        ((F[i-1][j-1])[0] + Dict[(s1[i-1],s2[j-1])],fromF),
                        ((X[i-1][j-1])[0] + Dict[(s1[i-1],s2[j-1])],fromX),
                        ((Y[i-1][j-1])[0] + Dict[(s1[i-1],s2[j-1])],fromY), 
                        )
            if (F[i][j])[0] >= best:
                best = (F[i][j])[0]
                optloc = (i,j)
    
    # Perform Traceback
    i = optloc[0] # the i postion of where the matrix should begin it's traceback
    j = optloc[1] # the j postion of where the matrix should begin it's traceback
    # Best aligned sequences
    align1 = np.array([])
    align2 = np.array([])
    
    current = F # The current matrix that we are investigating
    
    while current[i][j][0]!=0: # While the score is not equal to zero
        if current == F: # if in main matrix
            if current[i][j][1] == fromF: # if current cell is fromF
                a1 = s1[i-1] # aa to be added to first seq 
                a2 = s2[j-1] # aa to be added to second seq 
                i -= 1
                j -= 1
                current = F
            elif current[i][j][1] == fromX: # if current cell is fromX
                a1 = s1[i-1] # aa to be added to first seq 
                a2 = s2[j-1] # aa to be added to second seq 
                i -= 1
                j -= 1
                current = X
            elif current[i][j][1] == fromY:
                a1 = s1[i-1] # aa to be added to first seq 
                a2 = s2[j-1] # aa to be added to second seq 
                i -= 1
                j -= 1
                current = Y
        
        elif current == X:
            if current[i][j][1] == fromF:
                a1 = s1[i-1] # aa to be added to first seq 
                a2 = '-' # Perfrom a skip
                i -= 1
                current = F
            elif current[i][j][1] == fromX:
                a1 = s1[i-1]  # aa to be added to first seq 
                a2 = '-' 
                i -= 1
                current = X              
         
        elif current == Y: 
            if current[i][j][1] == fromF:
                a1 = '-'
                a2 = s2[j-1] # aa to be added to second seq 
                j -= 1
                current=F
            elif current[i][j][1] == fromY:
                a1 = '-'
                a2 = s2[j-1] # aa to be added to second seq 
                j -= 1
                current=Y
        align1 = np.append(a1, align1)
        align2 = np.append(a2, align2)
   
    
    # print(X)
    # print(Y)
    # print(F)
    
    # print(align1)
    # print(align2)
    # return the opt score and the best location
    return best, optloc, align1, align2

def sw_norm(x, y, gap, extention, d):
    """
    x is one seq
    y is the other seq
    gap is the gap penalty
    extention is the extention penalty
    text is the matrix
    output, is the score and location
    """
    minimum = min(len(x), len(y))
    Dict = d
    F = make_matrix(len(x) + 1, len(y) + 1) # Keeps track of alignment scores
    X = make_matrix(len(x) + 1, len(y) + 1) # Keeps sequence X gaps and extentions
    Y = make_matrix(len(x) + 1, len(y) + 1) # Keeps sequence Y gaps and extentions
    best = 0
    optloc = (0,0)
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            X[i][j]= max(0,
                        F[i-1][j] - gap, # new gap
                        X[i-1][j] - extention # contiue to extend
                        )
            Y[i][j]= max(0,
                        F[i][j-1] - gap, # New Gap
                        Y[i][j-1] - extention # Continue to extend
                        )
            F[i][j]= max(0,
                        F[i-1][j-1] + Dict[(x[i-1],y[j-1])],
                        X[i-1][j-1] + Dict[(x[i-1],y[j-1])],
                        Y[i-1][j-1] + Dict[(x[i-1],y[j-1])], 
                        )
            if F[i][j] >= best:
                best = F[i][j]
                optloc = (i-1,j-1)
            # print(X)
            # print(Y)
            # print(F)
    # return the opt score and the best location
    return best/minimum, optloc



def test_known(gap, extend, d):
    """
     Generates a data frame that has each sequence,
     whether it is a pos or neg control
     the score it recieved
     The local alignment
    """
    # Read in negative and pos score
    pos = readPosPairs()
    neg = readNegPairs()
    # Make empty array
    pair1 = np.array([])
    pair2 = np.array([])
    known = np.array([])
    score = np.array([])
    align1 = []
    align2 = []
    # Get score for known pos
    for x in range(len(pos[1])):
        # print(x)
        A = getSeq(pos[0][x])
        B = getSeq(pos[1][x])
        s = sw_new(A, B, gap, extend, d)
        pair1 = np.append(pair1, (pos[0][x]))
        pair2 = np.append(pair2, (pos[1][x]))
        score = np.append(score, s[0])
        known = np.append(known, 1)
        align1 += ([s[2]])
        align2 += ([s[3]])
    # get scores for known negative
    for x in range(len(neg[1])):
        # print(x)
        A = getSeq(neg[0][x])
        B = getSeq(neg[1][x])
        s = sw_new(A, B, gap, extend, d)
        pair1 = np.append(pair1, (neg[0][x]))
        pair2 = np.append(pair2, (neg[1][x]))
        score = np.append(score, s[0])
        known = np.append(known, 0)
        align1 += ([s[2]])
        align2 += ([s[3]])
    # print(align1)
    df = pd.DataFrame({'pair1': pair1,'pair2': pair2, 'score': score, 'known':known, 'alignment1':align1, 'alignment2':align2})
    return df

def test_known_norm(gap, extend, text):
    # Read in negative and pos score
    pos = readPosPairs()
    neg = readNegPairs()
    # Make empty array
    pair1 = np.array([])
    pair2 = np.array([])
    known = np.array([])
    score = np.array([])
    # Get score for known pos
    for x in range(len(pos[1])):
        # print(x)
        A = getSeq(pos[0][x])
        B = getSeq(pos[1][x])
        s = sw_norm(A, B, gap, extend, text)
        pair1 = np.append(pair1, (neg[0][x]))
        pair2 = np.append(pair2, (neg[1][x]))
        score = np.append(score, s[0])
        known = np.append(known, 1)
    # get scores for known negative
    for x in range(len(neg[1])):
        # print(x)
        A = getSeq(neg[0][x])
        B = getSeq(neg[1][x])
        s = sw_norm(A, B, gap, extend, text)
        pair1 = np.append(pair1, (neg[0][x]))
        pair2 = np.append(pair2, (neg[1][x]))
        score = np.append(score, s[0])
        known = np.append(known, 0)
    df = pd.DataFrame({'pair1': pair1,'pair2': pair2, 'score': score, 'known':known})
    return df

def score_70percent(test):
    """
    returns the score that contains 70% of all True_positives
    """
    known_pos = test[test.known==1] # based on the scores in test calc the score required to get 0.7 true pos called pos
    sorted_known_pos = known_pos.sort_values(by='score') # sort based on score
    total = len(sorted_known_pos.score) # get the postion of the score that has 70% of the positive scores above
    score = np.array(sorted_known_pos.score)
    score = score[int(total-(total*0.7))]
    return score

def false_pos(test, score):
    known_neg = test[test.known==0]
    total=len(known_neg.score)
    # print(known_neg)
    f =len(known_neg[known_neg.score>=score])
    return f/total

def score(test, m):
    known_pos = test[test.known==1] # based on the scores in test calc the score required to get 0.7 true pos called pos
    sorted_known_pos = known_pos.sort_values(by='score') # sort based on score
    
    known_neg = test[test.known==0]
    total_neg=len(known_neg.score)
    
    score = np.array(sorted_known_pos.score)
    # print(score)
    # print(known_neg.score)
    # y = np.array([])
    x = np.array([])
    for i in range(0,len(score)):
        true_p = 1-(i/len(score)) # percentage of true postive
        s = score[i] # threshold to be >= to be true
        f =len(known_neg[known_neg.score>=s])
        false_p = f/total_neg
        # y = np.append(y, true_p)
        x = np.append(x, false_p)
        
    df=pd.DataFrame({m:x})    
    return df


def readNegPairs():
	data = pd.read_csv('/Users/ijones1/Documents/HW3_skeleton/Negpairs.txt', header = None, sep=" ")
	return data

def readPosPairs():
	data = pd.read_csv('/Users/ijones1/Documents/HW3_skeleton/Pospairs.txt', header = None, sep=" ")
	return data

def get_scores(align1, align2, D, g, e, known):
    s = np.array([])
    for x in range(len(align1)):
        seq1=align1[x]
        seq2=align2[x]
        score = 0
        for y in range(len(seq1)):
            if seq1[y] == '-':
                if seq1[y-1] == '-':
                    score -= e
                else:
                    score -= g
            elif seq2[y] == '-':
                if seq2[y-1] == '-':
                    score -= e
                else:
                    score -= g
            else:
                score += D[(seq1[y],seq2[y])]
        s = np.append(s, score)
    df=pd.DataFrame({'known':known,'score':s}) 
    return df


def score_neg(df):
    known_pos = df[df.known==1] 
    total_pos=len(known_pos.score)
    
    known_neg = df[df.known==0]
    sorted_known_neg = known_neg.sort_values(by='score')
    
    score = np.array(sorted_known_neg.score)
    L = len(score)#obtain the length of negative controls
    final_s = 0
    for x in [0.0,0.1,0.2,0.3]:
        i = L-(x*L)
        if i == L:
            s = score[-1]+1
        else:
            s = score[int(i)]
        final_s += (len(known_pos[known_pos.score>=s])/len(known_pos.score))
        # print (x)
        # print(final_s)
    return final_s

def ga_init(text):
    n = 200 # This is the starting population
    pop = np.array([])
    for x in range(n):
        pop = np.append(pop, readMatrix(text) )
    return pop

def ga_mut(pop, aa):
    l = len(aa)-1
    for i in range(len(pop)):
        # print(i)
        # Cause three random mutations
        rand_1a = aa[random.randint(0, l)] # amino acid to be changed
        rand_1b = aa[random.randint(0, l)] # amino acid to be changed
        # print(rand_1a)
        # print(rand_1b)
        rand_2a = aa[random.randint(0,l)] # amino acid to be changed
        rand_2b = aa[random.randint(0,l)] # amino acid to be changed
        # print(ran_index2a)
        rand_3a = aa[random.randint(0,l)] # amino acid to be changed
        rand_3b = aa[random.randint(0,l)] # amino acid to be changed
        # print(ran_index3a)
        change = random.randint(-3,3) # the amount to change the three 
        # print(change)
        pop[i][(rand_1a,rand_1b)] = pop[i][(rand_1a,rand_1b)] + change
        pop[i][(rand_1b,rand_1a)] = pop[i][(rand_1b,rand_1a)] + change
        pop[i][(rand_2a,rand_2b)] = pop[i][(rand_2a,rand_2b)] + change
        pop[i][(rand_2b,rand_2a)] = pop[i][(rand_2b,rand_2a)] + change
        pop[i][(rand_3a,rand_3b)] = pop[i][(rand_3a,rand_3b)] + change
        pop[i][(rand_3b,rand_3a)] = pop[i][(rand_3b,rand_3a)] + change
    return pop

def ga_mut_test(pop, aa):
    l = len(aa)-1

    # print(i)
    # Cause three random mutations
    rand_1a = aa[random.randint(0, l)] # amino acid to be changed
    rand_1b = aa[random.randint(0, l)] # amino acid to be changed
    # print(rand_1a)
    # print(rand_1b)
    rand_2a = aa[random.randint(0,l)] # amino acid to be changed
    rand_2b = aa[random.randint(0,l)] # amino acid to be changed
    # print(ran_index2a)
    rand_3a = aa[random.randint(0,l)] # amino acid to be changed
    rand_3b = aa[random.randint(0,l)] # amino acid to be changed
    # print(ran_index3a)
    change = random.randint(-3,3) # the amount to change the three 
    # print(change)
    pop[(rand_1a,rand_1b)] = pop[(rand_1a,rand_1b)] + change
    pop[(rand_1b,rand_1a)] = pop[(rand_1b,rand_1a)] + change
    pop[(rand_2a,rand_2b)] = pop[(rand_2a,rand_2b)] + change
    pop[(rand_2b,rand_2a)] = pop[(rand_2b,rand_2a)] + change
    pop[(rand_3a,rand_3b)] = pop[(rand_3a,rand_3b)] + change
    pop[(rand_3b,rand_3a)] = pop[(rand_3b,rand_3a)] + change
    return pop


def ga_eval(pop, align1, align2, g, e, known): # Not working need to fix if time
    base = get_scores(align1, align2, readMatrix('PAM100'), g, e, known)
    baseline = score_neg(base)
    for i in range(len(pop)):
        # print(i)
        D = pop[i]
        # print(D)
        t = get_scores(align1, align2, D, g, e, known)
        # print(t)
        score = score_neg(t)
        print(score)
        if score < baseline:
            print('changed')
            pop[i] = base
            
            
def test(text, align1, align2, g, e, known):
    aa = get_AA()
    pop = ga_init(text)
    k = 0
    orig_m = readMatrix(text)
    while k < 20:
        for i in range(len(pop)):
            # print(i)
            t = get_scores(align1, align2, pop[i], g, e, known)
            # print(t)
            score = score_neg(t)
            # print(score)
            D_new = ga_mut_test(pop[i], aa)
            t_new = get_scores(align1, align2, D_new, g, e, known)
            # print(t_new)
            score_new = score_neg(t_new)
            # print(score_new)
            if score_new >= score:
                # print('change')
                pop[i] = D_new
            else:
                pop[i] = readMatrix(text)
        k += 1

    # recover the best matrix
    best_s = 0
    best_m = 0
    for i in range(len(pop)):
        x = pop[i]
        t = get_scores(align1, align2, x, g, e, known)    
        score = score_neg(t)
        if best_s < score:
            print('the best score is...')
            print(score)
            best_s = score
            best_m = x
    return orig_m, best_m