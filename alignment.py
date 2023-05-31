# encoding: utf8
from collections import Counter
import itertools
import numpy as np
import pandas as pd
import arbre
from constants import df


class Alignment():
    '''Class performing the sequence alignment'''
    def __init__(self,fasta_file):
        '''Initialising the class'''
        self.file=fasta_file
        self.all_pairs=None
        self.sequences=None
        self.score=None
        self.score_table=None
        self.distance_table_UPGMA=None
        self.distance_table_NJ=None
        self.tree_UPGMA=None
        self.tree_NJ=None
        self.alignment_names=[]
        self.final_alignment=None
        self.indices=[]
    
    def fasta_reader(self):
        '''Method reading th efasta file'''
        with open(self.file) as fasta_file: 
            sequences = dict() #Creating an empty dictionary to store the sequences
            lines = fasta_file.readlines() 
            sequence_id="" #Empty string for the sequence id
            for line in lines: #For line in file
                line=line.strip() #Remove the special charachter at the end of the line ('\n')
                if ">" in line: #If ">" is in line
                    sequence_id=line #Store the line in sequence_id
                    sequence="" #Create an empty str to store the sequence
                    names=sequence_id.split(' ')[0] #Spliting the sequence id and taking the first element
                    sequences[names]=sequence #Adding the values to the dictionary
                else:
                    sequence+="".join(line) #Add the line to the sequence variable
                    sequences[names]=sequence #Updating the value in the dictionary
        all_pairs=list(itertools.combinations(sequences,2)) #Creating a list of all possible pairs of the sequences
        self.sequences=sequences #Passing the dictionnary of sequences to the object
        self.all_pairs=all_pairs #Passing the list of all possible pairs to the object
    
    def nw_only_score(self,sequence1,sequence2,gap=-1):
        '''Method calculating the Needleman-Wunsch score of the pairwase alignment'''
        sequence1_length=len(sequence1) 
        sequence2_length=len(sequence2)
        score_matrix=np.zeros((sequence1_length+1,sequence2_length+1)) #Creating an empty matrix 
        score_matrix[:,0]=np.linspace(0,sequence1_length*gap,sequence1_length+1) #Initialising the first column with gap scores
        score_matrix[0,:]=np.linspace(0,sequence2_length*gap,sequence2_length+1) #Initialising the first row with gap scores
        current_value=np.zeros(3) #Creating a small matrix to store the possible values
        for i in range(sequence1_length):
            for j in range(sequence2_length):
                #Applying Needleman-Wunsch formula to claculate possible scores
                current_value[0]=score_matrix[i,j]+df.loc[sequence1[i],sequence2[j]] 
                current_value[1]=score_matrix[i,j+1]+gap
                current_value[2]=score_matrix[i+1,j]+gap
                tmax=np.max(current_value) #Choosing the maximal score
                score_matrix[i+1,j+1]=tmax #Writing the score in the score matrix
        score=score_matrix[i+1,j+1] #Getting the last final score of the alignment
        self.score=score #Passing the final score to the object


    def nw(self,sequence1,sequence2,gap=-1):
        '''Method calculating the Needleman-Wunsch score of the pairwase alignment with traceback'''
        sequence1_length=len(sequence1)
        sequence2_length=len(sequence2)
        score_matrix=np.zeros((sequence1_length+1,sequence2_length+1))
        score_matrix[:,0]=np.linspace(0,sequence1_length*gap,sequence1_length+1)
        score_matrix[0,:]=np.linspace(0,sequence2_length*gap,sequence2_length+1)
        current_value=np.zeros(3)
        traceback=np.zeros((sequence1_length+1,sequence2_length+1)) #Creating matrix to store the directions
        #Initialising the first column and row
        traceback[:,0]=3 
        traceback[0,:]=4
        for i in range(sequence1_length):
            for j in range(sequence2_length):
                current_value[0]=score_matrix[i,j]+df.loc[sequence1[i],sequence2[j]]
                current_value[1]=score_matrix[i,j+1]+gap
                current_value[2]=score_matrix[i+1,j]+gap
                tmax=np.max(current_value)
                score_matrix[i+1,j+1]=tmax
                if current_value[0]==tmax: #If maximal value is a match/mismatch score
                    traceback[i+1,j+1]+=2 #Put 2 in the cell, meaning diagonal direction
                if current_value[1]==tmax: #If maximal value is gap in sequence 2
                    traceback[i+1,j+1]+=3 #Put 3 in the cell, meaning vertical direction
                if current_value[2]==tmax: #If maximal value is gap in sequence 1
                    traceback[i+1,j+1]+=4  #Put 4 in the cell, meaning horizontal direction
        i=sequence1_length 
        j=sequence2_length
        alignment_1=[] #List for the first aligned sequence
        alignment_2=[] #List for the second aligned sequence
        while i>0 or j>0: 
            if traceback[i,j] in [2,5,6,9]: 
                alignment_1.append(sequence1[i-1]) #Append character with index  i-1  to the alignment 1
                alignment_2.append(sequence2[j-1]) #Append character with index  j-1  to the alignment 2
                #Update the indexes
                i-=1 
                j-=1
            elif traceback[i,j] in[3,5,7,9]:
                alignment_1.append(sequence1[i-1])#Append character with index  i-1  to the alignment 1
                alignment_2.append('-') #Append the gap to the alignment 2
                i-=1 #Update only the index of first sequence
            elif traceback[i,j] in [4,6,7,9]:
                alignment_1.append('-') #Append the gap to the alignment 1
                alignment_2.append(sequence2[j-1]) #Append character with index  j-1  to the alignment 2
                j-=1 #Update only the index of second sequence
        #Join all the elements of the list in reversed order
        alignment_1=''.join(alignment_1)[::-1] 
        alignment_2=''.join(alignment_2)[::-1]
        return alignment_1,alignment_2 

    def score_counter(self):
        '''Method allowing to calculate the score of all pairwase alignments'''
        counter=0 #Setting the counter to 0
        score_df=pd.DataFrame(index=self.sequences,columns=self.sequences) #Creating a dataframe to store the information (with sequence names)
        for pair in self.all_pairs: #For pair in all_pairs list
            counter+=1 #Add one to counter
            self.nw_only_score(self.sequences[pair[0]],self.sequences[pair[1]]) #Calculate the score of the pair
            print("alignment number "+str(counter) + " done") #Print a message to user
            score_df.loc[pair[0],pair[1]]=score_df.loc[pair[1],pair[0]]=self.score #Store the calculated value in the corresponding cell
        score_df=score_df.fillna(0) #Fill the NaN values with zeros
        self.score_table=score_df.to_numpy() #Transform into matrix and pass it to object
        return self.score_table #Return the score matrix
    
    def table_transformation(self):
        '''Method transforming the score matrix into UPGMA distance matrix'''
        for value in self.score_table: #For value in matrix
            #using the formula to calculate the distance
            self.distance_table_UPGMA=np.min(self.score_table)-((self.score_table)-np.max(self.score_table))
        return self.distance_table_UPGMA
    
    def build_UPGMA(self):
        '''Method allowing to build UPGMA tree'''
        tree=arbre.UPGMA(self.distance_table_UPGMA,self.sequences) #Create the UPGMA object and store it in tree value
        self.tree_UPGMA=tree.upgma() #Create the UPGMA tree and pass it to the object
        return self.tree_UPGMA.newick() #Return the tree in Newick format with branch lengths

    def freq_matrix_calc(self,seq_list):
        '''Method allowing to create frequency profile'''
        if not isinstance(seq_list,list): #If the passed seq_list is not a list
            if isinstance(seq_list,str): #and if it is a string
                #Transform it to the list format
                seq_list_str=seq_list 
                seq_list=list()
                seq_list.append(seq_list_str)
            else:
                #Transform it to the list format
                seq_list=list(seq_list)
        num_seq = len(seq_list) #Number of sequences in list
        seq_len = len(seq_list[0]) #Length of sequences
        freq_matrix = np.zeros((21, seq_len)) #Create matrix with 21 line (for 20 amino acids and gap) and 
        #as many columns as there are character in sequence
        amino_acids = "CSTAGPDEQNHRKMILVWYF-" #list of amino acids and gap
        for i in range(num_seq):
            for j in range(seq_len) :
                indice_aa_matrice = amino_acids.find(seq_list[i][j]) #Find the index of amino acid in amino acid list
                freq_matrix[indice_aa_matrice][j] += 1 #Add one to the value of the corresponding cell
        freq_matrix /= num_seq #Devide each value of the matrix by number of sequences
        return freq_matrix 
    
    def nw_profiles(self,alignA,alignB):
        '''Method for profile-to-profile alignment using Needleman-Wunsch algorithm'''
        SCORE_GAP_vector=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]) #Initialising the gap constant vector
        dt = np.dtype([('float', np.float32), ('list', 'O')]) #Creating array supporting multiple values in a cell
        matrix_glob_align = np.zeros((len(alignA[0])+1,len(alignB[0])+1),dtype=dt) #Initialising the matrix 
        #Initialising the first row and the column
        for i in range(1, len(alignA[0])+1):
            matrix_glob_align[i][0][0] = matrix_glob_align[i-1][0][0] + np.dot(SCORE_GAP_vector, alignA[:,i-1])
            matrix_glob_align[i][0][1] = (i-1,0)
        for i in range(1, len(alignB[0])+1):
            matrix_glob_align[0][i][0] = matrix_glob_align[0][i-1][0] + np.dot(SCORE_GAP_vector, alignB[:,i-1])
            matrix_glob_align[0][i][1] = (0,i-1)
        for i in range(1, len(alignA[0])+1):
            for j in range(1, len(alignB[0]) + 1):
                #Using the formula of the score for profile-to-ptofile alignment
                score_1 = matrix_glob_align[i-1][j][0] + np.dot(SCORE_GAP_vector,alignA[:,i-1]) 
                score_2 = matrix_glob_align[i][j-1][0] + np.dot(SCORE_GAP_vector, alignB[:,j-1])
                score_3 = matrix_glob_align[i-1][j-1][0] + np.dot(alignA[:,i-1],alignB[:,j-1])
                #Performing the traceback
                if max(score_1 , score_2 , score_3) == score_1:
                    matrix_glob_align[i][j][0] = score_1
                    matrix_glob_align[i][j][1] = (i-1,j)
                if max(score_1 , score_2 , score_3) == score_2:
                    matrix_glob_align[i][j][0] = score_2
                    matrix_glob_align[i][j][1] = (i,j-1)
                if max(score_1 , score_2 , score_3) == score_3: 
                    matrix_glob_align[i][j][0] = score_3
                    matrix_glob_align[i][j][1] = (i-1,j-1)
        best_score_x = matrix_glob_align.shape[0]-1 #Length of the first alignment
        best_score_y = matrix_glob_align.shape[1]-1 #Length of the second alignment
        list_coordonate=[] #Creating an empty list to store the coordinates of the traceback
        while best_score_x>0 or best_score_y>0:
            variable= matrix_glob_align[best_score_x,best_score_y] #Storing the coordinates for a given position in a variable
            list_coordonate.append(variable[1]) # Adding the coordinates to the list
            best_score_x= variable[1][0] #Updating the index of first sequence according to coordinates
            best_score_y= variable[1][1] #Updating the index of second sequence according to coordinates
        list_coordonate= list(reversed(list_coordonate)) #Reversing the list
        return list_coordonate

    def insert_gaps(self,alignA,alignB,list_coordonate):
        '''Method for inserting the gaps in already aligned sequences according to their profile alignment'''
        if not isinstance(alignA,list): #If the passed alignA is not a list
            alignA=list(alignA) #Transform it to the list format
        if not isinstance(alignB,list): #If the passed alignB is not a list
            if isinstance(alignB,str):  #and if it is a string
                #Transform it to the list format
                alignB_str=alignB
                alignB=list()
                alignB.append(alignB_str)
            else:
                #Transform it to the list format
                alignB=list(alignB)
        all_alignement=[] #Creating an empty list to store all the alignments
        for sequence in alignA: 
            align_sequences = '' #Empty string to store the new sequence
            current_index=0 #Setting the current index to 0
            for index in range(len(list_coordonate)):
                if current_index==list_coordonate[index][0]: #If current index is the same as the current coordinate
                    if current_index==len(sequence): #And if the current index equals the length of the sequence
                        align_sequences+= "-" #Insert gap
                    else:
                        align_sequences+= str(sequence[current_index]) #Insert the corresponding character
                        current_index+=1 #Update the current index
                else: 
                    align_sequences+= "-" #Insert the gap
            all_alignement.append(align_sequences) #Add the new sequences to the alignment list
        for sequence in alignB:
            align_sequences = ''#Empty string to store the new sequence
            current_index=0 #Setting the current index to 0
            for index in range(len(list_coordonate)): 
                if current_index==list_coordonate[index][1]: #If current index is the same as the current coordinate
                    align_sequences+= str(sequence[current_index]) #Insert the corresponding character
                    current_index+=1#Update the current index
                else: 
                    align_sequences+= "-" #Insert the gap
            all_alignement.append(align_sequences)  #Add the new sequences to the alignment list
        return all_alignement
    
    def multiple_alignment(self,tree):
        '''Method allowing to perform the multiple sequence alignment'''
        sequences=self.sequences #Takes the dictionnary of sequences from the object
        if not tree.left.is_leaf(): #if the left branch of the tree is not a leaf
            if tree.right is not None: #If the right branch has a value
                if tree.left.left.is_terminal(): #And if the left branch of the left branch is terminal
                    print(str(tree.left))
                    self.alignment_names.append(str(tree.left)) #Append the left branches value to list of alignment names
                    print(str(tree.right))
                    self.alignment_names.append(str(tree.right)) #Append the right branches value to list of alignment names
                    if tree.right.right.is_terminal(): #And if the right branch of the right branch is terminal
                        alignment=self.nw(sequences[str(tree.left)],sequences[str(tree.right)]) #perform the classic Needleman-Wunsch alignment
                        frequency_profile=self.freq_matrix_calc(alignment) #Calculate the frequency profile
                        return alignment,frequency_profile #Return the alignment and frequency profile
                else:
                    left_alignment,left_frequency_profile=self.multiple_alignment(tree.left) #Recursicely apply the function to the left branch
                    right_alignment,right_frequency_profile=self.multiple_alignment(tree.right) #Recursicely apply the function to the right branch
                    aligned_freq=self.nw_profiles(left_frequency_profile,right_frequency_profile) #Perform the profile-to-profile alignment
                    self.final_alignment=self.insert_gaps(left_alignment,right_alignment,aligned_freq) #Insert gaps in the aligned sequences
                    final_frequency_profile=self.freq_matrix_calc(self.final_alignment) #Calculate the frequency matrix 
                    return self.final_alignment,final_frequency_profile #Return the alignment and frequency profile
        else: #The case of the solitary sequence on the branch
            alignement=sequences[str(tree.val)] #The alignment takes the value of the sequence itself
            frequency_profile= self.freq_matrix_calc(sequences[str(tree.val)]) #Calculating the frequency profile
            print(str(tree.val))
            self.alignment_names.append(str(tree.val)) #Adding the name to the alignment names list
            return alignement, frequency_profile #Return the alignment and frequency profile
        
    def get_conserved_indices(self,alignement):
        '''This function get the conserved positions'''
        frequences = {}  #Create an empty dictionnary 
        for position in range(len(alignement[0])): #
            aa_join = [] # create an empty list to stock the column of alignement
            for sequence in alignement: 
                aa = sequence[position] # retrieve the amino acid at the given position 
                aa_join.append(aa)  #add to the list 
            frequences[position]= aa_join #add the list to the position 
        for position, elements in frequences.items(): 
            counts = Counter(elements) #count the number of different AA at the given position in alignement
            if len(counts) == 2: #If only to different amino acids at the position (or amino acid and gap)
                self.indices.append(position) #Add to the list
        return self.indices

    

    def build_NJ_matrix(self,sites):
        '''Method calculating the NJ distance matrix based on conserved positions'''
        dict_alignement = {} #Creating a dictionnary
        for i in range(len(self.alignment_names)):
            dict_alignement[self.alignment_names[i]]=self.final_alignment[i] #Adding to the dictionnary the sequence name and corresponding alignment
        score=pd.DataFrame(index=list(dict_alignement.keys()),columns=list(dict_alignement.keys())) #Dataframe containing all the alignment names
        nbr_seq = len(dict_alignement.keys()) #Number of sequences in alignment
        #To create all the possible sequence pairs
        for i in range(nbr_seq):
            for j in range(i+1,nbr_seq):
                count = 0 #Set counter
                for k in sites: #For position in index list
                    #If the character in the position in both sequences is different
                    if dict_alignement[list(dict_alignement.keys())[i]][k]!=dict_alignement[list(dict_alignement.keys())[j]][k]: 
                        count = count + 1 #Add 1 to the counter
                score.iloc[i,j]=score.iloc[j,i] = count #Add the number of different positions (counter) to the corresponding cell
        for i in range(nbr_seq):
            score.iloc[i,i]=0 #Fill the diagonal with zeros
        self.distance_table_NJ=score.to_numpy() #Transforming into array
        return self.distance_table_NJ

    def build_NJ(self):
        '''Method allowing to build NJ tree'''
        tree=arbre.NJ(self.distance_table_NJ,self.alignment_names) #Create the NJ object and store it in tree value
        self.tree_NJ=tree.nj() #Create the NJ tree and pass it to the object
        return self.tree_NJ.newick_no_len()  #Return the tree in Newick format without branch lengths

