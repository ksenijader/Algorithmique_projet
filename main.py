# encoding: utf8
'''Code allowing to execute the multiple sequence alignment'''
import sys
from alignment import Alignment

class Controller():
    '''Class executing the analysis'''
    def __init__(self,file):
        '''Initialising the class'''
        self.alignment=Alignment(file) #Creating variable of Alignment class with the file passed by the user
        self.alignment.fasta_reader() #Reading the file

    def output_terminal(self):
        '''Method allowing to display the output in terminal'''
        print("Starting analysis...")
        score=self.alignment.score_counter()#Getting the score matrix of all pairwase alignments
        print(score)
        transformation_table=self.alignment.table_transformation() #Transforming the score matrix to distance matrix
        print(transformation_table)
        upgma_tree=self.alignment.build_UPGMA() #building the UPGMA guide tree
        print(upgma_tree)
        self.alignment.multiple_alignment(self.alignment.tree_UPGMA) #Performing the multiple sequence alignment
        for i in range(len(self.alignment.alignment_names)): #For each sequence ID in list of sligned sequences
            print(str(self.alignment.alignment_names[i])+"\n"+str(self.alignment.final_alignment[i])) #Print the sequence ID and corresponding sequence
        indices=self.alignment.get_conserved_indices(self.alignment.final_alignment) #Find the indexes of conserved positions
        print(indices)
        NJ_distance=self.alignment.build_NJ_matrix(self.alignment.indices) #Create a distance matrix from conserved indexes positions
        print(NJ_distance)
        tree_NJ=self.alignment.build_NJ() #Build the Neighbor-joining tree 
        print(tree_NJ)

    def output_file(self,output):
        '''Method for writing an output in a file specified by user'''
        print("Starting analysis...")
        score=self.alignment.score_counter()#Getting the score matrix of all pairwase alignments
        transformation_table=self.alignment.table_transformation()#Transforming the score matrix to distance matrix
        print("Distance transformed")
        upgma_tree=self.alignment.build_UPGMA()#building the UPGMA guide tree
        print("Guide tree created")
        self.alignment.multiple_alignment(self.alignment.tree_UPGMA)#Performing the multiple sequence alignment
        print("Multiple alignment done")
        indices=self.alignment.get_conserved_indices(self.alignment.final_alignment)#Find the indexes of conserved positions
        print("Conserved positions found")
        NJ_distance=self.alignment.build_NJ_matrix(self.alignment.indices)#Create a distance matrix from conserved indexes positions
        print("Neighbor-joining distance calculated")
        tree_NJ=self.alignment.build_NJ()#Build the Neighbor-joining tree
        print("Neighbor-joining tree created")
        with open(output,'w') as output_file: #Open the output file in writing mode as output_file
            output_file.write("NW score matrix"+"\n")
            output_file.write(str(score)+"\n"+"\n") #Write the scoring matrix in the file
            output_file.write("UPGMA distance matrix"+"\n")
            output_file.write(str(transformation_table)+"\n"+"\n") #Write the transformed matrix in the file
            output_file.write("UPGMA tree in newick format:"+"\n")
            output_file.write(upgma_tree+"\n"+"\n") #Write the UPGMA tree in file
            output_file.write("Multiple sequence alignment:"+"\n")
            for i in range(len(self.alignment.alignment_names)): #For each sequence ID in list of sligned sequences
                #Write the sequence ID and corresponding sequence in file
                output_file.write(str(self.alignment.alignment_names[i])+"\n"+str(self.alignment.final_alignment[i])+"\n")
            output_file.write("\n"+"\n")
            output_file.write("Indexes of conserved positions:"+"\n")
            output_file.write(str(indices)+"\n"+"\n") #Write the indexes of conserved position in file
            output_file.write("Distance matrix from conserved positions:"+"\n")
            output_file.write(str(NJ_distance)+"\n"+"\n") #Write the distance matrix from conserved positions in file
            output_file.write("NJ tree in newick format:"+"\n")
            output_file.write(tree_NJ) #Write the NJ tree in file

argument=sys.argv #Creating a list to store the arguments passed by user

if __name__ == "__main__":
    # Input takes the value given by user after "-input"
    input_file=argument[argument.index("-input")+1]
    #Create the object of Controller class with the given imput
    C=Controller(input_file)
    if "-output" in argument: #If the argument "-output" was passed by the user
        # Output takes the value given by user after "-output"
        file_output=argument[argument.index("-output")+1]
        #Execute the method that writes the outpus in a file provided by the user
        C.output_file(file_output)
    else:
        #Execute the method diplaying the output in the terminal
        C.output_terminal()
