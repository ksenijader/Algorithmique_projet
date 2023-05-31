# encoding: utf8
import numpy as np
class BinaryTree:
    '''Class allowing to build binary tree'''
    def __init__(self,taxa):
        '''Initialising the Binary tree'''
        if len(taxa) == 0: #If length of taxa is 0
            self.val = None #Value is None
            self.nLeaf = 0 #Number of leafs is 0
            self.depth = 0 #Depth of tree is 0
            self.branch_length = 0 #Branch length is 0
            self.left=None #Left branch is None
            self.right=None #Right branch is None
        else:
            self.val = taxa[0] #Value is the taxa name
            self.branch_length = 0 #Branch length is 0
            self.left = BinaryTree(taxa[1]) #Value of left branch is second taxa
            self.right = BinaryTree(taxa[2]) #Value of right branch is third taxa
            if self.left.is_terminal() and self.right.is_terminal(): #If both left and right branches are terminal
                self.nLeaf = 1 #Number of leafes is 1
                self.depth = 0 #Depth of node is 0
            else:
                self.nLeaf = self.left.nLeaf + self.right.nLeaf #Number of leafs is the sum of leafes of left and right branches
                self.depth = 1 + max(self.left.depth,self.right.depth) #Depth is maximal depth of left and rigth branches + 1

    def is_terminal(self):
        '''Method checking if node is terminal'''
        return (self.val is None) #Returns True if self.val is None and false if value exists
    
    def is_leaf(self):
        '''Method checking if node is leaf'''
        return (self.left is None and self.right is None) #Returns True if self.left and self.right are None and false if values exist

    def join(self,rightTree):
        '''Method allowing to conect nodes of the tree'''
        tmp = BinaryTree([]) #The new node is also binary tree object
        tmp.val = "" #Value is an empty string
        tmp.branch_length = 0 #The branch length is 0
        tmp.left = self #the value of left branch is self value
        tmp.right = rightTree #The value of right branch is the right tree that we are joining to self
        tmp.nLeaf = self.nLeaf + rightTree.nLeaf #The number of leafes id the sum of right tree leafes and self leafes
        tmp.depth = 0 #Depth is 0
        return tmp 
    
    def newick(self):
        '''Method displaying tree in Newick format with branch lengths'''
        if self.left.is_terminal() and self.right.is_terminal(): #If left and right branches are terminal
            return f"{self.val}:{self.branch_length:.2f}" #Return the value of self (taxa) and branch length
        return f'({self.left.newick()},{self.right.newick()}):{self.branch_length:.2f}' #Apply the function recursively on the left and right branches
    def newick_no_len(self):
        '''Method displaying tree in Newick format'''
        if self.left.is_terminal() and self.right.is_terminal():#If left and right branches are terminal
            return f"{self.val}" #Return the value of self (taxa)
        return f'({self.left.newick_no_len()},{self.right.newick_no_len()})'#Apply the function recursively on the left and right branches

    def __str__(self):
        '''Initialising character string of the tree'''
        if self.is_terminal(): #If self is terminal
            return '' #Return an empty string
        else: 
            string_val = str(self.val) #The value of the character string is the taxa
            if self.left is not None: #if left branch exists
                string_val = string_val + str(self.left) #add the taxa of left branch to the current character string
            if self.right is not None:#if right branch exists
                string_val = string_val + str(self.right) #add the taxa of right branch to the current character string
            return string_val

class UPGMA:
    '''Class allowing to construt UPGMA tree'''
    def __init__(self,distance_matrix,taxa):
        '''Initisialising the class'''
        self.distance_matrix=distance_matrix 
        self.taxa=taxa

    def find_upgma(self,distance):
        '''Method allowing to find the minimal distance in matrix'''
        n = len(distance) #Number of taxa
        mindist = distance[0,1] #Initial minimal distance is the first cell of the matrix
        min_i,min_j = 0,1 #the indexes of the taxa corresponding to minimal distance
        for i in range(n-1): 
            for j in range(i+1,n):
                if distance[i,j] < mindist: #If the current distance is smaller than mindist value
                    mindist,min_i,min_j = distance[i,j],i,j #Update the current values
        return mindist,min_i,min_j 
    
    def update_upgma(self,distance,i,j,size_i,size_j):
        '''Method updating the distance matrix'''
        n = len(distance) #Number of taxa
        new = np.zeros((n+1,n+1)) #New matrix with one additional row and column 
        new[1:,1:] = distance #Fill the new matrix with old matrix value starting from the first row/column
        for m in range(1,n+1):
            new[0,m] = (size_i*new[i+1,m] + size_j*new[j+1,m])/(size_i+size_j) #Applying the UPGMA formula to update the values
            new[m,0]=new[0,m]
        new = np.delete(new,[i+1,j+1],axis=0) #Deleting unique values of i and j
        new = np.delete(new,[i+1,j+1],axis=1)#Deleting unique values of i and j
        return new
    
    def upgma(self):
        '''Method that creates the UPGMA tree'''
        distance=self.distance_matrix #Takes the passed matrix as the distance matrix
        names=self.taxa #Takes the passed taxa as names
        clades = [BinaryTree([name,[],[]]) for name in names] #Creates Binary tree for each name in names
        n = len(names)  #Number of taxa
        for c in range(n-1):
            mindist,i,j = self.find_upgma(distance) #Finds the minimal distance with find_upgma method
            clades[i].branch_length = mindist/2 - clades[i].depth #Applying UPGMA formula to compute branch lengths
            clades[j].branch_length = mindist/2 - clades[j].depth
            u = clades[i].join(clades[j]) #Joining the taxa in a node
            size_i=clades[i].nLeaf #Size of the i branch
            size_j=clades[j].nLeaf #Size of the j branch
            u.depth = mindist/2 #Depth of the node
            u.left.branch_length = u.depth - u.left.depth #Length of the left branch
            u.right.branch_length = u.depth - u.right.depth  #Length of the right branch
            distance = self.update_upgma(distance,i,j,size_i,size_j) #Updating the distance matrix
            clades = [u] + clades[:i]+clades[i+1:j]+clades[j+1:] #Adding created node to the clade
        return clades[0]
    
class NJ:
    '''Class allowing to build Neighbor-Joining tree'''
    def __init__(self,distance_matrix,taxa):
        '''Initisialising the class'''
        self.distance_matrix=distance_matrix
        self.taxa=taxa

    def find_nj(self,distance):
        '''Method allowing to find the minimal distance in matrix'''
        n = len(distance)
        mindist = distance[0,1]
        mini,minj = 0,1
        q=distance[0,1] #Initialisng the q value
        for i in range(n-1):
            for j in range(i+1,n):
                q_now=(n-2)*distance[i,j]-np.sum(distance[i,:])-np.sum(distance[j,:]) #Calculating the q value with NJ formula
                if q_now < q: #If current q is smaller than q value
                    mindist,mini,minj,q = distance[i,j],i,j,q_now #Update the values
        return mindist,mini,minj

    def update_nj(self,distance,i,j,mindist):
        '''Method updating the distance matrix'''
        n = len(distance)
        new = np.zeros((n+1,n+1))
        new[1:,1:] = distance
        for m in range(1, n + 1):
            new[0, m] = ((new[i + 1, m] + new[j + 1, m]) - mindist)/2 #Applying the NJ formaula to calculate the updated distances
            new[m,0]=new[0,m]
        new = np.delete(new,[i+1,j+1],axis=0)
        new = np.delete(new,[i+1,j+1],axis=1)
        return new

    def nj(self):
        '''Method that creates the NJ tree'''
        distance=self.distance_matrix
        names=self.taxa
        clades = [BinaryTree([name,[],[]]) for name in names]
        n = len(distance)
        for c in range(n-1):
            mindist,i,j = self.find_nj(distance)
            sum_i_k=np.sum(distance[i,:]) #Calculating the sum of distance of taxa i and the rest of the taxa
            sum_j_k=np.sum(distance[j,:]) #Calculating the sum of distance of taxa j and the rest of the taxa
            clades[i].branch_length = mindist/2 -(sum_i_k-sum_j_k)/(2*(n-2)) #Applying NJ formula to calculate branch lengths
            clades[j].branch_length = mindist-clades[i].branch_length
            u = clades[i].join(clades[j])
            u.depth = mindist/2
            u.left.branch_length = u.depth - u.left.depth
            u.right.branch_length = u.depth - u.right.depth
            distance = self.update_nj(distance,i,j,mindist)
            clades = [u] + clades[:i]+clades[i+1:j]+clades[j+1:]
        return clades[0]

    
