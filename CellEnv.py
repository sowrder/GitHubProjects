import simpy as sim
import os
import string as string
import numpy as num
from numpy import array as tensor #deal with multi dimensional values
import ChroDatabase_Search as s
from pathlib import Path 
from rdkit import Chem 
from rdkit.Chem import Draw 
import re 
import Proarray as pro 

os.chdir('/Users/ismailnashwan/Desktop/GeneDatabase_Project/GeneSeq/Temp_Seq')

#assume a concentration of basic functional proteins already exists, but requires maintenance 

#chatGPT provided functions for randomness and probabilistic regulation  

"""Weighted/randomizer function 
    Selects a random item from input_list based on the probabilities defined in prob_dict.

    Arguments:
    input_list: A list of either integers, floats (for ranges), or letters (for selection).
    prob_dict: A dictionary that defines the probabilities for items in input_list. 
               The dictionary should contain the keys corresponding to items in input_list, 
               with their respective probabilities as values.
               The dictionary may also contain a 'depend' key (True/False) to control probability dependence.

    Returns:
    A randomly selected item from input_list based on the given probabilities.
    """

def randomizer(input_list):
    # Check if the input is a list and contains the appropriate types
    if isinstance(input_list, list):
        # Case 1: List of 2 integers (range)
        if len(input_list) == 2 and all(isinstance(i, int) for i in input_list):
            min_val, max_val = sorted(input_list)
            random_int = int.from_bytes(os.urandom(4), byteorder='big') % (max_val - min_val + 1)
            return min_val + random_int
        
        # Case 2: List of capital letters
        elif all(isinstance(letter, str) and letter.isupper() for letter in input_list):
            random_index = os.urandom(1)[0] % len(input_list)
            return input_list[random_index]
        
        # Case 3: List of 2 floats (range)
        elif len(input_list) == 2 and all(isinstance(i, float) for i in input_list):
            min_val, max_val = sorted(input_list)
            random_int = int.from_bytes(os.urandom(8), byteorder='big')
            max_int = 2**64
            random_float = min_val + (random_int / max_int) * (max_val - min_val)
            return random_float
        
        else:
            raise ValueError("Input list is not in the correct format.")
    
    else:
        raise ValueError("Input must be a list.")

gene_type = ['v_low', 'low', 'mod', 'upreg', 'h_upreg', 'high_exp' ] #established according to known literature promoters of transciription 

pairs_DNA = {'A':'T',  'C':'G', 'U':'A', 'G':'C'} #Forming DNA units 
pairs_RNA = {'A': 'U', 'T': 'A', 'C': 'G', 'G': 'C'} #Forming RNA units 

def read(txt):
    bases_path = '/Users/ismailnashwan/Desktop/GeneDatabase_Project/GeneSeq/Temp_Seq'
    path = Path(f'{bases_path}/{txt}')
    return path.read_text() 


def transcribe(DNA_str):
    RNA_str = '' 
    # Iterate over each character in the DNA string
    for base in DNA_str.strip():  
        for key, pair in pairs_RNA.items(): #complimentary RNA formation 
            if base == key: #since our dictionaries are 1-1 
                RNA_str += pair  # Add the corresponding RNA base to the assignments string
    return RNA_str  # Return an RNA string after completing 

#get index from base position int
def mod(int):
    return int//6

#get base position int from index 
def revmod(int):
    return int*6

def revtranscribe(RNA_str):
    DNA_str = ''
    for base in RNA_str.strip(): 
       for key, pair in pairs_DNA.items(): #complimentary DNA formation 
            if base == key: #since our dictionaries are 1-1 
                DNA_str += pair
    return DNA_str  


amino_acids = {
     "F": ["UUU", "UUC"],  # Phenylalanine
    "L": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],  # Leucine
    "S": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],  # Serine
    "Y": ["UAU", "UAC"],  # Tyrosine
    "C": ["UGU", "UGC"],  # Cysteine
    "W": ["UGG"],  # Tryptophan
    "P": ["CCU", "CCC", "CCA", "CCG"],  # Proline
    "H": ["CAU", "CAC"],  # Histidine
    "Q": ["CAA", "CAG"],  # Glutamine
    "R": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],  # Arginine
    "I": ["AUU", "AUC", "AUA"],  # Isoleucine
    "M": ["AUG"],  # Methionine (start codon)
    "T": ["ACU", "ACC", "ACA", "ACG"],  # Threonine
    "N": ["AAU", "AAC"],  # Asparagine
    "K": ["AAA", "AAG"],  # Lysine 
    "V": ["GUU", "GUC", "GUA", "GUG"],  # Valine
    "A": ["GCU", "GCC", "GCA", "GCG"],  # Alanine
    "D": ["GAU", "GAC"],  # Aspartic acid
    "E": ["GAA", "GAG"],  # Glutamic acid
    "G": ["GGU", "GGC", "GGA", "GGG"],  # Glycine
    "M": ["AUG"],  # Start codon
    "stop": ["UAA", "UAG", "UGA"]  # Stop codons
}

#append function to check for mRNA attribute 


def insert_string_at_(original_str, index, element):
        return original_str[:index] + element + original_str[index:]
 # Pad rows with zeros



#used by ribosome function on processed RNA, mRNA 
def translate(mRNA_str, show=False): 
    Protein_str = ''
    START = 'AUG'
    STOP = False
    initiate = RNA_str.find(START)  # Find first occurrence of "AUG" 
    if initiate == -1: #returned if no AUG found
        return "Start codon (AUG) not found"
    RNA_str = insert_string_at_(RNA_str, initiate, '_')  # Add separator item
    separator = RNA_str.split('_')  # Split the string at '_'
    codable = separator[1]  # Get the part of the string after 'AUG'
    for base in range(0, len(codable), 3):  # Iterate in steps of 3
        matching_bases = codable[base:base+3]  # Get a matching codon
        if matching_bases in amino_acids['stop']:  # If it's a stop codon
            STOP = True
            break
        for amino_acid in amino_acids:
            if not show:
                if matching_bases in amino_acids[amino_acid]:  # Match codon to amino acid
                    Protein_str += amino_acid  # Add amino acid to protein string
            elif show:
                if matching_bases in amino_acids[amino_acid]:  # Match codon to amino acid
                    Protein_str += '-' + amino_acid  # Add amino acid to protein string
    if show:
        Protein_str += '-Stop'
    if STOP:
        print(f'Stop codon {matching_bases} found at iteration {base}')
        print(f'Translation complete, final string is \n{Protein_str}')
        return Protein_str
    



'''The Gene object
    The Gene object is a buildable object which is initially modified into a gene DNA sequence instead of the FASTA NCBI amino acid and UTR sequence.
    This is neccessary because a protein has several steps of processing which are ultimately key in it's existence. 

    Gene Object : is passed a string, a corresponding type representing upstream regulation, 
    and a dictionary representing degradability and a value within the corresponding range for the key set as seen above.

    seq : extracted via NCBI database using ChroDatabase_search.py alongside Ensamble  

    type : determines a modulation to the string sequence later appended into promoter upstream string objects 
    expression : uses type strings as keys to represent degradability of the gene sequence once transcribed 

    if a gene is highly expressed, it's either not appended or appended very slightly (for now, we'll avoid apended regulation)
'''



class Gene:
    def __init__(self, seq, type):
        if not isinstance(seq, str):
            raise ValueError('Please Enter a proper template sequence, or minus sequence according to NCBI')
        #self seq is defined with a promoter string and some number of emoty strings, this is set arbitrarily for simplicity 
        self._seq =  seq  #sequence from sequences file containing main promoter
        self._upstream = ''  #appendable upstream attribute based on type and other conditions 
        #self._expression = expression  #adjustable value reflecting relative access to translation into protein class 
        
        # Setting the upstream sequence based on the type, this will likely be implemented later 
        if type == 'v_low':
            self._upstream = '_' * 15000
        elif type == 'low':
            self._upstream = '_' * 7000
        elif type == 'mod':
            self._upstream = '_' * 5000
        elif type == 'upreg':
            self._upstream = '_' * 3000
        elif type == 'h_upreg':
            self._upstream = '_' * 2000
        elif type == 'high_exp':
            self._upstream = '_' * 1000 
        else:
            raise ValueError(f"Invalid gene type: {type}")

        self.seq = self._upstream + self._seq  # Combine upstream with seq

        #self.expression = expression if expression in expression_pairs else ValueError #may use random function to generate values...
    
    #if the gene is highly regulated through an upstream process, include other transcirption factors that regulate RNA pol 
    def Promoter(self): #define the promoter method to append sequences while building the gene object 
        pass
    
    
    @property
    def upstream(self):
        return self._upstream

    @upstream.setter
    def upstream(self, value):
        self._upstream = value
        # Every time upstream changes, update the full sequence
        self._seq = self._upstream + self._seq  # Update the seq with new upstream

    @property
    def seq(self):
        return self._seq

    @seq.setter
    def seq(self, value):
        self._seq = value



#transcription function converts any string according to desired RNA need 
        



'''The Protein class 
    The protein class defines a class with several attributes according to the assigned function 

    the protien object gains an attribute self.active which grants it access to attributes of the class as desired 
    this attribute is gained conditionally, usually through chaperoneing or some other folding inducing event. 

    Implementation of protiens - 
    The protien sequence is passed as a string which creates a padded raw array for the protien containing the alpha helix, beta sheet, and string groups 
    the fold method is passed onto the protien to initalize the inital domain according to multiple dictionaries 

    the first dictionary is through an index function in combination with unitProt (extremely nice and helpful)
    Unit pro provides the sequences (which will be generated from ribosomes in this case) as well as the domain active regions 
    these regions are defined by passing the region index to domain_index, this is for implementation purposes and defines which items 
    in the string comprise the protien active domains, which is quite useful. 



    defining states - 
    The second dictionary contains implemented information regarding the protien, specifically how to create an active protien object
    array from the defined components through domain_index.

    the protien name, such as CREB or Helicase, defines the access to the implemented states. 
    the values are a list of dictionaries containing states T0, T1, and T2

    the T keys are dictionaries coontaining the postiional information of states, where the different 
    states reflect different changes to the inital protien structure in 1T0, or 2T0 if the protien has multiple 
    resting states.

    T0 is defined as the resting state, it's a key which  has a list of a boolean variable at index 0 
    this defines the 0th state  gained through fold which initalizes the protien domain 

    the first dictionary argument, that is 1T0 or index 1,is the proper domain structure of the alpha 
    helicies, strings, and domains, based on their index along the assigned domains, and is a directionally
    orientated array that shows the directions of the different domain units 
    the direction appension is defined as follows, alongside an integer in a tuple which represents the unit position/s 
    
    |, -| : vert up or reverse sequence vert 
    ->,<- : orient units horizontally or in reverse horizontally 
    /, -/ : diagonal up right or down left 
    \\, -\\ : diagonal up left or downright 
    (*n, *n) : define counter and clockwise rotation with n, where n = 1 is a 180 turn in array according 
    to orientation, if aligned vertically, then unit is appended to a position in the reverse coordinate, with 
    it's direction being reflected as such 
    = : anti parallel alignment using second unit (evens beta)
    =- : anti parallel alignment using first unit  (odds beta)
    stack - use units to make a 2d array seperated by | 


    base dictionary format : 

    {... ,Protien_name: [{'T0': [False, {'1T0': 
    {a:[()] , b:[()] , s:[()]} }]}],

    {'T1': {'0T1': {a:[()] ,b:[()] ,s:[()]}},

    {'T2': {'0T2': {a: [()], b:[()], s:[()]}}}

    }}


'''




#defining TF2D through F 

TBP = 'MDQNNSLPPYAQGLASPQGAMTPGIPIFSPMMPYGTGLTPQPIQNTNSLSILEEQQRQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQAVAAAAVQQSTSQQATQGTSGQAPQLFHSQTLTTAPLPGTTPLYPSPMTPMTPITPATPASESSGIVPQLQNIVSTVNLGCKLDLKTIALRARNAEYNPKRFAAVIMRIREPRTTALIFSSGKMVCTGAKSEEQSRLAARKYARVVQKLGFPAKFLDFKIQNMVGSCDVKFPIRLEGLVLTHQQFSSYEPELFPGLIYRMIKPRIVLLIFVSGKVVLTGAKVRAEIYEAFENIYPILKGFRKTT'

TBP_ar = pro.pro_array(TBP)
TBP_ar = tensor(TBP_ar, dtype=object)

pro.filter(TBP_ar)

pro.domain_index(TBP_ar, 158, 337)

Protein_functions = ['chapp', 'kinase', 'phosphtase', 'receptor', 'ligand', 'transporter', 'polymerase', 'protease', 'DNA_repair', 'TF']
Enzyme = ['generator', 'slicer', 'processor']
TF = ['core', 'primary', 'secondary', 'tertiary']  #these ranks are somewhat arbitrarily decided depending on need 


class Protein:
   
    def __init__(self, units_aa, protien):
        self.parr = pro_array(units_aa)  #create raw array of protein using beta and alpha helicies 
        self.states = Tstates[f'{protien}'] #find the protien key in T states for implemented protiens 
        if self.states:
            pass
        else: 
            return NotImplemented 
        
    def fold(self):Mn'h[ig]

        if self.states[0]['T0'][0]: #check if zeroth variable is True 

            canon_arr =  tensor([[].[].[]]) #initalize x,y,z coordinates (shape 1,1, 3 )

            for aaseq, idx in pro.filter(self.parr): #first iterate over array tuple values 

                for struct, condition in enumerate(p åﬂ3.states[0]['T9'][1]['0T0']): #then iterate over canon state key 

                    for ui, rot in condition: #and finally iterate through the list values for each key 
                        range = list(ui):
                        domain_str = ''
                        domain_item = []
                        if idx in range and 'c' in rot: #zy plane, zcos, ysin 
                            domain_str += ''.join(aaseq)
                            deg = int(rot.split('c')[-1]) 
                            
                        r = len(domain_str)

                        elif  idx in range and '/' in rot: #xy plane, xcos, ysin 
                            pass
                        elif idx in range '_' in rot: #xz plane, xcos, zsin 
                            pass



        else: 
            return 'This protien is not ready for folding '

            
for aaseq, idx in pro.filter(TBP_ar):
    print(idx) 

class Psubunit(Protein):
    pass

        

expression_values = [list(range(start, end-1, -1)) for start, end in [(180, 151), (150, 121), (120, 91), (90, 61), (60, 31), (30, 0)]]
expression_pairs = [(key,exp) for key,exp in zip(gene_type, expression_values)] #list of tupules, keys in this case define degradabo;oty levels  

#RNA is processed in several stages for each type, this class is generated to allow that to occur

RNA_type = ['ribounit_n', 'splicounit_n', 'mRNA', 'tRNA', 'silencer'] 

'''The RNA class 
    The RNA class is comprised of strings that are transcribed from genes, or Diced in the case of small RNA, in general it's connected
    to the gene class however it is an output from it's processing and so it's attributes change, it's expression argument also controls it's existence 

    seq arg : the sequence used to generate the RNA class object 

    expression arg : a value containing information that establishes regulation parameters of the RNA class object after creation 

    type : a string describing the type of RNA, mRNA, rRNA, t-RNA, silencer RNA, etc.. each with their distinct attributes 

    stage : an initalization argument which describes which stage of processing the generated RNA class being generated is at 

'''

class RNA:
    def __init__(self, seq, expression, type, stage): 
        RNA = [seq[i] for i in range(len(seq))] 
        self.seq = tensor(RNA, dtype=str) #RNA after being passed by the transcribe function or other sources, this sequence contains U instead of T 
        self.type = type 
        self.stage = stage 
        self.expression = expression 
        if self.type == 'ribounit_n':
            pass


defined_types = ['NeuronM', 'NeuronS', 'Neuron_V1', 'Neuron_V2', 'microglia', 'astrocyte', 'olig_cyte', 'schwan'] #modelling nervous system, visual V1,V2 and motor peripheral cells 

defined_states = ["progenitor", "precurssor", "immature_precurssor", "mature_cell"] #arbitrarily set based on Dr.Cody Smith's class, also generally true 



class Cell: #when initalized, argument attributes determine access to methods, and are restricted based on input
    
    def __init__(self, type, state):
        self.type = type if type in defined_types else NotImplemented #define cell final state, neuron, liver, heartcell, etc.. 
        self.state = state if state in defined_states else NotImplemented #define current developmental state 
    #for each defined cell type, extract proper relevant sequence information. there are several, so focus is climited
        if self.type == 'NeuronM':
            pass
        if self.type == 'microglia':
            pass
        if self.type == 'astrocyte':
            pass
        if self.type == 'olig_cyte':
            pass
        if self.type == 'schwan':
            pass
    def transcription(self, gene):
        pass



