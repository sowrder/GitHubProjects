import simpy as sim
import os
import string as string
import numpy as num
from numpy import array as tensor #deal with multi dimensional values
import ChroDatabase_Search as s
import re 
import Proarray as pro 



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

    the protiens complicated structure is linearized through the combination of imaginary rotation and quarternions, where any 
    specefic sequence has with it 3 specific xyz coordinates that reflect the position of it in the protien array space 
    this space must be initalized in the exact way defined by the dictionary, conditions will be employed to regulate how this may happen. 

    In general, there are 3 main arguments that designate the plane for the class object 

    '_' which defines rotation along the k axis (perpendicular to xy)
    '/' defining rotation along the j axis (perppendicular to xz)
    'c' defining rotation along the i axis (perpendicular to yz)

    for any single rotation in space, then the string and corresponding n degrees are passed depending on dictionary domains 
    i.e '/45' , since this is a single rotation along one plane, then we simply use i and break down components with 
    rot_img(s1, n), then we sum the real component (in this case) to x, and the complex component to y, with the previous coordinates summed. 
    defining sqrt(2)/2 as a, our point would look like [[str, ua+x],[str,ua+y],[str, z]]
    notably, the points themselves represent the position, wheras drawing lines between the dots will reflect the sequence length between. 


    if '/45+_45', then we use rot_imag values as the arguments for rot_quart, since /45 assigns 
    x and y components, we also include a 0 z component that can be used to define the vector quarternion. 

    for every + passed, we pass the rotated vectors list, defining a condition based on the string used.

    To define a surface in plane or rotated, we do a very similar thing. imag_rot takes care of defining the direction of the surface 
    with the dir argument, where we pass the length of a unit L (based on some integer prioir to #) as s2 and the count of L as s1 
    these correspond to the horizontal and vertical component respectively, with a third reversed vertical component. 

    After all rotations are complete for a surface (that is, after returning the list of 9 coordinates)
    the vectors are summed in order with each other, where we define 
    p4 = the previous point coordinate 
    p1 = w0 + p4
    p2 = w1 + p1
    p3 = w3 + p2 

    after which the coordinates are returned in the order of 0-3 from 1-4 as a 4 item list 

    base dictionary format : 

    {... ,Protien_name: [{'T0': [False, {'1T0': 
    {a:[()] , b:[()] , s:[()]} }]}],

    {'T1': {'0T1': {a:[()] ,b:[()] ,s:[()]}},

    {'T2': {'0T2': {a: [()], b:[()], s:[()]}}}

    }}
'''

Tstates ={'TBP': {'T0':[False, 
    {'a':[(range(0,7),'_45'), (range(7,12),'c25'), (range(20,25), '/45+c25+_45')] ,
    'b':[(range(12,18), '4#!_45'), (range(30,41), '5#/45+_45')] ,
    's':[(range(19,20), '_45+/45')]}]}}


#defining TF2D through F 

TBP = 'MDQNNSLPPYAQGLASPQGAMTPGIPIFSPMMPYGTGLTPQPIQNTNSLSILEEQQRQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQAVAAAAVQQSTSQQATQGTSGQAPQLFHSQTLTTAPLPGTTPLYPSPMTPMTPITPATPASESSGIVPQLQNIVSTVNLGCKLDLKTIALRARNAEYNPKRFAAVIMRIREPRTTALIFSSGKMVCTGAKSEEQSRLAARKYARVVQKLGFPAKFLDFKIQNMVGSCDVKFPIRLEGLVLTHQQFSSYEPELFPGLIYRMIKPRIVLLIFVSGKVVLTGAKVRAEIYEAFENIYPILKGFRKTT'

TBP_ar = pro.pro_array(TBP)
test = pro.filter(TBP_ar)

pro.domain_index(TBP_ar, 158, 337)

revmod(33) 



Protein_functions = ['chapp', 'kinase', 'phosphtase', 'receptor', 'ligand', 'transporter', 'polymerase', 'protease', 'DNA_repair', 'TF']
Enzyme = ['generator', 'slicer', 'processor']
TF = ['core', 'primary', 'secondary', 'tertiary']  #these ranks are somewhat arbitrarily decided depending on need 


#string combiner for especially large strings 

def comb_strings(array, indices, target_length):
    concatenated_strings = []
    current_string = ""
    
    # Decode the byte strings in the 'seq' field
    selected_strings = [str(array[i]['seq'], 'utf-8') for i in indices if i < len(array)]  # Ensure index is valid
    
    for string in selected_strings:
        current_string += string
        
        # When the length of the current concatenated string reaches or exceeds the target length
        if len(current_string) >= target_length:
            # If the string is longer than required, truncate it to the target length
            concatenated_strings.append(current_string[:target_length])
            current_string = current_string[target_length:]  # Remaining part will be used for next iteration
    
    # If there's a remaining part in the current_string, it's discarded since it's incomplete
    return concatenated_strings



class Protein:
   
    def __init__(self, units_aa, protien):
        self.parr = pro.pro_array(units_aa)  #create raw array of protein using beta and alpha helicies 
        self.states = Tstates[f'{protien}'] #find the protien key in T states for implemented protiens 
        self.canon = []
        if self.states:
            pass
        else: 
            return NotImplemented 
        
#Later, this may implement actual physics by making each string item an object with different attributes that are 
#reflected through folding, turnining in that case would be calculated based on molecular interactions of sequences, too much work rn. 
    def fold(self):
        
        # Initialize a list to hold x, y, z coordinates
        canon_arr = [[], [], []]  # Holds the x, y, z coordinates
        Surface_arrays = [] #list of generated surfaces that contain the string sequences, seperate from the array 
        self.filtered = pro.filter(self.parr)
        # Check if the zeroth variable in states is True
        if self.states[0]['T0'][0]:
            # Iterate over array tuple values
            for aaseq, idx in pro.filter(self.parr):
                # Iterate over canon state key
                for struct, condition in enumerate(self.states[0]['T0'][1]['0T0']):
                    # Iterate through the list values for each key
                    for ui, rot in condition:
                        if '#' not in rot: #if the dictionary structure is simply a line (beta, helix, or string sequence) 
            
                            ui_range = list(ui)  # Convert `ui` to a list
                            matches = [seq.decode('utf-8') for seq, idx in self.parray if idx in ui_range] #get matching indecies from array in order 
                            single_unit= ''.join(matches)  # Combine all `aaseq` values in the range
                            uL = len(single_unit) 
                        
                            if '_' in rot and '+' not in rot: #if only a single rotation along xy 
                                deg = rot.split('_')[-1] 
                                vec = pro.rot_imag(uL, deg) #rotate scaled length 
                                if not canon_arr: #if this is empty 
                                    canon_arr[0, 0].append((single_unit,vec[0])) #x is real 
                                    canon_arr[0, 1].append((single_unit,vec[1]))  #y is imaginary 
                                    canon_arr[0, 2].append((single_unit,0)) #z is zero (k axis rotation) 
                                else: #if values exist 
                                    canon_arr[0, 0].append((single_unit, vec[0]+canon_arr[0, 0, -1])) 
                                    anon_arr[0, 2].append((single_unit, vec[1]+canon_arr[0, 1, -1])) 
                                    anon_arr[0, 1].append((single_unit, canon_arr[0, 0, -1])) #add previous vector values 
                        
                            if '/' in rot and '+' not in rot: #along j axis or xz plane 
                                deg = rot.split('/')[-1] 
                                vec = pro.rot_imag(uL, deg) 
                                if not canon_arr: #if this is empty 
                                    canon_arr[0, 0].append((single_unit,vec[0])) #x is real 
                                    canon_arr[0, 2].append((single_unit,vec[1])) #z is imaginary 
                                    canon_arr[0, 1].append((single_unit,0)) #y is zero 
                                else: #if values exist 
                                    canon_arr[0, 0].append((single_unit, vec[0]+canon_arr[0, 0, -1]))  
                                    anon_arr[0, 2].append((single_unit, vec[1]+canon_arr[0, 1, -1])) 
                                    anon_arr[0, 1].append((single_unit, canon_arr[0, 0, -1]))
                        
                            if 'c' in rot and '+' not in rot: #along the i axis or yz plane 
                                deg = rot.split('_')[-1] 
                                vec = pro.rot_imag(uL, deg) #rotate scaled length 
                                if not canon_arr: #if this is empty 
                                    canon_arr[0, 1].append((single_unit,vec[0])) #y is real 
                                    canon_arr[0, 2].append((single_unit,vec[1])) #z is imaginary 
                                    canon_arr[0, 0].append((single_unit,0)) #x is zero 
                                else: #if values exist 
                                    canon_arr[0, 1].append((single_unit, vec[0]+canon_arr[0, 0, -1])) 
                                    anon_arr[0, 2].append((single_unit, vec[1]+canon_arr[0, 1, -1])) 
                                    anon_arr[0, 0].append((single_unit, canon_arr[0, 0, -1])) 

                            elif '+' in rot: #if we have a string seperated by + signs meaning multiple rotations 
                                rotations = rot.split('+') #create a list of possible rotations 
                                plane = ''
                                if '_' in rotations[0]:
                                    deg = int(rotations[0],split('_')[-1])
                                    vec = pro.rot_imag(uL, deg)
                                    wi = [(vec[0], vec[1], 0)] #set xy values
                                if '/' in rotations[0]:
                                    deg = rotations[0],split('/')[-1]
                                    vec = pro.rot_imag(uL, deg)
                                    #the result is always tuple of 2 indecies 0 and 1  
                                    wi = [(vec[0], 0, vec[1])] #set xz values
                                if 'c' in rotations[0]:
                                    deg = rotations[0],split('c')[-1]
                                    vec = pro.rot_imag(uL, deg)
                                    wi = [(0, vec[0], vec[1])] #set zy values
                                #first generate the scalar components and assign appropriately, then call the function 
                            
                                wf = [] #create list to keep track of rotations and rotate previous item for each time 
                                for rots in rotations[1:]: #all items except the first 
                                    if '_' in rot and wf: #if the list is empty, then use the inital w list 
                                        deg = int(rot.split('_')[-1])
                                        w = pro.rot_quart('xy', n=deg, w= wi)
                                        wf.append((w.x, w.y, w.z)) #extract quarternion components into list, w is always zero so discard  
                                
                                    elif 'c' in rot and not wf:  #if not, then use the last element in the list 
                                            deg = int(rot.split('_')[-1])
                                            w = pro.rot_quart('yz', n=deg, w=wf[-1])
                                            wf.append((w.x, w.y, w.z))
                                
                                    elif '/' in rot and wf:
                                        deg = int(rot.split('/')[-1])
                                        w = pro.rot_quart('xz', n=deg, w= wi)
                                        wf.append((w.x, w.y, w.z)) #extract quarternion components into list, w is always zero so discard  
                                
                                    elif '_' in rot and not wf:
                                        deg = int(rot.split('_')[-1])
                                        w = pro.rot_quart('xy', n=deg, w=wf[-1])
                                        wf.append((w.x, w.y, w.z))
                                
                                    elif '/' in rot and wf:
                                        deg = int(rot.split('_')[-1])
                                        w = pro.rot_quart('xz', n=deg, w= wi)
                                        wf.append((w.x, w.y, w.z)) #extract quarternion components into list, w is always zero so discard  
                                
                                    elif '/' in rot and not wf:
                                        deg = int(rot.split('_')[-1])
                                        w = pro.rot_quart('xy', n=deg, w=wf[-1]) 
                                        wf.append((w.x, w.y, w.z))
                                #after rotation, add the final point components to the previous vector if it exists 
                                final_p = wf[-1]
                                if not canon_arr: #if this is empty i.e this is the first position 
                                    canon_arr[0, 0].append((single_unit, final_p[0])) 
                                    canon_arr[0, 1].append((single_unit, final_p[1])) 
                                    canon_arr[0, 2].append((single_unit, final_p[2])) 
                                else: #use the previous coordinates along the array 
                                    canon_arr[0, 0].append((single_unit, final_p[0]+canon_arr[0, 0, -1])) 
                                    canon_arr[0, 1].append((single_unit, final_p[1]+canon_arr[0, 1, -1])) 
                                    canon_arr[0, 2].append((single_unit, final_p[2]+canon_arr[0, 2, -1])) 
                                #with this then we have taken care of rotating any line any number of times along a different plane each time. 
                        elif '#' in rot:  
                            Sarr_x = []
                            Sarr_y = []
                            Sarr_z = []
                            #initalize vertical and horizontal scalars based on value prioir to #      
                            ul = rot.split('#')[0]
                            uL = int(ul)             
                            ui_range = list(ui)  
                            units = comb_strings(self.filtered, ui_range, uL) #get a string list of fixed string length 
                            uN = len(units) #the strings length is defined as the vertical component 
                            Sten = tensor(units, dtype=str)
                            Srep = Sten.reshape((uN, uL)) #string representation of the surface 
                            Surface_arrays.append(Srep) #add to list of surfaces 
                            def_str = 'S' + f'{len(Surface_arrays)}' #defines a surface iteration 

                            if '!' in rot: #if an exclamation point is in the string, this means the intalized direction is CW 
                                direction = 'CW'
                                rot.strip('!') #remove exclamation mark to avoid later splitting errors as this occurs right after # 
                            else:
                                direction = 'CC'
                            if '+' not in rot: #if this is only a single plane rotation of a surface 
                                rotation = rot.split('#')[1:] #any item after the length index (0) 
                                p_x = 0 
                                p_y = 0 
                                p_z = 0 
                                if '_' in rot: #in the xy plane 
                                    #surface dimensions are ultimately due to the scalar, so the total area is unchanged 
                                    deg = int(rot.split('_')[-1])
                                    #generate 3 points using the 2 scalars  
                                    S = pro.rot_imag(s1=uN, s2=uL, n=deg, dir=direction)
                                    if not canon_arr: 
                                        for tup in S: 
                                            p_x += tup[0] 
                                            p_y += tup[1]
                                            Sarr_x.append((def_str, p_x))
                                            Sarr_y.append((def_str, p_y)) #add defined components to array 
                                        for i in range(len(S)): #this is defined as 4 now  
                                            Sarr_z.append((def_str,0)) #add 4 zeros for each point 
                                        Sarr_x.append((def_str, 0))
                                        Sarr_y.append((def_str, 0)) #add the last 2 origin points to the x and y coordinate arrays 
                                        canon_arr[0, 0].append(Sarr_x)
                                        canon_arr[0, 1].append(Sarr_y)
                                        canon_arr[0, 2].append(Sarr_z)
                                    elif canon_arr: 
                                        #get the components from the last vector 
                                        p_x = canon_arr[0,0,-1]
                                        p_y = canon_arr[0,1,-1]
                                        p_z = canon_arr[0,2,-1]
                                        s4 = (p_x, p_y, p_z)
                                        #sum the vectors in order from S0 to S2 in order 
                                        s1 = tuple(sum(x) for x in zip(S[0], s4))
                                        s2 = tuple(sum(x) for x in zip(S[1], s1))
                                        s3 = tuple(sum(x) for x in zip(S[2], s2))
                                        pts = [s1,s2,s3,s4] 
                                        for idx, tup in enumerate(pts):
                                            for x,y,z in tup:
                                        #full the surface array then add it to the the canon array 
                                                Sarr_x.append((def_str, x))
                                                Sarr_y.append((def_str, y))
                                                Sarr_z.append((def_str, z))
                                        canon_arr[0, 0].append(Sarr_x)
                                        canon_arr[0, 1].append(Sarr_y)
                                        canon_arr[0, 2].append(Sarr_z)

                                elif '/' in rot: #in the xz plane  
                                    deg = int(rot.split('/')[-1])  
                                    S = pro.rot_imag(s1=uN, s2=uL, n=deg, dir=direction)
                                    if not canon_arr: 
                                        for tup in S: 
                                            p_x += tup[0] 
                                            p_z += tup[1]
                                            Sarr_x.append((def_str, p_x))
                                            Sarr_z.append((def_str, p_z)) 
                                        for i in range(len(S)):  
                                            Sarr_y.append((def_str,0)) 
                                        Sarr_x.append((def_str, 0))
                                        Sarr_z.append((def_str, 0)) 
                                        canon_arr[0, 0].append(Sarr_x)
                                        canon_arr[0, 1].append(Sarr_y)
                                        canon_arr[0, 2].append(Sarr_z)
                                    else:  
                                        p_x = canon_arr[0,0,-1]
                                        p_y = canon_arr[0,1,-1]
                                        p_z = canon_arr[0,2,-1]
                                        s4 = (p_x, p_y, p_z) 
                                        s1 = tuple(sum(x) for x in zip(S[0], s4))
                                        s2 = tuple(sum(x) for x in zip(S[1], s1))
                                        s3 = tuple(sum(x) for x in zip(S[2], s2))
                                        pts = [s1,s2,s3,s4] 
                                        for idx, tup in enumerate(pts):
                                            for x,y,z in tup:
                                                Sarr_x.append((def_str, x))
                                                Sarr_y.append((def_str, y))
                                                Sarr_z.append((def_str, z))
                                        canon_arr[0, 0].append(Sarr_x)
                                        canon_arr[0, 1].append(Sarr_y)
                                        canon_arr[0, 2].append(Sarr_z)

                                elif 'c' in rot: #in the yz plane 
                                    deg = int(rot.split('c')[-1])
                                    #generate 3 points using the 2 scalars  
                                    S = pro.rot_imag(s1=uN, s2=uL, n=deg, dir=direction)
                                    if not canon_arr: 
                                        for tup in S: 
                                            p_y += tup[0] 
                                            p_z += tup[1]
                                            Sarr_y.append((def_str, p_x))
                                            Sarr_z.append((def_str, p_y)) #add defined components to array 
                                        for i in range(len(S)): #this is defined as 4 now  
                                            Sarr_x.append((def_str,0)) #add 4 zeros for each point 
                                        Sarr_y.append((def_str, 0))
                                        Sarr_z.append((def_str, 0)) #add the last 2 origin points to the x and y coordinate arrays 
                                        canon_arr[0, 0].append(Sarr_x)
                                        canon_arr[0, 1].append(Sarr_y)
                                        canon_arr[0, 2].append(Sarr_z)
                                    else: 
                                        #get the components from the last vector 
                                        p_x = canon_arr[0,0,-1]
                                        p_y = canon_arr[0,1,-1]
                                        p_z = canon_arr[0,2,-1]
                                        s4 = (p_x, p_y, p_z)
                                        #sum the vectors in order from S0 to S2 in order 
                                        s1 = tuple(sum(x) for x in zip(S[0], s4))
                                        s2 = tuple(sum(x) for x in zip(S[1], s1))
                                        s3 = tuple(sum(x) for x in zip(S[2], s2))
                                        pts = [s1,s2,s3,s4] 
                                        for idx, tup in enumerate(pts):
                                            for x,y,z in tup:
                                                Sarr_x.append((def_str, x))
                                                Sarr_y.append((def_str, y))
                                                Sarr_z.append((def_str, z))
                                        canon_arr[0, 0].append(Sarr_x)
                                        canon_arr[0, 1].append(Sarr_y)
                                        canon_arr[0, 2].append(Sarr_z)

                            #The only difference is that we now have 3 points to pass instead of 1 
                            elif '+' in rot: #apply quarternion rotation similar to the single lines 
                                rotations = rot.split('+') #create a list of possible rotations 
                                plane = ''
                                if '_' in rotations[0]:  #alomg the xy plane  
                                    deg = int(rotations[0].split('_')[-1]) 
                                    vecs = pro.rot_imag(uL, uN, n = deg, dir=direction) 
                                    v1 = vecs[0]
                                    v2 = vecs[1]
                                    v3 = vecs[2]
                                    wi = [(v1[0], v1[1], 0),(v2[0],v2[1],0),(v3[0],v3[1],0)] #set xy values
                                
                                if '/' in rotations[0]: #along the xz plane 
                                    deg = rotations[0],split('/')[-1]
                                    vecs = pro.rot_imag(uL, uN, n = deg, dir=direction) 
                                    v1 = vecs[0]
                                    v2 = vecs[1]
                                    v3 = vecs[2]
                                    wi = [(v1[0], 0, v1[1]),(v2[0],0,v2[1]),(v3[0],0,v3[1])]
                                    #the result is always tuple of 2 indecies 0 and 1  
                                    wi = [(v1[0],v1, 0)]
                                if 'c' in rotations[0]: #along the yz plane 
                                    deg = rotations[0],split('c')[-1]
                                    vecs = pro.rot_imag(uL, uN, n = deg, dir=direction) 
                                    v1 = vecs[0]
                                    v2 = vecs[1]
                                    v3 = vecs[2]
                                    wi = [(0, v1[0], v1[1]),(0,v2[0],v2[1]),(0,v3[0],v3[1])] 
                                
                                #store each components rotation in a seperate list 
                                wf1 = []
                                wf2 = []
                                wf3 = []
                                #iterate over the possible rotations of the surface and apply each 
                                for rots in rotations: 
                                    #if all lists are empty, add components based on the rotation by summing the rotated components with each other 
                                    if '_' in rots and not wf1 and not wf2 and not wf3:  
                                        deg = rots.split('_')[-1]
                                        w1 = pro.rot_quart('xy', n=deg, w=wi[0])
                                        w2 = pro.rot_quart('xy', n=deg, w=wi[1])
                                        w3 = pro.rot_quart('xy', n=deg, w=wi[2])
                                        wf1.append(w1)
                                        wf2.append(tuple(sum(x) for x in zip(w1,w2)))
                                        wf3.append(tuple(sum(x) for x in zuo(w2, w3)))
                                    if '/' in rots and not wf and not wf2 and not wf3:
                                        deg = rots.split('_')[-1]
                                        w1 = pro.rot_quart('xz', n=deg, w=wi[0])
                                        w2 = pro.rot_quart('xz', n=deg, w=wi[1])
                                        w3 = pro.rot_quart('xz', n=deg, w=wi[2])
                                        wf1.append(w1)
                                        wf2.append(tuple(sum(x) for x in zip(w1,w2)))
                                        wf3.append(tuple(sum(x) for x in zuo(w2, w3)))
                                    if 'c' in rots and not wf and not wf2 and not wf3: 
                                        deg = rots.split('_')[-1]
                                        w1 = pro.rot_quart('yz', n=deg, w=wi[0])
                                        w2 = pro.rot_quart('yz', n=deg, w=wi[1])
                                        w3 = pro.rot_quart('yz', n=deg, w=wi[2])
                                        wf1.append(w1)
                                        wf2.append(tuple(sum(x) for x in zip(w1,w2)))
                                        wf3.append(tuple(sum(x) for x in zip(w2, w3)))
                                    #if the lists are not empty, perform rotation on the last rotated components iteratively 
                                    elif '_' in rots and wf1 and wf2 and wf3: 
                                        w1 = pro.rot_quart('xy', n=deg, w=wf1[-1])
                                        w2 = pro.rot_quart('xy', n=deg, w=wf2[-1])
                                        w3 = pro.rot_quart('xy', n=deg, w=wf3[-1])
                                        wf1.append(w1)
                                        wf2.append(tuple(sum(x) for x in zip(w1,w2)))
                                        wf3.append(tuple(sum(x) for x in zip(w2, w3)))
                                    elif '/' in rots and wf1 and wf2 and wf3: 
                                        w1 = pro.rot_quart('xz', n=deg, w=wf1[-1])
                                        w2 = pro.rot_quart('xz', n=deg, w=wf2[-1])
                                        w3 = pro.rot_quart('xz', n=deg, w=wf3[-1])
                                        wf1.append(w1)
                                        wf2.append(tuple(sum(x) for x in zip(w1,w2)))
                                        wf3.append(tuple(sum(x) for x in zip(w2, w3)))
                                    elif 'c' in rots and wf1 and wf2 and wf3: 
                                        w1 = pro.rot_quart('yz', n=deg, w=wf1[-1])
                                        w2 = pro.rot_quart('yz', n=deg, w=wf2[-1])
                                        w3 = pro.rot_quart('yz', n=deg, w=wf3[-1])
                                        wf1.append(w1)
                                        wf2.append(tuple(sum(x) for x in zip(w1,w2)))
                                        wf3.append(tuple(sum(x) for x in zip(w2, w3)))
                                #extract the rotated tuple values for the surfac points after filling wf1,wf2  and wf3 
                                s1 = wf1[-1]
                                s2 = wf2[-1]
                                s3 = wf3[-1]
                                s4 = (0,0,0)
                                #add the previous vector or the origin to every tupule to connect the structure 
                                if not canon_arr:
                                    pass
                                elif canon_arr: 
                                    s4 = (canon_arr[0,0,-1],canon_arr[0,1,-1],canon_arr[0,2,-1])
                                    pts = [s1,s2,s3,s4]
                                for s in pts[:3]: #exclude the last tuple as to not add it to itself  
                                    for x,y,z in s: 
                                        x += s4[0]
                                        y += s4[1]
                                        z += s4[2]
                                        Sarr_x.append((def_str,x))
                                        Sarr_y.append((def_str, y))                         
                                        Sarr_z.append((def_str, z))
                                canon_arr[0,0].append(Sarr_x)
                                canon_arr[0,1].append(Sarr_y)
                                canon_arr[0,2].append(Sarr_z)
        
        #set the protiens canon array and surfaces representation 
        self.T0 = canon_arr 
        self.surfaces = Surface_arrays 
        return self.T0, self.surfaces                      

test1 = Protein(TBP, 'TBP')
test2 = test1.fold
test2


                            



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



