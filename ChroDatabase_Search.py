import pandas as pnd # import pandas library for data manipulation
import os

'''Chromosome engine - use table keys to find specific info on protiens, create archive of all 24 chromosomes (22 + XY)'''

# Change the working directory to the correct folder
os.chdir('/Users/ismailnashwan/Desktop/GeneDatabase_Project')

# Function to generate the ordered sequence lists 
def Loci_adder(chromfile_label): 
    try:
        data = pnd.read_csv(chromfile_label, header=1)
        # Remove header and modify to get attributes 
        data.drop('Accession', axis=1, inplace=True)
        data.insert(0, 'Sequences', ['seq' + str(i+1) for i in range(len(data))])
        cols = [data[i] for i in data.columns] # 0: seq n, 1,2: start and stop, 3: plus/-, 4:ID, 5: desc
        return cols 
    except FileNotFoundError:
        print(f"The file {chromfile_label} was not found.")

# Define the chromosome loci contents 
class Chromosome:  
    def __init__(self, chro=''):  
        attrs = Loci_adder(chro) # Get CSV file to add to class 
        seq_list = attrs[0].tolist() # Store sequences list with indexes 
        self.seqn = [seq.strip('seq') for seq in seq_list]
        self.label = chro.strip('.csv') # Store chromosome label containing loci 
        self.Loci_list = attrs[6].tolist() # Define every single Loci as part of chromosome 
        self.Bases = (attrs[2] - attrs[1]).tolist()  # To list converts this to a list instead of []
        self.Proteins = (attrs[3]).tolist()
        self.ID = attrs[5].tolist()
        self.coding = []
        
        for index, relevance in enumerate(attrs[4]):
            if relevance == 'plus':
                self.coding.append('Within chromosome')
            elif relevance == 'minus':
                self.coding.append('Within partner chromosome')
            else:
                self.coding.append('Not Available ???')

    def __str__(self): 
        return f'{self.label}'


'''Chromosome archvie initialization'''

C1 = Chromosome('Chrom1.csv')
C2 = Chromosome('Chrom2.csv')
C3 = Chromosome('Chrom3.csv')
C4 = Chromosome('Chrom4.csv')
C7 = Chromosome('Chrom7.csv')
C8 = Chromosome('Chrom8.csv')
C11 = Chromosome('Chrom11.csv')
C12 = Chromosome('Chrom12.csv')
C14 = Chromosome('Chrom14.csv')
C16 = Chromosome('Chrom16.csv')
C17 = Chromosome('Chrom17.csv')
C19 = Chromosome('Chrom19.csv')
C21 = Chromosome('Chrom21.csv')
C22 = Chromosome('Chrom22.csv')

Chromosomes = {name: str(obj) for name, obj in globals().items() if isinstance(obj, Chromosome)}

class Loci: 
    def __init__(self, chro, seq):  
        if isinstance(chro, Chromosome): 
            if isinstance(seq, int):  # Ensure seq is an integer
                number = seq - 1 # Convert 1-based seq to 0-based index
                if number < 0 or number >= len(chro.Bases):  # Validate range
                    raise IndexError(f"{chro.label} does not possess this sequence.")
                # Extract info
                self.Origin = chro.label  
                self.seqn = seq  # Store 1-based sequence
                self.Bases = chro.Bases[number]
                self.desc = chro.Loci_list[number]
                self.ID = chro.ID[number]
                self.protein = chro.Proteins[number]
                self.coding = chro.coding[number]
            else: 
                raise TypeError(f'Please enter an integer number corresponding to sequence position along {chro.label}')
        else:
            raise TypeError("Input must be a Chromosome object.")

    def __str__(self):
        # Define the dictionary within the __str__ method
        dic = {
            'Sequence Number': self.seqn,
            'Length': self.Bases,
            'Protein': self.protein,
            'Description': self.desc
        }
        
        # Create the formatted string with newlines between key-value pairs
        formatted_string = "   ".join([f"{key}: {value}" for key, value in dic.items()])
        
        # Ensure newlines are properly returned
        return formatted_string
    
    def length(self):
        return f'This locus contains {self.Bases} bases'
    
    def codability(self):
        if self.coding == 'Within chromosome':
            return f'locus {self.seqn} is along the template strand'
        else:
            return f'locus {self.seqn} is along the complementary strand'
    
    def relation(self):
        return f'{self.desc}'
    
    def get_ID(self): 
        return f'{self.ID} is the NCBI id'

    def get_info(self):
        director = 'https://www.ncbi.nlm.nih.gov/gene/' 
        return f'{director}{self.ID}'

def Find_Loci(chrom, **kwargs):
    if not isinstance(chrom, Chromosome):
        raise TypeError("Input must be a Chromosome object.")
    
    # Parse keyword arguments
    seqn = kwargs.get('seqn', None)
    larger_than = kwargs.get('larger_than', None)
    smaller_than = kwargs.get('smaller_than', None)
    symbol = kwargs.get('symbol', None)
    ID = kwargs.get('ID', None)

    # Initialize matched list
    matched_list = []

    # If any key arguments, i.e seqn, ID, protein name, return result immediately for first match 
    if symbol is not None:
        #use an underscore to filter out n characters 
        #prefex determines number of filtered characters in symbol 
        if "_" in symbol:
            symbol_prefix, length = symbol.split('_')
            try:
                length = int(length)  # Convert length to an integer
            except ValueError:
                print(f"Invalid length value for symbol: {length}. It should be an integer.")
                return matched_list
        else:
            symbol_prefix = symbol 
            length = len(symbol)  # Default length to symbol length for match 
        
        # Match the first 'length' characters of the symbol
        for index, symb in enumerate(chrom.Proteins):
            if symb[:length] == symbol_prefix[:length]:  # Match the first 'length' characters
                matched_list.append(str(Loci(chrom, index + 1)))


    if ID is not None:
        for index, chrom_id in enumerate(chrom.ID):
            if chrom_id == ID:
                matched_list.append(str(Loci(chrom, index + 1)))
                return matched_list
    if seqn is not None:
        for index, loc_seqn in enumerate(chrom.seqn):
            if int(loc_seqn) == seqn:
                matched_list.append(str(Loci(chrom, index + 1)))
                return matched_list

    # If base range information is provided, iterate through base lengths 
    if larger_than is not None and smaller_than is not None:
        for seq, bases in enumerate(chrom.Bases):
            if larger_than < bases < smaller_than:
                matched_list.append(str(Loci(chrom, seq + 1)))

    elif larger_than is not None:
        for seq, bases in enumerate(chrom.Bases):
            if bases > larger_than:
                matched_list.append(str(Loci(chrom, seq + 1)))

    elif smaller_than is not None:
        for seq, bases in enumerate(chrom.Bases):
            if bases < smaller_than:
                matched_list.append(str(Loci(chrom, seq + 1)))

    # Return matched loci
    return matched_list



def seq_search():  # Sequence Search
    print("\nWelcome to the local search engine!")
    display = True
    while display:
        
        # Retrieve all Chromosome objects from the global namespace
        chromosomes = {name: obj for name, obj in globals().items() if isinstance(obj, Chromosome)}
        
        # Display the list of available chromosomes
        if chromosomes:
            chrom_list = ", ".join(f"{name} ({chrom.label})" for name, chrom in chromosomes.items())
            chrom_all = [chrom for _,chrom in chromosomes.items()] #get a list of all chromosomes available 
            print(f"Currently, your defined chromosomes are: {chrom_list}")
        else:
            print("No chromosomes have been defined yet.")
        
        chromosome_name = input('\nPlease select an available chromosome, "all" for a general archive search, or exit to shut down the engine: ') 
        if chromosome_name.lower() == 'exit':
            print('Exiting...')
            break 
 
        '''If all chromosomes are being analyzed'''
        
        if chromosome_name.lower() == 'all':
            print(f"\nYou've selected {chromosome_name}, Please insert attributes to search.")
            print("\nPossible filters to get info: \nseqn= sequence number\nID= NSBI ID\nsymbol=symb\nlarger/smaller_than or both")
            constraints = input('Please enter filtering items. If multiple arguments, please separate using +. If using ranges, please include the underscore\n')
            
            # Split the input if multiple arguments are included
            if '+' in constraints:
                input_list = [constraint for constraint in constraints.split('+')]
            else:  # If only one argument is used
                input_list = [constraints]
            
            # Convert the input list into keyword arguments
            kw_args = {}
            for item in input_list:
                # Split the item into key and value based on '='
                if '=' in item:
                    key, value = item.split('=')
                    key = key.strip()
                    value = value.strip()  # Remove any extra spaces in the argument
                    # Convert value to the correct type based on the key
                    if key == 'seqn' :
                        try:
                            kw_args[key] = int(value)  # Convert seqn to integer
                        except ValueError:
                            print(f"Invalid value for '{key}': {value}. It should be an integer.")
                            continue
                    elif key in ['larger_than', 'smaller_than']:
                        try:
                            kw_args[key] = float(value)  # Convert range arguments to float
                        except ValueError:
                            print(f"Invalid value for '{key}': {value}. It should be a float.")
                            continue
                    else:
                        kw_args[key] = value  # Keep other arguments as string (e.g., ID, symbol)
                else:
                    print(f"Invalid argument format for '{item}'. Please use 'key=value' format.")
            
            print(f'{kw_args} your input for all')  # Debug: Check the contents of the dictionary
            
            # If keyword arguments are provided, call the Find_Loci function 
            results = {}
            if kw_args: 
                for chrom in chrom_all:
                    loci = Find_Loci(chrom, **kw_args)  # Call the Find_Loci function for each chromosome
                    if loci:
                        results[chrom] = loci  # Add to results if matches are found
            if results:
                print('Your matches are as follows:')
                #iterate over chrom 
                for chrom in results: #iterate over each match chrom 
                    if len(results[chrom]) > 10:
                        print('\n' * 5 + '-' * 30)
                        print('\n' * 5 + f'For chromosome {chrom}:') #print chrom 
                        print('\n' * 5 + '-' * 30)
                    else:
                        print('\n' * 5 + f'For chromosome {chrom}:')
                    for locus in results[chrom]: #iterate over value list during chrom iteration  
                        print(f'\n{locus}') #print iteratable value 
                        
            resume = input('\nWould you like to continue or exit?\ninput "yes" for another search, or "exit" to exit\n')
            if resume == 'yes':
                pass
            elif resume == 'exit':
                print('Exiting...')
                break
            else:
                print('\nNo valid match found. Please try again or input "exit" to quit.')
                pass 
        else:
            '''If a specific chromosome is choosen'''
            chromosome = chromosomes.get(chromosome_name)  # Fetch the Chromosome object by name 
            
            if chromosome:
                print(f"\nYou've selected {chromosome_name}, Please insert attributes to search.")
                print("\nPossible filters to get info: \nseqn= sequence number\nID= NSBI ID\nsymbol=symb\nlarger/smaller_than or both")
                
                constraints = input('Please enter filtering items. If multiple arguments, please separate using +. If using ranges, please include the underscore\n')
                # Split the input if multiple arguments are included
                if '+' in constraints:
                    input_list = [constraint for constraint in constraints.split('+')]
                else:  # If only one argument is used
                    input_list = [constraints]
                
                # Convert the input list into keyword arguments
                kw_args = {}
                for item in input_list:
                    # Split the item into key and value based on '='
                    if '=' in item:
                        key, value = item.split('=')
                        key = key.strip()
                        value = value.strip()  # Remove any extra spaces in the argument
                        # Convert value to the correct type based on the key
                        if key == 'seqn' :
                            try:
                                kw_args[key] = int(value)  # Convert seqn to integer
                            except ValueError:
                                print(f"Invalid value for '{key}': {value}. It should be an integer.")
                                continue
                        elif key in ['larger_than', 'smaller_than']:
                            try:
                                kw_args[key] = float(value)  # Convert range arguments to float
                            except ValueError:
                                print(f"Invalid value for '{key}': {value}. It should be a float.")
                                continue
                        else:
                            kw_args[key] = value  # Keep other arguments as string (e.g., ID, symbol)
                    else:
                        print(f"Invalid argument format for '{item}'. Please use 'key=value' format.")
                
                print(f'{kw_args} your input for {chromosome}')  # Debug: Check the contents of the dictionary
                
                # If keyword arguments are provided, call the Find_Loci function
                if kw_args:
                    results = Find_Loci(chromosome, **kw_args)  # Call the Find_Loci function with keyword arguments
                    print('The following matches your filters:')
                    for result in results:
                        print(f'{result.label}')
                    resume = input('\nWould you like to continue or exit?\ninput "yes" for another search, or "exit" to exit\n\n')
                    if resume == 'yes':
                        pass
                    elif resume == 'exit':
                        print('Exiting...')
                        break
                else:
                    print('\nNo valid match found. Please try again or input "exit" to quit.')
            else:
                print(f'\n{chromosome} is not defined, please define alongside defined chromosomes')
if __name__ == "__main__":
    seq_search()
    
    