
import numpy as num
from numpy import array as tensor
import quaternion 


""""
numpy functions 

    create arrays 

    num.array()      # Create a tensor from a list or other array-like object
    num.zeros()      # Create a tensor filled with zeros
    num.ones()       # Create a tensor filled with ones
    num.empty()      # Create an uninitialized tensor
    num.full()       # Create a tensor filled with a specified value
    num.arange()     # Create a tensor with evenly spaced values
    num.linspace()   # Create a tensor with linearly spaced values
    num.random.rand()  # Create a tensor with random values (uniform distribution)
    num.random.randn() # Create a tensor with random values (normal distribution)

    num.eye()        # Create an identity matrix
    num.diag()       # Create a diagonal matrix or extract diagonal elements
    num.reshape()    # Change the shape of a tensor
    num.ravel()      # Flatten a tensor into a 1D array
    num.flatten()    # Create a copy of the tensor as a flattened 1D array
    num.expand_dims() # Add a new axis (increase dimensions)


    num.squeeze()    # Remove axes of size 1
    num.transpose()  # Transpose a tensor (swap axes)
    num.swapaxes()   # Swap two specific axes
    num.moveaxis()   # Move specific axes to new positions
    num.take()       # Extract elements based on indices
    num.put()        # Place elements into specific positions
    num.where()      # Find indices or conditionally select elements
    num.nonzero()    # Indices of non-zero elements
    num.extract()    # Extract elements based on a condition
    num.concatenate()  # Concatenate tensors along a specific axis
    num.stack()        # Stack tensors along a new axis
    num.hstack()       # Horizontally stack tensors
    num.vstack()       # Vertically stack tensors
    num.dstack()       # Stack tensors along the third axis

    num.split()        # Split a tensor into multiple sub-tensors
    num.hsplit()       # Split horizontally
    num.vsplit()       # Split vertically
    num.dsplit()       # Split along the third axis

    num.tile()         # Repeat a tensor along specified axes
    num.repeat()       # Repeat elements of a tensor

    Element operations - 
    num.add(), num.subtract(), num.multiply(), num.divide()
    num.power(), num.mod(), num.sqrt()

    num.sum(), num.prod(), num.mean(), num.std(), num.var()
    num.min(), num.max(), num.argmin(), num.argmax()

    properties 

    tensor.shape     # Shape of the tensor
    tensor.size      # Total number of elements
    tensor.ndim      # Number of dimensions
"""

def rangemod(n, k):
    return [(n//6, k//6)]

def create_grid(unit):
    length = (unit*2)+1 #including negative counterpart 
    orig = (0,0)
    #build the array axis by axis, as well as quadrant by quadrant
    c = (0,0) #origin 
    xp = [(x,0) for x in range(1, unit+1)] #posx
    xn = [(-x,0) for x in range(1, unit+1)][::-1] #negx
    yp = [(0,y) for y in range(1, unit+1)][::-1] #pos y vert
    yn = [(0,-y) for y in range(1, unit+1)] #neg y 
    yxp = [(x,y) for x in range(1, unit+1) for y in range(1, unit+1)] #pos xy 
    yxn = [(-x,-y) for x in range(1, unit+1) for y in range(1, unit+1)][::-1] #neg xy
    y_x = [(x,y) for x,_ in xn[::-1] for _,y in yp[::-1]] #neg x, y
    x_y = [(x,y) for x,_ in xp for _,y in yn] #negative y values #x, neg y

    #create positive quadrant (1 counter)
    yxp = tensor(yxp, dtype=[('x', 'i4'), ('y', 'i4')])
    yxp_d2 = yxp.reshape(unit,unit).T[::-1, :]

    #create _x quadrant (2 counter)
    y_x = tensor(y_x, dtype=[('x', 'i4'), ('y', 'i4')])
    y_xd2 = y_x.reshape(unit, unit).T
    y_xd2 = y_xd2[:, ::-1]
    y_xd2 = y_xd2[::-1, :] 

    #create negative quadrant (3 counter)
    yxn = tensor(yxn, dtype=[('x', 'i4'), ('y', 'i4')])
    yxn_d2 = yxn.reshape(unit, unit).T[::-1]

    #create _y quadrant (4 counter)
    x_y = tensor(x_y, dtype=[('x', 'i4'), ('y', 'i4')])
    x_yd2 = x_y.reshape(unit, unit) 
    x_yd2 = x_yd2.T 

    #initalize main template 
    main = num.zeros((length,length), dtype=object)

    #create axis and origin 
    main[unit, unit] = c  
    main[unit, unit+1:length] = xp 
    main[unit, 0:unit] = xn
    for j, tup in zip(range(unit+1, length), yn):  # Pair each `j` with a tuple
        if main[j, unit] == 0:  # Check if the current entry is 0
            main[j, unit] = tup  # Replace 0 with the tuple
    for j, tup in zip(range(0, unit), yp):
        if main[j, unit] == 0: 
            main[j, unit] = tup 

    #add quadrants in counter order, recall that the end is non inclusive while the start is 

    main[0:unit, unit+1:length] = yxp_d2 #quad 1 
    main[0:unit, 0:unit] = y_xd2 #quad 2 
    main[unit+1:length, 0:unit] = yxn_d2 #quad 3 
    main[unit+1:length, unit+1:length] = x_yd2 #quad 4 
    return main 

def pad_row(row, length):
    return row + [(0, 0)] * (length - len(row))


def filter(arr):
    # Create a mask to filter out rows where both 'seq' and 'idx' are (0, 0)
    mask = ~( (arr['seq'] == b'0') & (arr['idx'] == 0) )  # ~ inverts the mask, keeping only non-(0, 0) items
    filtered_arr = arr[mask]
    
    return filtered_arr



def pro_array(proseq, show=False):
    # Define amino acid sets
    helix_set = {'E', 'L', 'K', 'A', 'M', 'Q', 'S'}
    beta_set = {'V', 'I', 'T', 'G', 'P', 'F', 'Y', 'W'}
    amino_acids = set('ACDEFGHIKLMnumQRSTVWY')  # All amino acids
    string_set = {aa for aa in amino_acids if aa not in helix_set and aa not in beta_set}

    # Initialize dimensions and parameters
    helix_dim, beta_dim, string_dim = [], [], []
    search = 6  # Fixed chunk size
    i = 0  # Index to iterate over the sequence

    while i < len(proseq):
        chunk = proseq[i:i + search]

        # Count matches for each category
        helix_match = sum(1 for aa in chunk if aa in helix_set)
        beta_match = sum(1 for aa in chunk if aa in beta_set)
        string_match = sum(1 for aa in chunk if aa in string_set)

        # Show chunk comparison logic
        if show:
            print(f"Chunk: {chunk}")
            print(f"Helix Matches: {helix_match}, Beta Matches: {beta_match}, String Matches: {string_match}")

        # Check for G or P in the chunk
        if 'G' in chunk or 'P' in chunk:
            idx_gp = [idx for idx, aa in enumerate(chunk) if aa in {'G', 'P'}]
            for idx in idx_gp:
                left = chunk[:idx]
                right = chunk[idx + 1:]

                # Print the chunk split details for show=True
                if show:
                    print(f"  Found G/P at index {idx}. Left part: {left}, Right part: {right}")

                # Process the left part
                if len(left) >= 2:
                    helix_match_left = sum(1 for aa in left if aa in helix_set)
                    beta_match_left = sum(1 for aa in left if aa in beta_set)
                    if helix_match_left >= 2:
                        helix_dim.append((left, f'{i // search}'))
                        if show:
                            print(f"  Classified left part as Helix")
                    elif beta_match_left >= 3:
                        beta_dim.append((left, f'{i // search}'))
                        if show:
                            print(f"  Classified left part as Beta")
                    else:
                        string_dim.append((left + chunk[idx], f'{i // search}'))
                        if show:
                            print(f"  Classified left part as String")

                # Process the right part
                if len(right) >= 2:
                    helix_match_right = sum(1 for aa in right if aa in helix_set)
                    beta_match_right = sum(1 for aa in right if aa in beta_set)
                    if helix_match_right >= 2:
                        helix_dim.append((right, f'{i // search}'))
                        if show:
                            print(f"  Classified right part as Helix")
                    elif beta_match_right >= 3:
                        beta_dim.append((right, f'{i // search}'))
                        if show:
                            print(f"  Classified right part as Beta")
                    else:
                        string_dim.append((chunk[idx] + right, f'{i // search}'))
                        if show:
                            print(f"  Classified right part as String")

            # Move to the next chunk after this split
            i += search
            continue

        # Standard classification if no G or P is found
        if helix_match >= 4 or (helix_match == 3 and string_match + beta_match == 3):
            helix_dim.append((chunk, f'{i // search}'))
            if show:
                print(f"Classified as: Helix\n")
        elif beta_match >= 4:
            beta_dim.append((chunk, f'{i // search}'))
            if show:
                print(f"Classified as: Beta\n")
        else:
            string_dim.append((chunk, f'{i // search}'))
            if show:
                print(f"Classified as: String\n")

        # Move to the next chunk
        i += search

    # Find the maximum row length for padding
    max_len = max(len(helix_dim), len(beta_dim), len(string_dim))
    max_sqr = num.sqrt(max_len)

    if isinstance(max_sqr, float): #if not perfect square 
        perfect_sqr = num.ceil(max_sqr) #get next perfect square root  
        max_len = int(perfect_sqr**2) #adjust to next perfect square n 
    else:
        perfect_sqr = int(max_sqr)

    helix_dim_padded = pad_row(helix_dim, max_len)
    beta_dim_padded = pad_row(beta_dim, max_len)
    string_dim_padded = pad_row(string_dim, max_len)

    # Convert to a tensor-like array
    base_2 = tensor([helix_dim_padded, beta_dim_padded, string_dim_padded], dtype=[('seq', 'S10'), ('idx', 'i4')])
    return base_2.reshape(1, 3, int(perfect_sqr), int(perfect_sqr))

def domain_index(arrays, min_idx, max_idx, show=False):
    matched_components = []
    
    # Apply floor division by 6 to the range values
    min_group = min_idx // 6
    max_group = max_idx // 6
    
    # Iterate through the array, extracting sequences and indices
    for i in range(arrays.shape[1]):
        for j in range(arrays.shape[2]):
            subarray = arrays[0, i, j, :]
            
            for tup in subarray:
                seq, idx = tup
                
                # Check if the index falls within the adjusted range
                if min_group <= idx <= max_group:
                    matched_components.append((tup, i, j))  # Store the tuple along with the (i, j) indices
                    if show:
                        print(f"Matched tuple: {tup} (Index: {idx}) from dimensions ({i}, {j})")
    
    # Sort the matched components by the index
    matched_components.sort(key=lambda x: x[0][1])  # Sort by the index value (second element of the tuple)
    
    results = ''
    for (seq, idx), i, j in matched_components:
        results += '\n' + f"({seq}, {idx}) (Index: {idx}) from dimensions ({i}, {j})"
    
    print(results)


prion = 'MRKHLSWWWLATVCMLLFSHLSAVQTRGIKHRIKWNRKALPSTAQITEAQVAENRPGAFIKQGRKLDIDFGAEGNRYYEANYWQFPDGIHYNGCSEANVTKEAFVTGCINATQAANQGEFQKPDNKLHQQVLWRLVQELCSLKHCEFWLERGAGLRVTMHQPVLLCLLALIWLTVK'


prion_array = pro_array(prion)


'''
The domain_index function is used to find the units of interest required to form the protein 

            This construction occurs by subsequently adding each unit prioir as it's own vector, everytime we define
            the scalar for the position as the unit length of the structure, using complex rotation and quarternion rotation. 

            we will define the axis such that z is perpendicular to x and y. as such, we have 
            k the xy rotational axis, j, the xz rotational axis and i the zy rotational axis. 

            as such, the euler identities for each initial rotation will be defined accordingly 
            the general euler identity is e**num = cos(p)+sin(p)n for some imaginary axis n 


            / defines j rotation, where x is horizontal and z is vertical, so j defines the imaginary axis "z", q= -1j
            _ defines k rotation, where x is horizontal and y is vertical, so k defines the imaginary axis "y", q=1k 
            c defines i rotation, where y is horizontal and z is vertical, so i defines the imaginary axis 'x', q=1i 

            note that the "imaginary axis" variable simply designates the created quarternion for a line or surface point 
            that is, after inital rotation, a component carrying j is defined as a z axis component, k as a y component, and i as an x component  
            little confusing, but it's how quarternions are defined and the easiest way to assign components into quarternions. 

            mathematically, the inital rotation along an imaginary plane seperates the scalar value of the unit 
            into components representing it's position from origin according to the rotation. if + is passed, then 
            the coordinates are converted into a quarternion and multiplied by another quarternion defined with a unit 
            vector of 1 (rotating at 1 scaling of cos theta/2 + sintheta/2 axis) after which the inverse quarternion is 
            applied to reset the coordinates to the originally defined coordinate system, thus representing a 3 dimensional 
            position from an initially 2 dimensional vector. 

            once the final vector is extracted following all the applied rotations, it's components are summed to the previous
            vector components, and a specified position vector is created with 3 coordinates along the array. 

            For a surface, the exact same thing happens, except we define counter or clockwise rotation by using 
            the appropraite euler identity and a list [e**num, e**n2p] or [e**n-1, e**n-2] which is cc or cw 
            this defines the inital orientation relative to the previous vector, after which the appropraite rotations are 
            applied, where the scaled values are the horizontal vector = choosen lanegth and the vertical vectpr = count 
            assigned D[0] and D[1] respectively according to cc or cw. we need a third component to represent another position
            (rectangular) however this is simply the exact reverse of the original vertical component.
            After applying all the rotations and getting the reverse component, we sum the components with each other to 
            define 4 seperate points 
            point 4 is defined as the sum of all previous vectors and sets an exact connection with previous structure 
            point 1 is summed to 4 and is the vertically oriented point from 4 
            point 2 is summed to the original sum and desiginates the horizontal movement of that vertical position 
            point 3 is simply the same thing again where we sum the negative of the vertical component, which will 
            generate a horizontally shifted position from the original point 4 

            represent x,y,z points as 12 points total for the 4 points. 
'''

#if s1, return i = 0, 1 as real, imaginary 
#if s2m return i=0, 1, 2 as vert 1, horiz 1, and vert 1 

def rot_imag(s1, n=0, s2=None, dir='CC'): #rotation function using the lengths of the sequences and degrees provided 

    #first define the used euler identities, numpy already does this for us, so we only need 
    #to define the final result based on the plane. 
    i= 1j #define i as the imaginary unit 
    deg = num.radians(n) #convert to radians  
    p = num.radians(90)

    #each third return item designates the components as quarternion components for the rot function 
    
    vec = s1*num.exp(deg*i)
    v1 = vec.real 
    v2 = vec.imag 

    if not s2: #if only a single line is being rotated 
        return (v1, v2) #rotation complete, horizontal and vertical components returned  

    if dir == 'CC':
        direction_multiplier = 1  # Counterclockwise (positive direction)
    elif dir == 'CW':
        direction_multiplier = -1  # Clockwise (negative direction)
    else:
        raise ValueError("Direction must be either 'CC' or 'CW'")

    # Counterclockwise (CC) and Clockwise (CW) direction definitions
    CC = [num.exp(i * p), num.exp(2 * i * p)]
    CW = [num.exp(-i * p), num.exp(-2 * i * p)]

    # Select direction-dependent values
    if dir == 'CC':
        v1 = s1 * CC[0]  # set i for CC
        v2 = s2 * CC[1]  # set -1 for CC  
    elif dir == 'CW':
        v1 = s1 * CW[0]  # set -i for CW
        v2 = s2 * CW[1]  # set -1 for CW
    # Apply rotation to v1 and v2 and generate w1-w3 
    w1 = v1 * num.exp(direction_multiplier * i * deg)  
    w1r = w1.real
    w1i = w1.imag
    w2 = v2 * num.exp(direction_multiplier * i * deg) 
    w2r = w2.real
    w2i = w2.imag
    w3 = -w1
    w3r = w3.real
    w3i = w3.imag

    return [(w1r, w1i), (w2r, w2i), (w3r, w3i)]

rot_imag(3, 45) 

#define quarternion rotation for the generated components based on the plane 


def rot_quart(plane, n=0, w=[]): #tupule list of 2 or 3 
    deg = num.radians(n)/2 #half of any defined theta is used 
    i = tensor([1,0,0])
    j = tensor([0,-1,0]) #define counterclockwise correction 
    k = tensor([0,0,1])

    #define rotation quarternion components 
    if plane == 'xy':
        axis_vector = k
    elif plane == 'yz':
        axis_vector = j 
    elif plane == 'xz':
        axis_vector = j
    else:
        raise ValueError("Axis must be 'xy', 'yz', or 'xz'.")

    axis_vectpr = axis_vector / num.linalg.norm(axis_vector) #generate unit vector from vector object 

    #a unit quarternion has the form cos(n/2) + sin(n/2)*quart  where cos is the scalar part and sin is the rotational axis 
    #first generate the appropriate rotation quarternion components 

    qw = num.cos(deg) #scalar part of quarternion 
    qx, qy, qz = axis_vector * num.sin(deg) #axis part of quarternion  
    qA = quaternion.quaternion(qw, qx, qy, qz)

    qB = [] #create quarternions out of the vectors in w, either 1 or 3 coordinates 

    for x,y,z in w: #list of tupules  passed specifically in the form x,y,z 
        qb = quaternion.quaternion(0, x, y, z)
        qB.append(qb)
    qneg = qA.conjugate()

    rotated = []
    for quart in qB: #for each array present 
        rot = qA*quart*qneg 
        rotated.append(rot)

    return rotated #return the rotated vectors according to the passed plane and degrees provided 

