# Written by David Weber
# dsw7@sfu.ca

"""
This script determines bridge/aromatic chain superimposability for ALL PDB entries
(1) We import a structure
(2) We find Met-aromatic interactions 
(3) From Met-aromatic interactions we find 2-bridges
(4) We then find closely spaced Tyr/Trp residues
(5) We apply network theory approaches to determine superimposability

IMPORTANT: Here I only search for bridges of the form:
           A::MET::B
Where (A != PHE) & (B != PHE)
I do not however work this condition into the met_aromatic function.
Instead I filter out Phe containg pairs in the BP list at the end of the
try block. I do this as the met_aromatic function is now a literature
mathematical function, and modifying it would take away from functional
elegance. Moreover, reformatting the function could yield to all sorts
of downstream errors.
"""

"""
Geometry S counts (2,611 proteins) (total pairs = 12,685)
Type    Count       Proportion
NR   -> 10558   ->  0.83
PM   -> 848     ->  0.07
DS   -> 549     ->  0.04
IS   -> 730     ->  0.06


Geometry S counts (20,784 proteins) (total pairs = 119,038)
Type    Count       Proportion
NR   -> 99,166  ->  0.83
2C   -> 7,249   ->  0.06     // note I have renamed PM to 2C
DS   -> 5,606   ->  0.05
IS   -> 7,107   ->  0.06

"""

# ------------------------------------------------------------------------------
# dependencies

from sys                    import path as PATH
from os                     import path, getcwd, remove
PATH.append(path.join(getcwd(), 'lib'))

import networkx as nx
from traceback              import print_exc
from csv                    import DictWriter, writer, QUOTE_MINIMAL
from ma_lowlevel            import met_aromatic
from nn                     import get_nn
from PDB_filegetter         import PDBFile
from EC_classifier          import get_EC_classifier
from itertools              import chain
from get_organism_from_file import get_organism

# ------------------------------------------------------------------------------
# hard code + some trivial functions

CHAIN = 'A'
CUTOFF = 6.0
CUTOFF_FOR_NN = 7.4
ANGLE = 109.5
MODEL = 'cp'
BRIDGE_ORDER = 2
HEADER = ['AROMATIC', 'ARO POS', 'MET', 'MET POS', 'NORM', 'MET-THETA', 'MET-PHI']

# TODO: might not need this?
def write_to_csv(data, path, header):
    """ A function for writing to .csv... """
    with open(path, mode='w', newline='') as f:
        obj_header = DictWriter(f, fieldnames=header)
        obj_header.writeheader()
        obj_writer = writer(f, delimiter=',', quotechar='"', quoting=QUOTE_MINIMAL)
        for line in data:
            obj_writer.writerow(line)
 
# TODO: might not need this?           
def reduce_output(input_data, round_to=3):
    """ rounds output data from met_aromatic()
        to round_to number of decimal places
    """
    for entry in input_data:
        for idx in range(4, 7):
            entry[idx] = round(entry[idx], round_to)    
    return input_data

def bridge_compressor(met_list, pair_list):
    """
        [
         ('TRP340', 'MET390'),     
         ('PHE393', 'MET390')
        ] 
            ->
        'TRP340-PHE393'
    """
    
    f = [[p for p in pair_list if m in p] for m in met_list]  # isolate MET groups
    f = [list(chain(*p)) for p in f]  # flatten into one data structure
    f = [[i for i in p if 'MET' not in i] for p in f]  # strip out MET stuff  
    f = [tuple(i) for i in f]  
    return f
    
NR_CONDITION = 'NR : {} : {}'
PM_CONDITION = 'PM : {} : {}'
DS_CONDITION = 'DS : {} : {}'
IS_CONDITION = 'IS : {} : {}'
UK_CONDITION = 'UC'  # unknown condition

def compare_networks(chain_list, bridge_list):
    """
    Compares the relationship between a 
    set of disconnected graph components
    
    list(nx.connected_components(G1))
    
    Parameters
    ----------
    chain_list : a list of closely spaced aromatic chains
    bridge_list : a list A / B bridges where we have A::MET::B
            
    chain_list, bridge_list are of type list(nx.connected_components(G))
    where G is an nx.Graph() object

    Returns
    -------
    A list of any of NR_CONDITION, PM_CONDITION,
                     DS_CONDITION, IS_CONDITION,
                     UK_CONDITION
    """
    # beautiful, syntactic Python logic
    list_cond = []
    for aromatic_chain in chain_list:
        for bridge in bridge_list:    
            if aromatic_chain & bridge == set():                 # NR
                list_cond.append(NR_CONDITION.format(bridge, aromatic_chain))
            elif aromatic_chain - bridge == set():               # PM
                list_cond.append(PM_CONDITION.format(bridge, aromatic_chain))
            elif aromatic_chain & bridge == bridge:              # DS
                list_cond.append(DS_CONDITION.format(bridge, aromatic_chain))
            elif next(iter(aromatic_chain & bridge)) in bridge:  # IS
                list_cond.append(IS_CONDITION.format(bridge, aromatic_chain))
            else:
                list_cond.append(UK_CONDITION)
    return list_cond


# ------------------------------------------------------------------------------

# create a results file if one doesn't exist
if path.exists('results.txt') == False:
    f = open('results.txt', 'w')
    f.write('CODE : INTERACTION TYPE : {BRIDGE} : {CHAIN} : EC CLASSIFIER : ORGANISM \n')
    f.close()

# read in delimiter list
with open('delim.txt', mode='r') as f: CODES_ALL = f.readlines()

for count, CODE in enumerate(CODES_ALL):
    CODE = CODE.strip('\n')
    
    if '$' in CODE:  # for testing only
        continue
    
    print('-' * 25)
    print('{}. PDB CODE: {}'.format(count + 1, CODE))
    try:
        file_pdb = PDBFile(CODE)
        path_to_file = file_pdb.fetch_from_PDB()
        
        # first find bridging interactions using raw met aromatic algorithm
        data_retrieved = met_aromatic(CHAIN=CHAIN, 
                                      CUTOFF=CUTOFF, 
                                      ANGLE=ANGLE, 
                                      MODEL=MODEL,
                                      filepath=path_to_file)
        
        # then attempt to find EC classifier
        classifier_EC = get_EC_classifier(path_to_file)
        
        # then attempt to find organism
        list_org = get_organism(path_to_file)
        str_org = list_org[0] if len(list_org) == 1 else ' '.join(list_org)
        
        # then attempt to find all closely spaced Tyr/Trp residues
        NN = get_nn(path_to_file, cutoff=CUTOFF_FOR_NN)
        
        # delete file from pwd if everything was succesfully collected
        file_pdb.clear()
               
        # create met-aro string pairs
        pairs = [(i[0] + i[1], i[2] + i[3]) for i in data_retrieved]    
        
        # remove dupes
        pairs = list(set(pairs))        
        
        # pass to nx
        G_BRDG = nx.Graph()
        G_BRDG.add_edges_from(pairs)
                
        # we need to get met nodes to assess degree of each met
        mets = [m for m in G_BRDG.nodes() if 'MET' in m]
        
        # a degree of 2 means a 2-bridge
        mets = [m for m in mets if G_BRDG.degree(m) == BRIDGE_ORDER]
        
        # 2-bridges -> aromatic pairs
        BP = bridge_compressor(mets, pairs)
        
        # isolate to only (A != PHE) && (B != PHE) A::MET::B bridges
        BP = [i for i in BP if 'PHE' not in ''.join(i)]
        
    except Exception as exception:
        print_exc()
        if path.exists(path_to_file):
            remove(path_to_file)
        continue
            
    else:       
        # skip over any file that does not contain both a bridge and aromatic chain
        if NN == [] or BP == []:
            continue
        else:  # otherwise compare the networks G1 and G2 and append results to text file
            # match thesis / paper naming
            G1 = nx.Graph()
            G2 = nx.Graph()
            
            G1.add_edges_from(NN)
            G2.add_edges_from(BP)
            
            G1_disconnected = list(nx.connected_components(G1))  # chains
            G2_disconnected = list(nx.connected_components(G2))  # bridges   
    
            comparisons = compare_networks(G1_disconnected, G2_disconnected)
            
            with open('results.txt', 'a') as f:
                for line in comparisons:
                    f.write('{} : {} : {} : {} \n'.format(CODE, line, 
                                                          classifier_EC, str_org))
    
print('Done!')
        
        