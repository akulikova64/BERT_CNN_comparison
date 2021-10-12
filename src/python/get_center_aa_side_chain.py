# this code will find the center of amino acid side chain. 

from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one
import pandas as pd
import os
import sys
import csv

def process_residue(residue):
    '''
    Processes a single residue to determine the coordinates of the alpha-carbon
    and the sidechain center-of-mass. Also checks for missing atoms in a
    residue.
    '''
    output_dict = {}
    # Convert three letter amino acid to one letter
    try:
        output_dict['amino_acid'] = three_to_one(residue.resname)
    except KeyError:
        return 'Error'
    
    # Grab residue number AND any insertion site labeling (11A, 11B, etc.)
    output_dict['residue_number'] = str(residue.get_id()[1]) + residue.get_id()[2].strip()
    
    # Straightforward, grab the chain ID
    output_dict['chain'] = residue.get_full_id()[2]
    
    #Just checking correctness
    try:
        int(output_dict['residue_number'])
    except:
        return 'Error, residue number is... not an int?'
    
    # Coordinates of all sidechain atoms in this residue
    sidechain_coords = []
    atoms_seen = []
    for atom in residue:
        atoms_seen.append(atom.name)
        if atom.name == 'CA':
            # Save alpha-carbon coordinates separately
            output_dict['CA_coords'] = atom.get_coord()
            
            # If it's Glycine... call that the CB to
            if residue.resname == 'GLY':
                output_dict['CB_coords'] = atom.get_coord()
                
        if atom.name == 'CB':
            # Save beta-carbon coordinates
            output_dict['CB_coords'] = atom.get_coord()

        #Ignore the backbone and add the coordinates of all side-chain atoms
        if atom.name not in ['C', 'CA', 'O', 'N']:
            # Must be a sidechain atom...
            sidechain_coords.append(atom.get_coord())
            
    if 'CA_coords' not in output_dict or 'CB_coords' not in output_dict:
        return 'Error, this residue is missing a CA and/or a CB atom'
    
    for mainchain_atom in ['N', 'C', 'O']:
        # Warn about any missing mainchain atoms
        if mainchain_atom not in atoms_seen:
            return 'Error, this residue has a strange backbone'

    if len(sidechain_coords) == 0:
        ###Treat glycine separately. Normally CA is ignored in the side-chain atom calculations
        ###but for glycine it's all the information that we have.
        if output_dict['amino_acid'] == 'G':
            sidechain_coords.append(output_dict['CA_coords'])
        else:
            return 'Error'
        
    # Calculate side chain geometric center
    output_dict['SCcenter_coords'] = sum(sidechain_coords)/len(sidechain_coords)
    
    return output_dict

def get_residues_df(pdb_id, path):
  pdb_file = path + pdb_id 
  structure = PDBParser().get_structure(pdb_id, pdb_file)
  
  atoms = structure.get_atoms()

  ###Get residue coordinates
  temp_listy = []
  cols = ['residue_number', 'amino_acid', 'chain', 'CA_coords', 'CB_coords', 'SCcenter_coords']
  for residue in structure.get_residues():
      if is_aa(residue):
          temp_dict = process_residue(residue)
          if type(temp_dict) == str:
            continue
          temp_listy.append([temp_dict[col] for col in cols])
      else:
          #print('Problem with this residue: {}'.format(residue))
          qwe = 1
          
  residues_df = pd.DataFrame(temp_listy, columns=cols)
  return residues_df

#---------------main-----------------------------------------------------------

input_path = "../../data/pdb/"

fileList = os.listdir(input_path)

for file in fileList:
    pdb_id = file
    output_path = "../../output/side_chain_centers/SCC_" + file + ".csv" # make "side_chain_centers/" folder first.
    with open(output_path, "w", newline='\n', encoding='utf-8') as CSV_file:

        #writer = csv.writer(CSV_file)
        #writer.writerow(['col_num', 'residue_number', 'amino_acid', 'chain', 'CA_coords', 'CB_coords', 'SCcenter_coords', 'gene'])

    
        df = get_residues_df(pdb_id, input_path) # parameters: (pdb_id, path)
        gene = file[0:4]
        gene_col = [gene] * len(df)
        df['gene'] = gene_col
        #print(df[0:10])
        #sys.exit()

        df.to_csv(path_or_buf = output_path, sep = ',', header = True)


