#!/usr/bin/env python
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, default='device.pdb', help='name of file')
parser.add_argument('outfile', type=str, default='output.xyz', help='name of output file')
args = parser.parse_args()
filename = args.file
output = args.outfile


def read_xyz(filename, freq):
   """Read filename in XYZ format and return lists of atoms and coordinates.

   If number of coordinates do not agree with the statd number in
   the file it will raise a ValueError.
   """


#xyz file

   Atoms = []
   Coordinates = []

   xyz = open(filename)
   frame = 0
   while True:

       n_atoms = xyz.readline()

       if n_atoms == '':
           break
       else:
           n_atoms = int(n_atoms)
       title = xyz.readline()

       if frame%freq==0:
           atoms, coordinates = read_frame(xyz, n_atoms)
           Coordinates.append(coordinates)
           Atoms.append(atoms)

       else:
           read_frame(xyz, n_atoms)
       frame+=1

   return Atoms, Coordinates

def read_frame(xyz, n_atoms):

   atoms = []
   coordinates = []

   for i in range(n_atoms):
       atom,x,y,z,null = xyz.readline().split()
       coordinates.append([float(x), float(y), float(z)])
       atoms.append(atom)
   return atoms, coordinates
    
def find_atoms(atoms, coordinates, atom_type):
   
    type_coordinates = []
 
    i = 0
    for atom_name in atoms:
        if atom_type == atom_name:
            index = i
            type_coordinates.append(coordinates[index])
        i+=1
 
    return type_coordinates

def find_closest_atom(coords1, coords2):
    
    """ Calculates distances between OH (phenol oxygen) and all OT (water oxygen). """

    coords1 = np.array(coords1)
    coords2 = np.array(coords2)
    diff = coords2[:, np.newaxis] - coords1[np.newaxis, :]
    dist = np.einsum('ijk->ij', diff**2)**0.5
    index = np.argmin(dist)
    return index

def remove_atom(atom_list, atom_index):
    atom_list.pop(atom_index)
    return atom_list

def n_closest_waters(coordinates, atom, oxygens, n):

    """ Finds the n closest oxygen atoms to the OH (phenol group oxygen atom) as specified by 'n'."""

    waters = []
    for i in range(n):
        index = find_closest_atom(atom, oxygens)
        closest_oxygen = oxygens[index]
        if closest_oxygen in coordinates:
            oxygen_index = coordinates.index(closest_oxygen)
        OT = coordinates[oxygen_index]
        HT1 = coordinates[oxygen_index+1]
        HT2 = coordinates[oxygen_index+2]
        water = [OT, HT1, HT2]
        waters.append(water)
        oxygens = remove_atom(oxygens, index)
    return waters

def add_water(coordinates, atoms, waters):

    "Appends the n closest water molecules(OT and corresponding HT)."""

    for water in waters:
        atoms.append('OT')
        atoms.append('HT')
        atoms.append('HT')
        coordinates.append(water[0])
        coordinates.append(water[1])
        coordinates.append(water[2])
    return coordinates, atoms

def write_xyz(filename, atoms, coordinates, frame='no default set'):

    """ Writes a new xyz file containing the water molecules closest to the phenol oxygen.
    It will raise an error if the number of atoms is different than the number of xyz coordinates."""

    if len(atoms) != len(coordinates):
        raise ValueError('Number of atoms is  different than number of positions')

    xyz_file = open(filename,'a')
    xyz_file.write('{}\n'.format(len(atoms)))
    xyz_file.write('frame {}\n'.format(frame))
    for i in range(len(atoms)):
        xyz_file.write('{}\t{}\t{}\t{}\n'.format(atoms[i], coordinates[i][0],
                                  coordinates[i][1], coordinates[i][2]))
    xyz_file.close()

atoms, coordinates = read_xyz(filename, 1)
if os.path.isfile(output):
    os.remove(output)
for frame in range(len(atoms)):
    OT = find_atoms(atoms[frame], coordinates[frame], 'OT')
    OH = find_atoms(atoms[frame], coordinates[frame], 'OH')
    HT = find_atoms(atoms[frame], coordinates[frame], 'HT')
    HO = find_atoms(atoms[frame], coordinates[frame], 'HO')
    final_coordinates = [x for x in coordinates[frame] if x not in OT and x not in HT]
    final_atoms = [x for x in atoms[frame] if x != 'OT' and x not in 'HT']
    waters = n_closest_waters(coordinates[frame], OH, OT, 5)
    final_coordinates, final_atoms = add_water(final_coordinates,final_atoms,waters)
    write_xyz(output, final_atoms, final_coordinates, frame)

