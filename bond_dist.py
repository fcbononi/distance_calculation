#!/usr/bin/env python3
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
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)
    diff = coords2[:, np.newaxis] - coords1[np.newaxis, :]
    dist = np.einsum('ijk->ij', diff**2)**0.5
    index = np.argmin(dist)
    return dist
#
#def remove_atom(atom_list, atom_index):
#    atom_list.pop(atom_index)
#    return atom_list
#
#def n_closest_waters(coordinates, atom, oxygens, n):
#    hydroxys = []
#    for i in range(n):
#        index = find_closest_atom(atom, oxygens)
#        closest_oxygen = oxygens[index]
#        if closest_oxygen in coordinates:
#            oxygen_index = coordinates.index(closest_oxygen)
#        OH = coordinates[oxygen_index]
#        HO = coordinates[oxygen_index+1]
#        hydroxy = [OH, HO]
#        hydroxys.append(hydroxy)
#        oxygens = remove_atom(oxygens, index)
#    return hydroxys
#
#def add_hydroxy(coordinates, atoms, hydroxys):
#    for hydroxy in hydroxys:
#        atoms.append('OH')
#        atoms.append('HO')
#        coordinates.append(hydroxy[0])
#        coordinates.append(hydroxy[1])
#        coordinates.append(hydroxy[2])
#    return coordinates, atoms
#
#
atoms, coordinates = read_xyz(filename, 1)
if os.path.isfile(output):
    os.remove(output)
for frame in range(len(atoms)):
   bond_dist=[]
   timestep=[]
   OH = find_atoms(atoms[frame], coordinates[frame], 'OH')
   HO = find_atoms(atoms[frame], coordinates[frame], 'HO')
   index = find_closest_atom(OH, HO)
   # print OH, HO, index
   timestep.append(frame)
   bond_dist.append(index)
#    #print timestep, np.array(bond_dist)
   time_series=[timestep, bond_dist]
   with open(output, "a") as file:
      for x in zip(*time_series):
         file.write("{0}\t{1}\n".format(*x))

with open(output, 'r') as my_file:
	text = my_file.read()
	text = text.replace("[", "")
	text = text.replace("]", "")

with open(output, 'w') as my_file:
    my_file.write(text)

print("distance calculation...DONE")
