import MDAnalysis as mda
import pandas as pd


def build_universe(topology_path, *args):
    universe = mda.Universe(topology_path, args)
    return universe


# Selects residue atoms from an its integer residue number
def select_resnum_atoms(res_num, universe):
    resid_num = 'resid ' + str(res_num)
    atoms = universe.select_atoms(resid_num)
    return atoms


# Selects residue alpha carbon from its integer residue number
def select_resid_ca(res_num, universe):
    selection = 'resid ' + str(res_num) + ' and name CA'
    res_ca = universe.select_atoms(selection)
    return res_ca


def select_atom(atom_num, universe):
    pass


def save_resid_cm(res_num, universe):
    res_atoms = select_resnum_atoms(res_num, universe)
    coords = res_atoms.center_of_mass()
    return coords


def save_resid_cg(res_num, universe):
    res_atoms = select_resnum_atoms(res_num, universe)
    coords = res_atoms.center_of_geometry()
    return coords


# Writes atom group coordinates into numpy array
def save_atom_coords(atom_group):
    coords = atom_group.positions
    return coords


# Selects residue atoms and writes coordinates into numpy array
def save_resnum_coords(res_num, universe):
    res_atoms = select_resnum_atoms(res_num, universe)
    res_coords = save_atom_coords(res_atoms)
    return res_coords


def save_ca_coords(res_num, universe):
    ca = select_resid_ca(res_num, universe)
    ca_coords = save_atom_coords(ca)
    return ca_coords


def build_reslist_dict(res_pair_list):
    res_pairs = {}
    for respair in res_pair_list:
        key_name = str(respair[0]) + '-' + str(respair[1])
        res_pairs[key_name] = []
    res_keys = list(res_pairs.keys())
    return res_pairs, res_keys


def save_xyz_coords(atom_group):
    atom_list = []
    coords_array = []
    for atom in atom_group:
        atom_list.append(atom)
        coords_array.append(atom.positions)

    coords = dict.fromkeys(atom_list)
    for i in range(len(atom_list)):
        coords[atom_list[i]] = coords_array[i]

    xyz = pd.DataFrame(coords)
    return xyz
