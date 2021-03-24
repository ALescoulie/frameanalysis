from analysisclass.selection.selectionfuncs import *
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.distances import dist
import numpy as np


def dist_calc(coordinates_a, coordinates_b):
    group_dist = distance_array(coordinates_a, coordinates_b)
    return group_dist


def atom_dist(cooridnates_a, cooridnates_b):
    distance = dist(cooridnates_a, cooridnates_b)
    return distance


def dist_vec(atom_group_a, atom_group_b):
    distances = distance_array(atom_group_a, atom_group_b)
    dist_vector = distances.flatten()
    return dist_vector
