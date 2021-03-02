import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.distances import dist


def dist_calc(coordinates_a, coordinates_b):
    group_dist = distance_array(coordinates_a, coordinates_b)
    return group_dist


def atom_dist(cooridnates_a, cooridnates_b):
    distance = dist(cooridnates_a, cooridnates_b)
    return distance

