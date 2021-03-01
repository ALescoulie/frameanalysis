from MDAnalysis.analysis.base import AnalysisBase
from singleframe/framefuncs import *
from selections/selectionfuncs import *


class FrameAnalysis(AnalysisBase):
    def __init__(self, universe, res_pairs, func_set):
        super(FrameAnalysis, self).__init__(universe.trajectory)
        self._unv = universe        # Universe being analyzed
        self._rpl = res_pairs      # First list of int residue numbers
        self._fxn = func_set        # Bool values enabling functions

    def _prepare(self):
        """ Builds dictionary to store function outputs computed for inputted residue pairs.
            Keys are string of the respective residue numbers punctuated by a hyphen ex: '121-145'
            dist_calc outputs 2d numpy array with dimensions equal to atoms of each residue.
            Each frame will output an array of the distances of calculated for the given reside pairs"""
        # Distance array
        if self._fxn[0] is True:
            self.res_dists, self.res_keys = build_reslist_dict(self._rpl)

        # Distance between alpha carbons
        if self._fxn[1] is True:
            self.ca_dists, self.ca_keys = build_reslist_dict(self._rpl)

        # Distance between resid center of mass
        if self._fxn[2] is True:
            self.cm_dists, self.cm_keys = build_reslist_dict(self._rpl)

        # Distance between resid center of geometry
        if self._fxn[3] is True:
            self.cg_dists, self.cg_keys = build_reslist_dict(self._rpl)

    def _single_frame(self):
        # Only running functions specified in input file
        if self._fxn[0] is True:
            # Iterating through index pairs and returning distance array to dictionary
            for key in range(len(self.res_keys)):
                self.cord1 = save_resnum_coords(self._rpl[key][0], self._unv)
                self.cord2 = save_resnum_coords(self._rpl[key][1], self._unv)
                self.res_dists[self.res_keys[key]].dist_calc(self.cord1, self.cord2)

        if self._fxn[1] is True:
            # Iterating through index pairs and returning alpha carbon distance
            for key in range(len(self.ca_keys)):
                self.cord1 = save_ca_coords(self._rpl[key][0], self._unv)
                self.cord2 = save_ca_coords(self._rpl[key][1], self._unv)
                self.ca_dists[self.res_keys[key]].append(atom_dist(self.cord1, self.cord2))

        if self._fxn[2] is True:
            # Iterating through index pairs and returning resid center of mass distance
            for key in range(len(self.ca_keys)):
                self.cord1 = save_resid_cm(self._rpl[key][0], self._unv)
                self.cord2 = save_resid_cm(self._rpl[key][1], self._unv)
                self.ca_dists[self.res_keys[key]].append(atom_dist(self.cord1, self.cord2))

        if self._fxn[3] is True:
            # Iterating through index pairs and returning residue center of geometry distance
            for key in range(len(self.ca_keys)):
                self.cord1 = save_resid_cg(self._rpl[key][0], self._unv)
                self.cord2 = save_resid_cg(self._rpl[key][1], self._unv)
                self.ca_dists[self.res_keys[key]].append(atom_dist(self.cord1, self.cord2))
