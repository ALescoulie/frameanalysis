# Trajectory Frame Analysis
Trajectory analyzer build in [MDAnalysis](https://www.mdanalysis.org) for measuring residue pair distances over a simulation. The script accepts an input file in terminal to define its inputs, analysis settings, and residue pairs. 
## Dependencies
### 1. [MDAnalysis](https://www.mdanalysis.org/pages/installation_quick_start/)
Responsible for the selection of residue atoms and coordinates and analysis templates.

## Use
Create an input.in file based on the template in the repository. Residue pairs are separated by a space, and all other values are separated by an '=' sign. In '3:RUN_SETTINGS' booleans must be capitalized as is seen in python ex 'False'. frameanalysis.py references the functions each on their own files in the analysisclass folder, while res_dist.py contains all functions on a single script for greater portability.

To run program
```bash
python res_dist.py input.in
```
