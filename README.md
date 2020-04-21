# MCPB_Glycam_merge
Python script to merge an MCPB leap input script, an associated MCPB pdb file, and a Glycam pdb file.

# Requirements
* Python3 or greater

# Installation
Download and place MCPB_Glycam_merge.py in a directory included in your $PATH variable.

# Arguments
Execute `MCPB_Glycam_merge.py -h` for a brief description of arguments.

# Instructions
First model the metal center with MCPB from [AmberTools Leap](https://ambermd.org/AmberTools.ph).  Then, attach glycans to the pdb file using [Glycam Web](http://glycam.org/). Download "structure_AMBER.pdb" for use with `MCPB_Glycam_merge.py`.  

# Additional information
`MCPB_Glycam_merge.py` automatically reincludes residues that Glycam Web removes.  Residues and CONECT records are reindexed appropriately.  N-glycan sites are automatically renamed to NLN.
