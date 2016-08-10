# NC-membrane
Code to study the collective interactions of an NC with a cell membrane
##08/10/16: Minor bug fixes
1) ```interaction_parameter.cpp``` -- the bondlength for a non-interacting receptor-ligand pairs is set to 0 instead of this->L in the previous version
2) ```module_mcsmoves.f``` -- A error in using the current receptor locations has been fixed in ```Rosenbluth_Sampling_Bound_Antigens(antigen_no)```
