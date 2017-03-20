Here's the process:

1) Run prep_NBS_designmats.bash, with no arguments, in order to get the initial subject lists and corresponding design matrices. This will also generate separate lists/designmatrices by gender and diagnostic status. Calls the python script which outputs basic stats to command window. Some plots will also be saved.

2) Run NBSprep_make3Dconmat.m in matlab. This will take a subject list in order to generate 3D stacks of connection matrices.

3) (optional) Run COMprep_make3dmats.m to convert the connection matrices into communicability matrices.

4) Run NBS analyses either through the GUI (just type 'NBS' in matlab), through NBS2caret.m, or through NBS_parameterspace.m. Choose from other scripts to run other multivariate statistical analyses - e.g. elastic nets in 'multiple regression.'
