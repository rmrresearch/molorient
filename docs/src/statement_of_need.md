# Statement of Need

In computational chemistry, atoms are given in a structural file, such as:
Cartesian (.xyz), PDB, (.pdb), or MDL Molfile (.mol). All computational
chemistry packages require a structural file input to perform a simulation.
However, there is not a universal orientation that these structural files are
forced to follow, leading to an ambiguous orientation. This can become an issue
when analyzing and comparing results of a simulation. MolOrient attempts to fill
this need by offering a reproducible method of standardization utilizing
arbitrary-length fixed-precision mathematics to avoid floating-point errors. To
the best of our knowledge, there is currently no math package designed with
native support for Decimal objects, let alone a geometry standardization package
built around them.

MolOrient offers a reusable infrastructure to accommodate other orientation
conventions. MolOrient's default orientation convention is the Standard Nuclear
Orientation, with support for NWChem's orientation convention. NWChem's addition
was simple: a single file, nwchem_orientation.py was created, using entirely
Decimal arithmetic and reusing MolOrient functions. A small change to cli.py
was made to switch between MolOrient's default orientation and NWChem's in the
command line. A user can add an orientation convention in the same manner. Say a
user wishes to add support for GAMESS's orientation. They will create a file,
gamess_orientation.py, that reuses MolOrient's eigensolver, trigonometric
functions, and linear algebra functions to replicate GAMESS's orientation that
uses Decimal arithmetic. They will make 3 edits to cli.py: a `_MODES` entry, a
`--gamess` flag into `build_parser()`, and a branch in the mode selection in
`main()`. The user then runs `molorient water.xyz --gamess` and MolOrient will
standardize the molecule's orientation using GAMESS's convention. This can be
done for any orientation convention.
