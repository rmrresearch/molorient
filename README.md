# MolOrient

MolOrient is a molecular geometry standardization Python package utilizing Python's arbitrary-length fixed-precision objects, Decimals.

## Statement of Need
In computational chemistry, atoms are given in a structural file, such as: Cartestian (.xyz), PDB, (.pdb), or MDL Molfile (.mol). All computational chemistry packages require a structural file input to perform a simulation. However, there is not a universal orientation that these structural files are forced to follow, leading to an ambiguous orientation. This can become an issue when analyzing and comparing results of a simulation. MolOrient attempts to fill this need by offering a reproducible method of standardization utilizing arbitrary-length fixed-precision mathematics to avoid floating-point errors. To the best of our knowledge, there is currently no math package designed with native support for Decimal objects, let alone a geometry standardization package built around them.

MolOrient offers a reusable infrastructure to accommodate other orientation conventions. MolOrient's default orientation convention is the Standard Nuclear Orientation, with support for NWChem's orientation convention. NWChem's addition was simple: a single file, nwchem_orientation.py was created, using entirely Decimal arithmetic and reusing MolOrient functions. A small change to main.py was made to switch between MolOrient's default orientation and NWChem's in the command line. A user can add an orientation convention in the same manner. Say a user wishes to add support for GAMESS's orientation. They will create a file, gamess_orientation.py, that reuses MolOrient's eigensolver, trigonometric functions, and linear algebra functions to replicate GAMESS's orientation that uses Decimal arithmetic. They will then add a GAMESS option to `_MODES` on main.py. In the command line, the user can then add `gamess` after the .xyz file of choice and MolOrient will standardize the molecule's orientation using GAMESS's convention. This can be done for any orientation convention.

See [Background](./docs/background.md) for more details.

## Features
- Completely reproducible orientation standardization.
- Select number of significant figures.
- Every operation uses arbitrary-legnth fixed-precision.
- From scratch linear algebra, transcendental functions, and 3x3 symmetric eigensolver.
- Support for NWChem's orientation.

## Installation
MolOrient is on the Python Package Index, so installation is straightforward.

```Bash
pip install molorient
```

## Contributing

- [Contributor Guidelines](./docs/contributing.md)
- [Code of Conduct](./docs/code_of_conduct.md)

## Acknowledgments

Work at the Ames National Laboratory was supported by the U.S. Department of Energy Office of Science, Science Undergraduate Laboratory Internships (SULI) program under its contract with Iowa State University, Contract No. DE-AC02-07CH11358. John Lewis is grateful to the DOE for the assistantship and opportunity to participate in the SULI program.