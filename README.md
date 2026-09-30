# MolOrient

In computational chemistry *molecules* are sets of *atoms* where each atom is
often represented by a point in a three-dimensional Cartesian space. The
coordinates for each atom ultimately depend on the definition of the axis
system, including the distance and relative orientation with respect to the
molecule. For many problems, what matters is not the absolute Cartesian
coordinates of the atoms, but rather their relative positions, i.e., the
distances and angles among them. Since the relative positions of the atoms
are invariant to a translation and rotation of the axis system, this means that
Cartesian coordinates do **NOT** define a unique representation of the molecule.

The lack of uniqueness, complicates:

- using molecules as data descriptors such as in machine learning
- uniquely defining integration grids for density functional theory
- standardizing direction-dependent properties such as multipole moments,
  polarizabilities, etc.
- assigning symmetry labels
- reusing/restarting calculations (need to know how the intermediate quantities
  are oriented)
- comparing structures

One solution is to standardize the axis system. **MolOrient** is a Python
package for aligning molecules with a standardized axis systems.

## Features

- Completely reproducible orientation standardization. MolOrient relies on
  arbitrary-length fixed-precision arithmetic to ensure results contain a user-
  specified number of significant figures.
- Built-in support for Standard Nuclear Orientation, as well as NWChem's
  standard orientation.
- Easy to add other standards.

## Installation

MolOrient is on the Python Package Index, so installation is straightforward.

```Bash
pip install molorient
```

## Example Usage

### Standard Nuclear Orientation

MolOrient reads a molecule from an `.xyz` file and does not accept 
other chemical file formats. For example, `water.xyz`:

```text
3
Water
O	0.0000000	0.0000000	0.1168520
H	0.0000000	0.7546680	-0.4674100
H	0.0000000	-0.7546680	-0.4674100
```

To put it in Standard Nuclear orientation, run:

```Bash
molorient water.xyz
```

MolOrient will then prompt the desired number of significant figures to use, 
then writes the standardized geometry to `water_standardized.xyz` in the same
directory. With 6 significant figures, the result is:

```text
3
Standardized by molorient
H -0.754668 -0.467410 0.000000
H 0.754668 -0.467410 0.000000
O 0.000000 0.116852 0.000000
```

### NWChem Orientation

To use NWChem's orientation, pass `nwchem` as the second argument:

```Bash
molorient water.xyz nwchem
```

MolOrient will also prompt for the number of significant figures and write the 
new file as `water_standardized.xyz`. With 6 significant figures, the result is:

```text
3
Standardized geometry by molorient (NWChem orientation)
O 0.000000 0.000000 0.116852
H -0.754668 0.000000 -0.467410
H 0.754668 0.000000 -0.467410
```

## Resources

- TODO: link to documentation
- [Contributor Guidelines](.github/CONTRIBUTING.md)
- [Code of Conduct](.github/CODE_OF_CONDUCT.md)

## Acknowledgments

Work at the Ames National Laboratory was supported by the U.S. Department of
Energy Office of Science, Science Undergraduate Laboratory Internships (SULI)
program under its contract with Iowa State University, Contract No.
DE-AC02-07CH11358. John Lewis is grateful to the DOE for the assistantship and
for the opportunity to participate in the SULI program.
