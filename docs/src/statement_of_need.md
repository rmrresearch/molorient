# Statement of Need

In chemistry, the Born-Oppenheimer approximation enables modeling molecular
properties as a function of nuclear positions. In practice, the property usually
does not depend on the absolute nuclear positions (i.e., the x, y, and z
coordinates of each nucleus), but instead depends on the geometry of the
molecule. Here *geometry* refers to the relative orientation of the nuclei to
one another, e.g., distances between or angles among them.

For a variety of reasons, paramount of which is simplicity, the predominant
standard for describing a molecule remains a set of pairs where each pair
contains the atomic number (or atomic symbol) and Cartesian coordinates for one
nucleus. Since the numerical values of the Cartesian coordinates describe
offsets relative to a somewhat arbitrary axis system, the same geometry will
take on different Cartesian coordinates if a different axis system is adopted.

To establish whether the molecular geometries implied by two different sets of
Cartesian coordinates are identical it suffices to *align* the coordinates. To
align the coordinates: choose a reference axis system, reexpress the
coordinates in the reference system, and then compare Cartesian coordinates.
The comparison is usually done by taking the root-mean square deviation (RMSD)
of the coordinates. An RMSD of 0 means the coordinates (and hence the
geometries) are identical. For RMSD>0, the larger the RMSD, the larger the
difference in the geometries.

While many existing software packages provide such a feature, the alignment is
almost exclusively done with floating-point arithmetic. This means that not
only do the resulting coordinates suffer from floating-point round-off errors,
but the resulting coordinate values will also include a dependence on how the
alignment is implemented (e.g., since floating-point arithmetic is not
associative, the order of operations matters). This makes it difficult for
different implementations to obtain the same RMSD, and even the same
implementation may compute a non-zero RMSD for identical geometries.

The increasing interest in data science means there is an increasing interest
in creating tables of molecular properties. Conceptually these tables are key-
value arrays where the keys are the Cartesian coordinates of the molecule and
the value is the property. Assuming aligned molecular coordinates, the cost of
determining if the table contains a molecule is worst-case-scenario proportional
to the number of entries in the table times the number of atoms in the molecule.
For large tables and molecules, this can be prohibitive (N.B. databases such as
PubChem contain 100M+ molecules). While the cost can be brought down to
logarithmic in the number of table entries by sorting the molecules, the real
gain in lookup speed (and storage size) comes from being able to hash molecules.

This is the motivation for MolOrient, implement molecular alignment in a manner
that results in stable hashes. In particular, if two different sets of Cartesian
coordinates imply identical geometries, than their aligned coordinates should
be identical to within a user-specified tolerance. To accomplish this:

- MolOrient stores all Cartesian Coordinates as `Decimal` objects (from
  Python's `decimal` package).
- All linear algebra routines leverage `Decimal` objects.
- All operations maintain sufficient significant figures to ensure rounding
  does not affect the final result to within the user-specified tolerance.

By using the `decimal` package MolOrient is able to avoid the pitfalls of
floating-point arithmetic (`decimal` implements traditional fixed-precision
arithmetic). Moreover, if the same axis standard is leveraged by another
implementation (and that implementations also uses fixed-point arithmetic) the
results of MolOrient and the other implementation should be identical.

Finally, we want MolOrient to be a community resource. Whereas we are fine
adopting standard nuclear orientation for our standardization needs, we realize
that other communities may want to adopt other axis standards. Thus MolOrient
was designed to be modular. To illustrate this MolOrient also ships with
NWChem's axis standards. New standards can be added by creating python module
to store the standard, using the existing MolOrient arithmetic functions to
implement the standard, and registering the standard with the command line
interface.
