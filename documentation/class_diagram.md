```mermaid
classDiagram

class Molecule {
    +atoms: Atom[]
}

class Atom {
    +element: str
    +x: Decimal
    +y: Decimal
    +z: Decimal
    +mass: Decimal
    +charge: Decimal
}

Molecule "1" *-- "1..*" Atom