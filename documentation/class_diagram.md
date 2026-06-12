```mermaid
classDiagram

class Atom {
    +element: str
    +x: Decimal
    +y: Decimal
    +z: Decimal
    +mass: Decimal
    +charge: Decimal
}

class SquareMatrix {
    +elements: List[List[Decimal]]
    +assign(i: int, j: int, value: Decimal) void
    -add(other : SquareMatrix) SquareMatrix
    +multiply(other : SquareMatrix) SquareMatrix
    +multiply(other : Vector) Vector
    -transpose() SquareMatrix
    -scale(scalar : Decimal) SquareMatrix
}

class Vector {
    +elements: List[Decimal]
    +assign(i: int, value: Decimal) void
    -add(other : Vector) Vector
    -dot(other : Vector) Decimal
    +multiply(other : SquareMatrix) Vector
    -scale(scalar : Decimal) Vector
}

Vector ..> SquareMatrix
SquareMatrix ..> Vector
Atom ..> Vector
Atom ..> SquareMatrix