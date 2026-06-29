```mermaid
classDiagram

class Atom {
    +element: str
    +x: Decimal
    +y: Decimal
    +z: Decimal
    +charge: Decimal
}

class SquareMatrix {
    +elements: List[List[Decimal]]
    +assign(i: int, j: int, value: Decimal) void
    -add(other : SquareMatrix) SquareMatrix
    +multiply(other : SquareMatrix) SquareMatrix
    +multiply(other : Vector) Vector
    #transpose() SquareMatrix
    -scale(scalar : Decimal) SquareMatrix
    -inverse() SquareMatrix 
}

class Vector {
    +elements: List[Decimal]
    +assign(i: int, value: Decimal) void
    -add(other : Vector) Vector
    +dot(other : Vector) Decimal
    -scale(scalar : Decimal) Vector
    +cross(other : Vector) Vector
}

SquareMatrix ..> Vector
Atom ..> Vector
Atom ..> SquareMatrix