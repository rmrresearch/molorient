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
    -multiply(other : SquareMatrix) SquareMatrix
    -transpose() SquareMatrix
    -scale(scalar : Decimal) SquareMatrix
}