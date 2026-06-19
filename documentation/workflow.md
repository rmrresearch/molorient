```mermaid
flowchart LR

INPUT[Input Geometry] --> WRAPPER[Wrapper]

WRAPPER --> OUTPUT[Standardized Geometry]


subgraph WRAPPER
    DEC[Decimal Module Function\nInput: Input Geometry\nOutput: Decimal Geometry]
    TRANS[Origin Translation Function\nInput: Decimal Geometry\nOutput: Translated Geometry]
    ROT[Axis Standardization Function\nInput: Translated Geometry\nOutput: Rotated Geometry]
    SORT[Sort Atoms Function\nInput: Rotated Geometry\nOutput: Sorted Geometry]
    DEC --> TRANS --> ROT --> SORT
end

subgraph HELPERS[Supporting Functions]
    TRIG[Trig Helper Functions]
    DIAG[Diagonalization Functions]
    TRIG --> DIAG
end

DIAG --> ROT