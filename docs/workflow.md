```mermaid
flowchart LR

INPUT[Input Geometry] --> XYZ[XYZ Parser\nInput: Float Geometry\n Output: Decimal Geometry]
INPUT[Input Geometry] --> PREC[Set Precision]

XYZ --> ORIENTSYSTEM[ORIENT SYSTEM]

ORIENTSYSTEM --> OUTPUT[Standardized Geometry]

subgraph ORIENTSYSTEM[Orient System]
    TRANS[Origin Translation Function\nInput: Decimal Geometry\nOutput: Translated Geometry]
    ROT[Axis Standardization Function\nInput: Translated Geometry\nOutput: Rotated Geometry]
    SORT[Sort Atoms Function\nInput: Rotated Geometry\nOutput: Sorted Geometry]
    TRANS --> ROT --> SORT
end

subgraph HELPERS[Supporting Functions]
    TRIG[Trig Helper Functions]
    DIAG[Eigensolver]
    TRIG --> DIAG
end

DIAG --> ROT