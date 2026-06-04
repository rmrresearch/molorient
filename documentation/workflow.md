```mermaid
flowchart LR
    INPUT[Input Geometry] --> WRAPPER[Wrapper] 
    WRAPPER --> OUTPUT[Standardized Geometry] 
    subgraph WRAPPER
        DEC[Decimal Module Function Input: Input Geometry Output: Decimal Geometry] 
        TRANS[Origin Translation Function Input: Decimal Geometry Output: Translated Geometry] 
        ROT[Axis Standardization Function Input: Translated Geometry Output: Rotated Geometry] 
        SORT[Sort Atoms Function Input: Rotated Geometry Output: Sorted Geometry] 
        DEC --> TRANS --> ROT--> SORT
    
    end
