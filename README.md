# LinearAlgebra
A header only Linear Algebra library implemented in modern C++.  Constexpr compatible and tested in a compile time context.

## Contents
### Types
- **Vec**<Size, Type>
- **Pos**<Size, Type>
- **Ray**<Size, Type>
- **Mat**<Size, Size, Type>
- **Sphere3D**\<Type>
- **Triangle3D**\<Type>

### Methods
- **Vec**
    - getLength()
    - getLengthSquared()
    - getNormalised()
- **Ray**
    - getOrigin()
    - getDirection()
    - setOrigin()
    - setDirection() [*ensures that direction is normalised*]
    - getPointAlongRay(**Scalar**) -> **Pos**

### Free Functions
- **Vec**
    - dotProduct(**Vec**, **Vec**) -> **Scalar** 
    - dotProduct(**Mat**, **Vec**) -> **Vec**
    - dotProduct(**Vec**, **Mat**) -> **Vec** [*mathematically not true but allows for preference*]
    - getRotatedVec3(**Vec**, **Scalar**, **Scalar**, **Scalar**) -> **Vec**
    - getReflected(**Vec**, **Vec**) -> **Vec**
    - getRefracted(**Vec**, **Vec**, **Scalar**) -> **Vec**
- **Pos**
    - distance(**Pos**, **Pos**) -> **Scalar**
- **Mat**
    - dotProduct(**Mat**, **Mat**) -> **Mat**
    - transpose(**Mat**) -> **Mat**
    - deteminant(**Mat**) -> **Scalar**
    - getRotationMat3x3(**Scalar**, **Scalar**, **Scalar**) -> **Mat**
- **Sphere3D**
    - getNormalVec(**Pos**, **Sphere3D**) -> **Vec**
    - intersectionDist(**Ray**, **Sphere3D**) -> **Scalar**
- ***Misc***
    - toRadians(**Scalar**) -> **Scalar** 

### To Do
- Triangle Intersections
- Determinant for any size matrix