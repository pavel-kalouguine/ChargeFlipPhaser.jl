module Icosahedron

using SpaceGroups
import StaticArrays: SMatrix, SVector, @SMatrix, @SVector

export PI, PIh, PI_n, PIh_n, τ, epar;

# Golden mean
const τ=(sqrt(5)+1)/2

# Parallel space
const epar=(@SMatrix [1 τ 0 -1 τ 0; τ 0 1 τ 0 -1; 0 1 τ 0 -1 τ])/sqrt(1+τ^2)
# Perp space
const eper=(@SMatrix [-τ 1 0 τ 1 0; 1 0 -τ 1 0 τ; 0 -τ 1 0 τ 1])/sqrt(1+τ^2)

# Three-fold rotation
const r3=SMatrix{6,6}(
   [0 0 1 0 0 0;
    1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0 0 0 1;
    0 0 0 1 0 0;
    0 0 0 0 1 0])

# Five fold rotation
const r5=SMatrix{6,6}(
   [1 0 0 0 0 0;
    0 0 0 0 1 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 0 -1;
    0 0 0 -1 0 0])

# Central symmetry
const c=SMatrix{6,6}(
   [-1 0 0 0 0 0;
    0 -1 0 0 0 0;
    0 0 -1 0 0 0;
    0 0 0 -1 0 0;
    0 0 0 0 -1 0;
    0 0 0 0 0 -1])

# Translations
# Zero translation
const z=zeros(SVector{6, Rational{Int}})

# Translations for non-symmorphic icosahedral groups
const t1=SVector{6, Rational{Int}}([1//5, 0//1, -1//5, 0//1, 0//1, -1//5])
const t2=SVector{6, Rational{Int}}([0//1, 0//1, 0//1, 0//1, 1//2, -1//2])

# Rotations
const s3=SpaceGroupElement(r3,z)
const s5=SpaceGroupElement(r5,z)

# Central symmetry
const sc=SpaceGroupElement(c,z)

# Non-symmorphic operations
const s5_1=SpaceGroupElement(r5,t1)
const s5_2=SpaceGroupElement(r5,t2)

const PI=SpaceGroupQuotient([s3,s5])
const PIh=SpaceGroupQuotient([s3,s5,sc])
const PI_n=SpaceGroupQuotient([s3,s5_1])
const PIh_n=SpaceGroupQuotient([s3, s5_2, sc])

end
