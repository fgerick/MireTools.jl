module MireTools

using Reexport

using Arpack
@reexport using ArnoldiMethod
using CartesianSphericalHarmonics
using IncompleteLU
using LinearAlgebra
using LinearMaps
@reexport using Mire
using Statistics



include("eigen.jl")

include("polynomialtools.jl")
export truncpoly, truncvec, @fastfunc, @fastfuncr

include("exteriorfield.jl")
export exteriormagneticfields, bex, cleanexcess

include("basistools.jl")
export normbasis!

end
