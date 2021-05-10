module MireTools

using Reexport

using Arpack
@reexport using ArnoldiMethod
using CartesianSphericalHarmonics
using IncompleteLU
using LinearAlgebra
using LinearMaps
@reexport using Mire
using MultivariatePolynomials
using PyCall
using PyPlot
using Statistics
using TypedPolynomials



const cartopy = PyNULL()

function __init__()
    copy!(cartopy, pyimport("cartopy"))
end

export cartopy

include("eigen.jl")

include("polynomialtools.jl")
export truncpoly, truncvec, @fastfunc, @fastfuncr

include("exteriorfield.jl")
export exteriormagneticfields, bex, cleanexcess

include("basistools.jl")
export normbasis!

include("misc.jl")
export ekinmag

include("plotting.jl")

end
