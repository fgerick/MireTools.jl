using MireTools
using Test

@testset "MireTools.jl" begin
    p = HDProblem(5,Sphere(),[0,0,1.0],LebovitzBasis)
    assemble!(p)
    evals,evecs = MireTools.eigs(p.RHS,p.LHS; nev=80, which=LM())
     # Zhang et al., J. Fluid Mech. (2001), vol. 437, pp. 103–119. eq (4.5)
     zhang(m,N)=-2/(m+2)*(√(1+m*(m+2)/(N*(2N+2m+1)))-1)*im
     @test any(evals.≈zhang(1,1))
    # Write your tests here.
end
