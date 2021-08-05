function normbasis!(v,cmat)
    
    n = length(v)
    # T = coefficienttype(ptemp[1])
    u = zeros(coefficienttype(v[1][1]),n)
    # cm = T.(cmat)
    Threads.@threads for k=1:n
        u[k] = sqrt(Mire.inner_product(v[k],v[k],cmat))
        v[k]/=u[k]
    end
    return u
end

# function remove_factor!(u)
# 	u./= maximum(abs.(Mire.coefficients(u[1]+u[2]+u[3])))
# end

function normbasis!(P::T) where T<:MHDProblem
    bs1 = copy(P.bbasis.el)
    map(remove_factor!,bs1)
    nt = Threads.nthreads()
    T3 = promote_type(eltype(P.cmat),coefficienttype(P.bbasis.el[1][1]),coefficienttype(P.vbasis.el[1][1]))
    # ptemp = [zeros(Mire.Term{T3,Mire.Monomial{(x, y, z),3}},n_cache) for i=1:nt]

    normb = normbasis!(bs1,P.cmat)
    P.bbasis = typeof(P.bbasis)(P.N,P.V, bs1, P.bbasis.orthonorm)

    vs1 = copy(P.vbasis.el) 
    normu = normbasis!(vs1,P.cmat)
    P.vbasis = typeof(P.vbasis)(P.N,P.V, vs1, P.vbasis.orthonorm)
 
    return normu, normb
end