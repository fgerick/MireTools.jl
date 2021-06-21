function normbasis!(ptemp, v,cmat; n_cache=2*10^5)
    
    n = length(v)
    u = zeros(T,n)
    # cm = T.(cmat)
    Threads.@threads for k=1:n
        u[k] = sqrt(Mire._inner_product!(ptemp[Threads.threadid()],v[k],v[k],cmat))
        v[k]/=u[k]
    end
    return u
end

# function remove_factor!(u)
# 	u./= maximum(abs.(Mire.coefficients(u[1]+u[2]+u[3])))
# end

function normbasis!(P::T; n_cache=2*10^6) where T<:MHDProblem
    bs1 = copy(P.bbasis.el)
    map(remove_factor!,bs1)
    nt = Threads.nthreads()
    T3 = promote_type(T1,T2)
    ptemp = [zeros(Mire.Term{T3,Mire.Monomial{(x, y, z),3}},n_cache) for i=1:nt]

    normb = normbasis!(ptemp,bs1,P.cmat; n_cache)
    P.bbasis = typeof(P.bbasis)(P.N,P.V, bs1, P.bbasis.orthonorm)

    vs1 = copy(P.vbasis.el) 
    normu = normbasis!(ptemp,vs1,P.cmat; n_cache)
    P.vbasis = typeof(P.vbasis)(P.N,P.V, vs1, P.vbasis.orthonorm)
 
    return normu, normb
end