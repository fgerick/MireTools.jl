function normbasis!(v::Vector{Mire.vptype{T1}},cmat::Array{T2,3}; n_cache=2*10^5) where {T1,T2}
    T = promote_type(T1,T2)
    ptemp = zeros(Mire.Term{T,Mire.Monomial{(x, y, z),3}},n_cache);
    n = length(v)
    u = deepcopy(v)
    cm = T.(cmat)
    for k=1:n
        nrm=Mire.inner_product!(ptemp,u[k],u[k],cm)
        u[k]/=√nrm
    end
    return u
end

#remove some large factor of polynomial vector with possible large coefficients
# function remove_factor!(u)
# 	u./= maximum(abs.(Mire.coefficients(u[1]+u[2]+u[3])))
# end

function normbasis!(P::T; n_cache=2*10^6) where T<:MHDProblem
    bs1 = P.bbasis.el
    # map(remove_factor!,bs1)
    bs = normbasis!(bs1,P.cmat; n_cache)
    P.bbasis = typeof(P.bbasis)(P.N,P.V, bs, false)

    vs = normbasis!(P.vbasis.el,P.cmat; n_cache)
    P.vbasis = typeof(P.vbasis)(P.N,P.V, vs, false)
 
    return nothing
end