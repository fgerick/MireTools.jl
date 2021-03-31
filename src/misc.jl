function ekinmag(u,i,vbasis,bbasis,LHS)
    nLHS = size(LHS,1)
    nu = length(vbasis)
    nb = length(bbasis)
    α = u[1:nu,i]/mean(u[:,i])
    cα = conj.(α)
    β = u[nu+1:nu+nb,i]/mean(u[:,i])
    cβ = conj.(β)
    ekin =zero(eltype(α))
    @inbounds for j=1:nu
        for i=1:nu
            ekin+=α[i]*cα[j]*LHS[i,j]
        end
    end
    emag = zero(eltype(β))
    @inbounds for j=1:nb
        for i=1:nb
            emag+=β[i]*cβ[j]*LHS[nu+i,nu+j]
        end
    end
    ratio = abs(ekin/emag)
    return abs(ekin),abs(emag),ratio
end