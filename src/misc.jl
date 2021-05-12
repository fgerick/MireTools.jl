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

function ekinmag(u, nu, LHS)
    LHSu = LHS[1:nu,1:nu]
    LHSb = LHS[nu+1:end,nu+1:end]
    iu,ju,valu = findnz(LHSu)
    ib,jb,valb = findnz(LHSb)
    nev = size(u,2)

    ekin = zeros(eltype(u),nev)
    emag = zeros(eltype(u),nev)
     
    for iev in 1:nev
        phase=mean(u[:,iev])
        for (i,j,val) in zip(iu,ju,valu)
            ekin[iev]+=u[i,iev]*conj(u[j,iev])*val
        end
        ekin[iev]/=phase^2
        for (i,j,val) in zip(ib,jb,valb)
            emag[iev]+=u[i,iev]*conj(u[j,iev])*val
        end
        emag[iev]/=phase^2
    end 
    return abs.(ekin),abs.(emag),abs.(ekin/emag)
end



    

function ekinmag(u,i,nu,LHS)
    phase=mean(u[:,i])
    @views ekin = u[1:nu,i]'*LHS[1:nu,1:nu]*u[1:nu,i]/phase^2
    @views emag = u[nu+1:end,i]'*LHS[nu+1:end,nu+1:end]*u[nu+1:end,i]/phase^2
    return abs(ekin),abs(emag),abs(ekin/emag)
end
     