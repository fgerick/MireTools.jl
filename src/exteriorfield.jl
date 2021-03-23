
function gradIlm(l,m,x,y,z,r; kwargs...)
    rlmt = Mire.rlm(l,m,x,y,z; kwargs...)
    [Mire.∂(rlmt,ξ)/r^(2l+1)-(2l+1)*ξ*rlmt/r^(2l+3) for ξ in [x,y,z]]
end

function bex(Plmn::AbstractVector{Complex{T}},ls,ms,ns) where T
    R = Variable{:R}()
    b_ex = Plmn[1]*Mire.ve(ls[1],ns[1])*gradIlm(ls[1],ms[1],x,y,z,R; norm=Mire.Schmidt{T}())
    for i=2:length(Plmn)
        # if abs(Plmn[i]) > sqrt(eps())
        #
            b_ex += Plmn[i]*Mire.ve(ls[i],ns[i])*gradIlm(ls[i],ms[i],x,y,z,R; norm=Mire.Schmidt{T}())
        # end
    end
    b_ex
    return b_ex
end


function cleanexcess(p::RationalPoly)
    n = p.num
    d = p.den
    tn = terms(n)
    en = exponents.(tn)
    cn = coefficients(n)
    R = Variable{:R}()
    emax = exponent(monomial(tn[end]),R)
    tnew = map((c,(l,i,j,k))->c*x^i*y^j*z^k*R^(l-emax), cn, en)
    drest = R^(exponent(d,R)-emax)
    return polynomial(tnew)/drest
 end
 
function exteriormagneticfields(bbasis,evec,ls,ms,ns)
    nb = length(bbasis)
    n = length(evec)
    nstart = n-nb+1
    nend = nstart+length(ls)-1
    cs = @views evec[nstart:nend]
    return cleanexcess.(bex(cs,ls,ms,ns))
end


# function gausscoeffs(clmn::AbstractVector{Complex{T}},ls,ms,ns,normfac) where T
#     @assert length(clmn)==length(ls)
#     L = ls[end]
#     lsout = [l for l in 0:L for m in 0:l]
#     msout = [m for l in 0:L for m in 0:l]
#     glm = zeros(Complex{T},length(lsout))
#     hlm = zeros(Complex{T},length(lsout))
#     clmne = clmn.*Mire.ve.(ls,ns)./normfac
#     for (l,m,i) in zip(lsout,msout,1:length(glm))
#         glm[i] = sum(clmne[(ls.==l) .& (ms.==m)])
#         hlm[i] = sum(clmne[(ls.==l) .& (ms.==-m)])
#         if iseven(m)
#             glm[i]*=-1
#             hlm[i]*=-1
#         end
#     end
#     return glm,hlm,lsout,msout
# end

# function gauss2bext(glm::AbstractVector{Complex{T}},hlm::AbstractVector{Complex{T}},ls,ms) where T
#     b_ex = glm[1]*gradIlm(ls[1], ms[1],x,y,z,R; norm=Mire.Schmidt{T}())
#     # +
#     #        hlm[1]*gradIlm(ls[1], ms[1],x,y,z,R; norm=Schmidt{Float64}())
#     for i=2:length(glm)
#             b_ex += glm[i]*gradIlm(ls[i], ms[i],x,y,z,R; norm=Mire.Schmidt{T}())
#             if (ms[i]!=0)
#                 b_ex+=hlm[i]*gradIlm(ls[i],-ms[i],x,y,z,R; norm=Mire.Schmidt{T}())
#             end
#     end
#     return b_ex
# end
