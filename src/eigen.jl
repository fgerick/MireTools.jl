function eigs(A::AbstractMatrix{T},B::AbstractMatrix{T}; kwargs...) where T
    P = ilu(B,τ=zero(real(T)))
    LO = LinearMap{complex(T)}(x->P\(A*x),size(A,2))
    pschur,history = partialschur(LO; kwargs...)
    λ, u = partialeigen(pschur)
    return λ,u
end

function eigstarget(A,B,target; kwargs...)
    P = lu(A-target*B)
    LO = LinearMap{ComplexF64}((y,x)->ldiv!(y,P,B*x),size(A,2))
    evals,u = Arpack.eigs(LO; kwargs...)
    λ = 1 ./evals .+ target
    return λ,u
end


function track(
    p::MireProblem, 
    updateLHS!::Function,
    updateRHS!::Function, 
    c0::T,
    c1::T,
    σ0::Complex{T}; 
    verbose=false, 
    corrtol=0.99, 
    maxstepfac=10, 
    kwargs...) where T

    @assert c0 != c1

    k=0
    cout = [c0]
    λs, us = Complex{T}[], Vector{Complex{T}}[]
    target = σ0
    targets = [σ0]
    utarget = Complex{T}[]
    LHS = copy(p.LHS)
    RHS = copy(p.RHS) 

    direction = sign(c1-c0)
    
    dc = c0/maxstepfac*direction
    # dc_t = dc
    c_t = c0
    iter = 0
    c = c0
    c_old = c0

    while sign(c1 - (cout[end] + dc)) == direction

        if k != 0
            c = cout[end] + dc
        end

        if (c != c_t) || (sign(c - c_t) != direction)
            c = c_t
        end
        
        if k > 0
            updateLHS!(LHS, p, c, c_old)
            updateRHS!(RHS, p, c, c_old)
            c_old = c
        end
        λ,u = try
            eigstarget(RHS, LHS, target; v0 = utarget, kwargs...)
        catch
            @warn "arpack error, stopping at c = $c"
            return λs,us,cout
        end
        nev = length(λ)

        if k == 0
            imax = 1 # assuming that for the first value the target hits the right eigensolution!
            c_t = c0 + dc
        else
            # calculate correlations between eigenvectors and previous solution
            corrs = [abs(cor(u[:,i],utarget)) for i=1:nev]
            max_corr,imax = findmax(corrs)

            if abs(corrs[imax]) < corrtol #if no correlating eigenvector is found the parameter stepping is too high
                if verbose
                    @warn "Correlation is only $(abs(corrs[imax]))!, lowering step"
                    @show dc
                    # flush(stdout)
                    # flush(stderr)
                end
                dc /= 2
                if abs(dc) <10eps()
                    @warn "abort tracking, stuck"
                    return λs,us,cout 
                end
                c_t = cout[end] + dc
                continue
            else
                push!(cout,c)
                dc = 2*abs(cout[end]-cout[end-1])*direction
                if abs(dc) > abs(cout[end]/maxstepfac)
                    dc = cout[end]/maxstepfac*direction
                end
                c_t = cout[end] + dc 
            end
        end

        push!(λs,λ[imax])
        push!(us,u[:,imax])
        target = λ[imax]
        utarget= u[:,imax]

        if verbose
            ω=abs(target)
            @show ω
        end
        
        push!(targets,target)
        k+=1
        if verbose
            @show c
            flush(stdout)
            flush(stderr)
        end
    end

    return λs,us,cout

end


function tracking_lehnert(N,cmat,nu,les,σ0, LHS0::AbstractArray{T},RHS0::AbstractArray{T}; verbose=false, corrtol=0.99, maxdlefac=100, kwargs...) where T
    k=0
    lesout = [les[1]]
    λs,us = Complex{T}[],Vector{Complex{T}}[]
    target = σ0
    targets = [σ0]
    utarget = Complex{T}[]
    LHS = copy(LHS0)
    RHS = copy(RHS0)

    dle = les[1]/maxdlefac #diff(les) # les[2]-les[1]
    dlet = dle
    le_t = les[1]
    iter = 0
    mpeakold = 0


    while lesout[end] + dle <= les[end]

        if (k == 0 )
            Le = les[1]
        else
            Le = lesout[end]+dle #dle[k+1]
        end

        if Le>le_t
            Le = le_t
        end
        
		#only coriolis force depends on Le:
        if k>0
            view(RHS,1:nu,1:nu).*=lesout[end]/Le
        end

        λ,u = eigstarget(RHS,LHS, target; v0=utarget, kwargs...)
        nev = length(λ)
        if k == 0
            imax = 1
            le_t = les[1] + dle
        else
            corrs = [abs(cor(u[:,i],utarget)) for i=1:nev]
            max_corr,imax = findmax(corrs)
            if abs(corrs[imax]) < corrtol #if no correlating eigenvector is found the parameter stepping is too high
                if verbose
                    @warn "Correlation is only $(abs(corrs[imax]))!, lowering step"
                flush(stdout)
                flush(stderr)
                end
                dlet /= 2
                if dlet <10eps()
                    break
                end
                le_t = lesout[end] + dlet
                continue
            else
                push!(lesout,Le)
                dle = 2*abs(lesout[end]-lesout[end-1])
                if dle > lesout[end]/maxdlefac
                    dle = lesout[end]/maxdlefac
                end
                le_t = lesout[end] + dle #dle[k+1]
                dlet = dle #[k+1]
            end
        end

        push!(λs,λ[imax])
        push!(us,u[:,imax])
        target = λ[imax]
        utarget= u[:,imax]

        if verbose
            ω=abs(target)
            @show ω

        end
        push!(targets,target)
        k+=1
        if verbose
            @show Le
            flush(stdout)
            flush(stderr)
        end
    end
    return λs,us,lesout
end
