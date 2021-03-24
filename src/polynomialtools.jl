function truncpoly(u;atol=√(eps()))
    uo=zero(u)
    ti = terms(u)
    for t in ti
        c=coefficient(t)
        if abs(c)>atol
            uo+=t
        end
    end
    return uo
end

function truncvec(u;atol=√(eps()))
    uo=zero(u)
    for i=1:3
        ti = terms(u[i])
        for t in ti
            c=coefficient(t)
            if abs(c)>atol
                uo[i]+=t
            end
        end
    end
    return uo
end


macro fastfunc(polynomial)
	return quote
		local polystr = string($(esc(polynomial)))
		local vars = string(tuple(Mire.variables($(esc(polynomial)))...))
        # create expression for function definition
		eval(Meta.parse(vars*" -> @fastmath "*polystr))
	end
end

macro fastfuncr(polynomial)
	return quote
		local polystr = string($(esc(polynomial)))
		local vars = string(tuple(Mire.variables(($(esc(polynomial)).num))...))
        # create expression for function definition
		eval(Meta.parse(vars*" -> @fastmath "*polystr))
	end
end

#some functions that help us convert the datatypes of the polynomial coefficients:
function convert_polynom(S::Type{T},p) where T
	c = coefficients(p)
	m = monomials(p)
	return sum(S.(c).*m)
end

convert_polyvec(S::Type{T},u) where T = convert_polynom.(S,u)

#remove some large factor of polynomial vector with possible large coefficients
function remove_factor!(u)
	u./= maximum(abs.(coefficients(u[1]+u[2]+u[3])))
end

import Base.real
function real(p::AbstractPolynomialLike)
	c = coefficients(p)
	m = monomials(p)
	return dot(real.(c),m)
end
import Base.conj
function conj(p::AbstractPolynomialLike)
	c = coefficients(p)
	m = monomials(p)
	return dot(c,m)
end