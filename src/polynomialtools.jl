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
