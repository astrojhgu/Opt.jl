module Utils

struct Tolerance{T}
    rel_tol::Union{T, Missing}
    abs_tol::Union{T, Missing}
end

function abs_tol(x::T)::Tolerance{T} where {T<:AbstractFloat}
    Tolerance(missing, x)
end

function rel_tol(x::T)::Tolerance{T} where{T<:AbstractFloat}
    Tolerance(x, missing)
end

function accepted(tol::Tolerance{T}, x1::T, x2::T)::Bool where {T<:AbstractFloat}
    two=one(T)+one(T)
    if !ismissing(tol.rel_tol)
        abs(x2-x1)<abs(x1+x2)*two*tol.rel_tol
    elseif !ismissing(tol.abs_tol)
        abs(x2-x1)<tol.abs_tol
    else
        error("Either abs or rel tol must be set")
    end
end

function mysign(a::T, b::T)::T where {T<:AbstractFloat}
    abs(a)*sign(b)
end



end
