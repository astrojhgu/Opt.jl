module Opt

include("Utils.jl")
include("Linmin.jl")
include("OptResult.jl")
using LinearAlgebra
import .Linmin.linmin
import .Utils.accepted
import .Utils.Tolerance
import .OptResult.ResultType

const rel_tol=Utils.rel_tol
const abs_tol=Utils.abs_tol

function fmax(func::Function, p::U, ftol::Tolerance{T}, itmax::Int)::Tuple{U, ResultType} where {T<:AbstractFloat, U<:AbstractArray{T}}
    fmin(p->-func(p), p, ftol, itmax)
end

function fmin(func::Function, p::U, ftol::Tolerance{T}, itmax::Int)::Tuple{U, ResultType} where {T<:AbstractFloat, U<:AbstractArray{T}}
    two=one(T)+one(T)
    n=length(p)
    fret=func(p)
    xi=Array(Diagonal(ones(T, n)))
    pt=p
    for iter=1:itmax
        fp=fret
        ibig=0
        del=zero(T)
        for i in 1:n
            xit=xi[:,i]
            fptt=fret
            p, xit, fret=linmin(func, p, xit)
            if fptt-fret>del
                del=fptt-fret
                ibig=i
            end
        end

        if accepted(ftol, fp, fret)
            return (p, OptResult.Finished())
        end

        if iter==itmax
            return (p, OptResult.MaxIterReached())
        end

        ptt=two*p-pt
        xit=p-pt
        pt=p

        fptt=func(ptt)

        if fptt<fp
            t=two*(fp-two*fret+fptt)*(fp-fret-del)^2-del*(fp-fptt)^2
            if t<zero(T)
                p, xit, fret=linmin(func, p, xit)
                xi[:, ibig]=xi[:,n]
                xi[:, n]=xit
            end
        end
    end
    (p, OptResult.Finished())
end

end # module
