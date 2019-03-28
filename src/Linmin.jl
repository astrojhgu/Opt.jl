module Linmin

import .. Utils.mysign
import .. Utils.Tolerance
import .. Utils.accepted

function brent(func::Function, ax::T, bx::T, cx::T,  tol::T)::Tuple{T,T} where {T<:AbstractFloat}
    two=one(T)+one(T)
    itmax=100
    cgold=T(0.381_966_0)
    zeps=eps(T)*T(1e-3)
    a,b=if ax<cx
        (ax,cx) else (cx,ax) end
    d=zero(T)
    e=zero(T)
    v=w=x=bx
    fv=fw=fx=func(x)
    iter=0
    etemp=zero(T);
    p=zero(T)
    q=zero(T)
    r=zero(T)
    tol1=zero(T)
    tol2=zero(T)
    xm=zero(T)

    u=zero(T)
    while iter<itmax
        xm=(a+b)/two
        tol1=tol*abs(x)
        tol2=two*(tol1+zeps)

        if abs(x-xm)<=(tol2-(b-a)/two)
            return (x,fx)
        end

        if abs(e)>tol1
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=two*(q-r)
            if q>zero(T)
                p=-p
            else
                q=-q
            end
            etemp=e
            e=d
            if abs(p)>=abs(p*etemp/two) || p<=q*(a-x) || p>=q*(b-x)
                e= if x>=xm
                    (a-x) else (b-x) end
                d=cgold*e
            else
                d=p/q
                u=x+d
                if (u-a)<tol2 || (b-u)<tol2
                    d=mysign(tol1, xm-x)
                end
            end
        else
            e=if x>=xm
                a-x
            else
                b-x
            end
            d=cgold*e
        end
        u=if abs(d)>=tol1
            x+d
        else
            x+mysign(tol1, d)
        end

        fu=func(u)
        if fu<=fx
            if u>=x
                a=x
            else
                b=x
            end
            v,w,x=w,x,u
            fv,fw,fx=fw,fx,fu
        else
            if u<x
                a=u
            else
                b=u
            end
            if fu<=fw || w==x
                v,w=w,u
                fv,fw=fw,fu
            elseif fu<=fv || v==x || v==w
                v,fv=u,fu
            end
        end
        iter+=1
    end
    (x,fx)
end

function mnbrak(func::Function,
    ax::T, bx::T)::Tuple{T,T,T,T,T,T} where {T<:AbstractFloat}
    two=one(T)+one(T)
    tiny=eps(T)
    gold = T(1.618_034)
    glimit=T(100)
    ulim=u=r=q=fu=zero(T)
    fa=func(ax)
    fb=func(bx)
    if fb>fa
        ax,bx=bx,ax
        fa,fb=fb,fa
    end
    cx = bx + gold * (bx - ax);
    fc=func(cx)

    while fb>fc
        r = (bx - ax) * (fb - fc);
        q = (bx - cx) * (fb - fa);
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (two * mysign(max(abs(q - r), tiny), q - r))

        ulim = bx + glimit * (cx - bx)

        if (bx - u) * (u - cx) > zero(T)
            fu = func(u)
            if fu < fc
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                return (ax,bx,cx,fa,fb,fc)
            elseif fu > fb
                cx = u
                fc = fu
                return (ax,bx,cx,fa,fb,fc)
            end
            u = cx + gold * (cx - bx)
            fu = func(u)
        elseif (cx - u) * (u - ulim) > zero(T)
            fu = func(u)
            if fu < fc
                xx = cx + gold * (cx - bx)
                bx,cx,u=cx,u,xx
                fb,fc,fu=fc,fu,func(u)
            end
        elseif (u - ulim) * (ulim - cx) >= zero(T)
            u = ulim
            fu = func(u)
        else
            u = cx + gold * (cx - bx)
            fu = func(u)
        end
        ax,bx,cx=bx,cx,u
        fa,fb,fc=fb,fc,fu
    end
    (ax,bx,cx,fa,fb,fc)
end

function linmin(func::Function,
    p::AbstractArray{T},
    xi::AbstractArray{T},
    )::Tuple{AbstractArray{T}, AbstractArray{T}, T} where {T<:AbstractFloat}
    tol=sqrt(eps(T))
    fb=fa=bx=ax=fx=zero(T)
    xx=one(T)
    funcadapter(x)=func(p+x*xi)
    ax, xx, bx, fa, fx, fb=mnbrak(funcadapter, ax, xx)
    xmin, fret=brent(funcadapter, ax, xx, bx, tol)
    xi=xi*xmin
    p=p+xi
    (p, xi, fret)
end

end  # module Linmin
