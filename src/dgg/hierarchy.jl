#-----------------------------------------------------------------------------# Parent / Children / Siblings for DGGCell

function parent(o::DGGCell{P,A,T}) where {P,A,T}
    r = resolution(o)
    r == 0 && return nothing
    DGGCell{P,A,T}(dgg_base(o.index), dgg_digits(o.index, Val(A))[1:end-1])
end

function children(o::DGGCell{P,A,T}) where {P,A,T}
    r = resolution(o)
    maxd = _dgg_maxd(Val(A))
    r >= maxd && return DGGCell{P,A,T}[]
    base = dgg_base(o.index)
    ds = dgg_digits(o.index, Val(A))
    maxval = _dgg_maxval(Val(A))
    return [DGGCell{P,A,T}(base, [ds; d]) for d in 0:maxval]
end

function siblings(o::DGGCell{P,A,T}) where {P,A,T}
    r = resolution(o)
    r == 0 && return nothing
    return filter!(!=(o), children(parent(o)))
end

#-----------------------------------------------------------------------------# Pentagons
dgg_pentagons(::DGGGrid{P,A,:hex}, r::Integer) where {P,A} =
    [DGGCell{P,A,:hex}(b, fill(0, r)) for b in 0:11]
pentagons(g::DGGGrid{P,A,:hex}, r::Integer) where {P,A} = dgg_pentagons(g, r)
