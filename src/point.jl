struct Point
    x :: Float64
    y :: Float64
end

# Base extensions
Base.:+(p1::Point,p2::Point) = Point(p1.x+p2.x,p1.y+p2.y)
Base.:-(p1::Point,p2::Point) = Point(p1.x-p2.x,p1.y-p2.y)
Base.:*(a::Real,p::Point) = Point(p.x*a,p.y*a)
Base.:*(p::Point,a::Real) = Point(p.x*a,p.y*a)
Base.:/(p::Point,a::Real) = Point(p.x/a,p.y/a)

# Linear algebra
LinearAlgebra.dot(p1::Point,p2::Point) = p1.x*p2.x + p1.y*p2.y
squared_norm(p::Point) = dot(p,p)
LinearAlgebra.norm(p::Point) = sqrt(squared_norm(p))
@inline zcross(p1::Point,p2::Point) = p1.x*p2.y-p2.x*p1.y

# Transformations
function rotate(p::Point, θ::Real)
    sinθ,cosθ = sincos(θ)
    Point(p.x*cosθ-p.y*sinθ, p.x*sinθ+p.y*cosθ)
end

# Barycenter
center(p::AbstractVector{Point}) = sum(p) / length(p)

# Bounding box
_bounding_box(ps) = extrema(p.x for p in ps)..., extrema(p.y for p in ps)...
bounding_box(ps::AbstractArray{Point,1})    = _bounding_box(ps)
bounding_box(ps::NTuple{N,Point}) where {N} = _bounding_box(ps)

# Broadcast utility: treat Point as a scalar
Base.broadcastable(m::Point) = Ref(m)
