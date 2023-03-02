module Rayons

using LinearAlgebra

include("rectangles_from_png.jl")
include("unit_to_string.jl")
include("numericalParameters.jl")
include("material_properties.jl")
include("materialParameters.jl")
include("material_layers.jl")
include("point.jl")
include("readjson.jl")
include("abstractCelerityDomain.jl")
include("layerCelerityDomain.jl")
include("boxCelerityDomain.jl")
include("faceMesh.jl")
include("facemintime.jl")

include("ray.jl")
include("plotting.jl")
include("powersplit.jl")
include("locate_export.jl")
include("ray_tracing.jl")


end

   


















    


