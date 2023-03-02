
struct MaterialLayer
    height::Float64
    cₗ::Float64
    cₜ::Float64
    ρ::Float64   
end

function get_material_parameter_from_tag(material_parameters,tag)
    for mp ∈ material_parameters
        mp.name == tag && return mp
    end
    error("tag $tag not found")
end

function get_material_layers(layers,material_parameters)
    result=Vector{MaterialLayer}()
    h=0.0
    for l ∈ layers
        h += l.height
        tag = l.tag
        mp = get_material_parameter_from_tag(material_parameters,tag)
        push!(result,MaterialLayer(h,mp.cₗ,mp.cₜ,mp.ρ))
    end
    result
end

function material_layer_from_position(material_layers,width,x,y)
    for ml ∈ material_layers
        y < ml.height && return ml
    end

    error("y out of bounds $y")
end

function layer_mesh_parameters(material_layers,ny)
    height=material_layers[end].height
    minwidth=0.2
    h = height/ny
    width = ceil(minwidth/h)*h
    nx = Int(width/h)
    nx,ny,h,width,height
    # @show nx,ny,h,height,width
end


struct Layer
    tag::String
    height::Float64
end

function getlayers()
    layers = [
        Layer("fer",0.02),
        Layer("huile",0.03),
        Layer("carton",0.005),
        Layer("cuivre",0.02),
        Layer("carton",0.005),
        Layer("huile",0.03),
        Layer("acier",0.005)
        ]
end





   


















    


