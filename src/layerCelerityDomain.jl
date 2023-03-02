struct MLB
    hmin::Float64
    hmax::Float64
    cl::Float64
    ct::Float64
    ρ::Float64
end

struct LayerCelerityDomain <: AbstractCelerityDomain
    mlbs::Vector{MLB}
    width::Float64
end

function LayerCelerityDomain(material_layers::Vector{MaterialLayer})
    mlbs=Vector{MLB}()
    hmin=0.0
    for ml ∈ material_layers
        hmax=ml.height
        push!(mlbs,MLB(hmin,hmax,ml.cₗ,ml.cₜ,ml.ρ))
        hmin=hmax
    end
    LayerCelerityDomain(mlbs,0.2)
end

struct LayerCelerityState
    lcd::LayerCelerityDomain
    mlb_current::MLB
    mlb_previous::MLB
    mlb_zero::MLB
end

function initial_state(lcd::LayerCelerityDomain) 
    LayerCelerityState(lcd,lcd.mlbs[1],lcd.mlbs[2],MLB(0.0,0.0,0.0,0.0,0.0))
end
@inline outside_width(width,x) = x<0 || x > width 
@inline inside_mlb(mlb,x,y) = (y>=mlb.hmin) && (y<mlb.hmax)

function get_state(p::Point,state::LayerCelerityState)
    (x,y) = p.x,p.y
    outside_width(state.lcd.width,x) && return LayerCelerityState(state.lcd,state.mlb_zero,state.mlb_current,state.mlb_zero)
    inside_mlb(state.mlb_current,x,y) && return state
    if inside_mlb(state.mlb_previous,x,y)
        return LayerCelerityState(state.lcd,state.mlb_previous,state.mlb_current,state.mlb_zero)
    else
        for mlb ∈ state.lcd.mlbs
            if inside_mlb(mlb,x,y)
                return LayerCelerityState(state.lcd,mlb,state.mlb_current,state.mlb_zero)
            end
        end
    end
    return LayerCelerityState(state.lcd,state.mlb_zero,state.mlb_current,state.mlb_zero)
end  

cl_ct_rho(state::LayerCelerityState) = state.mlb_current.cl,state.mlb_current.ct,state.mlb_current.ρ
longitudinal_celerity(state::LayerCelerityState) = state.mlb_current.cl





    


