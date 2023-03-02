# using Infiltrator

struct BoxCelerityState
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    i::Int64
    j::Int64
    cl::Float64
    ct::Float64
    ρ::Float64
end
zero(::Type{BoxCelerityState}) = BoxCelerityState(0.0,0.0,0.0,0.0,0,0,0.0,0.0,0.0)

struct BoxCelerityDomain <: AbstractCelerityDomain
    cmg::CartesianMaterialGrid
    materials::Vector{MaterialParameters}
    function BoxCelerityDomain(cmg::CartesianMaterialGrid,mps_init::Vector{MaterialParameters})
        mps = Vector{MaterialParameters}()
        for material_name ∈ cmg.materials
            material_idx=findfirst(x->isequal(material_name,x.name),mps_init)
            @show material_name,material_idx
    
            if !isnothing(material_idx) 
                push!(mps,mps_init[material_idx])
            else
                push!(mps,MaterialParameters(material_name,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0))
            end
        end
        new(cmg,mps)
    end

end


@enum Neighbour Left Right Up Down

function show_error(l,p,u,cs)
    @show l,"+++ on ne devrait pas arriver là +++"
    @show p
    @show u
    @show cs
    error("")
end



@inline function  intersection_rectangle(pstart::Point, u::Point, cs::BoxCelerityState,bcd::BoxCelerityDomain,nx,ny)
    xm,xM,ym,yM = cs.xmin,cs.xmax,cs.ymin,cs.ymax
    # @show  xm,xM,ym,yM 
    i,j=cs.i,cs.j
    p = pstart

    if  u.x > 0   # on arrive sur un bord vertical à droite
        qvx=xM
        dx = qvx - p.x 
        sv = dx/u.x # maybe Inf but its OK because the interval condition is then false
        qvy = p.y + sv*u.y
        # @show "ixi",qvy
        if (ym <= qvy < yM)
            qv = Point(qvx,qvy)
            next_state = (i<nx) ? get_state_from_ij(i+1,j,bcd) : outside_state(cs) 
            return qv,next_state,u
        end
    else          # on arrive sur un bord vertical à gauche
        qvx=xm
        dx = qvx - p.x 
        sv = dx/u.x # maybe Inf but its OK because the interval condition is then false
        qvy = p.y + sv*u.y
        # @show "la",qvy
        if (ym <= qvy < yM)
            qv = Point(qvx,qvy)
            next_state = (i>1) ? get_state_from_ij(i-1,j,bcd) : outside_state(cs)
            return qv,next_state,Point(-u.x,-u.y)
        end
    end

    if u.y > 0  #on arrive sur un bord horizontal en haut
        qhy=yM
        dy = qhy - p.y
        sh  = dy/u.y # maybe Inf but its OK because the interval condition is then false
        qhx = p.x + sh*u.x
        # @show "ixi",qhx
        if (xm <= qhx < xM)
            qh = Point(qhx,qhy)
            next_state = (j<ny) ? get_state_from_ij(i,j+1,bcd) : outside_state(cs)
            return qh,next_state,Point(u.y,-u.x)
        end
    else             #on arrive sur un bord horizontal en bas
        qhy=ym
        dy = qhy - p.y
        sh  = dy/u.y # maybe Inf but its OK because the interval condition is then false
        qhx = p.x + sh*u.x
        # @show "la",qhx
        if (xm <= qhx < xM)
            qh = Point(qhx,qhy)
            next_state = (j>1) ? get_state_from_ij(i,j-1,bcd) : outside_state(cs)
            return qh,next_state,Point(-u.y,u.x)
        end
    end

    show_error(@__LINE__,pstart,u,cs)
    return pstart,outside_state(cs),Point(-u.y,u.x)
end

function check_state(line,p,u,cs)
    ok = state_inbounds(p,cs)
    !ok && show_error(line,p,u,cs)
    return ok
end


function advance_rectangle(pstart::Point, u::Point, celerity_state,bcd::BoxCelerityDomain,ϵ::Float64)
    next_celerity_state = celerity_state
    cl = longitudinal_celerity(celerity_state)
    cln = cl
    qn,uθ = pstart,pstart
    nx,ny = get_cell_sizes(bcd.cmg)
    while cln==cl 
        celerity_state = next_celerity_state
        qn,next_celerity_state,uθ = intersection_rectangle(pstart, u, celerity_state,bcd,nx,ny)
        cln = longitudinal_celerity(next_celerity_state)
    end
    return qn,qn,uθ,celerity_state,next_celerity_state
end


state_inbounds(p::Point,cs::BoxCelerityState) = cs.xmin<=p.x<cs.xmax && cs.ymin<=p.y<cs.ymax
@inline initial_state(bcd::BoxCelerityDomain) = zero(BoxCelerityState)
@inline outside_state(state::BoxCelerityState) = zero(BoxCelerityState)

@inline function get_state_from_ij(i,j,bcd::BoxCelerityDomain)
    cmg = bcd.cmg
    xnodes,ynodes = cmg.xnodes,cmg.ynodes
    @inbounds xmin,xmax = xnodes[i],xnodes[i+1]
    @inbounds ymin,ymax = ynodes[j],ynodes[j+1]
    @inbounds material_idx = cmg.material_grid[i,j]
    @inbounds mp = bcd.materials[material_idx]
    return BoxCelerityState(xmin,xmax,ymin,ymax,i,j,mp.cₗ,mp.cₜ,mp.ρ)
end


function get_state(p::Point,state::BoxCelerityState,bcd::BoxCelerityDomain)
    (x,y) = p.x,p.y
    # same box
    state_inbounds(p,state) && return state
    # outside domain
    out_of_bounds(x,y,bcd.cmg) && return outside_state(state)
    # new box
    i,j = get_cartesian_indexes_from_xy(x,y,bcd.cmg,state.i,state.j)

    return get_state_from_ij(i,j,bcd)
end  





function get_next_state(p::Point,state::BoxCelerityState,bcd::BoxCelerityDomain,neighbour::Neighbour)
    (x,y) = p.x,p.y
    cmg = bcd.cmg
    out_of_bounds(x,y,cmg) && return outside_state(state)  
    i,j=state.i,state.j
    neighbour==Right && return get_state_from_ij(i+1,j,bcd)
    neighbour==Left  && return get_state_from_ij(i-1,j,bcd)
    neighbour==Up    && return get_state_from_ij(i,j+1,bcd)
    neighbour==Down  && return get_state_from_ij(i,j-1,bcd)
    return get_state_from_ij(i,j,bcd)
end


cl_ct_rho(state::BoxCelerityState) = state.cl,state.ct,state.ρ
longitudinal_celerity(state::BoxCelerityState) = state.cl













    


