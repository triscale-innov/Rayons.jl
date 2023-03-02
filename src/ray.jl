# <image src= file:///Users/laurentplagne/Projects/Rayons.jl/min_traj.png>
import Base.@kwdef
using StaticArrays

struct Ray
    pstart::Point
    pend::Point
    c::Float64
    time::Float64
    power::Float64
end

Ray() = Ray(Point(0, 0), Point(0, 0), 0.0, -1.0,1.0)


# <image src= file:///Users/laurentplagne/Projects/Rayons.jl/min_traj.png>


"""
Linear search of the media change on a line begining at `pstart` by step `dp`
returns the point just before contact.
![mon image !](/assets/min_traj.png>)
"""
function last_point_before_contact(pstart, dp, celerity_state)
    i, maxsteps = 1, 1000
    cl = longitudinal_celerity(celerity_state)
    celerity_state = get_state(pstart + i * dp, celerity_state)

    while longitudinal_celerity(celerity_state) == cl
        i += 1
        i > maxsteps && error("maxstep reach in last_point_before_contact : i=$i,p=$(pstart+i*dp), c=$cl, $(longitudinal_celerity(celerity_state)) $(pstart + i * dp)")
        celerity_state = get_state(pstart + i * dp, celerity_state)
    end

    pstart + (i - 1) * dp   # last point before contact
end

# Bissection search for pend (and pnext) between ps and pe
# one assumes c=celerity(ps)!=celerity(pe)=c2
# return pend and pnext s.t. 
#          c==celerity(pend)!=celerity(p,next)
#           norm(pend-pnext)<ϵ
function bissection_to_pend(ps, pe, celerity_state, ϵ)
    i, maxsteps = 1, 100
    cl = longitudinal_celerity(celerity_state)
    while (squared_norm(ps - pe) > ϵ^2)
        pm = (ps + pe) / 2
        celerity_state = get_state(pm, celerity_state)
        cm = longitudinal_celerity(celerity_state)
        cm == cl ? (ps = pm) : (pe = pm)
        ((i+=1) > maxsteps) && error("maxstep reach in bissection_to_pend pm=$(pm),c$cm")
    end
    # pstart and pend are in the same speed domain 
    # pnext point is in the next speed domain
    ps, pe
end

@inline rotated_point(o, r, c, s, u, v) = o + r * (c * u + s * v)

#Compute angle (unit vector) between the media boundary normal and the incident 
# ray of direction u.
function bissection_angle(pend, pnext, ds, u, celerity_state, ϵ)
    pm = (pend + pnext) / 2 # pm almost on the surface
    r = ds / 20
    v = Point(-u.y, u.x) # u ⟂ v
    cθmax, cθmin = 1.0, -1.0
    cl = longitudinal_celerity(celerity_state)
    i, maxsteps = 1, 100

    # LP: could be accelerated by not imposing to look for a solution on a circle
    while (cθmax - cθmin) > ϵ
        cθmid = (cθmin + cθmax) / 2
        sθmid = sqrt(1.0-cθmid^2) 
        p = rotated_point(pm, r, cθmid, sθmid, u, v)
        celerity_state = get_state(p, celerity_state)
        cp = longitudinal_celerity(celerity_state)
        # cp, celerity_state = celerity_function(p, celerity_state)
        (cp != cl) ? (cθmax = cθmid) : (cθmin = cθmid)
        ((i+=1) > maxsteps) && error("maxstep reach in bissection_angle")
    end
    cθmid = (cθmin + cθmax) / 2
    sθmid = sqrt(1.0 - cθmid^2)
    umv(Point(0.0,1.0),Point(cθmid,sθmid))# π/2 - θmid : angle between u and the boundary normal.
end

# Define Ray segment between two media starting from pstart with direction u
# return pend,pnext and θ s.t
#          celerity(pstart)==celerity(pend)!=celerity(p,next)
#          norm(pend-pnext)<ϵ
#          θ = angle between the incident ray and the media boundary
function advance_ds(pstart::Point, u::Point, np::NumericalParameters, celerity_state,celerity_domain)
    ds, ϵ = unpack(np)
    # uncomment this to use non rectangular domain
    # dp = ds * u
    # ps = last_point_before_contact(pstart, dp, celerity_state)
    # pend, pnext = bissection_to_pend(ps, ps + dp, celerity_state, ϵ)
    # uθ = bissection_angle(pend, pnext, ds, u, celerity_state, ϵ)

    pend,pnext,uθ,celerity_state,next_celerity_state = advance_rectangle(pstart,u,celerity_state,celerity_domain,ϵ)
    # !ok && @show "something went wrong !!!",pstart,u,celerity_state
    return pend, pnext,uθ,celerity_state,next_celerity_state
end

function unit_vector(α)
    s, c = sincos(α)
    Point(c, s)
end

clipindex(i, nx) = max(min(i, nx), 1)

function ray_to_cell_mintimes!(fm::FaceMesh, pstart, pend, c, te, ray_idx, αorigin,porigin,power,store_history::Bool)
    l = norm(pend - pstart)
    l == 0.0 && return
    tb = te - l / c
    u = (pend - pstart) / l
    (nx, ny) = get_cell_sizes(fm)
    
    xaxis, yaxis = fm.xfaces, fm.yfaces
    
    xmin, xmax = first(xaxis), last(xaxis)
    ymin, ymax = first(yaxis), last(yaxis)
    
    # @show pstart.x,pstart.y,pend.x,pend.y
    # @show xmax,xmin,ymax,ymin
    
    # (pstart.x < xmin || pstart.x > xmax) && return 
    # (pstart.y < ymin || pstart.y > ymax) && return 
    
    # (pend.x < xmin || pend.x > xmax) && return 
    # (pend.y < ymin || pend.y > ymax) && return 

    hm1=inv(fm.h)
    
    i_from_x(x) = Int(floor(((x - xmin) * hm1))) + 1
    j_from_y(y) = Int(floor(((y - ymin) * hm1))) + 1
    
    xb, yb = pstart.x, pstart.y
    xe, ye = pend.x, pend.y
    
    # @show xb,xe,fm.h
    
    ib, ie = clipindex(i_from_x(xb), nx), clipindex(i_from_x(xe), nx)
    jb, je = clipindex(j_from_y(yb), ny), clipindex(j_from_y(ye), ny)
    
    coarse(i) = fine_to_coarse_index(i,fm.coarse_factor,fm.coarse_two_power)
    

    IB,IE,JB,JE = coarse.((ib,ie,jb,je))

    (IB,IE) = IE>IB ? (IB,IE) : (IE,IB)
    (JB,JE) = JE>JB ? (JB,JE) : (JE,JB)

    skip = true
    @inbounds for I ∈ IB:IE
        for J ∈ JB:JE
            if tb < fm.coarse_ray_min_times[I,J]
                skip = false
                break
            end
        end
    end

    # skip && @show IB:IE,JB:JE
    skip && return

    fαs=length(fm.αscore)/(2π)


    function setmin!(i, j, t, p::Point)
        use_coarse = true
        if use_coarse
            I,J = coarse.((i,j))
            @inbounds max_coarse_time = fm.coarse_ray_min_times[I,J]
            if t>max_coarse_time
                # @show "save some time"
                return nothing
            end
            # update_coarse_ray_min_time!(fm,I,J)
        end

        @inbounds mt = fm.ray_min_times[i, j].time
        if t < mt
            @inbounds fm.ray_min_times[i, j] = MinRay(t, p, u, ray_idx, αorigin,porigin)
            
            αidx=Int(floor(fαs*αorigin))+1
            αidx=min(αidx,length(fm.αscore))
            if fm.αscore[αidx]<1000
                fm.αscore[αidx]+=1
            end
            update_coarse_ray_min_time!(fm,I,J,i,j)
        end


        if store_history
            szs=size(fm.history)
            if (i<1 || i>szs[1] || j<1 || j>szs[2])
                @show i,j,szs,αorigin,porigin
            end
        
            hij = fm.history[i, j]
            nh = length(hij)
            # @show t
            time_ratio = min((t+1.e-12)/fm.tmax,1.0)
            # @show fm.tmax,(t+1.e-12),(t+1.e-12)/fm.tmax
            # @show time_ratio, min((t+1.e-16)/fm.tmax,1.0)
            history_index=Int(ceil(nh*time_ratio))
            @show history_index
            @inbounds hij[history_index] += power
        end

        return nothing
    end

    # first point
    ilast, jlast = ib, jb
    setmin!(ilast, jlast, tb, Point(xb, yb))

    txc = (te - tb) / (xe - xb)
    tyc = (te - tb) / (ye - yb)

    if abs(u.x) >= abs(u.y)
        a = u.y / u.x
        ri = ie > ib ? (ib+1:ie) : ib:-1:ie+1
        ishift = ie > ib ? 0 : 1
        for i ∈ ri
            @inbounds x = xaxis[i]
            y = yb + (x - xb) * a
            j = j_from_y(y)
            tx = tb + (x - xb) * txc
            setmin!(i - ishift, j, tx, Point(x, y))
            if (jlast < j)
                @inbounds y = yaxis[j]
                ty = tb + (y - yb) * tyc
                setmin!(ilast - ishift, j, ty, Point(x, y))
            elseif (jlast > j)
                @inbounds y = yaxis[jlast]
                ty = tb + (y - yb) * tyc
                setmin!(ilast - ishift, j, ty, Point(x, y))
            end
            ilast, jlast = i, j
        end
    else
        a = u.x / u.y
        tyc = (te - tb) / (ye - yb)
        rj = je > jb ? (jb+1:je) : jb:-1:je+1
        jshift = je > jb ? 0 : 1
        for j ∈ rj
            @inbounds y = yaxis[j]
            x = xb + (y - yb) * a
            i = i_from_x(x)
            ty = tb + (y - yb) * tyc
            setmin!(i, j - jshift, ty, Point(x, y))
            if (ilast < i)
                @inbounds x = xaxis[i]
                tx = tb + (x - xb) * txc
                setmin!(i, jlast - jshift, tx, Point(x, y))
            elseif (ilast > i)
                @inbounds x = xaxis[ilast]
                tx = tb + (x - xb) * txc
                setmin!(i, jlast - jshift, tx, Point(x, y))
            end

            ilast, jlast = i, j
        end
    end

    # last point
    ilast, jlast = ie, je
    setmin!(ilast, jlast, te, Point(xe, ye))
    return nothing
end

# function update_mintimes!(fm::FaceMesh, rays, αorigin)
#     for (ray_idx, ray) ∈ enumerate(rays)
#         pstart, pend, c, te = ray.pstart, ray.pend, ray.c, ray.time
#         ray_to_cell_mintimes!(fm, pstart, pend, c, te, ray_idx, αorigin)
#     end
# end

abstract type MaterialModel end
struct Fluid <: MaterialModel end
struct Solid <: MaterialModel end

# sinθr(θ₁, c₁, c₂) = (c₂ / c₁) * sin(θ₁)

# function maybeθr(α, θ, c₁, c₂)
#     c₂ == 0.0 && return nothing # out of domain
#     sθr = sinθr(θ, c₁, c₂)
#     abs(sθr) > 1.0 && return nothing # beyong limit angle
#     asin(sθr) # θᵣ
# end

# Compute refracted angle θr from the incident angle θ and the velocities c₁, c₂
function maybeuθr(uθ::Point, c₁, c₂)
    c₂ == 0.0 && return nothing # out of domain
    sθr = (c₂ / c₁) * uθ.y
    abs(sθr) > 1.0 && return nothing # beyong limit angle
    cθr = sqrt(1.0-sθr^2)
    Point(cθr,sθr)
end


maybe(x,f) = isnothing(x) ? nothing : f(x)

# αf(α,θ,θᵣ) = π + α - θ - θᵣ 
# αr(α,θ,θᵣ) = α + θᵣ - θ

cosapb(u::Point,v::Point) = u.x*v.x - u.y*v.y #cos a cos b − sin a sin b
cosamb(u::Point,v::Point) = u.x*v.x + u.y*v.y #cos a cos b + sin a sin b
sinapb(u::Point,v::Point) = u.y*v.x + u.x*v.y #sin a cos b + sin b cos a 
sinamb(u::Point,v::Point) = u.y*v.x - u.x*v.y #sin a cos b - sin b cos a
upv(u::Point,v::Point) = Point(cosapb(u,v),sinapb(u,v))  
umv(u::Point,v::Point) = Point(cosamb(u,v),sinamb(u,v)) 
double_u(u::Point)=upv(u,u)
coordinates(u::Point) = (u.x,u.y)
uαr(uα::Point,uθ::Point,uθᵣ::Point) = upv(uα,umv(uθᵣ,uθ)) # α + θᵣ - θ
uαf(uα::Point,uθ::Point,uθᵣ::Point) = -1.0*umv(uα,upv(uθ,uθᵣ)) # π + α - θ - θᵣ 

r_coeffs(uα, uθ, c₁, c₂, ρ₁, ρ₂) = (c₂*ρ₂*uα.x - c₁*ρ₁*uθ.x)/(c₂*ρ₂*uα.x + c₁*ρ₁*uθ.x)
t_coeffs(uα, uθ, c₁, c₂, ρ₁, ρ₂) = 2*c₁*ρ₁*uα.x/(c₂*ρ₂*uα.x + c₁*ρ₁*uθ.x)


# function amatrix(uθd1,uθs1,uθd2,uθs2,zd1,zd2,zs1,zs2,T::DataType)
#     res =@SMatrix [
#         -uθd1.x  -uθd2.x -uθs1.y +uθs2.y ;
#         -uθd1.y  +uθd2.y +uθs1.x +uθs2.x ;
#         -uθd1.x  -uθd2.x -uθs1.y uθs2.y ;
#         -uθd1.x  -uθd2.x -uθs1.y uθs2.y 
#     ]



# amplitude coefficients (from http://www.fast.u-psud.fr/~martin/acoustique/support/r%C3%A9flection-r%C3%A9fraction.pdf)
function rt_coeffs_fluid_fluid(uθ, uθr, c₁, c₂, ρ₁, ρ₂)
    ci = uθ.x
    ct = uθr.x
    z₁ = c₁*ρ₁
    z₂ = c₂*ρ₂ 
    a = z₂*ci
    b = z₁*ct
    r = (a-b)/(a+b)
    t = 2*z₁*ci/(a+b)
    r,t
end

# power coefficients (from http://www.fast.u-psud.fr/~martin/acoustique/support/r%C3%A9flection-r%C3%A9fraction.pdf)
function rtpower_coeffs_fluid_fluid(uθ, uθr, c₁, c₂, ρ₁, ρ₂)
    ci , ct = uθ.x  , uθr.x
    z₁ , z₂ = c₁*ρ₁ , c₂*ρ₂
    a  , b  = z₂*ci , z₁*ct
    denom = inv(a+b)^2
    r = denom * (a-b)^2
    t = denom * 4*a*b
    r,t
end


# maybeαf(α, θ, c₁, c₂) = maybe(maybeθr(α, θ, c₁, c₂), θᵣ ->  αf(α,θ,θᵣ))
# maybeαr(α, θ, c₁, c₂) = maybe(maybeθr(α, θ, c₁, c₂), θᵣ ->  αr(α,θ,θᵣ))

maybeuαf(uα, uθ, c₁, c₂) = maybe(maybeuθr(uθ, c₁, c₂), uθᵣ ->  uαf(uα,uθ,uθᵣ))
maybeuαr(uα, uθ, c₁, c₂) = maybe(maybeuθr(uθ, c₁, c₂), uθᵣ ->  uαr(uα,uθ,uθᵣ))


# function maybeuθrrt(uθ::Point, c₁, c₂, ρ₁, ρ₂)
#     c₂ == 0.0 && return nothing # out of domain
#     sθr = (c₂ / c₁) * uθ.y
#     abs(sθr) > 1.0 && return nothing # beyong limit angle
#     cθr = sqrt(1.0-sθr^2)
#     uθ=Point(cθr,sθr)
#     r,t=r_coeffs(uα, uθ, c₁, c₂, ρ₁, ρ₂)
#     uθ,r,t
# end

# maybeuαrrt(uα, uθ, c₁, c₂, ρ₁, ρ₂) = maybe(maybeuθrrt(uθ, c₁, c₂, ρ₁, ρ₂), (uθᵣ,r,t) ->  (uαr(uα,uθ,uθᵣ),r,t))



# function fluid_fluid_reflection_transmission(iang, d1, d2, c1, c2)
#     sint =(c2/c1)*sin(iang)
     
#     abs(sθr) > 1.0 && return nothing # beyong limit angle
# end

# cost = sqrt(1-sint.^2).*(sint <= 1) + i.*sqrt(sint.^2 -1).*(sint > 1);
# R = ((c2*d2).*cos(iang) - (c1*d1).*cost)./((c2*d2).*cos(iang) + (c1*d1).*cost);
# T = ((2*d1*c1).*cos(iang))./((c2*d2).*cos(iang) + (c1*d1).*cost);

# longitudina

function normalize(u::Point) 
    nu=norm(u)
    u = u ./ nu
    u
end  


function recursive_generate_rays_fluid!(reduce_ray!, state,celerity_domain, idx, time,u::Point, pstart::Point, np::NumericalParameters, depth, power)
    c₁,_,ρ₁ = cl_ct_rho(state)
    u=normalize(u)
    # c₁ == 0.0 && return 
    # @show cl₁,ρ₁

    pend, pnext,uθ,state,next_state = advance_ds(pstart, u, np, state,celerity_domain)
    # !ok && return
    # @show c₁,pstart,pend,pnext
    # next_state = get_state(pnext,state)
    c₂,_,ρ₂ = cl_ct_rho(next_state)

    Δt = norm(pend - pstart) / c₁   
    time += Δt
    
    reduce_ray!(pstart, pend, c₁, time, idx,power)
    
    if depth <= np.maxdepth && time <= np.tmax  && power>np.power_threshold  # 100 DB attenuation cutoff

        uθr = maybeuθr(uθ, c₁, c₂) # refracted angle θr from incident angle θ      
        r,t = isnothing(uθr) ? (1.0,0.0) : rtpower_coeffs_fluid_fluid(uθr, uθ, c₁, c₂, ρ₁, ρ₂)
        # @assert r+t ≈ 1.0 
        # @show r,t,r+t
        # reflection αf(α,θ,θ) =  π + α - θ - θᵣ 
        recursive_generate_rays_fluid!(reduce_ray!, state,celerity_domain, 2idx, time,uαf(u,uθ,uθ), pend, np, depth+1, r*power) 
        # refraction α + θᵣ - θ
        !isnothing(uθr) && recursive_generate_rays_fluid!(reduce_ray!, next_state,celerity_domain, 2idx+1, time,uαr(u,uθ,uθr), pnext, np, depth+1, t*power)
    end
end

#---------------------------- solid (two speeds) -----------------------------------------
function recursive_generate_rays_solid!(reduce_ray!, state,celerity_domain, is_longitudinal, idx, time,
    u, pstart::Point, np::NumericalParameters, depth, power)
    cl₁,ct₁,ρ₁ = cl_ct_rho(state)
    # @show cl₁,ct₁,ρ₁

    (c₁ , ĉ₁) = is_longitudinal ? (cl₁,ct₁) : (ct₁,cl₁)
    
    ok,pend, pnext,uθ,state,next_state = advance_ds(pstart, u, np, state,celerity_domain)
    !ok && return
    # next_state = get_state(pnext,state)
    cl₂,ct₂,ρ₂ = cl_ct_rho(next_state)
    (c₂ , ĉ₂) = is_longitudinal ? (cl₂,ct₂) : (ct₂,cl₂)

    Δt = norm(pend - pstart) / c₁
    time += Δt
    reduce_ray!(pstart, pend, c₁, time, idx,power)
    # other_celerity = celerity_function == celerity ? transversal_celerity : celerity

    if depth <= np.maxdepth && time <= np.tmax
        # reflection with the same  celerity function c₁
        recursive_generate_rays_solid!(reduce_ray!, state,celerity_domain, is_longitudinal, 4idx-2, time,uαf(u,uθ,uθ), pend, np, depth+1,power)
        # refraction with the same  celerity function c₂
        uθrc₂= maybeuθr(uθ, c₁, c₂)   
        !isnothing(uθrc₂) && recursive_generate_rays_solid!(reduce_ray!, next_state,celerity_domain, is_longitudinal , 4idx-1, time,uαr(u,uθ,uθrc₂), pnext, np, depth+1,power)
        # maybe(maybeuαr(u, uθ, c₁, c₂), uαr1 -> recursive_generate_rays_solid!(reduce_ray!, state, is_longitudinal , 4idx-1, time,uαr1, pnext, np, depth+1, maxtime,power))
        # reflection with the other celerity function ĉ₁
        uθrĉ₁ = maybeuθr(uθ, c₁, ĉ₁) 
        !isnothing(uθrĉ₁) && recursive_generate_rays_solid!(reduce_ray!, state,celerity_domain, !is_longitudinal, 4idx+0, time,uαf(u,uθ,uθrĉ₁), pend, np, depth+1,power)
        # maybe(maybeuαf(u, uθ, c₁, c₃), uαf2 -> recursive_generate_rays_solid!(reduce_ray!, state,!(is_longitudinal), 4idx+0, time,uαf2, pend, np, depth+1, maxtime,power))
        # refraction with the other celerity function ĉ₂
        uθrĉ₂ = maybeuθr(uθ, c₁, ĉ₂) 
        !isnothing(uθrĉ₂) && recursive_generate_rays_solid!(reduce_ray!, next_state,celerity_domain, !is_longitudinal,4idx+1,time,uαr(u,uθ,uθrĉ₂), pnext, np, depth+1,power)

        # c₄, cs4 = other_celerity(pnext, cs)
        # maybe(maybeuαr(u, uθ, c₁, c₄), uαr2 -> recursive_generate_rays_solid!(reduce_ray!, c₄, cs4, 4idx+1, time,uαr2, pnext, np, depth+1, maxtime,power))
    end
end

function recursive_generate_rays!(M::Type{Fluid}, reduce_ray!, state,celerity_domain, incident_idx, time,u, pstart::Point, np::NumericalParameters, depth, power)
    recursive_generate_rays_fluid!(reduce_ray!, state,celerity_domain, incident_idx, time,u,pstart::Point, np::NumericalParameters, depth,power)
end

function recursive_generate_rays!(M::Type{Solid}, reduce_ray!, state,celerity_domain, incident_idx, time,u, pstart::Point, np::NumericalParameters, depth, power)
    recursive_generate_rays_solid!(reduce_ray!, state,celerity_domain,true, incident_idx, time,u, pstart::Point, np::NumericalParameters, depth,power)
end


nrays(::Type{Fluid}, maxdepth) = sum(2^i for i ∈ 0:maxdepth)
nrays(::Type{Solid}, maxdepth) = sum(4^i for i ∈ 0:maxdepth)

function create_rays(M, porigin::Point, αorigin, np::NumericalParameters, celerity_domain)
    pstart = porigin
    istate = initial_state(celerity_domain)
    state = get_state(pstart,istate,celerity_domain)
    cl,_,_ = cl_ct_rho(state)
    rays = [Ray() for i ∈ 1:nrays(M, np.maxdepth)]
    accumulate_ray!(pstart, pend, cl, time, ray_idx,power) = (rays[ray_idx] = Ray(pstart, pend, cl, time,power))
    recursive_generate_rays!(M, accumulate_ray!, state,celerity_domain, Int128(1), 0.0, unit_vector(αorigin), pstart, np, 1,1.0)
    rays
end



function generate_rays!(M, fm::FaceMesh, porigin::Point, αorigin, np::NumericalParameters, celerity_domain::AbstractCelerityDomain)
    pstart = porigin
    istate = initial_state(celerity_domain)
    state = get_state(pstart,istate,celerity_domain)
    cl,_,_ = cl_ct_rho(state)
    accumulate_ray!(pstart, pend, cl, time, ray_idx,power) = ray_to_cell_mintimes!(fm, pstart, pend, cl, time, ray_idx, αorigin,porigin,power,true)
    recursive_generate_rays!(M, accumulate_ray!, state,celerity_domain, Int128(1), 0.0,unit_vector(αorigin), pstart, np, 1,1.0)
end

function generate_rays_fast!(M, fm::FaceMesh, porigin::Point, αorigin, np::NumericalParameters, celerity_domain::AbstractCelerityDomain)
    pstart = porigin
    istate = initial_state(celerity_domain)
    state = get_state(pstart,istate,celerity_domain)
    cl,_,_ = cl_ct_rho(state)
    accumulate_ray!(pstart, pend, cl, time, ray_idx,power) = ray_to_cell_mintimes!(fm, pstart, pend, cl, time, ray_idx, αorigin,porigin,power,false)
    recursive_generate_rays!(M, accumulate_ray!, state,celerity_domain, Int128(1), 0.0,unit_vector(αorigin), pstart, np, 1,1.0)
end

ascendent_idx(::Type{Fluid}, idx) = idx ÷ 2
ascendent_idx(::Type{Solid}, idx) = (idx + 2) ÷ 4

function get_ray_trajectory_from_ntree!(M, ray_idx, rays::Vector{Ray}, tinit)

    firstray = rays[ray_idx]
    pstart, pend = firstray.pstart, firstray.pend
    u = (pend - pstart) / norm(pend - pstart)
    te = firstray.time
    tb = te - norm(pend - pstart) / firstray.c
    dt = (tinit - tb)
    firstpoint = pstart + dt * u * firstray.c

    xs = [firstpoint.x,pstart.x]
    ys = [firstpoint.y,pstart.y]

    idx = ascendent_idx(M, ray_idx)
    while (idx != 0)
        ray = rays[idx]
        # @show idx,ray.pend,ray.pstart 
        push!(xs, ray.pend.x)
        push!(ys, ray.pend.y)
        push!(xs, ray.pstart.x)
        push!(ys, ray.pstart.y)
        idx = ascendent_idx(M, idx)
        # @show idx
    end

    xs, ys
end

function get_ray_from_ray_idx(M, ray_idx, porigin::Point, αorigin, np::NumericalParameters, celerity_domain)
    ray_indexes=Vector{Int128}()
    current_ray_idx=ray_idx
    while current_ray_idx!=1
        push!(ray_indexes,current_ray_idx)
        current_ray_idx = ascendent_idx(M, current_ray_idx)
    end
    push!(ray_indexes,1)
    reverse!(ray_indexes)
    
    # for (i,ray_idx) ∈ enumerate(ray_indexes)
    #     @show i,ray_idx
    # end

 
    depth=1
    ray_trajectory=Vector{Ray}(undef,length(ray_indexes))
    function build_ray_trajectory!(pstart, pend, cl, time, ray_idx,power)
        if depth <=length(ray_indexes)
            if ray_idx==ray_indexes[depth]
                # @show depth,ray_idx
                ray_trajectory[depth]=Ray(pstart, pend, cl, time,power)
                depth+=1
            end
        end
    end


    pstart = porigin
    istate = initial_state(celerity_domain)
    state = get_state(pstart,istate,celerity_domain)
    cl,_,_ = cl_ct_rho(state)


    recursive_generate_rays!(M, build_ray_trajectory!, state,celerity_domain, Int128(1), 0.0, unit_vector(αorigin), pstart, np, 1,1.0)
    
    
    # for (i,ray_idx) ∈ enumerate(ray_trajectory)
    #     @show i,ray_idx
    # end
    return ray_trajectory

end


function get_ray_trajectories(M, fms,source_point, observation_points, np, lcd)
    observation_trajectories = Dict{String,Tuple{Vector{Float64},Vector{Float64}}}()
    source_name = source_point[1]
    @show source_name
    for op ∈ observation_points
        op_name,op_pos = get_name_and_position_of_source(op)
        # @show op_name
        isx, isy = get_cell_indexes(fms[source_name],op_pos.x,op_pos.y) 
        mr = fms[source_name].ray_min_times[isx, isy]
        if mr.idx!=-1
            rα = get_ray_from_ray_idx(M, mr.idx, mr.porigin , mr.αorigin, np, lcd)
            xs = Vector{Float64}()
            ys = Vector{Float64}()
            for i in 1:length(rα)
                ray = rα[i]
                push!(xs, ray.pstart.x)
                push!(ys, ray.pstart.y)
                push!(xs, ray.pend.x)
                push!(ys, ray.pend.y)
            end
            observation_trajectories[op_name] = (xs,ys)
        else
            println("no ray have reached $op")
        end
    end
    # @show observation_points
    observation_trajectories
end

# function get_ray_trajectories(M, fms, observation_points, np, lcd, isx, isy)
#     observation_trajectories = Dict{String,Tuple{Vector{Float64},Vector{Float64}}}()
#     for op ∈ observation_points
#         name = op[1]
#         mr = fms[name].ray_min_times[isx, isy]
#         # @show name,mr.porigin,isx,isy
#         rα = create_rays(M, mr.porigin , mr.αorigin, np, lcd)
#         # @show rα
#         # @show name,length(rα),mr
#         if mr.idx!=-1
#             observation_trajectories[name] = get_ray_trajectory_from_ntree!(M, mr.idx, rα, mr.time)
#         else
#             println("no ray have reached $op")
#         end
#     end
#     # @show observation_points
#     observation_trajectories
# end


































