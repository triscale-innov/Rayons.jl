# Rayons.jl

This project simulates optical rays interacting with a 2D geometry composed of a set of rectangles with different optical velocities.

Minimal travel times from a given source to a set of sensors are computed.

The ray interacting with a (plane) boundary is split into **two sub-rays** (**reflected** and **refracted**) via the following recursive function :

```julia
function recursive_generate_rays_fluid!(reduce_ray!, state,celerity_domain, idx, time,
                        u::Point, pstart::Point, np::NumericalParameters, depth, power)
    c₁,_,ρ₁ = cl_ct_rho(state)
    u=normalize(u)
    pend, pnext,uθ,state,next_state = advance_ds(pstart, u, np, state,celerity_domain)
    c₂,_,ρ₂ = cl_ct_rho(next_state)

    Δt = norm(pend - pstart) / c₁   
    time += Δt
    
    reduce_ray!(pstart, pend, c₁, time, idx,power)
    
    if depth <= np.maxdepth && time <= np.tmax  && power>np.power_threshold

        uθr = maybeuθr(uθ, c₁, c₂) # refracted angle θr from incident angle θ      
        r,t = isnothing(uθr) ? (1.0,0.0) : rtpower_coeffs_fluid_fluid(uθr, uθ, c₁, c₂, ρ₁, ρ₂)
        # reflection αf(α,θ,θ) =  π + α - θ - θᵣ 
        recursive_generate_rays_fluid!(reduce_ray!, state,celerity_domain, 2idx, time,
                                                uαf(u,uθ,uθ), pend, np, depth+1, r*power) 
        # refraction α + θᵣ - θ
        !isnothing(uθr) && 
        recursive_generate_rays_fluid!(reduce_ray!, next_state,celerity_domain, 2idx+1, time,
                                                uαr(u,uθ,uθr), pnext, np, depth+1, t*power)
    end
end
```
Note how different **reduction** functions (`reduce_ray!`) can be passed to the recursive function. In particular, we can pass a function computing the minimal arrival time of a given ray into each cell of a 2D Cartesian grid.

The first example is a simple geometry where the ray source is placed to the left in a fast material (v=1000). A layer of slow material (v=100) is placed in the middle and we place 3 sensors on the right.

The input file `test0.json` describe this geometry :

```json
{
  "material_properties": [
    {"name": "Fast","density": 7900,"longitudinal_velocity": 1000,"transversal_velocity": 1000},
    {"name": "Slow","density": 1000,"longitudinal_velocity": 100,"transversal_velocity": 100} 
  ],
  "name": "Square1",
  "bounding_box": {"sw": [-1.0001,-1.0001],"ne": [1.0001,1.0001]},
  "sensors": [ 
    ["right_up",0.999,0.5],["right_middle",0.999,0.0],["right_down",0.999,-0.5]
  ],
  "rects": [
    {"name": "ExtSquare","material": "Fast","sw": [-1.0,-1.0],"ne": [1.0,1.0]},
    {"name": "L1","material": "Slow","sw": [-0.3,-1.0],"ne": [0.3,1.0]}
  ],
  "sources": [ [-1.0,0.0]]
}
```
The following command launch the computation

```julia
julia> include("go_test0.jl")
````
Producing the following contour plot of arrival times as well as the shortest path from the source to the sensor placed at the opposite side of the maze :

![](docs/test0_min_traj_source_1.png)

This case is simple and only refracted rays are significant.
The following animation illustrates the arrival times through the maze :

![](docs/test0_source_1.gif)


A second example corresponds to a wall with two holes :

```julia
julia> include("go_test1.jl")
````

![](docs/two_holes_min_traj_source_1.png)

In this case both reflected and refracted rays must be considered.

he following animation illustrates the arrival times through the maze :

![](docs/two_holes_source_1.gif)


A last example corresponding to a maze can be launched via the following script

```julia
julia> include("go_maze.jl")
````
Producing the following contour plot of arrival times as well as the shortest path from the source to the sensor placed at the opposite side of the maze :

![](docs/min_traj_source_1.png)

This case is super simple because the velocity of the walls is super small. In such case, only reflected rays are significant.
The following animation illustrates the arrival times through the maze :

![](docs/maze_source_1.gif)


