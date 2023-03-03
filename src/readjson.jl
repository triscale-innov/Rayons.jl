using Rayons
 using JSON
using CairoMakie.GeometryBasics: Polygon,Point2f
using Colors


struct Rectangle
    material::String
    xmin::Float64
    ymin::Float64
    xmax::Float64
    ymax::Float64
    area::Float64
end

function getRectangles(jdata)
    jrectangles=jdata["rects"]

    rectangles=Vector{Rectangle}()
    for jr ∈ jrectangles
        xmin,ymin = jr["sw"][1],jr["sw"][2]
        xmax,ymax = jr["ne"][1],jr["ne"][2]
        material = jr["material"]
        area = (ymax-ymin)*(xmax-xmin)
        push!(rectangles,Rectangle(material,xmin,ymin,xmax,ymax,area))
    end
    sort(rectangles,by= r-> r.area)
end

struct CartesianMaterialGrid
    xnodes::Vector{Float64}
    ynodes::Vector{Float64}
    material_grid::Matrix{Int}
    materials::Vector{String}
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    ncx::Int
    ncy::Int
    rectangles::Vector{Rectangle}
end

function get_material_from_xy(x,y,rectangles)

    for r ∈ rectangles
        ( x >= r.xmin ) && (x<=r.xmax) && ( y >= r.ymin ) && (y<=r.ymax) && return r.material
    end
    return "unknown"
end


function CartesianMaterialGrid(jdata)
    rectangles=getRectangles(jdata)
    materials = unique([r.material for r ∈ rectangles])
    materials = vcat(materials,"unknown")

    material_indexes = Dict{String,Int}()
    for (i,m) ∈ enumerate(materials)
        material_indexes[m] = i
    end

    xs = sort(unique(vcat([r.xmin for r ∈ rectangles],[r.xmax for r ∈ rectangles])))
    ys = sort(unique(vcat([r.ymin for r ∈ rectangles],[r.ymax for r ∈ rectangles])))

    nx,ny=length.((xs,ys))
    material_grid=zeros(Int,nx-1,ny-1)

    for i ∈ 1:nx-1
        xm = (xs[i]+xs[i+1])/2
        for j ∈ 1:ny-1
            ym = (ys[j]+ys[j+1])/2
            mm = get_material_from_xy(xm,ym,rectangles)
            material_grid[i,j] = material_indexes[mm]
        end
    end

    xmin,xmax=first(xs),last(xs)
    ymin,ymax=first(ys),last(ys)

    ncx=length(xs)-1
    ncy=length(ys)-1

    CartesianMaterialGrid(xs,ys,material_grid,materials,xmin,xmax,ymin,ymax,ncx,ncy,rectangles)
end

get_cell_sizes(cmg::CartesianMaterialGrid) = cmg.ncx,cmg.ncy

@inline valid_index(i,nx) = (i>0) && (i<nx)
@inline inside_xbounds(x,i,xnodes) = (@inbounds xnodes[i]<= x) && (x < @inbounds xnodes[i+1])


function getgridindex(x,xnodes,io)
    nx = length(xnodes)
    # @show io
    valid_index(io  ,nx) && inside_xbounds(x,io  ,xnodes) && return io
    valid_index(io+1,nx) && inside_xbounds(x,io+1,xnodes) && return io+1
    valid_index(io-1,nx) && inside_xbounds(x,io-1,xnodes) && return io-1

    for (i,xi) ∈ enumerate(xnodes)
        # @show x,xi,i
        (xi>x) && return i-1
    end
    return 0
end


out_of_bounds(x,y,cmg) = (x < cmg.xmin) || (x >= cmg.xmax) || (y < cmg.ymin) || (y >= cmg.ymax)

function get_cartesian_indexes_from_xy(x,y,cmg::CartesianMaterialGrid,io,jo)
    i = getgridindex(x,cmg.xnodes,io) 
    j = getgridindex(y,cmg.ynodes,jo)
    return i,j
end


function color_map()
    cm = Dict(
        "Water"=>colorant"white",
        "Steel"=>colorant"steelblue"
    )
end

function get_material_color(material_name) 
    cm = color_map()
    mn = lowercase(material_name)
    c = haskey(cm,mn) ? cm[mn] : colorant"lightgreen"
    c
end

function draw_material_grid(cmg::CartesianMaterialGrid)
    (xs,ys) = (cmg.xnodes,cmg.ynodes)

    polygons = Vector{Polygon}()
    colors_indexes = Vector{Int}()

    nx,ny=length.((xs,ys))

    for i ∈ 1:nx-1
        xm,xM = xs[i],xs[i+1]
        for j ∈ 1:ny-1
            ym,yM = ys[j],ys[j+1]
            push!(polygons,Polygon([Point2f(p) for p ∈ ((xm,ym),(xM,ym),(xM,yM),(xm,yM))]))
            push!(colors_indexes,cmg.material_grid[i,j])
        end
    end

    mcolors = [get_material_color(cmg.materials[i]) for i ∈ colors_indexes]

    @show mcolors

    poly(polygons,color = mcolors,strokewidth = 0.5, colormap = :Hiroshige,
        figure = (; resolution = (1200, 1200), fontsize = 22),
        axis = (; aspect = DataAspect(), title = "Cartesian Material Grid"))
end


