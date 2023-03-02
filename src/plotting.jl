using CairoMakie
using PrettyTables


function plot_vcells!(ax,xf,yf)
    xys_b = ((xfi,yf[1]) for xfi ∈ xf)
    xys_e = ((xfi,yf[end]) for xfi ∈ xf)
    xys_n = ((NaN,NaN) for xfi ∈ xf) #NaN separates segments in Makie
    xys = collect(Iterators.flatten(zip(xys_b,xys_e,xys_n))) # [(),(),()]
    lines!(ax,xys)
end  
function plot_hcells!(ax,xf,yf)
    xys_b = ((xf[1],yfi) for yfi ∈ yf)
    xys_e = ((xf[end],yfi) for yfi ∈ yf)
    xys_n = ((NaN,NaN) for yfi ∈ yf) #NaN separates segments in Makie
    xys = collect(Iterators.flatten(zip(xys_b,xys_e,xys_n))) # [(),(),()]
    lines!(ax,xys)
end    
function plot_cells!(ax,xf,yf)
    plot_vcells!(ax,xf,yf)
    plot_hcells!(ax,xf,yf)
end

function plot_cells!(ax,fm::FaceMesh)
    plot_cells!(ax,fm.xfaces,fm.yfaces)
end

function save_to_file(outdir,basefilename,fig)
    extensions=(".png",)
    # extensions=(".png",".svg")
    # out_dir="out"
    mkpath(outdir)
    for e ∈ extensions
        filename=joinpath(outdir,basefilename*e)
        print("saving $filename ")
        save(filename,fig)
    end
    println("")
end


function plot_face!(ax,pts,face)
    colors=[v==Inf ? :red : :black for v ∈ min_times(face)]
    # texts=[string(round(10*v,digits=3)) for v ∈ min_times(face)]

    cpts=collect(pts)
    values=collect(min_times(face))
    # @show cpts
    # @show values
    # cpts_ends=[(x+0.1v,y+0.1v) for ((x,y),v) ∈ zip(cpts,values)]
    # cpts_NaN=[(NaN,NaN) for p ∈ cpts]
    # xys = collect(Iterators.flatten(zip(cpts,cpts_ends,cpts_NaN))) # [(),(),()]
    scatter!(ax,cpts,color = colors)
    # @show cpts_ends
    # Makie.lines!(ax,xys,color=:orange,linewidth=5)
    # Makie.text!(ax,texts,position=cpts)

end

plot_vface!(ax,fm::FaceMesh,i,j) = plot_face!(ax,get_vface_points(fm,i,j),fm.vfaces[i,j])
plot_hface!(ax,fm::FaceMesh,i,j) = plot_face!(ax,get_hface_points(fm,i,j),fm.hfaces[i,j])

function plot_vfaces!(ax,fm::FaceMesh)
    vnx,vny=size(fm.vfaces)
    for i ∈ 1:vnx, j ∈ 1:vny
        plot_vface!(ax,fm,i,j)
    end
end


function plot_hfaces!(ax,fm::FaceMesh)
    vnx,vny=size(fm.hfaces)
    for i ∈ 1:vnx, j ∈ 1:vny
        plot_hface!(ax,fm,i,j)
    end
end

function plot_trajectory!(ax,fm::FaceMesh,faceid::FaceId,l)

    xs=Vector{Float64}()
    ys=Vector{Float64}()
    zs=Vector{Float64}()

    while true
        face = get_face(fm,faceid)
        xmin,ymin = fm.xfaces[faceid.i],fm.yfaces[faceid.j]
        d=fm.lfaces[l]
        (x,y) = faceid.vertical ? (xmin,ymin+d) : (xmin+d,ymin)
        minpath = face.minpaths[l]
        z = minpath.time
        push!(xs,x)
        push!(ys,y)
        push!(zs,z)
        pred = minpath.pred
        # @show minpath.time
        if pred.index==-1
            # println("Final Time : ",minpath.time)
            break
        end
        faceid=pred.faceid
        l = pred.index
    end
    println("Final Time : ",zs[1])


    lines!(ax,xs,ys,linewidth=2)
end


function plot_contour!(ax,fm::FaceMesh)

    
    # vfaces
    (nxp,ny) = size(fm.vfaces)
    lfaces = fm.lfaces
    m = length(lfaces)
    # xs=zeros(nxp*m,ny)
    # ys=zeros(nxp*m,ny)
    zs=zeros(nxp*m,ny)
    for i ∈ 1:nxp
        for j ∈ 1:ny
            points=get_vface_points(fm,i,j)
            values=min_times(fm.vfaces[i,j])
            for (l,p) ∈ enumerate(points)
                (x,y) = p
                # xs[(i-1)*m+l,j]=x
                # ys[(i-1)*m+l,j]=y
                zs[(i-1)*m+l,j]=fm.vfaces[i,j].minpaths[l].time
            end
        end
    end

    xs = Vector{Float64}()
    for i ∈ 1:nxp
        xmin=fm.xfaces[i]
        for l ∈ 1:m
            push!(xs,xmin+lfaces[l])
        end
    end

    ys = Vector{Float64}()
    for j ∈ 1:ny
        push!(ys,fm.yfaces[j])
    end


    # @show typeof(xs),xs
    # @show typeof(xs),ys
    # @show typeof(xs),zs

   surface!(ax,xs,ys,zs)
end

function plot_materials!(ax,fm::FaceMesh)
    (nx,ny) = get_cell_sizes(fm)
    xcs=[(fm.xfaces[i]+fm.xfaces[i+1])/2 for i ∈ 1:nx]
    ycs=[(fm.yfaces[j]+fm.yfaces[j+1])/2 for j ∈ 1:ny]
    heatmap!(ax,xcs,ycs,fm.speed_l)
end



function plot_ray!(ax,ray)
    xs=[ray.pstart.x,ray.pend.x]
    ys=[ray.pstart.y,ray.pend.y]
   scatter!(ax,xs,ys,linewidth=5)
   lines!(ax,xs,ys,linewidth=5)
end

function replace_infbynan(x) 
     (isinf(x)||x==0.0) ? NaN : x
end


function extrema_skip_NaNs(rmt)
    minMT=typemax(eltype(rmt))
    maxMT=typemin(eltype(rmt))
    for mt ∈ rmt
        if !(isnan(mt))
            minMT=min(minMT,mt)
            maxMT=max(maxMT,mt)
        end
    end
    minMT,maxMT
    # @show minMT,maxMT
end



function plot_mintimes!(ax,fm::FaceMesh,nlevels,maxlevel,tmax)
    (nx,ny) = get_cell_sizes(fm)
    xcs=[(fm.xfaces[i]+fm.xfaces[i+1])/2 for i ∈ 1:nx]
    ycs=[(fm.yfaces[j]+fm.yfaces[j+1])/2 for j ∈ 1:ny]
    min_times = getproperty.(fm.ray_min_times,:time)
    rmt = map(replace_infbynan,min_times)

    cm=[RGBAf(c.r, c.g, c.b, 1.0) for c in to_colormap(:lightrainbow)]

    minMT,maxMT = extrema_skip_NaNs(rmt)
    maxMT=min(maxMT,tmax)

    levels= collect(range(minMT, maxMT, length=nlevels))
    tightlimits!(ax) 

    mtl = map(x->x>levels[maxlevel] ? NaN : x,rmt)
   

    hm=contourf!(ax,xcs,ycs,mtl,levels=levels,colormap=cm)
    contour!(ax,xcs,ycs,mtl,levels=levels,color=:lavender,colormap=cm)

    # hm=contourf!(ax,xcs,ycs,rmt,levels=20,colormap=cm)
    # contour!(ax,xcs,ycs,rmt,levels=19,color=:lavender,colormap=cm)

    # hm=contourf!(ax,xcs,ycs,rmt,levels = 0.0:0.1:0.9,mode=:relative,colormap=cm)
    # contour!(ax,xcs,ycs,rmt,levels=19,color=:lavender,colormap=cm)

    return hm
end


function plot_ray_trajectory!(ax,xs,ys)
   lines!(ax,xs,ys,linewidth=7,color=:black)
    lines!(ax,xs,ys,linewidth=7,color=:white,linestyle = :dash)
    
    lines!(ax,xs,ys,linewidth=5)
   scatter!(ax,xs,ys,markersize=10,color=:black)
   scatter!(ax,xs,ys,markersize=8,color=:red)
end


function plot_rays!(ax,rays)
  
    xs=Vector{Float64}()
    ys=Vector{Float64}()
    ts=Vector{Float64}()
    for ray ∈ rays
        ray.time==-1 && continue
        push!(xs,ray.pstart.x)
        push!(xs,ray.pend.x)
        push!(xs,NaN)

        push!(ys,ray.pstart.y)
        push!(ys,ray.pend.y)
        push!(ys,NaN)

        push!(ts,ray.time)
        push!(ts,ray.time)
        push!(ts,ray.time)
    end

   scatter!(ax,xs,ys,linewidth=10,color=:white)
   scatter!(ax,xs,ys,linewidth=4,color=:orange)
   lines!(ax,xs,ys,linewidth=2,color=:red)
   scatter!(ax,xs,ys,linewidth=20,markersize=15,color=:orange,strokecolor=:white,strokewidth=2.0)
   scatter!(ax,xs,ys,linewidth=20,markersize=10,color=:orange,strokecolor=:black,strokewidth=2.0)
end

function plot_rays(outdir,fm::FaceMesh,rays,withcontour)

    fig = Figure(resolution=(24000,21000))
    ax = Axis(fig[1, 1],rightspinevisible = false, xlabel="x",ylabel="y")
    ax.aspect = DataAspect()

    plot_cells!(ax,fm.xfaces,fm.yfaces)
    plot_rays!(ax,rays)

    Colorbar(fig[2, 1], hm,label = "Arrival Time (s)")
    save_to_file(outdir,"arrival",fig)

    fig
end

function plot_mintimes(fm::FaceMesh,tmax,withcontour=true)

    fig = Figure(resolution=(1200,1050))
    # ga = fig[1, 1] = GridLayout()
    # gb = fig[2, 1] = GridLayout()
    ax = Makie.Axis(fig[1, 1],rightspinevisible = false, xlabel="x",ylabel="y")
    ax.aspect = DataAspect()
    # xlims!(ax,fm.xfaces[begin],1.05*fm.xfaces[end])
    xlims!(ax,fm.xfaces[begin],1.05*fm.xfaces[end])
    nlevels=40
    maxlevel=nlevels
    hm = plot_mintimes!(ax,fm::FaceMesh,nlevels,maxlevel,tmax)
    # Colorbar(fig[2, 1], hm)
    Colorbar(fig[2, 1], hm,label = "Arrival Time (s)", vertical = false,
    flipaxis = false, ticksize=15, tickalign = 1, width = Relative(3.5/4))
    fig,ax
end

function plot_animated_mintimes(outdir,source_name,source_center,fms,cmg::CartesianMaterialGrid,tmax)
    @show source_name,source_center
    fm = fms[source_name]
    fig = Figure(resolution=(1200,1050))
    ax = Makie.Axis(fig[1, 1],rightspinevisible = false, xlabel="x",ylabel="y")
    ax.aspect = DataAspect()
    xlims!(ax,fm.xfaces[begin],1.05*fm.xfaces[end])
    nlevels=40

    (nx,ny) = get_cell_sizes(fm)
    xcs=[(fm.xfaces[i]+fm.xfaces[i+1])/2 for i ∈ 1:nx]
    ycs=[(fm.yfaces[j]+fm.yfaces[j+1])/2 for j ∈ 1:ny]
    min_times = getproperty.(fm.ray_min_times,:time)
    rmt = map(replace_infbynan,min_times)

    cm=[RGBAf(c.r, c.g, c.b, 1.0) for c in to_colormap(:lightrainbow)]

    minMT,maxMT = extrema_skip_NaNs(rmt)
    maxMT = min(maxMT,tmax)

    levels= collect(range(minMT, maxMT, length=nlevels))
    tightlimits!(ax) 

    mtl = Observable(copy(rmt))

    hm=contourf!(ax,xcs,ycs,mtl,levels=levels,colormap=cm)
    contour!(ax,xcs,ycs,mtl,levels=levels,color=:lavender,colormap=cm)
    Colorbar(fig[2, 1], hm,label = "Arrival Time (s)", vertical = false,
    flipaxis = false, ticksize=15, tickalign = 1, width = Relative(3.5/4))
    draw_rectangles(ax,cmg)
    text!(ax,"S",position=(source_center.x,source_center.y),align=(:right,:top),glowwidth=1.0,offset = (-20, 0),
        space = :data,
        fontsize=50)

    @show nlevels
    filename=joinpath(outdir,source_name*".gif")
    record(fig,filename,2:nlevels;framerate = 10) do ml
        mtl[] = map(x->x>levels[ml] ? NaN : x,rmt)
    end
    nothing
end


function plot_opmintimes(outdir,observation_points,op_mintimes)
    names = first.(observation_points)
    xs = getproperty.(last.(observation_points),:x)
    ys = getproperty.(last.(observation_points),:y)

    table=hcat(names,xs.* 1.e3,ys.* 1.e3,op_mintimes .* 1.e6)

    nr=length(observation_points)

    table=table[sortperm(table[:, 4]), :]
    headers=(["obs point ","x","y","Arrival Time"],
             ["name","[mm]","[mm]","[μs]"])
    pretty_table(table,header=headers)
    fig = Figure(resolution = (800, 600))
    ax = Axis(fig,xticks=(1:6,table[:,1]),
        xlabel="Observation Points",ylabel="Arrival Times (μs)")
    barplot!(ax,1:nr,Vector{Float32}(table[:,4]),color = 1:nr,colormap = (:Spectral_10, 0.85))
    fig[1,1] = ax
    save_to_file(outdir,"omintimes",fig)
    fig
end

function draw_rectangles(ax,cmg::CartesianMaterialGrid,transparency_ratio=0.1)
    rectangles=cmg.rectangles
    polygons = Vector{Makie.Polygon}()
    for r ∈ rectangles
        xm,ym,xM,yM=r.xmin,r.ymin,r.xmax,r.ymax
        push!(polygons,Makie.Polygon([Point2f(p) for p ∈ ((xm,ym),(xM,ym),(xM,yM),(xm,yM))]))
    end
    pcolors =[get_material_color(r.material) for r ∈ rectangles]
    pcolors = [(mc,transparency_ratio) for mc ∈ pcolors]
    poly!(ax,polygons,color = pcolors,strokewidth = 0.5)
end

function draw_rectangles_contour(ax,cmg::CartesianMaterialGrid,line_color=:white)
    rectangles=cmg.rectangles
    xs = Vector{Float64}()
    ys = Vector{Float64}()
    for r ∈ rectangles
        xm,ym,xM,yM=r.xmin,r.ymin,r.xmax,r.ymax
        append!(xs,[xm,xM,xM,xm,xm,NaN])
        append!(ys,[ym,ym,yM,yM,ym,NaN])
    end
    lines!(ax,xs,ys,linewidth = 0.05,color=line_color)
end

function plot_all(outdir,source_name,source_center,fms,observation_points,observation_trajectories,cmg::CartesianMaterialGrid,tmax)
    @show source_name,source_center
    fm = fms[source_name]
    fig, ax = plot_mintimes(fm,tmax)

    for (name,(xs,ys)) ∈ observation_trajectories
        plot_ray_trajectory!(ax,xs,ys)
    end

    draw_rectangles(ax,cmg)

    text!(ax,"S",position=(source_center.x,source_center.y),align=(:right,:top),glowwidth=1.0,offset = (-20, 0),
        space = :data,
        fontsize=50)
    # ysn = [(name, p.y) for (name, p) ∈ observation_points if (name != "bottom_middle" && name !="source")]
    # ticknames, tickys = first.(ysn), last.(ysn)
    # ax2 = Axis(fig[1, 1], yticks = (tickys, ticknames))

    # ax2.yaxisposition = :right
    # ax2.yticklabelalign = (:left, :center)
    # ax2.xticklabelsvisible = false
    # ax2.xticklabelsvisible = false
    # ax2.xlabelvisible = false
    # ax2.aspect = DataAspect()
    # xlims!(ax2, fm.xfaces[1], fm.xfaces[end])
    # ylims!(ax2, fm.yfaces[1], fm.yfaces[end])

    save_to_file(outdir,"min_traj_"*source_name, fig)
    fig
end

function animate_facemesh(outdir,source_name,fm::FaceMesh,cmg::CartesianMaterialGrid)
    
    nx,ny = get_cell_sizes(fm)
    nframe = length(fm.history[1,1])

    xcs=[(fm.xfaces[i]+fm.xfaces[i+1])/2 for i ∈ 1:nx]
    ycs=[(fm.yfaces[j]+fm.yfaces[j+1])/2 for j ∈ 1:ny]
    
  
    fig = Figure(resolution=(500,400))
    ax = Axis(fig[1, 1],rightspinevisible = false, xlabel="x",ylabel="y")

    frame = rand(nx,ny)
    draw_rectangles(ax,cmg)
    hm =heatmap!(ax,xcs,ycs,frame,colorrange=(0.0,0.5),colormap=:viridis, interpolate=true)
    record(fig, joinpath(outdir,"animate_"*source_name*".gif"), 1:nframe, framerate = 10) do frame_index
        hm[3] = [fm.history[i,j][frame_index] for i ∈ 1:nx, j ∈ 1:ny]
        draw_rectangles_contour(ax,cmg)
    end

    fig,ax
end

function animate_mintimes(outdir,source_name,fm::FaceMesh,cmg::CartesianMaterialGrid)
    
    nx,ny = get_cell_sizes(fm)
    nframe = length(fm.history[1,1])

    xcs=[(fm.xfaces[i]+fm.xfaces[i+1])/2 for i ∈ 1:nx]
    ycs=[(fm.yfaces[j]+fm.yfaces[j+1])/2 for j ∈ 1:ny]
    
  
    fig = Figure(resolution=(800,800))
    ax = Makie.Axis(fig[1, 1],rightspinevisible = false, xlabel="x",ylabel="y")

    min_times=[NaN for i ∈ 1:nx, j ∈ 1:ny]
    threshold=0.001

    function update!(frame_index)
        for i ∈ 1:nx, j ∈ 1:ny
            old_value = min_times[i,j]
            intensity = fm.history[i,j][frame_index]
            if (isnan(old_value) && (intensity>threshold))
                min_times[i,j] = frame_index
            end
        end
    end

    frame = rand(nx,ny)
    # draw_rectangles(ax,cmg)

    cm=[RGBAf(c.r, c.g, c.b, 0.99) for c in to_colormap(:lightrainbow)]
    update!(1)
    # hm = Makie.heatmap!(ax,xcs,ycs,frame,colorrange=(0.0,180),colormap=cm)
    hm =heatmap!(ax,xcs,ycs,frame,colorrange=(0.0,nframe/2),colormap=cm)
    # hm = Makie.contour!(ax,xcs,ycs,frame)

    record(fig, joinpath(outdir,s"animate_arrival"*source_name*".gif"), 1:nframe, framerate = 10) do frame_index
        update!(frame_index)
        hm[3] = min_times
        draw_rectangles_contour(ax,cmg,:black)
    end

    fig,ax
end


















    


