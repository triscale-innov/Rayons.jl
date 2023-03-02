using JLD2
using JSON
using Dates


function jld2save(jld2_filename::String, d::Dict)
    jldopen(jld2_filename, "w") do f
        for key ∈ keys(d)
            f[key] = d[key]
        end
    end
end

function export_forward(outdir, source_point::Tuple{String,Point,Float64},
    observation_points::Vector{Tuple{String,Point,Point}},
    face_mesh::FaceMesh)

    forward = Dict{String,Any}()

    sname, spos, radius = source_point
    source = []
    push!(source, spos.x, spos.y)
    forward["source"] = source

    arrivals = Dict{String,Float64}()
    for op ∈ observation_points
        opname, op_sw, op_ne = op
        oppos = (op_sw + op_ne)/2
        isx, isy = get_cell_indexes(face_mesh, oppos.x, oppos.y)
        arrival_time = face_mesh.ray_min_times[isx, isy].time
        #fixme : check for Inf
        arrivals[opname] = arrival_time
    end
    forward["arrival times"] = arrivals
    @show source
    @show arrivals

    forward_filename = "forward_" * sname * ".jld2"
    jld2save(joinpath(outdir,forward_filename), forward)
end

function export_all_forward(outdir,
    source_points::Vector{Tuple{String,Point,Float64}},
    observation_points::Vector{Tuple{String,Point,Point}},
    face_meshes::Dict{String,FaceMesh})

    for s ∈ source_points
        export_forward(outdir,s, observation_points, face_meshes[first(s)])
    end
end

function push_point!(bounding_box,s::String,x,y)
    result = []
    push!(result, x)
    push!(result, y)
    bounding_box[s]=result
end

function extract_bounding_box(face_meshes::Dict{String,FaceMesh})
    ffm = last(first(face_meshes))
    xm, xM, ym, yM = first(ffm.xfaces), last(ffm.xfaces), first(ffm.yfaces), last(ffm.yfaces)
    bounding_box = Dict{String,Any}()
    push_point!(bounding_box,"sw",xm,ym)
    push_point!(bounding_box,"ne",xM,yM)
    hx = ffm.xfaces[2] - ffm.xfaces[1]
    hy = ffm.yfaces[2] - ffm.yfaces[1]
    @assert hx ≈ hy

    bounding_box, hx
end


function extract_sensor(observation_point::Tuple{String,Point,Point}, fm::FaceMesh,threshold=0.0)
    opname, op_sw, op_ne = observation_point
    oppos = (op_sw+op_ne)/2 
    sensor = Dict{Any,Any}()
    sensor["sensor"] = [opname, oppos.x, oppos.y]
    nx, ny = get_cell_sizes(fm::FaceMesh)

    if threshold==0.0
        @show opname, nx, ny
        times = [fm.ray_min_times[i, j].time for i ∈ 1:nx, j ∈ 1:ny]
        sensor["times"] = times
    else
        @show threshold,opname, nx, ny
        threshold_times=zeros(nx,ny)
        for i ∈ 1:nx, j ∈ 1:ny
            h = fm.history[i,j]
            nh = length(h)
            for l ∈ 1:nh
                if h[l]>threshold
                    threshold_times[i,j] = (l/nh)*fm.tmax
                    break
                end
            end
        end
        sensor["times"] = threshold_times
    end
    sensor
end


function export_backward(outdir,
    observation_points::Vector{Tuple{String,Point,Point}}, face_meshes::Dict{String,FaceMesh})

    backward = Dict{String,Any}()
    bounding_box, h = extract_bounding_box(face_meshes)
    backward["bounding_box"] = bounding_box
    @show bounding_box, h
    backward["h"] = h
    sensors = Vector{Dict{Any,Any}}()

    minmax_intensity = typemax(Float64)
    for op ∈ observation_points
        opname, _ = op
        fm = face_meshes[opname]
        nx, ny = get_cell_sizes(fm::FaceMesh)
        for i ∈ 1:nx, j ∈ 1:ny
            h = fm.history[i,j]
            max_intensity=0.0
            for l ∈ 1:length(h)
                max_intensity=max(h[l],max_intensity)
            end
            minmax_intensity=min(max_intensity,minmax_intensity)
        end
    end
    threshold = minmax_intensity
    # threshold = 0.0
    @show minmax_intensity,threshold

    for op ∈ observation_points
        opname, _ = op
        fm = face_meshes[opname]
        sensor = extract_sensor(op, fm,threshold)
        push!(sensors, sensor)
    end
    backward["sensors"] = sensors
    @show (backward["bounding_box"]["ne"] .- backward["bounding_box"]["sw"]) ./ backward["h"]

    jld2save(joinpath(outdir,"backward_area.jld2"), backward)
end


function export_data(outdir,json_fname::String,
    source_points::Vector{Tuple{String,Point,Float64}},
    observation_points::Vector{Tuple{String,Point,Point}}, 
    face_meshes::Dict{String,FaceMesh},
    cartesian_material_grid::CartesianMaterialGrid,
    numerical_parameters::NumericalParameters)

    export_all_forward(outdir,source_points,observation_points,face_meshes)
    export_backward(outdir,observation_points,face_meshes)

    open(joinpath(outdir,"parameters.json"), "w") do f
        write(f, json(json_fname, 2))
        write(f, json(numerical_parameters, 2))
    end

    alldata=Dict{String,Any}()

    alldata["face_meshes"]=face_meshes
    alldata["cmg"]=cartesian_material_grid
    alldata["source_points"]=source_points
    alldata["observation_points"]=observation_points
    alldata["numerical_parameters"]=numerical_parameters
    alldata["json_fname"]=json_fname
    jld2save(joinpath(outdir,"alldata.jld2"),alldata)
end

function import_data(outdir)

    alldata= JLD2.load(joinpath(outdir,"alldata.jld2"))

    face_meshes=alldata["face_meshes"]
    cartesian_material_grid=alldata["cmg"]
    source_points=alldata["source_points"]
    observation_points=alldata["observation_points"]
    numerical_parameters=alldata["numerical_parameters"]
    json_fname=alldata["json_fname"]

    face_meshes,cartesian_material_grid,source_points,observation_points,numerical_parameters,json_fname
end






