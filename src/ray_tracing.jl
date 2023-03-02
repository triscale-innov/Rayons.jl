using BenchmarkTools
using FLoops
using CairoMakie
# using ProgressMeter

function nl(lmin,h)
    n = Int(ceil(lmin/h))
    l = n*h
    n,l
end

function get_face_mesh(cmg::CartesianMaterialGrid,numerical_parameters::NumericalParameters,celerity_domain::AbstractCelerityDomain) 
    lxmin,lymin = cmg.xmax-cmg.xmin, cmg.ymax-cmg.ymin
    nmin = numerical_parameters.nmin
    h = min(lxmin,lymin)/nmin
    nx,lx = nl(lxmin,h)
    ny,ly = nl(lymin,h)

    nx,ny,h,lx,ly
    FaceMesh(cmg.xmin, cmg.ymin,nx,ny,h, 20, celerity_domain,numerical_parameters.tmax,numerical_parameters.history_size)
end

struct Problem
    material_parameters::Vector{MaterialParameters}
    cartesian_material_grid::CartesianMaterialGrid
    celerity_domain::BoxCelerityDomain
    observation_points::Vector{Tuple{String,Point,Point}}
    source_points::Vector{Tuple{String,Point,Float64}}
    face_mesh::FaceMesh
end

function read_observation_points(jdata)
    ops=Vector{Tuple{String,Point,Point}}()
    for sensor ∈ jdata["sensors"]
        name = sensor[1]
        x,y = sensor[2],sensor[3]
        push!(ops,(name,Point(x,y),Point(x,y)))
    end
    ops
end

function read_sources_points(jdata)
    source_points=Vector{Tuple{String,Point,Float64}}()
    for (i,source) ∈ enumerate(jdata["sources"])
        name = "source_"*string(i)
        x,y = source[1],source[2]
        push!(source_points,(name,Point(x,y),0.002))
    end
    source_points
end

function init_problem_from_json(fname::String,numerical_parameters::NumericalParameters)
    jdata=load(fname)
    material_parameters=read_material_parameters(jdata)
    cmg=CartesianMaterialGrid(jdata)
    celerity_domain = BoxCelerityDomain(cmg,material_parameters)
    face_mesh = get_face_mesh(cmg,numerical_parameters,celerity_domain)
    observation_points = read_observation_points(jdata)
    sources_points = read_sources_points(jdata)
    Problem(material_parameters,cmg,celerity_domain,observation_points,sources_points,face_mesh)
end

function test_trajectories(outdir,fname)
    numerical_parameters = NumericalParameters(ds=2.e-3, ϵ=1.e-10,nmin=50,maxdepth=4,nα=10000,model=Solid)
    p=init_problem_from_json(fname,numerical_parameters)
    model = numerical_parameters.model
    op = p.observation_points[2]
    α=π/4
    generate_rays!(model,p.face_mesh,op[2],α, numerical_parameters, p.celerity_domain)
    rays=create_rays(model,op[2], α, numerical_parameters, p.celerity_domain)
    plot_rays(outdir,p.face_mesh,rays,false)
end
    
function create_outdir(json_fname::String,np::NumericalParameters)
    bn=replace(basename(json_fname),".json"=>"")
    np=create_out_dir_from_np(np)  
    outdir = joinpath("out",bn,np)
    mkpath(outdir)
    println("creating output directory : ",outdir)
    return outdir
end


function load_observation_points_from_forward(fn)
    forward=JLD2.load(fn)

    ops=Vector{Tuple{String,Point,Point}}()
    for o in forward["observation_points"]
        @show o
        xbl,ybl=o["bottom_left"]
        xtr,ytr=o["top_right"]
        push!(ops,(o["name"],Point(xbl,ybl),Point(xtr,ytr)))
    end
   ops
end

function get_name_and_position_of_source(op::Tuple{String,Point,Float64})
    return op[1],op[2]
end

function get_name_and_position_of_source(op::Tuple{String,Point,Point})
    return op[1],(op[2]+op[3])/2
    # x = op[2].x+rand()*(op[3].x-op[2].x)
    # y = op[2].y+rand()*(op[3].y-op[2].y)
    # return op[1],Point(x,y)
end

# function Rayons(fname::String,numerical_parameters::NumericalParameters,forward_name="")

#     println("entering Rayons2D")
#     outdir=create_outdir(fname,numerical_parameters)
#     test_trajectories(outdir,fname)

   
#     problem =  init_problem_from_json(fname,numerical_parameters)
#     celerity_domain,empty_face_mesh,cmg = problem.celerity_domain,problem.face_mesh,problem.cartesian_material_grid
#     model = numerical_parameters.model

#     face_meshes = Dict{String,FaceMesh}()

#     # return face_meshes,cmg


#     source_points = problem.source_points


#     observation_points = isfile(forward_name) ?  load_observation_points_from_forward(forward_name) : problem.observation_points
#     @show isfile(forward_name)

#     all_sources = vcat(source_points,observation_points)

#     print("[")
#     tray = @elapsed begin

#         Threads.@threads for op ∈ all_sources
#             op_name,op_pos = get_name_and_position_of_source(op)
#             print(op_name,"...")
#             face_meshes[op_name] = deepcopy(empty_face_mesh) # one FaceMesh per observation point
#             for α ∈ range(0,stop=2pi,length=numerical_parameters.nα)
#                 generate_rays!(model,face_meshes[op_name], op_pos, α, numerical_parameters, celerity_domain)
#             end
#         end
#     end
#     println("]")

#     nr = Rayons.nrays(model,numerical_parameters.maxdepth)
#     @show tray,(length(observation_points)*numerical_parameters.nα*nr)/(tray*1.e6) 

#     for source_point in source_points
#         source_name,source_center,radius = source_point
#         isx, isy = get_cell_indexes(empty_face_mesh,source_center.x,source_center.y) 
#         observation_trajectories = get_ray_trajectories(model,face_meshes,observation_points,numerical_parameters,celerity_domain,isx,isy)
#         plot_all(outdir,source_name,source_center,face_meshes,observation_points,observation_trajectories,cmg)
#         @show source_name
#         animate_facemesh(outdir,source_name,face_meshes[source_name],cmg)
#         animate_mintimes(outdir,source_name,face_meshes[source_name],cmg)
#     end

#     export_data(outdir,fname,source_points,observation_points,face_meshes,cmg,numerical_parameters)

#     # op_mintimes = [face_meshes[op[1]].ray_min_times[isx, isy].time for op ∈ observation_points]
#     # plot_opmintimes(outdir,observation_points, op_mintimes)
#     face_meshes,cmg
# end


function get_randα(αproba,maxαproba)
    α,y,ys=0.0,1.0,0.0
    fαs=length(αproba)/(2π)
    i=0
    while y>ys && i<100
        α,y = rand()*2π,rand()*maxαproba
        αidx=Int(floor(fαs*α))+1
        ys=αproba[αidx]
        # @show ys,y
        i+=1
    end

    i>30 && @show i
    # @show i,α
    return α
end

function get_randα(hinv)
    x=rand()
    n=length(hinv)
    idx=Int(ceil(x*n))
    
    ϵ = 1.e-14
    δα=2π/n

    idxmin = max(idx-1,1)
    idxmax = min(idx+1,n)

    hidxmin = max(min(hinv[idxmin],hinv[idx]-1),0)
    hidxmax = min(max(hinv[idxmax],hinv[idx]+1),n)

    αmin = hidxmin*δα+ϵ
    αmax = hidxmax*δα-ϵ
    α=αmin + (αmax-αmin)*rand()

    return α
end

# build_inv_proba(αproba,maxαproba)


function update_nsi!(nsi,α)
    n=length(nsi)
    fαs=n/(2π)
    αidx=Int(floor(fαs*α))+1
    αidx=max(min(αidx,n),1)
    nsi[αidx]+=1
    return nothing
end

# function update_αproba!(αproba,αscore,nsi)
#     αproba .= 0.0
#     w = 4
#     n = length(αscore)
#     for i ∈ 1:n

# end

function update_αproba!(αproba,αscore,nsi)
    αproba .= 0.0
    w = 1
    n = length(αscore)
    for i ∈ 1:n
        αsi = αscore[i]
        if  nsi[i]>0.0 
            αsi /= nsi[i]
        end
        for j=max(1,i-w):min(i+w,n)
            αproba[j] += αsi*exp(-(i-j)^2/w^2)
        end
    end
    return nothing
end
function update_distribution!(αproba,h,hinv)
    n = length(αproba)
    sα = sum(αproba)
    αproba .*= n*inv(sα)

    h[1]=αproba[1]
    for i ∈ 2:n
        h[i] = max(min(αproba[i]+h[i-1],n),1)
    end
    # @show h[n]

    for i ∈ 1:n
        j=1
        while (i-0.5)>=h[j]
            j+=1
        end
        hinv[i]=j
    end
    return nothing
end

function ray_tracing!(face_mesh,op_name,op_pos,problem,numerical_parameters)
    model = numerical_parameters.model

    celerity_domain= problem.celerity_domain
    # face_mesh = deepcopy(problem.face_mesh)

    # op_name,op_pos = get_name_and_position_of_source(op)
    # print(op_name,"...")

    αscore = face_mesh.αscore
    αscore .= 0.1

    αproba=ones(length(αscore))
    h=ones(length(αscore))
    hinv=ones(length(αscore))
    update_distribution!(αproba,h,hinv)
    nsi=zeros(length(αscore))

    nαbatch = min(numerical_parameters.nα,1000)
    nbatch = Int(round(numerical_parameters.nα/nαbatch))

    # @show nαbatch,nbatch

    # nbatch = 100
    αbatch = zeros(nαbatch)
    for batch_idx ∈ 1:nbatch
        # @show op_name,batch_idx 
    
        αscore .= 0.1
        for αidx ∈ 1:nαbatch
            α=get_randα(hinv)
            αbatch[αidx]=α
            update_nsi!(nsi,α)

            generate_rays_fast!(model,face_mesh, op_pos, α, numerical_parameters, celerity_domain)
        end

        # maxαproba = maximum(αproba)    
        # minαproba = minimum(αproba) 
        # # @show maxαproba,minαproba
        # if Threads.threadid()==1
        #     fig = Figure(resolution=(1200,1050))
        #     slabel="op ="*op_name*" batch="*string(batch_idx)*" /"*string(nbatch)*"\n"
        #     ax0 = Axis(fig[1, 1],title=slabel)
        #     ax1 = Axis(fig[1, 1][1, 1],title="nsi")
        #     ax2 = Axis(fig[1, 1][1, 2],title="αbatch")
        #     ax3 = Axis(fig[1, 1][2, 1],title="αproba")
        #     ax4 = Axis(fig[1, 1][2, 2],title="hinv")
        #     # @show op_name,batch_idx,maxαproba,minαproba
        #     barplot!(ax1,nsi,bins=100)
        #     hist!(ax2,αbatch,bins=100)
        #     barplot!(ax3,αproba)
        #     barplot!(ax4,hinv)
        #     display(fig)
        # end

        update_αproba!(αproba,αscore,nsi)
        update_distribution!(αproba,h,hinv)
        nsi.=0.0
    end #for batch_idx

    return face_mesh
end


function rayons_fast(fname::String,numerical_parameters::NumericalParameters,forward_name="")

    println("entering Rayons2D")
    outdir=create_outdir(fname,numerical_parameters)
    # test_trajectories(outdir,fname)

   
    problem =  init_problem_from_json(fname,numerical_parameters)
    celerity_domain,empty_face_mesh,cmg = problem.celerity_domain,problem.face_mesh,problem.cartesian_material_grid
    model = numerical_parameters.model

    source_points = problem.source_points
    observation_points = isfile(forward_name) ?  load_observation_points_from_forward(forward_name) : problem.observation_points
    @show isfile(forward_name)

    all_sources = vcat(source_points,observation_points)

    # p=barplot(rand(100),label="dummy")
    # display(p)
    # all_sources = (source_points[1],)

    face_meshes = Dict{String,FaceMesh}()
    for (i,op) ∈ enumerate(all_sources)
        op_name,op_pos = get_name_and_position_of_source(op)
        face_meshes[op_name] = deepcopy(problem.face_mesh)
    end
    # progress = Progress(length(all_sources), 0.1)
    tray = @elapsed begin
        @floop for op ∈ all_sources
        # for op ∈ all_sources
            op_name,op_pos = get_name_and_position_of_source(op)
            fm = face_meshes[op_name]
            ray_tracing!(fm,op_name,op_pos,problem,numerical_parameters)
            # next!(progress)
        end
    end

            

    # map!_args = [op,face_meshes[op_name]] 

    # tray = @elapsed begin
    #     fms=ThreadsX.map!(op->ray_tracing!(face_meshes[op[1]],op,problem,numerical_parameters),all_sources)  
    #     # fms=Base.map(op->ray_tracing(op,problem,numerical_parameters),all_sources)  
    # end
    # nr = Rayons.nrays(model,numerical_parameters.maxdepth)
    # @show tray,(length(observation_points)*numerical_parameters.nα*nr)/(tray*1.e6) 



    post_traitement = true
    if post_traitement
        
        export_data(outdir,fname,source_points,observation_points,face_meshes,cmg,numerical_parameters)

        for source_point in source_points
            source_name,source_center,radius = source_point
            observation_trajectories = get_ray_trajectories(model,face_meshes,source_point,observation_points,numerical_parameters,celerity_domain)
            # observation_trajectories = get_ray_trajectories(model,face_meshes,observation_points,numerical_parameters,celerity_domain,isx,isy)
            plot_all(outdir,source_name,source_center,face_meshes,observation_points,observation_trajectories,cmg,numerical_parameters.tmax)
            plot_animated_mintimes(outdir,source_name,source_center,face_meshes,cmg,numerical_parameters.tmax)
            # animate_mintimes(outdir,source_name,face_meshes[source_name],cmg)
        end

    end

    face_meshes,cmg
end



