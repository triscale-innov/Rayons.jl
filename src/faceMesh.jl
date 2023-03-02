struct FaceId
    i::Int
    j::Int
    vertical::Bool
end

FaceId() = FaceId(-1,-1,true)

struct Predecessor
    faceid::FaceId
    index::Int
end

Predecessor() = Predecessor(FaceId(),-1)

struct MinPath
    time::Float64
    pred::Predecessor
end


MinPath() = MinPath(Inf,Predecessor()) 

struct Faces
    minpaths::Vector{MinPath} # [m] 
    Faces(m) = new([MinPath() for k ∈ 1:m])
end

min_times(f::Faces) = ( mp.time for mp ∈ f.minpaths)

struct MinRay
    time::Float64
    p::Point
    u::Point
    idx::Int128
    αorigin::Float64
    porigin::Point
end


MinRay() =  MinRay(Inf,Point(0,0),Point(0,0),-1,0.0,Point(0,0))

# fine_to_coarse_index(i,c) = (i+c-1) ÷ c
fine_to_coarse_index(i,c,cp) = (i+c-1) >> cp

struct FaceMesh
    vfaces::Matrix{Faces} #[nx+1,ny]
    hfaces::Matrix{Faces} #[nx,ny+1]
    h::Float64
    hh::Float64
    xfaces::Vector{Float64} #[nx+1]
    yfaces::Vector{Float64} #[ny+1]
    lfaces::Vector{Float64} #[m]
    speed_l::Matrix{Float64} #[nx,ny]
    ray_min_times::Matrix{MinRay} #[nx,ny]
    coarse_ray_min_times::Matrix{Float64} #[mx,my]
    coarse_two_power::Int
    coarse_factor::Int
    coarse_counter::Matrix{Int} #[mx,my]
    coarse_ijmax::Matrix{Tuple{Int,Int}} #[mx,my]
    # coarse_jmaxtime::Matrix{Int} #[mx,my]
    tmax::Float64
    history::Matrix{Vector{Float64}} #[nx,ny]
    αscore::Vector{Float64}

    function FaceMesh(xmin,ymin,nx,ny,h,m,celerity_domain::AbstractCelerityDomain,tmax,history_size)
        xfaces = [xmin+i*h for i ∈ 0:nx]
        yfaces = [ymin+j*h for j ∈ 0:ny]
        hh = h/m
        lfaces = [hh/2 + k*hh for k ∈ 0:(m-1)]
        vfaces = [Faces(m) for i ∈ 0:nx, j ∈ 1:ny]
        hfaces = [Faces(m) for i ∈ 1:nx, j ∈ 0:ny]
      
        state = initial_state(celerity_domain)
        speed_l=zeros(nx,ny)
        for i ∈ 1:nx , j ∈ 1:ny
            x=(xfaces[i]+xfaces[i+1])/2
            y=(yfaces[j]+yfaces[j+1])/2
            state =  get_state(Point(x,y),state,celerity_domain)
            c,_,_ = cl_ct_rho(state)
            # c,cs=celerity(Point(x,y),cs)
            speed_l[i,j]=c
        end

        ray_min_times = [MinRay() for i ∈ 1:nx , j ∈ 1:ny]
        history = [zeros(history_size) for i ∈ 1:nx , j ∈ 1:ny]
        αscore_size = 100
        αscore=ones(αscore_size)
        coarse_two_power = 5
        coarse_factor=2^coarse_two_power
        (mx,my) = fine_to_coarse_index.((nx,ny),coarse_factor,coarse_two_power)
        coarse_ray_min_times = [Inf for I ∈ 1:mx, J ∈ 1:my]
        coarse_counter = [0 for I ∈ 1:mx, J ∈ 1:my]
        coarse_ijmax =  [(coarse_factor*(I-1)+1,coarse_factor*(J-1)+1) for I ∈ 1:mx, J ∈ 1:my]

        new(vfaces,hfaces,h,hh,xfaces,yfaces,lfaces,speed_l,ray_min_times,
            coarse_ray_min_times,coarse_two_power,coarse_factor,coarse_counter,coarse_ijmax,
            tmax,history,αscore)
    end
end


function update_coarse_ray_min_time!(fm::FaceMesh,I,J,inew,jnew)
    (nx, ny) = get_cell_sizes(fm)

    @inbounds (imax,jmax) = fm.coarse_ijmax[I,J]
    if (inew==imax && jnew==jmax) 
        max_coarse_time = typemin(Float64)
      
        c = fm.coarse_factor
        ib,ie = c*(I-1)+1,min(c*I,nx)
        jb,je = c*(J-1)+1,min(c*J,ny)
        for i ∈ ib:ie
            for j ∈ jb:je
                @inbounds tij=fm.ray_min_times[i, j].time
                if max_coarse_time<tij
                    max_coarse_time=tij
                    imax,jmax=i,j
                end
            end
        end
        @inbounds fm.coarse_ray_min_times[I,J]=max_coarse_time
        @inbounds fm.coarse_ijmax[I,J]=(imax,jmax)
    end
end

function sum_min_times(fm::FaceMesh)
    rmt = fm.ray_min_times
    s=0.0
    m=0.0
    for mr ∈ rmt
        s+=mr.time
        m = max(m,mr.time) 
        s==Inf && break
    end
    s,m
end

function get_cell_sizes(fm::FaceMesh)
    (nxp,ny) = size(fm.vfaces)
    nx=nxp-1
    (nx,ny)
end

function get_cell_index(xfaces,x)
    xmin,xmax = first(xfaces),last(xfaces)
    if ( (x>=xmin) && x<= last(xmax) )
        for i ∈ 2:length(xfaces)
            xfaces[i]> x && return i-1
        end
    else
        error("coordinate out of bound : x=$x ∉ [$xmin, $xmax]")
    end
    return 0
end

function get_cell_indexes(fm::FaceMesh,x,y)
    i=get_cell_index(fm.xfaces,x)
    j=get_cell_index(fm.yfaces,y)
    i,j
end

function get_hface_points(fm::FaceMesh,i,j)
    xmin=fm.xfaces[i]
    ymin=fm.yfaces[j]
    ((xmin+lf,ymin) for lf ∈ fm.lfaces)
end


function get_vface_points(fm::FaceMesh,i,j)
    xmin=fm.xfaces[i]
    ymin=fm.yfaces[j]
    ((xmin,ymin+lf) for lf ∈ fm.lfaces)
end


function get_cell_center(fm::FaceMesh,is,js)
    xc=(fm.xfaces[is]+fm.xfaces[is+1])/2
    yc=(fm.yfaces[js]+fm.xfaces[js+1])/2
    xc,yc
end
    
function set_face_values!(face,facepoints,point_function)
    for (l,(x,y)) ∈ enumerate(facepoints)
        face.minpaths[l]=MinPath(point_function(x,y),face.minpaths[l].pred)
    end
end

function set_vface_values!(fm::FaceMesh,iface,jface,point_function)
    face=fm.vfaces[iface,jface]
    set_face_values!(face,get_vface_points(fm,iface,jface),point_function)
end
function set_hface_values!(fm::FaceMesh,iface,jface,point_function)
    face=fm.hfaces[iface,jface]
    set_face_values!(face,get_hface_points(fm,iface,jface),point_function)
end

function get_face(fm::FaceMesh,faceid::FaceId)
    faceid.vertical ? fm.vfaces[faceid.i,faceid.j] : fm.hfaces[faceid.i,faceid.j]
end

function get_face_point(fm::FaceMesh,faceid::FaceId)
    faceid.vertical ? get_vface_points(fm,faceid.i,faceid.j) : get_hface_points(fm,faceid.i,faceid.j)
end




   


















    


