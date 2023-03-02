using DataStructures

function distance(pa,pb)
    (xa,ya) = pa
    (xb,yb) = pb
    # return 1
    sqrt((xa-xb)^2+(ya-yb)^2)
end

function initialize_source_cell(fm::FaceMesh,is,js)
    xc,yc=get_cell_center(fm,is,js)
    min_time(x,y) = distance((xc,yc),(x,y))/fm.speed_l[is,js]
    for shift ∈ (0,1)
        set_vface_values!(fm,is+shift,js,min_time)
        set_hface_values!(fm,is,js+shift,min_time)
    end
end

function face_to_face_scatter!(fm::FaceMesh,fin::FaceId,fout::FaceId,speed)

    xmin_i,ymin_i = fm.xfaces[fin.i],fm.yfaces[fin.j]
    xmin_j,ymin_j = fm.xfaces[fout.i],fm.yfaces[fout.j]
    lf = fm.lfaces
    m = length(lf)
    cxi = fin.vertical ? 0.0 : 1.0
    cyi = fin.vertical ? 1.0 : 0.0
    cxj = fout.vertical ? 0.0 : 1.0
    cyj = fout.vertical ? 1.0 : 0.0

    face_i  = get_face(fm,fin)
    face_j  = get_face(fm,fout)

    mpj=face_j.minpaths
    mpi=face_i.minpaths

    speed⁻¹=inv(speed)

    minfacetime=Inf

    @inbounds for lj=1:m
        dxj=lf[lj]
        xj=xmin_j+cxj*dxj
        yj=ymin_j+cyj*dxj
        @simd for li ∈ 1:m
            tmin_i = mpi[li].time
            dxi=lf[li]
            xi=xmin_i+cxi*dxi
            yi=ymin_i+cyi*dxi

            new_tmin_j = distance((xi,yi),(xj,yj))speed⁻¹ + tmin_i
            # new_tmin_j = sqrt((xi-xj)^2+(yi-yj)^2)/speed + tmin_i
            mpj[lj].time > new_tmin_j && (mpj[lj]=MinPath(new_tmin_j,Predecessor(fin,li)))
            minfacetime=min(minfacetime,mpj[lj].time)
        end
    end
    minfacetime
end

function face_to_face_scatter_old!(fm::FaceMesh,fin::FaceId,fout::FaceId,speed)
    face_i  = get_face(fm,fin)
    face_j  = get_face(fm,fout)
    for (li,pti) ∈ enumerate(pts_i)
        tmin_i = face_i.minpaths[li].time
        for (lj,ptj) ∈ enumerate(pts_j)
            tmin_j = face_j.minpaths[lj].time
            new_tmin_j = distance(pti,ptj)*speed + tmin_i
            (tmin_j>new_tmin_j) && (face_j.minpaths[lj]=MinPath(new_tmin_j,Predecessor(fin,li)))
        end
    end
end


function cell_face_scatter!(fm::FaceMesh,i,j)
    faces_indexes = ( FaceId(i,j,true), FaceId(i+1,j,true), FaceId(i,j,false), FaceId(i,j+1,false) )
    speed=fm.speed_l[i,j]
    mincelltime=Inf
    for fin ∈ faces_indexes
        # @show fin
        for fout ∈ faces_indexes
            # @show fout
            if fout != fin
                minfacetime=face_to_face_scatter!(fm,fin,fout,speed)
                mincelltime=min(mincelltime,minfacetime)
            end
        end
    end
    mincelltime
end

function cell_scatter!(fm::FaceMesh,i,j)
    (nxp1,ny) = size(fm.vfaces) 
    nx = nxp1-1
    mincelltime=Inf
    for (is,js) ∈ ((i-1,j),(i+1,j),(i,j-1),(i,j+1))
        (is<1 || is>nx || js<1 || js >ny) && continue  
        # @show is,js
        mincelltime=min(mincelltime,cell_face_scatter!(fm,is,js))
    end
    mincelltime
end


struct CellTask
    i::Int
    j::Int
    mintime::Float64
end

Base.isless(a::CellTask,b::CellTask) = a.mintime < b.mintime

function string_ss(ss) 
    if length(ss)==0
        return ""
    else
        return mapreduce(c->"("*string(c.i)*","*string(c.j)*") ",*,ss) 
    end
end

function insert_ss!(celltasks,todocells,ct)
    ij=(ct.i,ct.j)
    insert!(todocells,ij)
    insert!(celltasks,ct)
end




function add_cell!(celltasks,todocells,donecells,i,j,nx,ny,speed,mincelltime,h)
    if ((i>0) && (i<=nx) && (j>0) && (j<=ny)) 
        if (!haskey(donecells,(i,j)) && !haskey(todocells,(i,j)))
            # @show "insert",i,j,donecells
            # ct=CellTask(i,j,speed[i,j])
            ct=CellTask(i,j,mincelltime+h/speed[i,j])
            # @show string_ss(celltasks)
            insert_ss!(celltasks,todocells,ct)
            # @show string_ss(celltasks)
        end
    end
end



function add_neighbours(celltasks,todocells,donecells,i,j,nx,ny,speed,mincelltime,h)
    ij=(i,j)
    insert!(donecells,ij)
    # @show donecells
    for (inext,jnext) ∈ ((i-1,j),(i+1,j),(i,j-1),(i,j+1))
        add_cell!(celltasks,todocells,donecells,inext,jnext,nx,ny,speed,mincelltime,h)
    end
end

function robust_cell_scatters(fm::FaceMesh,ic,jc)

    (nx,ny) = get_cell_sizes(fm)
    @show "ncells=",nx,ny,nx*ny
    speed = fm.speed_l

    donecells=SortedSet{Tuple{Int,Int}}()
    todocells=SortedSet{Tuple{Int,Int}}()
    celltasks=SortedSet{CellTask}()

    h=fm.xfaces[2]-fm.xfaces[1]
    add_neighbours(celltasks,todocells,donecells,ic,jc,nx,ny,speed,0,h)

    treated=1
    while (!isempty(celltasks))
        # @show length(celltasks)
        # @show celltasks
        cell=pop!(celltasks)
        treated+=1
        (i,j)=cell.i,cell.j
        # @show i,j
        # mincelltime=cell_face_scatter!(fm,i,j)
        mincelltime=cell_scatter!(fm,i,j)
        add_neighbours(celltasks,todocells,donecells,i,j,nx,ny,speed,mincelltime,h)
    end
    @show treated
end

function range_cells_scatter!(fm,ir,jr)
    (nxp1,ny) = size(fm.vfaces) 
    nx = nxp1-1
    for i ∈ ir
        (i<1 || i>nx) && continue 
        for j ∈ jr
            (j<1 || j>ny) && continue
            cell_scatter!(fm,i,j)
        end
    end
end

function all_cells_scatter!(fm::FaceMesh,ic,jc)
    (nxp1,ny) = size(fm.vfaces) 
    nx = nxp1-1
    for d ∈ 0:max(nx,ny)

        imin = ic-d
        imax = ic+d
        jmin = jc-d
        jmax = jc+d

        range_cells_scatter!(fm,imin:imax,jmin:jmin)
        range_cells_scatter!(fm,imin:imax,jmax:jmax)
        range_cells_scatter!(fm,imin:imin,jmin:jmax)
        range_cells_scatter!(fm,imax:imax,jmin:jmax)

    end
end








    


