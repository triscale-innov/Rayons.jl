using Images

struct Rect
    xmin::Int
    ymin::Int
    xmax::Int
    ymax::Int
end

# Create a json input file (rectangles, materials, source, sensors) from a maze png image
function maze_rectangles_from_png(filename)
    img = Images.load(filename)
    nx,ny=size(img)
    f(i,j) = img[i,j] == RGBA(1.0,1.0,1.0,1.0) ? RGBA(1.0,1.0,1.0,1.0) : RGBA(0.0,0.0,0.0,0.0)
    f2(i,j) = img[i,j] == RGBA(1.0,1.0,1.0,1.0) 

    bool_matrix = [f2(i,j) for i ∈ 1:nx , j ∈ 1: ny]

    bool_matrix = transpose(bool_matrix)
    nx,ny = ny,nx

    in_rect(x,y,r::Rect) = (r.xmin <= x <= r.xmax) && (r.ymin <= y <= r.ymax)

    rectangles = Vector{Rect}()

    @. bool_matrix = !bool_matrix 

    for i ∈ 1:nx
        for j ∈ 1:ny
            !bool_matrix[i,j] && continue
            in_rectangles(i,j) = any(r->in_rect(i,j,r), rectangles)

            in_rectangles(i,j) && continue
            
            imin = i
            jmin = j

            imax = imin
            while imax <nx && bool_matrix[imax,jmin] 
                imax +=1
            end
            imax <= nx && (imax=imax-1)
            
            line_ok(j) = all(i-> bool_matrix[i,j],imin:imax)
            jj = jmin
            while (jj <= ny) && line_ok(jj) 
                jj +=1
            end

            jmax = min(jj,ny)
            push!(rectangles,Rect(imin-1,jmin-1,imax+1,jmax+1))
        end
    end
    rectangles

    bw_in_rectangles(i,j) = any(r->in_rect(i,j,r), rectangles) ?  RGBA(1.0,1.0,1.0,1.0) : RGBA(0.0,0.0,0.0,0.0)
    im3 = [bw_in_rectangles(i,j) for i ∈ 1:nx , j ∈ 1:ny]

    save_to_file("out","bw.png",img)

    water = Dict(
        "name"=>"Water", 
        "density"=>1000,
        "longitudinal_velocity"=> 1.0,
        "transversal_velocity"=> 1.0
    )

    steel = Dict(
        "name"=>"Steel", 
        "density"=>1000,
        "longitudinal_velocity"=> 3100,
        "transversal_velocity"=> 3100
    )

   rect_dict(imin,jmin,imax,jmax) = Dict("sw"=>[imin,jmin],"ne"=>[imax,jmax])

   a=Dict(
    "material_properties"=>[steel,water],
    "name"=>"maze",
    "bounding_box"=>rect_dict(0,0,nx+1,ny+1),
    "sensors"=>[["sortie",nx-100,ny-1]],
    "sources"=>[[100,10]])



    rds = [Dict("name"=>"r",
        "material"=>"Water",
        "sw"=>[r.xmin,r.ymin],"ne"=>[r.xmax,r.ymax]) for r ∈ rectangles]

    bkg = [Dict("name"=>"background",
    "material"=>"Steel",
    "sw"=>[1,1],"ne"=>[nx,ny])]

    a["rects"]=vcat(rds,bkg)

   open("maze.json","w") do f 
    write(f, JSON.json(a))
   end

end
