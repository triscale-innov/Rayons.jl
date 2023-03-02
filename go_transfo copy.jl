using Rayons
using Random
using CairoMakie



function go()

    Rayons.rectangles_from_png("maze.png")

    fontsize_theme = Theme(fontsize = 30)
    set_theme!(fontsize_theme)

    Random.seed!(1234)

    numerical_parameters = Rayons.NumericalParameters(
        ds=2.e-3, # unused for rectangular based geometry
        ϵ=1.e-15, # unused for rectangular based geometry
        nmin=50, # cartesian grid size (for min time accumulation)
        maxdepth=200, # max trajectory tree depth : probably unused (ray stop because of power decay or max time reached)
        nα=20000, # number of initial angles for rays 
        tmax=4, # max time allowed for rays
        power_threshold=1.e-3, #min power allowed for rays
        history_size=1, # unused
        model=Rayons.Fluid) # Fluid or Solid
    
        # fname=joinpath(pwd(),"data","US_GeoGmsh_2DZX_1E05Hz.json")
        # fname = "test1.json"
        fname = "foo.json"
        # fname="US_ShellType_NoPB_2DXY_5E04Hz_Gmsh.json"
        # fname="US_ShellType_NoPB_2DXY_5E04Hz_Gmsh.json"
        # fname="US_ShellType_NoPB_2DXY_5E04Hz_Gmsh.json"
        # fname=joinpath(pwd(),"data","US_GeoGmsh_2DXY_1E05Hz.json") # Json describing the problem
   
    fms,cmg=Rayons.rayons_fast(fname,numerical_parameters); # 
end

fms,cmg=go();
nothing


