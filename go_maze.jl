using Rayons
using Random
using CairoMakie

function go_maze()

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
        fname = "maze.json"
    fms,cmg=Rayons.rayons_fast(fname,numerical_parameters); # 
end

fms,cmg=go_maze();
nothing


