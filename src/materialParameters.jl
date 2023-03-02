using PrettyTables
using JSON


load(fname) = read(fname, String) |> JSON.parse


struct MaterialParameters
    name::String
    E::Float64
    ν::Float64
    ρ::Float64
    cₗ::Float64
    cₜ::Float64
    λ::Float64
    μ::Float64
    K::Float64
end

# function Base.show(io::IO, ::MIME"text/plain", mp::MaterialParameters)
function Base.show(io::IO, mp::MaterialParameters)
        @show mp.name
    @show mp.E
    @show mp.ν
    @show mp.ρ
    @show mp.cₗ
    @show mp.cₜ
    @show mp.λ
    @show mp.μ
    @show mp.K
end

cl(μ,λ,ρ) = sqrt((2μ+λ)/ρ)
ct(μ,ρ) = sqrt(μ/ρ)
λ_ρcl(ρ,cₗ,μ)=ρ*(cₗ^2) - 2μ

μ_ctρ(ct,ρ) = ρ*ct^2
λ_clctρ(cl,ct,ρ) = ρ*cl^2 -2μ_ctρ(ct,ρ)

##https://fr.wikipedia.org/wiki/Coefficient_de_Lamé
#--------------------------  3D ------------------------------
function materialParametersFromEν3D(name,E,ν,ρ)
    λ = E*ν/((1+ν)*(1-2ν))
    μ = E/(2*(1+ν))
    K = λ+(2/3)μ
    cₗ= cl(μ,λ,ρ) 
    cₜ= ct(μ,ρ)
    return MaterialParameters(name,E,ν,ρ,cₗ,cₜ,λ,μ,K)
end

function materialParametersFromλμ3D(name,λ,μ,ρ)
    # K = λ+(2/3)μ
    K = λ+(2/3)μ
    ν = λ/(2(λ+μ))
    E = μ*(3λ+2μ)/(λ+μ)
    cₗ= cl(μ,λ,ρ) 
    cₜ= ct(μ,ρ)
    return MaterialParameters(name,E,ν,ρ,cₗ,cₜ,λ,μ,K)
end

function materialParametersFromclctρ3D(name,cl,ct,ρ)
    λ = λ_clctρ(cl,ct,ρ)
    μ = μ_ctρ(ct,ρ)
    K = λ+(2/3)μ
    ν = λ/(2(λ+μ))
    E = μ*(3λ+2μ)/(λ+μ)
    return MaterialParameters(name,E,ν,ρ,cl,ct,λ,μ,K)
end

function liquidParametersFromρμc3D(name,ρ,μ,cₗ)
    λ = λ_ρcl(ρ,cₗ,μ)
    return materialParametersFromλμ3D(name,λ,μ,ρ)
end

#--------------------------  2D ------------------------------

function materialParametersFromEν2D(name,E,ν,ρ)
    λ = E*ν/((1+ν)*(1-ν))
    μ = E/(2*(1+ν))
    K = λ + μ
    cₗ= cl(μ,λ,ρ)
    cₜ= ct(μ,ρ)
    return MaterialParameters(name,E,ν,ρ,cₗ,cₜ,λ,μ,K)
end

function materialParametersFromλμ2D(name,λ,μ,ρ)
    K = λ + μ
    ν = λ/(λ+2μ)
    E = 4μ*(λ+μ)/(λ+2μ)
    cₗ= cl(μ,λ,ρ)
    cₜ= ct(μ,ρ)
    return MaterialParameters(name,E,ν,ρ,cₗ,cₜ,λ,μ,K)
end

function liquidParametersFromρμc2D(name,ρ,μ,cₗ)
    λ = λ_ρcl(ρ,cₗ,μ)
    return materialParametersFromλμ2D(name,λ,μ,ρ)
end

function materialParametersFromclctρ2D(name,cl,ct,ρ)
    λ = λ_clctρ(cl,ct,ρ)
    μ = μ_ctρ(ct,ρ)
    K = λ+μ
    ν = λ/(2(λ+μ))
    E = μ*(3λ+2μ)/(λ+μ)
    return MaterialParameters(name,E,ν,ρ,cl,ct,λ,μ,K)
end

function dump_material(mps::Vector{MaterialParameters})
    fns = fieldnames(MaterialParameters)
    vc = [getproperty.(mps,field) for field ∈ fns]
    hvc=hcat(vc...)
    headers=[string(field) for field ∈ fns]
    pretty_table(hvc,headers,formatters=ft_printf("%5.2e"))
end

function dump_material_tex(mps::Vector{MaterialParameters})
    fns = fieldnames(MaterialParameters)
    vc = [getproperty.(mps,field) for field ∈ fns]
    hvc=hcat(vc...)
    headers=[string(field) for field ∈ fns]
    pretty_table(hvc,backend = Val(:latex),headers,formatters=ft_printf("%5.2e"))
end

function read_material_parameters(jdata)
    jmaterial = jdata["material_properties"]
    mps=Vector{MaterialParameters}()
    for mp ∈ jmaterial
        name=mp["name"]
        cl=mp["longitudinal_velocity"]
        ct=mp["transversal_velocity"]
        ρ=mp["density"]
        push!(mps,materialParametersFromclctρ2D(name,cl,ct,ρ))
    end
    mps
end


function getTestMaterialParameters()
    materialParams=[
    #
    # materialParametersFromEν3D("acier",210.e9,0.27,7850), #ν=0.27 => cₗ = 5781.7
    # materialParametersFromEν3D("carton",11.e9,0.3,600),
    # materialParametersFromEν3D("cuivre",128.e9,0.33,8920),
    # materialParametersFromEν3D("fer",208.e9,0.21,7860),
    # liquidParametersFromρμc2D("huile",900.,0.001,1000)

    #Andreas parameters
    #materialParametersFromctclρ3D(name,ct,cl,ρ)
    materialParametersFromclctρ2D("acier",5900,3230,7700), #ν=0.27 => cₗ = 5781.7
    materialParametersFromclctρ2D("carton",2300,1100,600),
    materialParametersFromclctρ2D("cuivre",3600,2100,8920),
    # materialParametersFromclctρ2D("cuivre",10000,2100,8920),
    materialParametersFromclctρ2D("fer",5900,3230,7700),
    materialParametersFromclctρ2D("huile",1400,0,1000)
    # liquidParametersFromρμc2D("huile",1000.,0.001,1400)


    # materialParametersFromEν2D("acier",210.e9,0.27,7850), #ν=0.27 => cₗ = 5781.7
    # materialParametersFromEν2D("carton",11.e9,0.3,600),
    # materialParametersFromEν2D("cuivre",128.e9,0.33,8920),
    # materialParametersFromEν2D("fer",208.e9,0.21,7860),
    # liquidParametersFromρμc2D("huile",900.,0.001,1000)
    ]

    dump_material(materialParams)

    return materialParams
end

function getFluidifiedTestMaterialParameters()
    ct_fluid = 0.1
    materialParams=[
        materialParametersFromclctρ3D("acier",5900,3230,7700),
        materialParametersFromclctρ3D("carton",2300,ct_fluid,600),
        materialParametersFromclctρ3D("cuivre",3600,ct_fluid,8920),
        materialParametersFromclctρ3D("fer",5900,ct_fluid,7700),
        materialParametersFromclctρ3D("huile",1400,ct_fluid,1000)
    ]
    dump_material(materialParams)

    return materialParams
end

function getHomogeneousFluidTestMaterialParameters()
    ct_fluid = 0.1
    cl_fluid = 1400.0
    rho_fluid = 2000.0
    materialParams=[
        materialParametersFromclctρ3D("acier",5900,3230,7700),
        materialParametersFromclctρ3D("carton",cl_fluid,ct_fluid,rho_fluid),
        materialParametersFromclctρ3D("cuivre",cl_fluid,ct_fluid,rho_fluid),
        materialParametersFromclctρ3D("fer",cl_fluid,ct_fluid,rho_fluid),
        materialParametersFromclctρ3D("huile",cl_fluid,ct_fluid,rho_fluid)
    ]
    dump_material(materialParams)

    return materialParams
end


