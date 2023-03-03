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





