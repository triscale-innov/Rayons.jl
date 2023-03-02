function get_material_closures(model,materialParams)

    labels = get_face_labeling(model)
    label(model, name) = Gridap.Geometry.get_tag_from_name(labels, name)
    tags = get_face_tag(labels, #=dimension=# 2)

    function get_mps(model,tag)
        for mps in materialParams
            if (tag == label(model,mps.name))
                return mps
            end
        end
        @show tag
        error("unknown material")
        nothing
    end
    get_material_properties(s::Symbol) = getproperty.(get_mps.(Ref(model), tags),s)

    velocities = get_material_properties(:cₗ)
    cₛ²(x, v) = x * v^2
    cₛ²(x) = cₛ²∘(x, velocities)

    cₛ⁻²(x, v) = x /( v^2 )
    cₛ⁻²(x) = cₛ⁻²∘(x, velocities)

    densities = get_material_properties(:ρ)
    ρₛ(x, density) = x * density
    ρₛ(x) = ρₛ∘(x, densities)


    ρₛ⁻¹(x, density) = x / density
    ρₛ⁻¹(x) = ρₛ⁻¹∘(x, densities)


    lambdas = get_material_properties(:λ)
    λₛ(x,lambda) = x * lambda
    λₛ(x) = λₛ∘(x, lambdas)

    mus = get_material_properties(:μ)
    μₛ(x,mu) = x * mu
    μₛ(x) = μₛ∘(x, mus)

    bulk_moduli = get_material_properties(:K)
    Kₛ(x,bm) = x * bm
    Kₛ(x) = Kₛ∘(x, bulk_moduli)
    pression(u) = -Kₛ(divergence(u))
    # pression(u) = -λₛ(tr(ε(u)))
    # pression(u) = -(divergence(u))

    function get_λμ(model, tag)
        for mps ∈ materialParams
            if (tag == label(model,mps.name))
                return (mps.λ,mps.μ)
            end
        end
        @show tag
        error("unknown material")
        nothing
    end

    function σₛ(ε, tag)
        (λ, μ) = get_λμ(Ref(model),tag)
        λ*tr(ε)*one(ε) + 2*μ*ε
        # λₚ=2μ*λ/(λ+2μ)
        # λₚ*tr(ε)*one(ε) + 2*μ*ε
    end
    σₛ(u) = σₛ∘(ε(u), tags)

    trσₛ(u) = tr(σₛ(u))

    return cₛ²,cₛ⁻²,ρₛ,ρₛ⁻¹,λₛ,μₛ,Kₛ,pression,σₛ,trσₛ
end

 