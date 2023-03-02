Base.@kwdef struct NumericalParameters
    ds::Float64
    ϵ::Float64
    nmin::Int #min cell number
    maxdepth::Int
    nα::Int
    tmax::Float64 = 1.0
    power_threshold::Float64 = 0.0
    history_size::Int = 100
    model::Type
end

function create_out_dir_from_np(np::NumericalParameters)
    snh="nh_"*unit_to_string(np.nmin,"")
    smd="md_"*unit_to_string(np.maxdepth,"")
    sr="nr_"*unit_to_string(np.nα,"")
    stm="tm_"*unit_to_string(np.tmax,"")
    spt="pt_"*unit_to_string(np.power_threshold,"")
    return snh*smd*sr*stm*spt
end

unpack(np::NumericalParameters) = np.ds, np.ϵ