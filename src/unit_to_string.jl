using Printf

function unit_to_string(val::Integer,unit) 
    sval=string(val)
    sval=replace(sval,"."=>"p")
    return sval*unit
end

function unit_to_string(val::Float64,unit) 
    sval=@sprintf("%.2E",val)
    sval=replace(sval,"."=>"p")
    sval=replace(sval,"p00" => "")
    return sval*unit
end

function unit_from_string(sval,unit)
    @show sval,unit
    sval=replace(sval,unit=>"")
    sval=replace(sval,"p"=>".")
    return parse(Float64,sval)
end

