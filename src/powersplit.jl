
import Base.@kwdef
using StaticArrays


function build_rt_a_matrix_solid_solid(
    uθd1, uθd2, uθs1, uθs2, 
    ρ1,ρ2,cd1,cs1,cd2,cs2)

    zd1,zs1 = ρ1*cd1,ρ1*cs1
    zd2,zs2 = ρ2*cd2,ρ2*cs2

    cθd1,sθd1 = coordinates(uθd1)
    cθs1,sθs1 = coordinates(uθs1)
    cθd2,sθd2 = coordinates(uθd2)
    cθs2,sθs2 = coordinates(uθs2)

    c2θd1,s2θd1 = coordinates(double_u(uθd1))
    c2θs1,s2θs1 = coordinates(double_u(uθs1))
    c2θd2,s2θd2 = coordinates(double_u(uθd2))
    c2θs2,s2θs2 = coordinates(double_u(uθs2))

    @SMatrix [
        -cθd1  -cθd2 -sθs1 +sθs2 ;
        -sθd1  +sθd2 +cθs1 +cθs2 ;
        -zd1*c2θs1  +zd2*c2θs2 -zs1*s2θs1 -zs2*s2θs2 ;
        -zs1*cs1*s2θd1/cd1  -zs2*cs2*s2θd2/cd2 zs1*c2θs1 -zs2*c2θs2 ;
    ]

    # @SMatrix [
    #     -cθd1  -cθd2 -sθs1 +sθs2 ;
    #     -sθd1  +sθd2 +cθs1 +cθs2 ;
    #     -zd1*c2θd1  +zd2*c2θd2 -zs1*s2θs1 -zs2*s2θs2 ;
    #     -zd1*cs1*s2θd1/cd1  -zd2*cs2*s2θd2/cd2 zs1*c2θs1 -zs2*c2θs2 ;
    # ]
end

function  build_rt_b_vector_solid_solid(uθdi,uθsi,ρ1,cd1,cs1)

    cθdi,sθdi = coordinates(uθdi)
    c2θdi,s2θdi = coordinates(double_u(uθdi))
    c2θsi,s2θsi = coordinates(double_u(uθsi))
    zd1,zs1 = ρ1*cd1,ρ1*cs1

    @SVector [ -cθdi , -sθdi , -zd1*c2θsi , -zs1*cs1*s2θdi/cd1 ]
    # @SVector [ -cθdi , -sθdi , -zs1*c2θdi , -zs1*cs1*s2θdi/cd1 ]
end

function  build_rt_c_vector_solid_solid(uθsi,ρ1,cd1,cs1)
    cθsi,sθsi = uθsi.x,uθsi.y
    u2θsi = upv(uθsi,uθsi) 
    c2θsi,s2θsi = u2θsi.x,u2θsi.y
    zd1,zs1 = ρ1*cd1,ρ1*cs1

    @SVector [ sθsi , cθsi , -zs1*s2θsi , -zs1*c2θsi ]
end

complex_cos_from_sin(s) = (s < 1) ? sqrt(1 - s^2) : im*sqrt(s^2-1)


function solid_solid(iang, d1, d2, cp1, cs1, cp2, cs2)
    # % [ tp, ts, rp, rs] = solid_solid(iangd, d1, d2, cp1, cs1, cp2, cs2, type)
    # % returns the transmitted P-wave (tp), transmitted SV-wave (ts), reflected
    # % P-wave (rp), and reflected SV-wave (rs) transmission/reflection 
    # % coefficients (based on velocity ratios) for two solids in welded contact.
    # % The inputs are the 
    # incident angle(s),iangd,(in degrees), 
    # (d1, d2), the  densities of the two media, 
    # (cp1, cs1), the compressional and shear wave speeds of the first medium, 
    # (cp2, cs2) the compressional and shear wave speeds of the second medium. 
    # type is a string ('P' or 'S') indicating the type of incident wave in medium one. 
    # If cs1 =0 and type = 'P' the function returns the coefficients for a fluid-solid interface with rs = 0. 
    # The wave speed cs2 cannot be set equal to zero.
    

    sinp1 = sin(iang)
    sins1 = (cs1/cp1)*sin(iang)
    sinp2 = (cp2/cp1)*sin(iang)
    sins2 = (cs2/cp1)*sin(iang)

    # take into account cosines may be imaginary beyond
    # critical angles 
    (cosp1,coss1,cosp2,coss2) = complex_cos_from_sin.((sins1,sins1,sinp2,sins2))

    # %double angle functions
    sin2s1 = 2*sins1*coss1
    sin2s2 = 2*sins2*coss2
    sin2p1 = 2*sinp1*cosp1
    sin2p2 = 2*sinp2*cosp2
    
    cos2s1 = 1 - 2*sins1^2
    cos2p1 = 1 - 2*sinp1^2
    cos2s2 = 1 - 2*sins2^2
    cos2p2 = 1 - 2*sinp2^2
    
    # % form up terms
    den1 = (cs1/cp1)*sin2s1*sinp1 +cos2s1*coss1
    den2 = (cs1/cp1)*sin2p1*sins1 + cos2s1*cosp1
    l1 = (cs1/cp2)*sin2s1*sinp2 +(d2/d1)*cos2s2*coss1
    m1 = (cs1/cs2)*sin2s1*coss2 - (d2/d1)*sin2s2*coss1
    l2 = (cs1/cp2)*cos2s1*sinp2 -(d2/d1)*cos2s2*sins1
    m2 = (cs1/cs2)*cos2s1*coss2 +(d2/d1)*sin2s2*sins1
    l3 = (cp1/cp2)*cos2s1*cosp2 +(d2*cs2^2/(d1*cp2^2))*sin2p2*sinp1
    m3 = -(cp1/cs2)*cos2s1*sins2 +(d2/d1)*cos2s2*sinp1
    l4 = -(cs1^2/(cp1*cp2))*sin2p1*cosp2 +(d2*cs2^2/(d1*cp2^2))*sin2p2*cosp1
    m4 = (cs1^2/(cp1*cs2))*sin2p1*sins2 +(d2/d1)*cos2s2*cosp1
    
    den = (l2/den1 + l4/den2)*(m1/den1 + m3/den2) - (l1/den1 +l3/den2)*(m2/den1 +m4/den2)
    
    # %calculate potential ratios
    atai = -2*(m2/den1 + m4/den2)/den;
    btai = 2*(l2/den1 + l4/den2)/den;
    arai = ((l2/den1 + l4/den2)*(m1/den1 - m3/den2) - (l1/den1 - l3/den2)*(m2/den1 + m4/den2))/den
    brai = (2*(l2/den1)*(m4/den2) - 2*(m2/den1)*(l4/den2))/den
    atbi = 2*(m1/den1 + m3/den2)/den
    btbi = -2*(l1/den1 + l3/den2)/den
    brbi = ((l4/den2 - l2/den1)*(m1/den1 + m3/den2) - (l1/den1 + l3/den2)*(m4/den2 - m2/den1))/den
    arbi = (2*(l1/den1)*(m3/den2) - 2*(m1/den1)*(l3/den2))/den
    
    # %calculate velocity ratios
    tp = (cp1/cp2)*atai;
    ts = (cp1/cs2)*btai;
    rp= arai
    rs = cs1 > 0 ? (cp1/cs1)*brai : 0.

    tp, ts, rp, rs
end



function test()

    cd1,cd2=rand(),rand()
    cs1,cs2=rand()*cd1,rand()*cd2
    ρ1,ρ2=rand(),rand()

    # θsi = π/2 *rand() 
    θdi = π/2 *rand()
    
    # uθsi = Point(cos(θsi),sin(θsi))
    
    # longitudinal incident ray 
    uθdi = Point(cos(θdi),sin(θdi))
    muθd1 = maybeuθr(uθdi, cd1, cd1) # longitudinal  reflection angle
    muθs1 = maybeuθr(uθdi, cd1, cs1) # transversal   reflection angle
    muθd2 = maybeuθr(uθdi, cd1, cd2) # longitudinal  refraction angle
    muθs2 = maybeuθr(uθdi, cd1, cs2) # transversal   refraction angle

    allvalid = !isnothing(muθd1) && !isnothing(muθs1) && !isnothing(muθd2) && !isnothing(muθs2)

    if allvalid
        uθd1::Point = muθd1
        uθs1::Point = muθs1
        uθd2::Point = muθd2
        uθs2::Point = muθs2

        uθsi::Point = muθs1

        a=build_rt_a_matrix_solid_solid(uθd1,uθd2,uθs1,uθs2,ρ1,ρ2,cd1,cs1,cd2,cs2)

        b=build_rt_b_vector_solid_solid(uθdi,uθsi,ρ1,cd1,cs1)
        c=build_rt_c_vector_solid_solid(uθsi,ρ1,cd1,cs1)

        println("")
        println("a")
        show(stdout, MIME"text/plain"(),a)
        println("")
        println("b")
        show(stdout, MIME"text/plain"(),b)
        println("")

        a = hcat(a[1,:],a[2,:],-im .* a[3,:], -im .* a[4,:])'
        b = @SVector [ b[1] , b[2] , -im*b[3] , -im*b[4] ]

        println("a")
        show(stdout, MIME"text/plain"(),a)
        println("")
        
        a1 = hcat(b,a[:,2],a[:,3],a[:,4])
        a2 = hcat(a[:,1],b,a[:,3],a[:,4])
        a3 = hcat(a[:,1],a[:,2],b,a[:,4])
        a4 = hcat(a[:,1],a[:,2],a[:,3],b)

        deta = det(a)
        Rd=det(a1)/det(a)
        Td=det(a2)/det(a)
        Rs=det(a3)/det(a)
        Ts=det(a4)/det(a)

        @show Rd,Td,Rs,Ts
        @show Rd+Td+Rs+Ts
        @show Rd^2+Td^2+Rs^2+Ts^2

        # Rd,Td,Rs,Ts = map(x-> x<0 ? 0. : x,(Rd,Td,Rs,Ts))


        zd1,zs1 = ρ1*cd1,ρ1*cs1
        zd2,zs2 = ρ2*cd2,ρ2*cs2
        cθd1,sθd1 = uθd1.x,uθd1.y
        cθs1,sθs1 = uθs1.x,uθs1.y
        cθd2,sθd2 = uθd2.x,uθd2.y
        cθs2,sθs2 = uθs2.x,uθs2.y
 
        # α=1 (longitudinal incident)
        # Rd_power = Rd                         # j=1 β=1
        # Td_power = Td*(zd2/zd1)*(cθd2/cθd1)   # j=2 β=1
        # Rs_power = Rs*(zs1/zd1)*(cθs1/cθd1)   # j=1 β=2
        # Ts_power = Ts*(zs2/zd1)*(cθs2/cθd1)   # j=2 β=2

        # # Rd_power,Td_power,Rs_power,Ts_power = map(x-> x<0 ? 0. : x,(Rd_power,Td_power,Rs_power,Ts_power))

        # @show Rd_power,Td_power,Rs_power,Ts_power
        # @show Rd_power+Td_power+Rs_power+Ts_power

        Rd_power = (Rd^2)                         # j=1 β=1
        Td_power = (Td^2)*(zd2/zd1)*(cθd2/cθd1)   # j=2 β=1
        Rs_power = (Rs^2)*(zs1/zd1)*(cθs1/cθd1)   # j=1 β=2
        Ts_power = (Ts^2)*(zs2/zd1)*(cθs2/cθd1)   # j=2 β=2

        # Rd_power,Td_power,Rs_power,Ts_power = map(x-> x<0 ? 0. : x,(Rd_power,Td_power,Rs_power,Ts_power))

        @show Rd_power,Td_power,Rs_power,Ts_power
        @show Rd_power+Td_power+Rs_power+Ts_power

      
        td, ts, rd, rs = solid_solid(θdi,ρ1,ρ2, cd1, cs1, cd2, cs2)
        @show td,ts,rd,rs
        @show  td+ts+rd+rs
        @show td^2,ts^2,rd^2,rs^2

          # α=1 (longitudinal incident)
          rd_power = rd^2                         # j=1 β=1
          td_power = td^2*(zd2/zd1)*(cθd2/cθd1)   # j=2 β=1
          rs_power = rs^2*(zs1/zd1)*(cθs1/cθd1)   # j=1 β=2
          ts_power = ts^2*(zs2/zd1)*(cθs2/cθd1)   # j=2 β=2
    
          @show rd_power,td_power,rs_power,ts_power
          @show rd_power+td_power+rs_power+ts_power

                # α=1 (longitudinal incident)
            rd_power = rd                         # j=1 β=1
            td_power = td*(zd2/zd1)   # j=2 β=1
            rs_power = rs*(zs1/zd1)   # j=1 β=2
            ts_power = ts*(zs2/zd1)   # j=2 β=2
    
            @show rd_power,td_power,rs_power,ts_power
            @show rd_power+td_power+rs_power+ts_power

              # α=1 (longitudinal incident)
              rd_power = rd                         # j=1 β=1
              td_power = td*(cθd2/cθd1)   # j=2 β=1
              rs_power = rs*(cθs1/cθd1)   # j=1 β=2
              ts_power = ts*(cθs2/cθd1)   # j=2 β=2
        
              @show rd_power,td_power,rs_power,ts_power
              @show rd_power+td_power+rs_power+ts_power

        
      




        return a,b


    else
        @show muθd1,muθd2,muθs1,muθs2
    end
    nothing
end

# a=build_rt_a_matrix_solid_solid(uθd1,uθd2,uθs1,uθs2,ρ1,ρ2,cd1,cs1,cd2,cs2)



# test()










































