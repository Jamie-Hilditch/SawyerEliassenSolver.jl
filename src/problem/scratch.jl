struct Scratch 
    rtmp1 :: Array{Float64,2}
    rtmp2 :: Array{Float64,2}
    rtmp3 :: Array{Float64,2}
    ctmp :: Array{ComplexF64,2}

    function Scratch(NX,NZ)
        SX = NX ÷ 2 + 1
        new(
            zeros(Float64,NX,NZ),
            zeros(Float64,NX,NZ),
            zeros(Float64,NX,NZ),
            zeros(ComplexF64,SX,NZ)
        )
    end
end
