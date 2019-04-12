include("GetFree.jl")
include("Fit.jl")


#-2.05 + 2.43T + 0.055T^2 = a2
#0.65 + -1.083T + 0.058T^2 = a3
#GetFree(a, b, a2_1, a2_2, a2_3, a3_1, a3_2, a3_3, a4, bx)
#GetFree(6, 0.5, -2.05, 2.43, 0.055, 0.65, -1.083, 0.058, 0.1, 3.5)
function forRange()
    val_range = -1.8:-0.02:-3
    for v in val_range
        GetFree(5.25, 0.5, v, 2.43, 0.055, 0.65, -1.083, 0.058, 0.1, 3.5)
    end
end

GetFree(6, 0.5, -2.05, 2.43, 0.055, 0.65, -1.083, 0.058, 0.1, 3.5)
