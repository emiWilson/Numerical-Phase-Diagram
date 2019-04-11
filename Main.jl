include("GetFree.jl")

#rho_bar = 0.55 #~average density - to convert real density to Gabiels ficticious density
#T0 = 2.6 #temperature scaling param tao (in PFC) = T (real) / T0
#=
a_range = 5.5:0.5:8.5

for a in a_range
    GetFree(a, 0.5)
end
=#

rho_bar = 0.55 #~average density - to convert real density to Gabiels ficticious density
T0 = 2.6 #temperature scaling param tao (in PFC) = T (real) / T0

data = readdlm("figures/DensTempData5.5.txt", ',')
T = data[:,1]
D = data[:,2]

pre_density = [0.906, 1.004, 0, 0.855,
                              0.025, 0.729,
                              0.946, 1.040, 0.0598, 0.644,
                              0.966, 1.055, 0.1433,0.5087]
                              #1.0583, 1.162]

pre_temperature = [0.75, 0.75,0.75, 0.75,
                      1,1,
                      1.15, 1.15, 1.15, 1.15,
                      1.35, 1.35, 1.35, 1.35]
                      #2.75,2.75]
function scalePRE(rho_bar_in, T0_in)
        for i in 1:length(pre_temperature)
                #convert real density into PFC's ficticious density
                pre_density[i] =(pre_density[i] - rho_bar_in)/ rho_bar_in
                #scale real temperature to PFC temperature
                pre_temperature[i] =pre_temperature[i]*T0_in
        end

        return pre_density, pre_temperature
end

scaled = scalePRE(rho_bar, T0)
sD = scaled[1]
sT = scaled[2]
function getChiSq()
        #returns index of nearest point to position d, t on plot
        function findnearest(d,t)
                distance = Array{Float64}(undef, length(D))
                for i in 1:length(D)
                        distance[i] = (D[i] - d)^2 + (T[i] - t)^2
                end
                #findmin(A, dims) -> (minval, index)
                min = findmin(distance)
                println(min[2])
                return min[2]
        end

        nearestIndex = Array{Int}(undef, length(sD))
        for i in 1:length(sD)
                nearestIndex[i] = findnearest(sD[i], sT[i])
        end

        TT = Array{Float64}(undef, length(sD))
        DD = Array{Float64}(undef, length(sD))

        for i in 1:length(sD)
                j = nearestIndex[i]

                TT[i] = T[j]
                DD[i] = D[j]
        end

        #ADD CHI SQ CALC HERE!!!
end
#make phase diagram plot
scatter(D, T, label = "PFC")
scatter!(sD, sT, label = "PRE")
scatter!(DD, TT, color = "red", label = "near")
