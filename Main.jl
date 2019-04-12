include("GetFree.jl")
include("Fit.jl")

#rho_bar = 0.55 #~average density - to convert real density to Gabiels ficticious density
#T0 = 2.6 #temperature scaling param tao (in PFC) = T (real) / T0


GetFree(5.5, 0.5)
#=
data = readdlm("figures/DensTempData5.5.txt", ',')
T = data[:,1]
D = data[:,2]



fit = fitPRE()
rho_bar = fit[2]
T0 = fit[3]
PRE_data = scalePRE(rho_bar, T0)

#make phase diagram plot
scatter(D, T, label = "PFC")
scatter!(PRE_data)
#scatter!(ChiSq[5], ChiSq[4], label = "PRE")
#scatter!(ChiSq[3], ChiSq[2], color = "red", label = "near")
=#
