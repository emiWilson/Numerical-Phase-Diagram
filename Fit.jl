
function scalePRE(rho_bar_in, T0_in)

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


        for i in 1:length(pre_temperature)
                #convert real density into PFC's ficticious density
                pre_density[i] =(pre_density[i] - rho_bar_in)/ rho_bar_in
                #scale real temperature to PFC temperature
                pre_temperature[i] = pre_temperature[i]*T0_in
        end

        return pre_density, pre_temperature
end

#good default input values are rho_bar = 0.55 and T0 = 2.6
function getChiSq(rho_bar, T0)
        scaled = scalePRE(rho_bar, T0)
        sD = scaled[1]
        sT = scaled[2]
        #returns index of nearest point to position d, t on plot
        #for each point in PRE paper find nearest point on my curve
        #return the distance between my point and the PRE point, will give accuracy of fit
        function findNearest(d,t)
                distance = Array{Float64}(undef, length(D))
                for i in 1:length(D)
                        distance[i] = (D[i] - d)^2 + (T[i] - t)^2
                end
                #findmin(A, dims) -> (minval, index)
                min = findmin(distance)
                return min
        end

        nearestIndex = Array{Int}(undef, length(sD))
        ChiSum = 0
        for i in 1:length(sD)
                nn = findNearest(sD[i], sT[i])
                nearestIndex[i] = nn[2]
                ChiSum = ChiSum + nn[1]
        end

        TT = Array{Float64}(undef, length(sD))
        DD = Array{Float64}(undef, length(sD))

        for i in 1:length(sD)
                j = nearestIndex[i]
                TT[i] = T[j]
                DD[i] = D[j]
        end
        #returns the value of Chi Squares sums which is an indication of the goodness
        #of fit, which
        return ChiSum
end

#find chi squares for an array of rho, T values(scaling factors for the PRE paper)
function testChi(p_start, p_step, p_end, T_start, T_step, T_end)
        #rho_list = 0.45:0.1:0.55
        #T_list = 2:0.1:3

        #list of densities/ temp scaling parameters to try
        rho_list = p_start:p_step:p_end
        T_list = T_start:T_step:T_end
        leng = length(rho_list)*length(T_list)

        oldChi = 100000

        #empty arrays to store min values in

        ChiArr = Array{Float64}(undef, (length(rho_list)*length(T_list)))
        pArr = Array{Float64}(undef, (length(rho_list)*length(T_list)))
        tArr = Array{Float64}(undef, (length(rho_list)*length(T_list)))

        let index = 1
                for i in 1:length(rho_list)
                        for j in 1:length(T_list)
                                p = rho_list[i]
                                t = T_list[j]
                                ChiArr[index] = getChiSq(p,t)

                                pArr[index] = p
                                tArr[index] = t
                                index = index + 1

                        end
                end
        end

        min = findmin(ChiArr)
        val = min[1]
        index = min[2]
        #(0.5747588264462813, 0.55, 2.6)
        return val, pArr[index], tArr[index]
end

#maybe would be good to make it recursive... but don't need insane accuracy, yet
#maybe have some sort of convergence test too
function fitPRE(minDen = 0.4, maxDen = 0.5, minT = 2, maxT = 4)
        #first interation
        numSteps = 20
        DenStep = (maxDen - minDen)/numSteps
        TStep = (maxT - minT)/numSteps
        tester1 = testChi(minDen, DenStep, maxDen, minT, TStep, maxT)

        it1dens = tester1[2]
        it1temp = tester1[3]
        minDens2 = it1dens - DenStep
        maxDens2 = it1dens + DenStep

        minT2 = it1temp - TStep
        maxT2 = it1temp + TStep
        DenStep2 = (maxDens2 - minDens2)/numSteps
        TStep2 = (maxT2 - minT2)/numSteps

        tester2 = testChi(minDens2, DenStep2, maxDens2, minT2, TStep2, maxT2)
        return tester2
end
