using Plots
#using Interact
#using Maxima
using Optim
#import Maxima: expand
using DelimitedFiles

using Dates; import Dates;

#all other parameters will be read in through readParams() file
ko = 2/sqrt(3) #radius of bragg ring

#get number from input line
function getNum(str)
        str = split(str, "=")[1]
        str = split(str, " ")[1]#just incase there was a space before the =
        return parse(Float64, str)
end

#get coefficients of input function (in string form)
#should be of form [num1] + [num2]T + [num3]T^2
function getFcn(str)
        str = split(str, "=")[1]
        str = split(str, " ")#just incase there was a space before the =

        #get rid of the T's
        num1 = str[1]
        num2 = split(str[3], "T")[1]
        num3 = split(str[3], "T")[1]

        num1 = parse(Float64, num1)

        num2 = parse(Float64, num2)
        num3 = parse(Float64, num3)
        return num1, num2, num3

end
f = open("Params.txt")
lines = readlines(f)

#look into function scope, no need to brute force this one
T_start = getNum(lines[2])
T_end = getNum(lines[3])
T_step = getNum(lines[4])

no_start = getNum(lines[7])
no_end = getNum(lines[8])
no_step = getNum(lines[9])


b = getNum(lines[12])
a1 = getNum(lines[13])

coeff_a2 = getFcn(lines[14])
a2(dB) = eval(coeff_a2[1] + coeff_a2[2]*dB + coeff_a2[3]*dB^2)
a2_printer = string(coeff_a2[1], " + ", coeff_a2[2], "T  +  ", coeff_a2[3], "T^2")

coeff_a3 = getFcn(lines[15])
a3(dB) = eval(coeff_a3[1] + coeff_a3[2]*dB + coeff_a3[3]*dB^2)
a3_printer = string(coeff_a3[1], " + ", coeff_a3[2], "T  +  ", coeff_a3[3], "T^2")

a4 = getNum(lines[16])
Bx = getNum(lines[17])

Pi = pi #different notation between maple and julia

T_range = T_start:T_step:T_end
T_numPts = length(T_range)

no_range = no_start:no_step:no_end
dB_range = T_start:T_step:T_end

length_dB = length(dB_range)
length_no = length(no_range)

F_id(no, phi, dB) = dB*(no+1)*log((no+1)*(1-b)/(1-(no+1)*b))-a1*no^2
"""
(1/((-1 + 1.05*b)^10)*(-1+b)^7)*(((1.628894626*no+1.628894626)*dB*b^17+(-26.91554453*no-26.91554453)*dB*b^16+(209.2852571*no+209.2852571)*dB*b^15+(-1017.040752*no-1017.040752)*dB*b^14+(3459.556021*no+3459.556021)*dB*b^13+(-8741.620595*no-8741.620595)*dB*b^12+(16990.43405*no+16990.43405)*dB*b^11+(-25945.73321*no-25945.73321)*dB*b^10+(31515.71121*no+31515.71121)*dB*b^9+(-30624.04845*no-30624.04845)*dB*b^8+(23805.23311*no+23805.23311)*dB*b^7+(-14719.12582*no-14719.12582)*dB*b^6+(7150.569638*no+7150.569638)*dB*b^5+(-2672.023812*no-2672.023812)*dB*b^4+(741.7025000*no+741.7025000)*dB*b^3+(-144.1125000*no-144.1125000)*dB*b^2+(17.50000000*no+17.50000000)*dB*b+(-1.000000000*no-1.000000000)*dB)*log((-1+b)/(-20+21*b))+((4.959206238*no+4.959206238)*dB-1.628894626*a1*no^2)*b^17+((-0.131250000062622e-2*no^6+0.787485714042191e-4*no^5+.500000000250000*no^9-no^10-1.90481230788464*10^(-9)*no^2+9.42857099062167*10^(-8)*no^3-3.28142852623614*10^(-6)*no^4-83.4963075217481*no-.112500001056244*no^8+0.149999995274879e-1*no^7-81.8674129009336)*dB+26.91554453*a1*no^2)*b^16+((635.967687600006+.814484209499955*no^4-.814448371399817*no^3+.814447328599890*no^2-6.38095238599999*no^9+11.76190475*no^10+.829187504799740*no^6-.815332495299827*no^5+2.07321429499995*no^8-.982678566699853*no^7+661.254337600006*no)*dB-209.2852571*a1*no^2)*b^15+((-10.7431498599998*no^4+11.2859185699992*no^3-11.8288776699993*no^2-3087.61467100002+38.3008314800000*no^9-66.42705967*no^10-9.73760214299894*no^6+10.2048248599991*no^5-15.4373299499997*no^8+10.0330130499992*no^7-3271.61327900003*no)*dB+1017.040752*a1*no^2)*b^14+((65.8924576200033*no^4-72.8274923800031*no^3+80.1704271399849*no^2+10492.8446700000-144.162761400000*no^9+239.8122952*no^10+53.5271152400019*no^6-59.3803742899980*no^5+66.4907881900006*no^8-50.7784342900007*no^7+11325.8868200000*no)*dB-3459.556021*a1*no^2)*b^13+((-26488.2202500000-249.244473299997*no^4+290.246701899986*no^3-336.350646200037*no^2+382.257315700000*no^9-619.7640643*no^10-181.960762399996*no^6+213.056684299988*no^5-191.748968099999*no^8+161.884258099999*no^7-29114.7341300000*no)*dB+8741.620595*a1*no^2)*b^12+((51434.3274799999+650.128329499978*no^4-798.776621000023*no^3+976.906300000133*no^2-759.054417600000*no^9+1213.759476*no^10+427.931620499963*no^6-527.282252899943*no^5+400.212368600001*no^8-361.055338099982*no^7+57549.4341900001*no)*dB-16990.43405*a1*no^2)*b^11+((-78469.7245699999-1238.48964799992*no^4+1607.50613299995*no^3-2080.64706699997*no^2+1168.95409400000*no^9-1857.422330*no^10-738.515187599984*no^6+953.884142899941*no^5-633.247290000001*no^8+596.661436699995*no^7-89345.0519099999*no)*dB+25945.73321*a1*no^2)*b^10+((95224.9393400000+1779.84580499990*no^4-2443.05982899988*no^3+3357.01660999985*no^2-1423.81677300000*no^9+2257.145383*no^10+968.043999999973*no^6-1304.04659999996*no^5+781.318204299997*no^8-756.112686699986*no^7+110295.345200000*no)*dB-31515.71121*a1*no^2)*b^9+((-92442.8732899999-1964.81577100010*no^4+2854.13241900004*no^3-4178.18634300021*no^2+1384.45416900000*no^9-2193.392741*no^10-983.656514300030*no^6+1373.24985200000*no^5-763.570719500003*no^8+749.846793800020*no^7-108888.178600000*no)*dB+30624.04845*a1*no^2)*b^8+((71791.0562399999+1679.97377599993*no^4-2583.06301399994*no^3+4044.46631899994*no^2-1076.28584600000*no^9+1705.007279*no^10+783.095149999986*no^6-1124.64083999996*no^5+594.478431899998*no^8-587.732458100003*no^7+85969.7993400000*no)*dB-23805.23311*a1*no^2)*b^7+((-1112.83670799998*no^4+1809.99789099997*no^3-3044.90525199998*no^2-44347.2638399999+665.483370999999*no^9-1054.231083*no^10-489.536627600013*no^6+717.502247199999*no^5-367.663091400003*no^8+364.333729500008*no^7-53973.7538599999*no)*dB+14719.12582*a1*no^2)*b^6+((566.719807600013*no^4-973.310398600005*no^3+1768.33975899997*no^2+21523.4655100000-323.292649500000*no^9+512.1467719*no^10+238.751038600020*no^6-354.389470000033*no^5+178.611187100002*no^8-177.074746700009*no^7+26616.1013300000*no)*dB-7150.569638*a1*no^2)*b^5+((-218.085879100001*no^4+394.114401899996*no^3-777.978156200003*no^2-8035.23293799999+120.807949800000*no^9-191.3789304*no^10-89.2911896200029*no^6+133.364638100000*no^5-66.7434023300002*no^8+66.1692650000009*no^7-10093.1667500000*no)*dB+2672.023812*a1*no^2)*b^4+((2228.30340900000+61.4753640499994*no^4-116.200421099997*no^3+250.988750099998*no^2-33.5339670000000*no^9+53.12311610*no^10+24.7855194299998*no^6-37.0888580499986*no^5+18.5266867599999*no^8-18.3673173799996*no^7+2842.39340900000*no)*dB-741.7025000*a1*no^2)*b^3+((-432.546217300006-12.0094052399989*no^4+23.5187509099991*no^3-56.0562499999996*no^2+6.51563642900000*no^9-10.32180027*no^10-4.81581654799989*no^6+7.20635031899957*no^5-3.59972785200000*no^8+3.56876243799999*no^7-560.158717200007*no)*dB+144.1125000*a1*no^2)*b^2+((52.4753147600004+1.45833700099992*no^4-2.91666677599992*no^3+7.75000000999996*no^2-.791212679500000*no^9+1.253406226*no^10+.584798609999980*no^6-.875088077199947*no^5+.437125423300000*no^8-.433365201900003*no^7+68.9753147600004*no)*dB-17.50000000*a1*no^2)*b+(0.452121531400000e-1*no^9-3.99573227400001*no+0.247637258600004e-1*no^7+0.500050330000007e-1*no^5-0.7162321291e-1*no^10+.166666672300000*no^3-0.249785955700000e-1*no^8-.500000000000001*no^2-0.334170633800008e-1*no^6-0.833335433399992e-1*no^4-2.99573227400001)*dB+1.000000000*a1*no^2)
"""
#dB*(no+1)*log((no+1)*(1-b)/(1-(no+1)*b))-a1*no^2
F_corr(no, phi, dB) = -(1/16)*dB*Bx*exp(-dB)*(8*no^2*exp(-2/3)+3*phi^2)
F_multi(no, phi, dB) = 8*sqrt(3)*(27*sqrt(3)*a4*no^2*phi^2*(1/128)+
        9*sqrt(3)*a4*no*phi^3*(1/128)+135*sqrt(3)*a4*phi^4*(1/4096)+
        9*sqrt(3)*a3(dB)*no*phi^2*(1/64)+3*sqrt(3)*a3(dB)*phi^3*(1/128)+
        9*sqrt(3)*a2(dB)*phi^2*(1/128))*(1/9)

F(no, phi, dB) = F_id(no, phi, dB) + F_corr(no, phi,dB) + F_multi(no, phi, dB)

function FF(no, dB)
        sol = optimize(phi -> F(no, phi, dB), 0.0, 10.0)
        return sol.minimum
end

function FF_ϕ(no, dB)
        sol = optimize(phi -> F(no, phi, dB), 0.0, 10.0)
        return sol.minimizer
end

function make_free()

        #Find equillibrium amplitude for given density, temperature


        let db_index = 1
                for db in dB_range
                        n_index = 1
                        for n in no_range
                                F_arr[db_index, n_index] = FF(n, db)
                                n_index += 1
                        end
                        db_index += 1
                end
        end


        for ii in 1:length(dB_range)
                str = string("in/data" ,ii, ".txt")
                for jj in 1:length(no_range)
                        data[jj, 1] = F_arr[ii, jj]
                        data[jj, 2] = F_arr[ii, jj]
                end

                writedlm(str, data, ',')
        end

end


function conv_hull()
        n1 = Array{Float64}(undef, 1, T_numPts)
        n2 = Array{Float64}(undef, 1, T_numPts)
        for a in 1:T_numPts
                filename = string("out/convex" ,a, "_coexist_lines.txt")
                f = open(filename)
                line = 1
                X = readline(f)
                pieces = split(X, ',')

                data = ["", "", "", ""]
                let index = 1
                        for str in pieces
                                for ch in str
                                        if ch == '}' || ch == '{'
                                                #do nothing
                                        else
                                                data[index] = string(data[index], ch)

                                        end
                                end
                                index = index + 1
                        end
                end

                n1[a] = parse(Float64, data[1]) + 0.25
                n2[a] = parse(Float64, data[3]) + 0.25
                close(f)
        end

        density = vcat(transpose(n1), transpose(n2))
        temp = vcat(T_range,T_range)

        #Data from PRE (1993)
        #"Pertubation weighted-density approximation: The phase diagram
        #        of a Lennard-Jones system"
        # - L.Mederos, G. Navasues, P. Tarazona and E. Chacon
        pre_density = [0.906, 1.004,
                        0.946, 1.040,
                        0.966, 1.055]

        pre_temperature = [0.75, 0.75,
                        1.15, 1.15,
                        1.35, 1.35]

        scatter(density, temp/3.3)
        scatter!(pre_density, pre_temperature)
        date = Dates.now()
        savefig("figures/fig$date.png")

        #write params to file
        toWrite = string("Tempertature field is : ", T_start, " to ", T_end, " with ",
                length(dB_range), " points with spacing of ", T_step, "\n",
                "Density field is : ", no_start, " to ", no_end, " with ",
                length(no_range), " points with spacing of ", no_step, "\n",
                "b = ", b, "\na1 = ", a1, "\na2 = ", a2_printer, "\na3 = ", a3_printer,
                "\na4 = ", a4,"\nBx = ", Bx)
        open("figures/data$date.txt", "w") do f
                write(f, toWrite)
        end


end

function main()


        global data = Array{Float64}(undef, length(no_range), 2)
        println("\n---------------INPUT PARAMETERS---------------")
        println("Tempertature field is : ", T_start, " to ", T_end, " with ",
                length(dB_range), " points with spacing of ", T_step)
        println("Density field is : ", no_start, " to ", no_end, " with ",
                length(no_range), " points with spacing of ", no_step)
        println("b = ", b, "\na1 = ", a1, "\na2 = ", a2_printer, "\na3 = ", a3_printer,
                "\na4 = ", a4,"\nBx = ", Bx)
        println("----------------------------------------------\n")
        global F_arr = Array{Float64}(undef, length(dB_range), length(no_range))
        make_free()

        for i in 1:length_dB
                run(`./run $length_no 2 in/data$i.txt out/convex$i $no_start 0 $no_step 0 `)
        end

        conv_hull()


end

main()

"""
coeff_a2 = getFcn(lines[14])
a2(dB) = eval(coeff_a2[1] + coeff_a2[2]*dB + coeff_a2[3]*dB^2) #getNum(lines[14])#
a2_printer = string(coeff_a2[1], " + ", coeff_a2[2], "T  +  ", coeff_a2[3], "T^2")

coeff_a3 = getFcn(lines[15])
a3(dB) = eval(coeff_a3[1] + coeff_a3[2]*dB + coeff_a3[3]*dB^2) #getNum(lines[15])
a3_printer = string(coeff_a3[1], " + ", coeff_a3[2], "T  +  ", coeff_a3[3], "T^2")
"""


"""

"""
