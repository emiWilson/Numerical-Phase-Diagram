using Plots
#using Interact
#using Maxima
using Optim
#import Maxima: expand
using DelimitedFiles

using Dates; import Dates;
using TaylorSeries;

#all other parameters will be read in through readParams() file
ko = 2/sqrt(3) #radius of bragg ring
TaylorAbout = 0.05
TaylorExpOrder = 10
rho_bar = 0.7 #~average density - to convert real density to Gabiels ficticious density
T0 = 2.3

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
        num3 = split(str[5], "T")[1]

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

F_id_helper(no, phi, dB) = dB*(no+1)*log((no+1)*(1-b)/(1-(no+1)*b))-a1*no^2

affine(a) = a + Taylor1(typeof(a), TaylorExpOrder) #expand to 10th order about a
t = affine(TaylorAbout) #expand t about t = 0.05 (as gabriel does)
#evaluate the taylor expansion at no
F_id(no, phi, dB) = evaluate(dB*(t+1)*log((t+1)*(1-b)/(1-(t+1)*b))-a1*t^2, no)

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

function FF_phi(no, dB)
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
function plot_Free(dB)
        Free_eng = Array{Float64}(undef, length(no_range), 1)
        let index = 1
                for n in no_range
                        Free_eng[index] = FF(n, dB)
                        index = index + 1
                end

        end

        scatter(no_range, Free_eng)
        date = Dates.now()
        savefig("FreeEng")
end

function plot_Phi(dB)
        scatter(no_range, FF_phi.(no_range, dB))
        date = Dates.now()
        savefig("PhivsDensityPlot")
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

                n1[a] = parse(Float64, data[1])
                n2[a] = parse(Float64, data[3])
                close(f)
        end

        density = vcat(transpose(n1), transpose(n2))
        temp = vcat(T_range,T_range)

        #Data from PRE (1993)
        #"Pertubation weighted-density approximation: The phase diagram
        #        of a Lennard-Jones system"
        # - L.Mederos, G. Navasues, P. Tarazona and E. Chacon
        pre_density = [0.906, 1.004, 0, 0.855,
                              0.025, 0.729,
                              0.946, 1.040, 0.0598, 0.644,
                              0.966, 1.055, 0.1433,0.5087,
                              1.0583, 1.162]
        for i in 1:length(pre_density) #convert real density into gabriel's ficticious density
                pre_density[i] =(pre_density[i] - rho_bar )/ rho_bar
        end

        pre_temperature = [0.75, 0.75,0.75, 0.75,
                              1,1,
                              1.15, 1.15, 1.15, 1.15,
                              1.35, 1.35, 1.35, 1.35,
                              2.75,2.75]
        for i in 1:length(pre_temperature) #convert real density into gabriel's ficticious density
                pre_temperature[i] =pre_temperature[i]*T0
        end

        scatter(density, temp)
        scatter!(pre_density, pre_temperature)
        date = Dates.now()
        savefig("figures/fig$date.png")

        #write params to file
        toWrite = string("Tempertature field is : ", T_start, " to ", T_end, " with ",
                length(dB_range), " points with spacing of ", T_step, "\n",
                "Density field is : ", no_start, " to ", no_end, " with ",
                length(no_range), " points with spacing of ", no_step, "\n",
                "b = ", b, "\na1 = ", a1, "\na2 = ", a2_printer, "\na3 = ", a3_printer,
                "\na4 = ", a4,"\nBx = ", Bx,
                "\nTaylor expand ideal part about n0 = ", TaylorAbout, " to order ",
                TaylorExpOrder, "\n T0 is ", T0, "\n Avg density is ", rho_bar)

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
