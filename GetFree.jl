"""
Using the free energy in chapter 8 of Gaby's thesis.
Use Maple (or possibly Nathan's Maxima package in julia) to simplify the free
energy (carry out the integrals).
Then use Mathew Seymours convex hull code (in C++) to find the common tangent
for each dB value.
Syntax for running MS's code in command line is:
    ./run [columns] [rows] [file] [output prefix] [xmin] [ymin] [xstep] [ystep]
Then plot in Julia.

Right now just make the phase diagram by outputting the free energy for a given
temperature then run in MS's code then make figure back in julia
"""


#Made largely using Nathan's code
#using Plots
#using Interact
#using Maxima
using Optim
#import Maxima: expand
using DelimitedFiles

#Parameters
const b = 0.5
const a1 = 7
const a2(dB) = 5.45#0.4958 + 11.6438*dB + 1.8*dB^2
const a3(dB) = -2#-0.220867 - 4.09696*dB + 1.8*dB^2
const a4 = 0.1
const Bx = 3.5
const ko = 2*(2*pi)/sqrt(3) #ko = 2*q/sqrt(3) and q = 2*pi



function make_text_file()
        F_id(no, phi, dB) = dB*(no+1)*log((no+1)*(1-b)/(1-(no+1)*b))-a1*no^2
        F_corr(no, phi, dB) = -dB*Bx*no^2*exp(-(1/2)*ko^2)/(2*exp(dB))-
            3*dB*Bx*no*phi*exp(-2/3)*exp(2*sqrt(3)*ko*(1/3))*exp(-(1/2)*ko^2)/(4*exp(dB))-
            3*dB*Bx*phi*no*exp(-(1/2)*ko^2)/(4*exp(dB))-
            9*dB*Bx*phi^2*exp(-2/3)*exp(2*sqrt(3)*ko*(1/3))*exp(-(1/2)*ko^2)/(8*exp(dB))
        F_multi(no, phi, dB) = (2/sqrt(3))*
                (9*phi^2*sqrt(3)*no^2*a4/32 + 3*phi^3*sqrt(3)*a4*no/32 +
                45*phi^4*sqrt(3)*a4/1024 + 2*phi^2*sqrt(3)*a3(dB)*no/16 +
                phi^3*sqrt(3)*3/32 + 3*phi^2*sqrt(3)*a2(dB)/32)

        F(no, phi, dB) = F_id(no, phi, dB) + F_corr(no, phi,dB) + F_multi(no, phi, dB)

        """ Seventh order taylor exp on log term
        dB*(no+1)*(no/(1-b)+(-b/((-1+b)*(1-b))+1/((2*(-1+b))*(1-b)))*no^2+(b^2/((-1+b)^2*(1-b))-b/(3*(-1+b)^2*(1-b))-(2*b-1)/(3*(-1+b)^2*(1-b)))*no^3+(-b^3/((-1+b)^3*(1-b))+b^2/(4*(-1+b)^3*(1-b))+(2*b-1)*b/(4*(-1+b)^3*(1-b))+(3*b^2-3*b+1)/(4*(-1+b)^3*(1-b)))*no^4+(b^4/((-1+b)^4*(1-b))-b^3/(5*(-1+b)^4*(1-b))-(2*b-1)*b^2/(5*(-1+b)^4*(1-b))-(3*b^2-3*b+1)*b/(5*(-1+b)^4*(1-b))-(4*b^3-6*b^2+4*b-1)/(5*(-1+b)^4*(1-b)))*no^5+(-b^5/((-1+b)^5*(1-b))+b^4/(6*(-1+b)^5*(1-b))+(2*b-1)*b^3/(6*(-1+b)^5*(1-b))+(3*b^2-3*b+1)*b^2/(6*(-1+b)^5*(1-b))+(4*b^3-6*b^2+4*b-1)*b/(6*(-1+b)^5*(1-b))+(5*b^4-10*b^3+10*b^2-5*b+1)/(6*(-1+b)^5*(1-b)))*no^6)-a1*no^2-dB*Bx*no^2*exp(-(1/2)*ko^2)/(2*exp(dB))-3*dB*Bx*no*phi*exp(-2/3)*exp(2*sqrt(3)*ko*(1/3))*exp(-(1/2)*ko^2)/(4*exp(dB))-3*dB*Bx*phi*no*exp(-(1/2)*ko^2)/(4*exp(dB))-9*dB*Bx*phi^2*exp(-2/3)*exp(2*sqrt(3)*ko*(1/3))*exp(-(1/2)*ko^2)/(8*exp(dB))+2*sqrt(3)*(9*phi^2*sqrt(3)*a4*no^2*(1/32)+3*phi^3*sqrt(3)*a4*no*(1/32)+45*phi^4*sqrt(3)*a4*(1/1024)+3*phi^2*sqrt(3)*a3*no*(1/16)+(1/32)*phi^3*sqrt(3)*a3+3*phi^2*sqrt(3)*a2*(1/32))*(1/3)
        """

        #Find equillibrium amplitude for given density, temperature
        function FF(no, dB)
                sol = optimize(phi -> F(no, phi, dB), 0.0, 10.0)
                return sol.minimum
        end

        no_range = -0.9999:0.001:1
        dB_range = 1:0.01:4
        F_arr = Array{Float64}(undef, length(dB_range), length(no_range))


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

        data = Array{Float64}(undef, length(no_range), 2)
        for ii in 1:length(dB_range)
                str = string("data" ,ii, ".txt")
                str = string("in/", str)
                for jj in 1:length(no_range)
                        data[jj, 1] = F_arr[ii, jj]
                        data[jj, 2] = F_arr[ii, jj]
                end
                writedlm(str, data, ',')
        end

        print("Made free energy text files")
end
