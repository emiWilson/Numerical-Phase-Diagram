using Plots
#include("GetFree.jl")
function getNum(str)
        str = split(str, "=")[1]
        str = split(str, " ")[1]#just incase there was a space before the =
        return parse(Float64, str)
end

function readParams()
        f = open("Params.txt")
        lines = readlines(f)

        global T_start = getNum(lines[2])
        global T_end = getNum(lines[3])
        global T_step = getNum(lines[4])


        global no_start = getNum(lines[7])
        global no_end = getNum(lines[8])
        global no_step = getNum(lines[9])


        global b = getNum(lines[12])
        global a1 = getNum(lines[13])
        global a2 = getNum(lines[14])#will get slightly more complicated when want a function, but should actually be ok
        global a3 = getNum(lines[15])
        global a4 = getNum(lines[16])
        global Bx = getNum(lines[17])




end
function conv_hull()
        readParams()
        T_range = T_start:T_step:T_end
        T_numPts = length(T_range)

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


        scatter(transpose(n1), T_range)
        scatter!(transpose(n2), T_range)
        savefig("fig.png")

end
conv_hull()
