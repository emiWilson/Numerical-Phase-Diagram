using Plots
include("GetFree.jl")
function Read_Common_Tangent_Data()
        num_dB = 301
        n1 = Array{Float64}(undef, 1, num_dB)
        n2 = Array{Float64}(undef, 1, num_dB)
        for a in 1:num_dB
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
        dB= 1:0.01:4
end
