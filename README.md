# Numerical-Phase-Diagram

$alias julia="/Applications/Julia-1.1.app/Contents/Resources/julia/bin/julia"
$ julia GetFree.jl
$ make
$ for i in {1 .. [T_numPts]}; do ./run [no_numPts] 2 in/data${i}.txt out/covex${i} no_atart 0 no_step 0
$ julia Plotting.jl
