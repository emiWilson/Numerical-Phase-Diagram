# Numerical-Phase-Diagram

```bash
$ alias julia="/Applications/Julia-1.1.app/Contents/Resources/julia/bin/julia" \\
$ julia GetFree.jl \\
$ make \\
$ for i in {1 .. [T_numPts]}; do ./run [no_numPts] 2 in/data${i}.txt out/covex${i} no_atart 0 no_step 0 \\
$ julia Plotting.jl
```

Use the free energy featured in Gabriel Kocher's paper equation 55 to construct a phase diagram for the two dimensional hexagonal lattice. Code from Nathan Smith (maxima package - julia) and Mathew Smith (covex hull - c++) is used. 

Task list
- [x] Do the basic Automation for eqn 55 phase diagram code
- [] Figure out why my free energy is getting cut off
- [] Read value from text file into bash (may help automize more)
- [] Delete contents of /in and /out files before each new run
- [] timestamp image output in name and save in a folder labelled "figures"
- [] get rid of unecessary files
- [] change all dB -> T variables
- [] get Nathan's Maxima code to work
- [] supress GKS file popup
- [] write nice little intro here
- [] comment the shit out of your code :) ... I guess it's pretty self explanatory but just to be safe. 
