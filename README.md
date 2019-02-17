# FTDWave - Julia 1.1.0

Have a look in the Animate branch for the working code - i'm having trouble merging it to master

solutions to the scalar 2d wave eqn with the finite difference method


This package is written in a programming language called julia, see more here: julialang.org

You can chose any function to animate in the 'Main.jl' file, along with a lot of other options. 
Make sure to set the filepath to something sensible. As default no source is selected, you can change this in the ftd_propagate args

The medium is defined in the refrindex function in main
The starting field is defined in the f(x,y) function in main
The starting time derivative of the field is also defined in main and is called g(x,y) 

The simulation is defined in the FTDWave.jl file, using some utility functions from the other modules in this package and also
some external packages.

Don't hesitate to send me a message if you run in to trouble or have any questionsp

Please submit any aesthetic results to Aesthetic function graphposting on FB :)

Add dependencies as follows:
> using Pkg 
> Pkg.add("depends1")
> Pkg.add("depends2")
> ...

List of depends:

DSP
FFTW
ProgressMeter
Makie
Plots

