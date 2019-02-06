using LinearAlgebra
using Plots
using ProgressMeter
using ImageMagick
using FFTW
FFTW.set_num_threads(4)

function pwd_path!(LOAD_PATH)
    for path in LOAD_PATH
        if path == pwd()
            return LOAD_PATH
        end
    end
    push!(LOAD_PATH, pwd())
end

pwd_path!(LOAD_PATH)

include("includes.jl")
##############################################################################
# DEFINE CONSTANTS
##############################################################################
filepath = homedir()*"/Videos/MathAnimations/"

const c₀ = 3*10^8

##############################################################################
# DECLARE FUNCTIONS
##############################################################################


# Medium definition
function n(x,y)
    Nx, Ny = [length(x), length(y)]
    n = ones(Nx, Ny)
    d1 = 3
    d2 = 3

    for  ix in 1:Nx, iy in 1:Ny
        if -d1/2  <= x[ix] <= d1/2
            if -d2/2  <= y[iy] <= d2/2
                n[ix,iy] = sqrt(2.4)
            end
        end
    end
    return n
end


# Wave Equation
function Wave!(u, δxx, δyy,  t, factor, S=0)
    u[:,:,t+1] = 2*u[:,:,t] .+ factor.*(⊗(δxx, u[:,:,t], 1) .+ ⊗(δyy, u[:,:,t], 2)+S[:,:,t+1]).- u[:,:,t-1]
    return u
end

###############################################################################
# Main Start
###############################################################################
    Nx = Ny = 100
    xmin, xmax = -5, 5
    ymin = xmin ; ymax = xmax


    Δx = (xmax-xmin)/(Nx-1) ; Δy = (ymax-ymin)/(Ny-1)
     x = xmin:Δx:xmax    ; y = ymin:Δy:ymax

    Δt = Δx/(c₀*2)
    Nt = Int(4*(xmax-xmin)/(c₀*Δt))
    tmin = 0 ; tmax = tmin + Δt*Nt
    t = tmin:Δt:tmax

    f₀ = c₀/(4*Δx)
    ω₀ = 2*pi*f₀

    fs = 1/Δt
    fgrid = fs*(0:(Nt))/(Nt)

    ###############################################################################
    # CALCULATE DERIVATIVE OPERATORS
    ###############################################################################
    δx⁻ = δ⁻(Nx, Δx)
    δx⁺ = δ⁺(Nx, Δx)
    δx = δ(Nx, Δx)
    δxx = δδ(Nx, Δx)

    δy⁻ = δ⁻(Ny, Δy)
    δy⁺ = δ⁺(Ny, Δy)
    δy = δ(Nx, Δy)

    δyy = δδ(Ny, Δy)
    δt⁻ = δ⁻(Nx, Δt)


    ###############################################################################
    # FIELD DEFINITION
    ###############################################################################
    f(x,y) = 0
    g(x,y) = 0
    source(x,y,t) = sin(ω₀*t)

    ###############################################################################
    # FIELD CALCULATION
    ###############################################################################
    # compute the source

    S= [xi == 2 ? source(x[xi], y[yi], t[ti]) : 0.0 for xi in 1:Nx, yi in 1:Ny, ti = 1:(Nt+1)]
    G = [g(x[xi],y[yi]) for xi in 1:Nx, yi in 1:Ny]
#    G[2,:] = ω₀
    # compute the medium
    refrindex = n(x, y)
    factor = (Δt.^2).*(c₀./refrindex).^2

    factorbg = (Δt.*c₀).^2

    # calculate the initial field
    u = zeros(Nx,Ny,Nt+2)
    U =  [f(x[xi],y[yi]) for xi in 1:Nx, yi in 1:Ny]
    u[:,:,2] = U

    # calculate the field at t = -1
    u[:,:,1] = u[:,:,2] .- 2*Δt.*G #setup the initial conditions
    #initialise background
    bg = copy(u)
    ###############################################################################
    # CALCULATE AND GENERATE VISUALISATION
    ###############################################################################
    #setup plots and progress meter
    pyplot()
    prog = Progress(Nt+1, 1)
    #surface(x,y,factor)
    # propagate and plot

    anim = Animation()
    for i in 2:(Nt+2)
        #if u != bg
            cameraAngle = (i-2)/Nt*2pi
            plot = surface(x,y, u[:,:,i], xlims=(xmin,xmax), ylims=(ymin,ymax) , zlims=(-0.5,0.5))
            surface!(plot, camera=(15*cos(cameraAngle), 40))
            frame(anim, plot)
        #end
        i == Nt+2 && break
        Wave!(u, δxx, δyy, i, factor, S)
        Wave!(bg, δxx, δyy, i, factorbg, S)


        next!(prog)
    end
    ###############################################################################
    # COMPILE GIF
    ###############################################################################
    gif(anim, filepath*"animation$Nt.gif", fps = 30)

function compile_gif(x, y, field, source, xval)
    pyplot(leg=false, ticks=nothing)
    Nx, Ny, Nt = size(field)
    prog = Progress(Nt-1, 1)
    global zs = zeros(0)
    global anim = Animation()

    anim = @animate for i in 1:(Nt-1)
        # create a plot with 3 subplots and a custom layout
        l = @layout [a; b{0.2h}]
        surf = surface(
            x, y, field[: , : ,i],
            xlims=(xmin,xmax), ylims=(ymin,ymax) , zlims=(-0.5,0.5)

         );
        # induce a slight oscillating camera angle sweep, in degrees (azimuth, altitude)
        plot!(surf, camera=(15*cos(i), 40));

        # add a tracking line
        fixed_x = xval*ones(Nx)
        response = field[2, :, i] - source[2, :, i]
        #response_sum = sum(response)
        plot!(surf, fixed_x, y, response, line = (:black, 5, 0.2));
        # add to and show the tracked values over time
        global zs = vcat(zs, sum(response))
        bar = plot(1:i, zs, palette = cgrad(:blues).colors);
        #assemble the plot
        plot(surf, bar, layout=l)

        # increment the progress bar
        next!(prog)
    end
    gif(anim, filepath*"animation$Nt.gif", fps = 30)

end
###############################################################################
# PROPAGATE
###############################################################################



#compile_gif(x, y, field, source, xval)
