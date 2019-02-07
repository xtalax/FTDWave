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
function Wave!(u, t, factor, δxx, δyy, s)
    u[:,:,t] = @views factor.*(∇(u[:,:,t-1], δxx, δyy) .+ s[:,:,t-1]) .+ 2.0.*u[:,:,t-1] .- u[:,:,t-2]

end

###############################################################################
# Main Start
###############################################################################


function EM_Propagate(xmin, xmax, Nx, ymin = xmin, ymax = xmax, Ny = Nx)


    ###############################################################################
    # Main Start
    ###############################################################################

        xmin, xmax = -5, 5
        ymin = xmin ; ymax = xmax

        Δx = (xmax-xmin)/(Nx-1) ; Δy = (ymax-ymin)/(Ny-1)
         x = xmin:Δx:xmax    ; y = ymin:Δy:ymax

        Δt = Δx/(c₀*2)
        global Nt = Int(2*(xmax-xmin)/(c₀*Δt))
        tmin = 0 ; tmax = tmin + Δt*Nt
        t = tmin:Δt:tmax

        k = 16

        f₀ = c₀/(k*Δx)
        ω₀ = 2*pi*f₀

        fs = 1/Δt
        fgrid = fs*(0:(Nt))/(Nt)

        w = rect.(fgrid/(2*f₀))
        wt= circshift(imag.(ifft(w)*Nt)/100, Nt/4)
        wt = unitize(wt)

        δxx = cu(δδ(Nx, Δx))
        δyy = cu(δδ(Ny, Δy))
        ###############################################################################
        # FIELD DEFINITION
        ###############################################################################
        f(x,y) = 0
        g(x,y) = 0
        source(x,y,t) = sin(ω₀*t)
        u = cu(zeros(Nx,Ny,Nt))
        u[:,:,2] =  cu([f(x[xi],y[yi]) for xi in 1:Nx, yi in 1:Ny])
        G = cu([g(x[xi],y[yi]) for xi in 1:Nx, yi in 1:Ny])


        S= zeros(Nx,Ny,Nt-1)

        S = cu([xi == 2 ? wt[ti] : 0.0 for xi in 1:Nx, yi in 1:Ny, ti in 1:(Nt-1)])

        ###############################################################################
        # MEDIUM INITIALISATION
        ###############################################################################
        # compute the medium
        refrindex = n(x, y)
        factor = cu((Δt.^2).*(c₀./refrindex).^2)

        factorbg = (Δt.*c₀).^2

        # callate the field at t = -1
        u[:,:,1] = u[:,:,2] .- 2 .*Δt.*G#setup the initial conditions
        #initialise background
        bg = copy(u)


    ###############################################################################
    # SIMULATE FIELD
    ###############################################################################
    #setup plots and progress meter
    prog = Progress(Nt+1, 1)
    # propagate and plot
    println("------- Simulating Field -------")
    for i in 3:Nt
        Wave!(u, i, factor, δxx, δyy, S)
        Wave!(bg, i, factorbg, δxx, δyy, S)

        next!(prog)
    end
    println("Done!")

    #gif(anim, filepath*"animation$Nt.gif", fps = 30)
    return x, y, t, bg, u, S, xmin+Δx*2
end

###############################################################################
# COMPILE GIF
###############################################################################
function compile_gif_simple(x,y,u)
    Nx, Ny, Nt = size(u)
    xmin = minimum(x) ;ymin = minimum(y)
    xmax = maximum(x) ;ymax = maximum(y)
    println("---------Compiling Animation-----------")
    pyplot()
    prog = Progress(Nt, 1)
    global anim = Animation()
    for i in 1:Nt
        cameraAngle = 2*pi*i/Nt

        plot = surface(x,y, u[:,:,i], xlims=(xmin,xmax), ylims=(ymin,ymax) , zlims=(-0.5,0.5))
        surface!(plot, camera=(15*cos(cameraAngle), 30))
        frame(anim, plot)

        next!(prog)
    end
    gif(anim, filepath*"animation$Nt.gif", fps = 30)
end

function compile_gif(x, y, bg, field, source, xval)
    pyplot(leg=false, ticks=nothing)
    Nx, Ny, Nt = size(field)
    prog = Progress(Nt-1, 1)
    global zs = zeros(0)
    global anim = Animation()

    anim = @animate for i in 1:(Nt-1)

        if bg[:,:,i] ≈ field[:,:,i]
            continue
            prog = Progress(Nt-1-i, 1)
        end

        # create a plot with 3 subplots and a stom layout
        l = @layout [a; b{0.2h}]
        surf = surface(
            x, y, field[: , : ,i],
            xlims=(xmin,xmax), ylims=(ymin,ymax) , zlims=(-0.5,0.5)

         );
        # induce a slight oscillating camera angle sweep, in degrees (azimuth, altitude)
        plot!(surf, camera=(15*cos(i), 40));

        # add a tracking line
        fixed_x = xval*ones(Nx)
        response = field[2, :, i] - bg[2, :, i]
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

xmin, xmax = -5, 5
ymin = xmin ; ymax = xmax


@time x, y, t, background, field, source, xval = EM_Propagate(xmin, xmax, 400)

#response = field.-background

#surface(x,y,response[:,:,end-10])

#compile_gif(x, y, background, field, xval)
@time compile_gif_simple(x,y,field)
