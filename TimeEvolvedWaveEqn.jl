using LinearAlgebra
using ProgressMeter
using ImageMagick
using FFTW
using Makie
using Rotations
#using Plots


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

function n(x,y, origin_prime, θ=0.0)
    Rθ =[ cos(θ) -sin(θ) ; sin(θ) cos(θ)]
    Nx = length(x)
    Ny = length(y)

    n = ones(Nx, Ny)
    d1 = 0.305
    d2 = 0.4


    for yi in 1:Ny, xi in 1:Nx
        r = [x[xi], y[yi]]
        r_prime = Rθ*(r.-origin_prime)
        if -d1/2  <= r_prime[1] <= d1/2
            if -d2/2  <= r_prime[2] <= d2/2

                n[xi,yi] = sqrt(2.3)

            end
        end
    end
    return n
end


# implements a reradiating boundary condition at all edges - rABC :https://personalpages.manchester.ac.uk/staff/fumie.costen/tmp/HuygensABC.pdf p 379
function rABC!(u,t)
    u[1, :, t] = u[2, :, t-1]
    u[end, :, t] = u[end-1, :, t-1]
    u[:, 1, t] = u[:, 2, t-1]
    u[:, end, t] = u[:, end-1, t-1]
end

function HeugensABC!(u, t, Δx, Δy, Δt)
    A⁺ = (1 - Δx/(Δt*c₀))
    A⁻ = 1/A⁺

    if Δx == Δy
        u[1,:,t] =     A⁺.*(u[1,:,t].-u[1,:,t-1])+u[1,:,t]
        u[end,:,t] = A⁻.*(u[end,:,t].-u[end,:,t-1])+u[end,:,t]
        u[:, 1,t] =    A⁺.*(u[1,:,t].-u[1,:,t-1])-u[1,:,t]
        u[:,end,t] = A⁻.*(u[:,end,t].-u[:,end,t-1])+u[:,end,t]
    else
        B⁺ = (1 - Δy/(Δt*c₀))
        B⁻ = 1/B⁺

        u[1,:,t] =     A⁺.*(u[1,:,t].-u[1,:,t-1])+u[1,:,t]
        u[end,:,t] = A⁻.*(u[end,:,t].-u[end,:,t-1])+u[end,:,t]
        u[:, 1,t] =    B⁺.*(u[1,:,t].-u[1,:,t-1])+u[1,:,t]
        u[:,end,t] = B⁻.*(u[:,end,t].-u[:,end,t-1])+u[:,end,t]
    end
end

# Wave Equation
function Wave!(u, t, Δx, Δy, Δt, factor, δxx, δyy, s)
    u[:,:,t] = factor.*(∇(u[:,:,t-1], δxx, δyy) .+ s[:,:,t-1]) .+ 2.0.*u[:,:,t-1] .- u[:,:,t-2]
    HeugensABC!(u, t, Δx, Δy, Δt)
end

###############################################################################
# Main Start
###############################################################################


function EM_Propagate(xmin, xmax, Nx, theta, ymin = xmin, ymax = xmax, Ny = Nx)


    ###############################################################################
    # Main Start
    ###############################################################################

        xmin, xmax = -5, 5
        ymin = xmin ; ymax = xmax

        Δx = (xmax-xmin)/(Nx-1) ; Δy = (ymax-ymin)/(Ny-1)
         x = xmin:Δx:xmax    ; y = ymin:Δy:ymax

        Δt = Δx/(c₀*2)
        global Nt = Int(floor(3*(xmax-xmin)/(c₀*Δt)))
        tmin = 0 ; tmax = tmin + Δt*Nt
        t = tmin:Δt:tmax

        k = 8

        f₀ = c₀/(k*Δx)
        ω₀ = 2*pi*f₀
        k₀ = c₀/ω₀

        fs = 1/Δt
        fgrid = fs*(0:(Nt))/(Nt)

        w = rect.(fgrid/(2*f₀))
        W = ifft(w)
        wt= circshift((real.(w)+imag.(w)*Nt)/100, floor(Nt/4))
        wt = unitize(wt)

        δxx = δδ(Nx, Δx)
        δyy = δδ(Ny, Δy)
        δ⁻t = δ⁻(Nt, Δt)
        HuygensABC⁺ = δ⁻t*(Δt-Δx/c₀)
        HeugensABC⁻ =  δ⁻t.*1/(Δt-Δx/c₀)
        ###############################################################################
        # FIELD DEFINITION
        ###############################################################################
        f(x,y) = 0
        g(x,y) = 0
        source(x,y,t) = sin(ω₀*t)
        u = zeros(Nx,Ny,Nt)
        u[:,:,2] =  [f(x[xi],y[yi]) for xi in 1:Nx, yi in 1:Ny]
        G = [g(x[xi],y[yi]) for xi in 1:Nx, yi in 1:Ny]

        S= zeros(Nx,Ny,Nt-1)

     S = [xi == 3 ? wt[ti] : 0.0 for xi in 1:Nx, yi in 1:Ny, ti in 1:(Nt-1)]
        #S = [xi == 3 ? sin(ω₀*t[ti]) : 0.0 for xi in 1:Nx, yi in 1:Ny, ti in 1:(Nt-1)]
        ###############################################################################
        # MEDIUM INITIALISATION
        ###############################################################################
        # compute the medium
        refrindex = n(x, y,[0,0], pi/6)
        factor = ((Δt.^2).*(c₀./refrindex).^2)

        factorbg = (Δt.*c₀).^2

        # callate the field at t = -1
        u[:,:,1] = u[:,:,2] .- 2*Δt*G .+ S[:,:,1]#setup the initial conditions
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
        Wave!(u, i, Δx, Δy, Δt, factor, δxx, δyy, S)
        Wave!(bg, i, Δx, Δy, Δt, factorbg, δxx, δyy, S)

        next!(prog)
    end
    println("Done!")
    return x, y, t, bg, u, S, xmin+Δx*2
end

###############################################################################
# COMPILE GIF
###############################################################################


function makie_animation(x,y,t,u)
    prog = Progress(Nt, 1)
    println("---------Animating Data-----------")
    scene = Scene(resolution = (1024,768), limits = FRect(x[1],x[end],y[1],y[end]),axis= frame);

    surf = surface!(scene, x, y, u[:,:,end],)[end]
    center!(scene)
    camera(scene)
    record(scene, filepath*"EM_animation_Makie_$(Nt)$([2,3,50]).mp4", 1:Nt) do i
        surf[3] = u[:,:,i]
        next!(prog)
    end
end
###############################################################################
# PROPAGATE
###############################################################################

xmin, xmax = -2, 2
ymin = xmin ; ymax = xmax


@time x, y, t, background, field, source, xval = EM_Propagate(xmin, xmax, 200, 2*pi)

#response = field.-background

#surface(x,y,response[:,:,end-10])

#compile_gif(x, y, background, field, xval)
@time makie_animation(x,y,t,field)
