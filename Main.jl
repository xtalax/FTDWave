# !~/bin/julia
    pwd() ∉ LOAD_PATH && push!(LOAD_PATH, pwd())


    using EZAnimate
    using DiscreteAxis
    using FTDWave
    using SignalProcessing
    using LinearAlgebra

    #my modules

    ##############################################################################
    # DEFINE CONSTANTS
    ##############################################################################

    const c₀ = 3*10^8

    ##############################################################################
    # DEFINE MEDIUM
    ##############################################################################


    # Medium definition
    function refrindex(x,y)
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
    function domain(y::LinearAxis)
        maximum(y.pts) - minimum(y.pts)
    end
    function domain(ymin,ymax)
        ymax-ymin
    end

#refrindex - the function determining the refractive index - because of julia's
# multiple dispatch feature the correct function will be called from the one above and one below
#


    function refrindex(x,y, origin_prime= [0.0,0.0], θ=0.0)
        Rθ =[ cos(θ) -sin(θ) ; sin(θ) cos(θ)]
        Nx = length(x)
        Ny = length(y)

        n = ones(Nx, Ny)
        d1 = 3
        d2 = 3


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




    ###############################################################################
    # SET UP AXES
    ###############################################################################
    # ********IMPORTANT**********
    # You might break it if the axies are not the same length - I do have a build in which this is fixed
    # and a lot of other optimisations have been implemented, but in that version it is not possible
    # to animate the field and generate aesthetic graphs
    filepath = "/path/to/your/animations/folder/here"

    N = 75

    xmin, xmax = -5.0, 5.0
    ymin, ymax = -5.0, 5.0

    time_multiplier = 2.0 # 1.0 would stop after a wave had propagated from xmin to xmax, 2.0 lets it be reflected back

    multiplier = domain(ymin,ymax)/ domain(xmin,xmax)

    x = LinearAxis(xmin,xmax,N)
    y = LinearAxis(ymin,ymax,round(Int,(N*multiplier)))
    println("Δx is $(x.Δ), Δy is $(y.Δ)")


    Δt = minimum([x.Δ,x.Δ])/(c₀*√2) # don't mess with this, this is the longest possible step which leads to a stable sim
    Nt = round(Int, time_multiplier*(maximum(x.pts)-minimum(x.pts))/(c₀*Δt))
    tmin = 0.0 ; tmax = tmin + Δt*Nt
    t = LinearAxis(tmin, tmax, Δt)
    println("Δt is $(t.Δ)")# if you increase Δt you may destabilise the vaccum and doom us all - you have been warned


    ###############################################################################
    # FIELD DEFINITION - mess with this to your heart's desire
    ###############################################################################

    # f(x,y) is the starting function
    f(x,y) = exp(-x^2-y^2)

    # g(x,y) is the starting time derivative
    g(x,y) = 0.0

    #this is the main function, please note that the boundary_conditions variable
    #does nothing as of yet. the definition of this function is found in the FTDWave file
    @time t, background, field, source = ftd_propagate(
                                        x, y, t;
                                        boundary_conditions = "Dirichlet",
                                        initialfield = f,
                                        initialderivative = g,
                                        refrindex = refrindex,
                                        source_type = "none", # change this to "sinusoid" or "sinc" to feed in such waveforms from the minumum of x
                                        )

    # plotting response will give an animation only of the response from any target (refractive index)
    # that you program in to the region
    response = field.-background


    # this uses EZanimate, an animation utilities module that I made
    @time makie_animation3D(x,y,t, field, filepath) # Feed in axies and a spacetime to animate it, order the dims like x, y, t
