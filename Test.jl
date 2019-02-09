# !~/bin/julia
    pwd() ∉ LOAD_PATH && push!(LOAD_PATH, pwd())


    using EZAnimate
    using DiscreteAxis
    using FTDWave
    using SignalProcessing
    using FFTW
    using LinearAlgebra
    using DSP
    FFTW.set_num_threads(4)
    import Plots.plot

    #using Traceur
    #using Plots
    #my modules



    ##############################################################################
    # DEFINE CONSTANTS
    ##############################################################################
    filepath = homedir()*"/Videos/MathAnimations/"

    const c₀ = 3*10^8

    ##############################################################################
    # DECLARE FUNCTIONS
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

    N = 150

    xmin, xmax = -8.0, 8.0
    ymin = -5.0 ; ymax = 5.0

    #multiplier = domain(ymin,ymax)/ domain(xmin,xmax)

    x = LinearAxis(xmin,xmax,N)
    y = LinearAxis(ymin,ymax,round(Int,(N)))
    println("Step size is $(x.Δ)")

    ###############################################################################
    # FIELD DEFINITION
    ###############################################################################

    #f(x,y) = 0.0
    f(x,y) = 0.0
    g(x,y) = 0.0
    #s(x,y,t) = sin(ω₀*t)




    @time t, background, field, source = ftd_propagate(
                                        x, y;
                                        boundary_conditions = "HeugensABC",
                                        initialfield = f,
                                        initialderivative = g,
                                        refrindex = refrindex,
                                        source_type = "sinc",
                                        time_multiplier = 1.99
                                        )

    #response = field.-background

    #surface(x,y,response[:,:,end-10])

    #compile_gif(x, y, background, field, xval)
    #spectrum = fft(field)/(t.N)
    #@time makie_animation3D(x,y,t, field,filepath)
    tmp=(field-background)
    flattenedsource = sum(source[2,:,:], dims = 1)./domain(y)
    out = sum(tmp[2,:,:], dims =1)./domain(y)
    #reflection = xcorr(reshape(out, (length(out))),reshape(flattenedsource,(length(out))))
    #reflection_dB = amp2db.(reflection)
    r = t.pts.*c₀./2
    spectrum = (fft(reshape(out,length(out))))
    plot(r,[real.(spectrum),imag.(spectrum),abs.(spectrum)])
