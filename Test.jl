

    function pwd_path!(LOAD_PATH)
        for path in LOAD_PATH
            if path == pwd()
                return LOAD_PATH
            end
        end
        push!(LOAD_PATH, pwd())
    end

    pwd_path!(LOAD_PATH)


    using EZAnimate
    using DiscreteAxis
    using FTDWave
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
    function n(x,y)
        n = 1
        d1 = 3
        d2 = 3
            if -d1/2  <= x <= d1/2
                if -d2/2  <= y <= d2/2
                    n = sqrt(2.4)
                end
            end

        return n
    end



    ###############################################################################
    # PROPAGATE
    ###############################################################################

    N = 150

    xmin, xmax = -5, 5
    ymin = xmin ; ymax = xmax
    x = LinearAxis(xmin,xmax,N)
    y = LinearAxis(ymin,ymax,N)
    println("Step size is $(x.Δ)")

    ###############################################################################
    # FIELD DEFINITION
    ###############################################################################

    f(x,y) = 0
    g(x,y) = 0
    #s(x,y,t) = sin(ω₀*t)




    @time t, background, field, source = ftd_propagate(
                                        x, y;
                                        boundary_conditions = "HeugensABC",
                                        initialfield = f,
                                        initialderivative = g,
                                        refrindex = n,
                                        source_type = "sinc"
                                        )

    #response = field.-background

    #surface(x,y,response[:,:,end-10])

    #compile_gif(x, y, background, field, xval)
    @time makie_animation3D(x,y,t,abs.(field), filepath)
