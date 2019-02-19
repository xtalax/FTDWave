module FTDWave
export ftd_propagate

    using LinearAlgebra
    using FFTW
    using ProgressMeter
    using PDEUtils
    using DSP
    using SignalProcessing

    using DiscreteAxis

    const c₀ = 3.0*10^8
    FFTW.set_num_threads(4)

    function source(x,y,t; type, start)
        if type == "impulse"
            println("x $(x.N) y = $(t.N) t = $(t.N)")
            out = zeros(x.N,y.N,t.N)
            out[2,:,2] .= 1.0
            return out
        elseif type == "sinc"
            k = 10.0
            f₀ = c₀/(k*x.Δ)
            wt = sinc.(f₀.*(t .- maximum(t)/5)) .* gaussian.(f₀.*(t .- maximum(t)/5)/5)

            return [xi==start ? wt[ti] : 0.0 for xi in x.i, Y in y, ti in t.i]
        elseif type == "sinusoid"
            k = 10.0
            f₀ = c₀/(k*x.Δ)
            ω₀ = 2*pi*f₀
            return [xi==start ? sin.(ω₀*T) : 0.0 for xi in x.i, Y in y, T in t]

        else
            return zeros(x.N,y.N,t.N)
        end
    end

    # this implements a kind of absorbing boundary condition that doesn't actually work that well
    function HeugensABC!(u, t, x, y)
        A⁺ = (1.0 - x.Δ/(t.Δ*c₀))
        A⁻ = 1/A⁺

        if x.Δ == y.Δ
            @views u[1,:,t] =     A⁺.*(u[1,:,t].-u[1,:,t-1])+u[1,:,t]
            @views u[end,:,t] = A⁻.*(u[end,:,t].-u[end,:,t-1])+u[end,:,t]
            @views u[:, 1,t] =    A⁺.*(u[1,:,t].-u[1,:,t-1])-u[1,:,t]
            @views u[:,end,t] = A⁻.*(u[:,end,t].-u[:,end,t-1])+u[:,end,t]
        else
            B⁺ = (1.0 - y.Δ/(t.Δ*c₀))
            B⁻ = 1/B⁺

            u[1,:,t] =     A⁺.*(u[1,:,t].-u[1,:,t-1])+u[1,:,t]
            u[end,:,t] = A⁻.*(u[end,:,t].-u[end,:,t-1])+u[end,:,t]
            u[:, 1,t] =    B⁺.*(u[1,:,t].-u[1,:,t-1])+u[1,:,t]
            u[:,end,t] = B⁻.*(u[:,end,t].-u[:,end,t-1])+u[:,end,t]
        end
    end

    # Wave Equation
    function Wave!(u, i, x, y, t, factor, δxx, δyy, s; bcs)
        @views u[:,:,i] = factor.*(∇(u[:,:,i-1], δxx, δyy) .+ s[:,:,i-1])  .+ 2.0.*u[:,:,i-1] .- u[:,:,i-2]
    end

    ###############################################################################
    # Main Start
    ###############################################################################


    function ftd_propagate(x, y, t;
                            refrindex,
                            initialfield,
                            initialderivative,
                            source_type,
                            boundary_conditions = "Dirichlet",
                            source_index = 2,
                            )


        ###############################################################################
        # Main Start
        ###############################################################################
            #y.N = Int64(floor((x[end]-x[1])/(y[end]-y[1])*x.N))



            #These are magic matricies which do the 2nd partial wrt space, they have Dirichlet baked in to them
            # they have 2 along the diagonal, and -1 one off the diagonal
            δxx = (δδ(x.N, x.Δ))
            δyy = (δδ(y.N, y.Δ))

            u = zeros(Float64,x.N,y.N,t.N) # make spacetime
            u[:,:,2] =  [initialfield(X,Y) for X in x, Y in y]
            G = [initialderivative(X,Y) for X in x, Y in y]
            n = refrindex(x.pts, y.pts, [0,0], 0)

            S= source(x,y,t,type = source_type, start = source_index)

            ###############################################################################
            # MEDIUM INITIALISATION
            ###############################################################################
            # compute the medium

            factor = ((t.Δ.*c₀./n).^2)

            factorbg = (t.Δ.*c₀).^2

            # callate the field at t = -1
            @views u[:,:,1] = u[:,:,2] .- 2.0.*t.Δ.*G .+ S[:,:,1]#setup the initial conditions
            #initialise background
            bg = copy(u)
        ###############################################################################
        # SIMULATE FIELD
        ###############################################################################
        #setup plots and progress meter
        prog = Progress(t.N+1, 1)
        # propagate and plot
        println("------- Simulating Field -------")
        for i in 3:t.N
            Wave!(u, i, x, y, t, factor, δxx, δyy, S; bcs = boundary_conditions)
            Wave!(bg, i, x, y, t, factorbg, δxx, δyy, S; bcs = boundary_conditions)
            #=for j in x.i, k in y.i
                isnan.(@view u[j,k,i]) && println("NaNs have infected the sim at timestep $(i)!")
            end=#
            next!(prog)

        end
        println("\n------------Done!----------------")
        return t, bg, u, S
    end
end
