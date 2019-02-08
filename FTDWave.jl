module FTDWave
export ftd_propagate

    using LinearAlgebra
    using ProgressMeter
    using PDEUtils
    using SignalProcessing
    using FFTW
    using DiscreteAxis
    const c₀ = 3*10^8
    FFTW.set_num_threads(4)

    function source(x,y,t; type)
        if type == "impulse"
            out = zeros(Complex{Float64}, length(t))
            out[2] = 1
        elseif type == "sinc"
            k = 4
            f₀ = c₀/(k*x.Δ)
            ω₀ = 2*pi*f₀
            fs = 1/t.Δ
            fgrid = fs*(0:(t.N))/(t.N)
            w = rect.(fgrid/(2*f₀))
            wt= circshift(ifft(w)*(t.N), round(Int,t.N/4))
            out = unitize(wt)
            return wt
        else
            return sin.(ω₀*t.pts)
        end
    end

    # Wave Equation
    function HeugensABC!(u, t, x, y)
        A⁺ = (1 - x.Δ/(t.Δ*c₀))
        A⁻ = 1/A⁺

        if x.Δ == y.Δ
            @views u[1,:,t] =     A⁺.*(u[1,:,t].-u[1,:,t-1])+u[1,:,t]
            @views u[end,:,t] = A⁻.*(u[end,:,t].-u[end,:,t-1])+u[end,:,t]
            @views u[:, 1,t] =    A⁺.*(u[1,:,t].-u[1,:,t-1])-u[1,:,t]
            @views u[:,end,t] = A⁻.*(u[:,end,t].-u[:,end,t-1])+u[:,end,t]
        else
            B⁺ = (1 - y.Δ/(t.Δ*c₀))
            B⁻ = 1/B⁺

            u[1,:,t] =     A⁺.*(u[1,:,t].-u[1,:,t-1])+u[1,:,t]
            u[end,:,t] = A⁻.*(u[end,:,t].-u[end,:,t-1])+u[end,:,t]
            u[:, 1,t] =    B⁺.*(u[1,:,t].-u[1,:,t-1])+u[1,:,t]
            u[:,end,t] = B⁻.*(u[:,end,t].-u[:,end,t-1])+u[:,end,t]
        end
    end

    # Wave Equation
    function Wave!(u, i, x, y, t, factor, δxx, δyy, s; bcs)
        @views u[:,:,i] = factor.*(∇(u[:,:,i-1], δxx, δyy) .+ s[i-1]) .+ 2.0.*u[:,:,i-1] .- u[:,:,i-2]
        #if boundary_conditions == "HeugensABC"
            #HeugensABC!(u, i, x.Δ, y.Δ, t.Δ)
        #end
    end

    ###############################################################################
    # Main Start
    ###############################################################################


    function ftd_propagate(x, y,;
                            refrindex,
                            initialfield,
                            initialderivative,
                            source_type,
                            boundary_conditions = "Dirichlet" )


        ###############################################################################
        # Main Start
        ###############################################################################
            #y.N = Int64(floor((x[end]-x[1])/(y[end]-y[1])*x.N))

            Δt = x.Δ/(c₀*2)
            Nt = round(Int, 2*(maximum(x.pts)-minimum(x.pts))/(c₀*Δt))
            tmin = 0 ; tmax = tmin + Δt*Nt
            t = LinearAxis(tmin, tmax, Δt)



            δxx = Complex.(δδ(x.N, x.Δ))
            δyy = Complex.(δδ(y.N, y.Δ))

            u = zeros(Complex{Float64},x.N,y.N,t.N)
            u[:,:,2] =  [initialfield(X,Y) for X in x, Y in y]
            G = [initialderivative(X,Y) for X in x, Y in y]
            n =  [refrindex(X,Y) for X in x, Y in y]

            S= source(x,y,t,type = source_type)


        #    S = [xi == 2 ? sin(ω₀*t[ti]) : 0.0 for xi in x.i, yi in y.i, ti in 1:(t.N-1)]
            ###############################################################################
            # MEDIUM INITIALISATION
            ###############################################################################
            # compute the medium

            factor = Complex.(((t.Δ.^2).*(c₀./n).^2))

            factorbg = (t.Δ.*c₀).^2

            # callate the field at t = -1
            @views u[:,:,1] = u[:,:,2] .- 2*t.Δ*G .+ S[1]#setup the initial conditions
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

            next!(prog)
        end
        println("Done!")
        return t, bg, u, S
    end
end
