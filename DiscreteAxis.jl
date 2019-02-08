module DiscreteAxis
export LinearAxis, LogAxis, FuncAxis, CoordinateSystem, DAxis

    using StaticArrays
    import Base.iterate
    import Base.convert
    import Base.length
    abstract type DAxis end

    struct LinearAxis <: DAxis
        pts::AbstractVector{Number}
        Δ::Number
        N::Int
        i
        function LinearAxis(min,max,N::Int)
            points = range(min, stop=max, length=N)
            Δ = abs(points[1]-points[2])

            new(points,Δ,N,1:N)
        end
        function LinearAxis(min,max,Δ::Float64)
            points = min:Δ:max
            N = length(points)
            new(points,Δ,N,1:N)
        end
    end
    iterate(x::LinearAxis, n::Int) = iterate(x.pts,n)
    iterate(x::LinearAxis) = iterate(x.pts)
    length(x::LinearAxis) = x.N
#=
    struct LogAxis <: DAxis
        pts::SVector{Number}
        Δ::SVector{Number}
        N::Int
        base::String
        function LogAxis(min, max, N::Int, base::String)
            pts = logspace(min,max,N)
            Δ = SVector[points[i+1]-points[i] for i in 1:(N-1)]
            new(points,Δ,N,base)
        end
    end

    struct FuncAxis <: Axis
        pts::SVector{Number}
        Δ::SVector{Number}
        N::Int
        func::function
        inverse::function

        function FuncAxis(min, max, N::Int, func, inverse)
            SVector(func.(range(inverse(min),stop=inverse(max),N)))
            Δ = SVector[points[i+1]-points[i] for i in 1:(N-1)]
            new(points,Δ,N,base)
        end
    end

    function convert(T, x::DAxis)
        convert(T, x.pts)
    end

    struct CoordinateSystem
        axes::SVector{DAxis}
        metric::SVector{Char}
        CoordinateSystem(axes, metric) = new(axes,metric)
    end
    =#
end
