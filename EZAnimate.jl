module EZAnimate
export makie_animation3D, plots_gif3D
###########################################################
###########################################################
    using Makie
    using ProgressMeter
    # it works, dw makie is confusing
    function makie_animation3D(x,y,t,u, filepath)
        prog = Progress(t.N, 1)
        println("---------Animating Data-----------")
        scene = Scene(resolution = (500, 500));
        limits = FRect3D(Vec3f0(x.pts[1], y.pts[1], x.pts[1]) , Vec3f0(x.pts[end]-x.pts[1], y.pts[end]-y.pts[1], y.pts[end]-y.pts[1]))
        surf = Makie.surface!(scene, x.pts, y.pts, u[:,:,end], limits = limits)[end]
        scene
        record(scene, filepath*"EM_animation_Makie_Complex$(t.N)$([2,3,50]).mp4", 1:t.N) do i
            surf[3] = u[:,:,i]
            next!(prog)
        end
    end

    using Plots

    function plots_gif3D(x,y,u)

        Nx, Ny, t.N = size(u)
        xmin = minimum(x) ;ymin = minimum(y)
        xmax = maximum(x) ;ymax = maximum(y)
        println("---------Compiling Animation-----------")
        pyplot()
        prog = Progress(t.N, 1)
        global anim = Animation()
        #frames = Dict()
        for i in 1:t.N
            cameraAngle = 2*pi*i/t.N
            plot = surface(x,y, u[:,:,i], xlims=(xmin,xmax), ylims=(ymin,ymax) , zlims=(-0.5,0.5))

            surface!(plot, camera=(15*cos(cameraAngle), 30))
            #push!(a, i=>plot)
            frame(anim,plot)
            next!(prog)

        end
        #=@sync for key in sort(collect(keys(frames)))
            frame(anim,frames[key])
        end=#

        gif(anim, filepath*"animation$t.N.gif", fps = 30)
    end

    function plotsgif3D(x, y, bg, field, source, xval)
        pyplot(leg=false, ticks=nothing)
        Nx, Ny, t.N = size(field)
        prog = Progress(t.N-1, 1)
        global zs = zeros(0)
        global anim = Animation()

        anim = @animate for i in 1:(t.N-1)

            if bg[:,:,i] ≈ field[:,:,i]
                continue
                prog = Progress(t.N-1-i, 1)
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
        gif(anim, filepath*"animation$t.N.gif", fps = 30)

    end

end
