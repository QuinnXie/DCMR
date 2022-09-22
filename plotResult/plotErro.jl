function plotErro(soln::DOSolution, alg::Collocation, meshRefi::MeshRefi, solnState::SolnState=Final())
    nsegment = 16
    ts = soln.solnError.errorTime
    errorQuad = soln.solnError.errorQuad
    errorFunc = soln.solnError.errorFunc
    
    N = alg.numInterval
    nstate = size(errorQuad, 2)
    # tnsegment = LinRange(0.0, ts[end], nsegment*(length(ts)-1)+1)
    tnsegment = mapreduce(x->x, vcat, (ts[1]:(ts[2]-ts[1])/nsegment:ts[2])[1:end-1])
    for i in 2:(size(ts,1)-1)
        append!(tnsegment, mapreduce(x->x, vcat, (ts[i]:(ts[i+1]-ts[i])/nsegment:ts[i+1])[1:end-1]))
    end
    append!(tnsegment,ts[end])
    errorDynamics = zeros(length(tnsegment), nstate);

    # errorDynamic
    for i in 1:N
        error = errorFunc[i](tnsegment[(nsegment*(i-1)+1):(nsegment*i)])
        for j in 1:nsegment
            errorDynamics[nsegment*(i-1)+j,:] = error[j][:]       
        end        
    end

    pltErrorDynamics = Array{Plots.Plot{Plots.GRBackend},1}(undef, nstate)
    pltErrorQuad = Array{Plots.Plot{Plots.GRBackend},1}(undef, nstate)

    [pltErrorDynamics[k] = plot(
        tnsegment[1:nsegment:end], errorDynamics[1:nsegment:end, k], 
        seriestype=:scatter, markersize=2, 
        label=nothing, legend = :bottomright, 
        ) for k in 1:nstate]
    [plot!(pltErrorDynamics[k], tnsegment, errorDynamics[:, k], label=nothing) for k in 1:nstate]

    plot!(pltErrorDynamics[1], title = "Residuals: dx/dt - f(t,x,u)")
    ylabel!(pltErrorDynamics[1], "position dynamics error (m/s)")
    ylabel!(pltErrorDynamics[2], "angle dynamics error (rad/s)")
    xlabel!(pltErrorDynamics[2], "time (s)")

    [pltErrorQuad[k] = plot(ts[1:end],([errorQuad[1, k]; errorQuad[:, k]]), label=nothing, linetype=:steppre, yaxis=:log) for k in 1:nstate]

    title!(pltErrorQuad[1],"Absolute local error")
    ylabel!(pltErrorQuad[1],"position error")
    ylabel!(pltErrorQuad[2],"pole angle")
    xlabel!(pltErrorQuad[2],"time (s)")
    
    pltErro = plot(
        pltErrorDynamics[1],
        pltErrorQuad[1],
        pltErrorDynamics[2],
        pltErrorQuad[2],
        layout=grid(2, 2),
        linewidth=2,
        size=(700, 700),
        dpi = 1000
    )
    # pltError = plot(
    #     pltErrorDynamics[1],
    #     pltErrorQuad[1],
    #     pltErrorDynamics[2],
    #     pltErrorQuad[2],
    #     pltErrorDynamics[3],
    #     pltErrorQuad[3],
    #     pltErrorDynamics[4],
    #     pltErrorQuad[4],
    #     layout=grid(2, 4),
    #     linewidth=2,
    #     size=(700, 700),
    # )

    algType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\$algType Error")
    savefig("Cart_Pole_v4_plot\\pltError\\$algType Error.svg")
    display(pltErro)

    return pltErro;
end