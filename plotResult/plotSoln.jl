function plotSoln(soln::DOSolution, alg::Collocation, meshRefi::MeshRefi, solnState::SolnState = Final())
    ts = cumsum([0; soln.solnTime])[1:end-1]
    q = soln.solnState
    u = soln.solnInput
    f_q = soln.solnFunc.FuncState
    f_u = soln.solnFunc.FuncInput
    N = alg.numInterval

    nstate = size(q, 2)
    ninput = size(u, 2)

    # plot initialization
    pltSolnState = Array{Plots.Plot{Plots.GRBackend},1}(undef, nstate)
    pltSolnInput = Array{Plots.Plot{Plots.GRBackend},1}(undef, ninput)

    if (typeof(alg) <: Trapezoidal{Continuous}) || (typeof(alg) <: Trapezoidal{Discontinuous})
        [pltSolnState[k] = plot(ts, q[:, k],seriestype=:scatter, markersize=2,label = "$(typeof(solnState))") for k in 1:nstate]
    elseif (typeof(alg) <: HermiteSimpson{Continuous}) || (typeof(alg) <: HermiteSimpson{Discontinuous})
        [pltSolnState[k] = plot(ts, q[1:2:end, k],seriestype=:scatter, markersize=2,label = "$(typeof(solnState))",) for k in 1:nstate]
    else
        @error "Unexpected method and iput rule '$(alg)'"
    end

    if (typeof(alg) <: Trapezoidal{Continuous})
        # title!(pltSolnState[1], "Trapezoidal Collocation with Continuous Input")
        title!(pltSolnState[1], "TRPuC")
        pltSolnInput[1] = plot(ts, u[:, 1], seriestype=:scatter, markersize=2, label = "$(typeof(solnState))")  
    elseif (typeof(alg) <: Trapezoidal{Discontinuous})
        title!(pltSolnState[1], "TRPuD")
        # title!(pltSolnState[1], "Trapezoidal Collocation with Discontinuous Input")
        pltSolnInput[1] = plot(ts[1:N], u[:, 1], seriestype=:scatter, markersize=2,label = "$(typeof(solnState))")
    elseif (typeof(alg) <: HermiteSimpson{Continuous})
        title!(pltSolnState[1], "HSSuC")
        # title!(pltSolnState[1], "Hermite-Simpson Collocation with Continuous Input")
        pltSolnInput[1] = plot(ts, u[1:2:end, 1], seriestype=:scatter, markersize=2,label = "$(typeof(solnState))")  
    elseif (typeof(alg) <: HermiteSimpson{Discontinuous})
        title!(pltSolnState[1], "HSSuD")
        # title!(pltSolnState[1], "Hermite-Simpson Collocation with Discontinuous Input")
        pltSolnInput[1] = plot(ts[1:N], u[:, 1],seriestype=:scatter, markersize=2,label = "$(typeof(solnState))")
    else
        @error "Unexpected method and iput rule '$(alg)'"
    end

    [plot!(pltSolnState[k], f_q[i, k], extrema(ts[(i):(i+1)])..., label=nothing, linecolor="blue") for i in 1:N for k in 1:nstate]
    [plot!(pltSolnInput[1], f_u[i], extrema(ts[(i):(i+1)])..., label=nothing, linecolor="blue") for i in 1:N] # for k in 1:ninput]

    ylabel!(pltSolnState[1],"Cart position (m)")
    ylabel!(pltSolnState[2],"Pole angle (rad)")
    ylabel!(pltSolnState[3],"Cart velocity (m/s)")
    ylabel!(pltSolnState[4],"Cart anglular velocity (rad/s)")
    ylabel!(pltSolnInput[1],"Motor force (N)")
    xlabel!(pltSolnInput[1],"Time (s)")

    # pltSoln = plot(
    #     pltSolnState[1], 
    #     pltSolnState[2],
    #     pltSolnState[3],
    #     pltSolnState[4],
    #     pltSolnInput[1],
    #     layout=grid(5,1),
    #     linewidth=2,
    #     size=(700, 700),
    # )
    pltSoln = plot(
        pltSolnState[1], 
        pltSolnState[2],
        pltSolnInput[1],
        layout=grid(3, 1),
        linewidth=2,
        size=(700, 700), dpi = 1000)

    algType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\$algType Trajectory")
    # savefig("a1.png")
    display(pltSoln)

    return pltSolnState, pltSolnInput, pltSoln;
end

# the knots are not shown due to dense points in Trapezoidal Collocation final results
function plotSoln(soln::DOSolution, alg::Trapezoidal, meshRefi::MeshRefi, solnState::Final)
    ts = cumsum([0; soln.solnTime])[1:end-1]
    q = soln.solnState
    u = soln.solnInput
    f_q = soln.solnFunc.FuncState
    f_u = soln.solnFunc.FuncInput
    N = alg.numInterval

    nstate = size(q, 2)
    ninput = size(u, 2)

    # plot initialization
    pltSolnState = Array{Plots.Plot{Plots.GRBackend},1}(undef, nstate)
    pltSolnInput = Array{Plots.Plot{Plots.GRBackend},1}(undef, ninput)

    
    [pltSolnState[i] = plot() for i in 1:nstate]
    [pltSolnInput[1] = plot() for i in 1:ninput]

    if (typeof(alg) <: Trapezoidal{Continuous}) || (typeof(alg) <: Trapezoidal{Discontinuous})
        [pltSolnState[k] = plot(ts, q[:, k],seriestype=:scatter, markersize=2,label = "$(typeof(solnState))",
            ) for k in 1:nstate]
    elseif (typeof(alg) <: HermiteSimpson{Continuous}) || (typeof(alg) <: HermiteSimpson{Discontinuous})
        [pltSolnState[k] = plot(ts, q[1:2:end, k],seriestype=:scatter, markersize=2,label = "$(typeof(solnState))",
            ) for k in 1:nstate]
    else
        @error "Unexpected method and iput rule '$(alg)'"
    end

    if (typeof(alg) <: Trapezoidal{Continuous})
        title!(pltSolnState[1], "TRPuC")
        pltSolnInput[1] = plot(ts, u[:, 1], seriestype=:scatter, markersize=2,label = "$(typeof(solnState))")  
    elseif (typeof(alg) <: Trapezoidal{Discontinuous})
        title!(pltSolnState[1], "TRPuD")
        pltSolnInput[1] = plot(ts[1:N], u[:, 1], seriestype=:scatter, markersize=2,label = "$(typeof(solnState))")
    elseif (typeof(alg) <: HermiteSimpson{Continuous})
        title!(pltSolnState[1], "HSSuC")
        pltSolnInput[1] = plot(ts, u[1:2:end, 1], seriestype=:scatter, markersize=2,label = "$(typeof(solnState))")  
    elseif (typeof(alg) <: HermiteSimpson{Discontinuous})
        title!(pltSolnState[1], "HSSuD")
        pltSolnInput[1] = plot(ts[1:N], u[:, 1],seriestype=:scatter, markersize=2,label = "$(typeof(solnState))")
    else
        @error "Unexpected method and iput rule '$(alg)'"
    end
    
    if (typeof(alg) <: Trapezoidal{Continuous})
        title!(pltSolnState[1], "TRPuC")
    elseif (typeof(alg) <: Trapezoidal{Discontinuous})
        title!(pltSolnState[1], "TRPuD")
    elseif (typeof(alg) <: HermiteSimpson{Continuous})
        title!(pltSolnState[1], "HSSuC")
    elseif (typeof(alg) <: HermiteSimpson{Discontinuous})
        title!(pltSolnState[1], "HSSuD")
    else
        @error "Unexpected method and iput rule '$(alg)'"
    end

    [plot!(pltSolnState[k], f_q[i, k], extrema(ts[(i):(i+1)])..., label=nothing, linecolor="blue") for i in 1:N for k in 1:nstate]
    [plot!(pltSolnInput[1], f_u[i], extrema(ts[(i):(i+1)])..., label=nothing, linecolor="blue") for i in 1:N] # for k in 1:ninput]

    ylabel!(pltSolnState[1],"Cart position (m)")
    ylabel!(pltSolnState[2],"Pole angle (rad)")
    ylabel!(pltSolnState[3],"Cart velocity (m/s)")
    ylabel!(pltSolnState[4],"Cart anglular velocity (rad/s)")
    ylabel!(pltSolnInput[1],"Motor force (N)")
    xlabel!(pltSolnInput[1],"Time (s)")

    pltSoln = plot(
        pltSolnState[1], 
        pltSolnState[2],
        pltSolnInput[1],
        layout=grid(3, 1),
        linewidth=2,
        size=(700, 700), dpi = 1000)

    algType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\$algType Trajectory")
    display(pltSoln)

    return pltSoln;
end

# function plotSoln(soln::DOSolution, alg::Collocation, meshRefi::MeshRefi, solnState::Initial())
#     ts = cumsum([0; soln.solnTime])[1:end-1]
#     q = soln.solnState
#     u = soln.solnInput
#     f_q = soln.solnFunc.FuncState
#     f_u = soln.solnFunc.FuncInput
#     N = alg.numInterval

#     nstate = size(q, 2)
#     ninput = size(u, 2)

#     # plot initialization
#     pltSolnState = Array{Plots.Plot{Plots.GRBackend},1}(undef, nstate)
#     pltSolnInput = Array{Plots.Plot{Plots.GRBackend},1}(undef, ninput)

#     if (typeof(alg) <: Trapezoidal{Continuous}) || (typeof(alg) <: Trapezoidal{Discontinuous})
#         [pltSolnState[k] = plot(ts, q[:, k],seriestype=:scatter, markersize=2,label=nothing,) for k in 1:nstate]
#     elseif (typeof(alg) <: HermiteSimpson{Continuous}) || (typeof(alg) <: HermiteSimpson{Discontinuous})
#         [pltSolnState[k] = plot(ts, q[1:2:end, k],seriestype=:scatter, markersize=2,label=nothing,) for k in 1:nstate]
#     else
#         @error "Unexpected method and iput rule '$(alg)'"
#     end

#     if (typeof(alg) <: Trapezoidal{Continuous})
#         title!(pltSolnState[1], "TRPuC")
#         pltSolnInput[1] = plot(ts, u[:, 1], seriestype=:scatter, markersize=2,label=nothing)  
#     elseif (typeof(alg) <: Trapezoidal{Discontinuous})
#         title!(pltSolnState[1], "TRPuD")
#         pltSolnInput[1] = plot(ts[1:N], u[:, 1], seriestype=:scatter, markersize=2,label=nothing)
#     elseif (typeof(alg) <: HermiteSimpson{Continuous})
#         title!(pltSolnState[1], "HSSuC")
#         pltSolnInput[1] = plot(ts, u[1:2:end, 1], seriestype=:scatter, markersize=2,label=nothing)  
#     elseif (typeof(alg) <: HermiteSimpson{Discontinuous})
#         title!(pltSolnState[1], "HSSuD")
#         pltSolnInput[1] = plot(ts[1:N], u[:, 1],seriestype=:scatter, markersize=2,label=nothing)
#     else
#         @error "Unexpected method and iput rule '$(alg)'"
#     end

#     [plot!(pltSolnState[k], f_q[i, k], extrema(ts[(i):(i+1)])..., label=nothing, linecolor="blue") for i in 1:N for k in 1:nstate]
#     [plot!(pltSolnInput[1], f_u[i], extrema(ts[(i):(i+1)])..., label=nothing, linecolor="blue") for i in 1:N] # for k in 1:ninput]

#     ylabel!(pltSolnState[1],"Cart position (m)")
#     ylabel!(pltSolnState[2],"Cart angle (rad)")
#     ylabel!(pltSolnState[3],"Cart velocity (m/s)")
#     ylabel!(pltSolnState[4],"Cart anglular velocity (rad/s)")
#     ylabel!(pltSolnInput[1],"Motor force (N)")

#     return pltSolnState, pltSolnInput;
# end

function plotSoln(pltSolnState, pltSolnInput, soln::DOSolution, alg::Collocation, meshRefi::MeshRefi, solnState::Final)
    ts = cumsum([0; soln.solnTime])[1:end-1]
    q = soln.solnState
    u = soln.solnInput
    f_q = soln.solnFunc.FuncState
    f_u = soln.solnFunc.FuncInput
    N = alg.numInterval

    nstate = size(q, 2)
    ninput = size(u, 2)

    if (typeof(alg) <: Trapezoidal{Continuous}) || (typeof(alg) <: Trapezoidal{Discontinuous})
        [pltSolnState[k] = plot!(pltSolnState[k], ts, q[:, k], seriestype=:scatter, markersize=2, markercolor = "green",label = "$(typeof(solnState))", legend = :bottomright) for k in 1:nstate]
    elseif (typeof(alg) <: HermiteSimpson{Continuous}) || (typeof(alg) <: HermiteSimpson{Discontinuous})
        [pltSolnState[k] = plot!(pltSolnState[k], ts, q[1:2:end, k], seriestype=:scatter, markersize=2, markercolor = "green", label = "$(typeof(solnState))", legend = :bottomright) for k in 1:nstate]
    else
        @error "Unexpected method and iput rule '$(alg)'"
    end

    if (typeof(alg) <: Trapezoidal{Continuous})
        # title!(pltSolnState[1], "TRPuC")
        pltSolnInput[1] = plot!(pltSolnInput[1], ts, u[:, 1], seriestype=:scatter, markersize=2, markercolor = "green", label = "$(typeof(solnState))", legend = :bottomright)  
    elseif (typeof(alg) <: Trapezoidal{Discontinuous})
        title!(pltSolnState[1], "TRPuD")
        pltSolnInput[1] = plot!(pltSolnInput[1], ts[1:N], u[:, 1], seriestype=:scatter, markersize=2, markercolor = "green", label = "$(typeof(solnState))", legend = :bottomright)
    elseif (typeof(alg) <: HermiteSimpson{Continuous})
        title!(pltSolnState[1], "HSSuC")
        pltSolnInput[1] = plot!(pltSolnInput[1], ts, u[1:2:end, 1], seriestype=:scatter, markersize=2, markercolor = "green", label = "$(typeof(solnState))", legend = :bottomright)  
    elseif (typeof(alg) <: HermiteSimpson{Discontinuous})
        title!(pltSolnState[1], "HSSuD")
        pltSolnInput[1] = plot!(pltSolnInput[1], ts[1:N], u[:, 1],seriestype=:scatter, markersize=2, markercolor = "green", label = "$(typeof(solnState))", legend = :bottomright)
    else
        @error "Unexpected method and iput rule '$(alg)'"
    end

    [plot!(pltSolnState[k], f_q[i, k], extrema(ts[(i):(i+1)])..., label=nothing, linecolor="green") for i in 1:N for k in 1:nstate]
    [plot!(pltSolnInput[1], f_u[i], extrema(ts[(i):(i+1)])..., label=nothing, linecolor="green") for i in 1:N] # for k in 1:ninput]

    pltSoln = plot(
        pltSolnState[1], 
        pltSolnState[2],
        pltSolnInput[1],
        layout=grid(3, 1),
        linewidth=2,
        size=(700, 700), dpi = 1000)

    algType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\$algType TrajectoryContrast")
    display(pltSoln)

    return pltSoln;
end
