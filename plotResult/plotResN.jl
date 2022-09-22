function plotResN(soln::DOSolution, alg::Collocation, meshRefi::MeshRefi, solnState::SolnState=Final())
    MRType = "$(typeof(meshRefi))"
    pltResN = plot(soln.solnError.numInterv, soln.solnError.errorTNIR, label = "$MRType", dpi = 1000)
    pltResN = plot!(pltResN, [4,5000], [1.0e-6,1.0e-6], linecolor="red", linestyle=:dash, label=nothing)

    yticks!(pltResN,[1e1, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8], 
    ["10", "1", "0.1", "0.01", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"])
    # plot!(pltStep, xticks = ([10,100,1000], ["10¹", "\\pi", "2\\pi"]))
    ylims!(pltResN, (1e-8,100))
    xaxis!(pltResN, xaxis=:log)
    yaxis!(pltResN, yaxis=:log)
    title!(pltResN,"Time normalized integrated residuals")
    ylabel!(pltResN,"Residuals")
    xlabel!(pltResN,"Number of intervals")

    algSolnType = "$(typeof(alg))$(typeof(solnState))$(typeof(meshRefi))"
    savefig("Cart_Pole_v4_plot\\$algSolnType Residual")
    display(pltResN)

    return pltResN;
end

function plotResN(prob_alg, solnState::SolnState=Final())
    solnInit, algInit, soln, alg, meshRefi= prob_alg;

    if typeof(meshRefi) <: MRMinimax
        MRType = "MRMinimax"
    elseif typeof(meshRefi) <: MREquidistributed
        MRType = "MREquidistributed"
    elseif typeof(meshRefi) <: MRBisectionMulti
        MRType = "MRBisectionMulti"
    elseif typeof(meshRefi) <: MRBisectionAll
        MRType = "MRBisectionAll"
    elseif typeof(meshRefi) <: MRBisectionOne
        MRType = "MRBisectionOne"
    else
        MRType = "$(typeof(meshRefi))"
    end

    # MRType = "$(typeof(meshRefi))"
    pltResN = plot(soln.solnError.numInterv, soln.solnError.errorTNIR, label = "$MRType", dpi = 1000)

    yticks!(pltResN,[1e1, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8], 
    ["10", "1", "0.1", "0.01", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"])
    ylims!(pltResN, (1e-8,50))
    xaxis!(pltResN, xaxis=:log)
    yaxis!(pltResN, yaxis=:log)
    ylabel!(pltResN,"Residuals")
    xlabel!(pltResN,"Number of intervals")

    if typeof(alg) <: Trapezoidal{Continuous}
        algType = "TRPuC"
        pltResN = plot!(pltResN, [4,300], [1.0e-6,1.0e-6], linecolor="red", linestyle=:dash, label=nothing)
    elseif typeof(alg) <: Trapezoidal{Discontinuous}
        algType = "TRPuD"
        pltResN = plot!(pltResN, [4,300], [1.0e-6,1.0e-6], linecolor="red", linestyle=:dash, label=nothing)
    elseif typeof(alg) <: HermiteSimpson{Continuous}
        algType = "HSSuC"
        pltResN = plot!(pltResN, [4,300], [1.0e-6,1.0e-6], linecolor="red", linestyle=:dash, label=nothing)
    elseif typeof(alg) <: HermiteSimpson{Discontinuous}
        algType = "HSSuD"
        pltResN = plot!(pltResN, [4,300], [1.0e-6,1.0e-6], linecolor="red", linestyle=:dash, label=nothing)
    end
    title!(pltResN,"($(algType)) Integrated Residuals")
    # title!(pltResN,"($(algType)) Time normalized integrated residuals")

    algSolnType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\pltResN\\$algSolnType Residual")
    display(pltResN)

    return pltResN;
end

function plotResN(pltResN, soln::DOSolution, alg::Collocation, meshRefi::MeshRefi, solnState::SolnState=Final())
    if typeof(meshRefi) <: MRMinimax
        MRType = "MRMinimax"
    elseif typeof(meshRefi) <: MREquidistributed
        MRType = "MREquidistributed"
    elseif typeof(meshRefi) <: MRBisectionMulti
        MRType = "MRBisectionMulti"
    elseif typeof(meshRefi) <: MRBisectionAll
        MRType = "MRBisectionAll"
    elseif typeof(meshRefi) <: MRBisectionOne
        MRType = "MRBisectionOne"
    else
        MRType = "$(typeof(meshRefi))"
    end

    pltResN = plot!(pltResN, soln.solnError.numInterv, soln.solnError.errorTNIR, label = "$MRType", dpi = 1000)

    algSolnType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\pltResN\\$algSolnType Residual")
    display(pltResN)

    return pltResN;
end

function plotResN(prob_alg, pltResN)
    _, _, solnFina, algFina, meshRefi = prob_alg;
    pltResN = plotResN(pltResN, solnFina, algFina, meshRefi)
    return pltResN
end