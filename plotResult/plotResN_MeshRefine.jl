function plotResN_MeshRefine(prob_alg, solnState::SolnState=Final())
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

    if typeof(alg) <: Trapezoidal{Continuous}
        algType = "TRPuC"
     elseif typeof(alg) <: Trapezoidal{Discontinuous}
        algType = "TRPuD"
     elseif typeof(alg) <: HermiteSimpson{Continuous}
        algType = "HSSuC"
    elseif typeof(alg) <: HermiteSimpson{Discontinuous}
        algType = "HSSuD"
    end

    # MRType = "$(typeof(meshRefi))"
    pltResN = plot(soln.solnError.numInterv, soln.solnError.errorTNIR, label = "$(algType)", dpi = 1000)

    yticks!(pltResN,[1e1, 1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8], 
    ["10", "1", "0.1", "0.01", "1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8"])
    ylims!(pltResN, (1e-8,50))
    xaxis!(pltResN, xaxis=:log)
    yaxis!(pltResN, yaxis=:log)
    ylabel!(pltResN,"Residuals")
    xlabel!(pltResN,"Number of intervals")

    pltResN = plot!(pltResN, [4,300], [1.0e-6,1.0e-6], linecolor="red", linestyle=:dash, label=nothing)
    
    title!(pltResN,"Integrated Residuals ($MRType)")
    # title!(pltResN,"($(algType)) Time normalized integrated residuals")

    algSolnType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\pltResN\\$algSolnType Residual")
    display(pltResN)

    return pltResN;
end

function plotResN_MeshRefine(prob_alg, pltResN)
    _, _, soln, alg, meshRefi = prob_alg;
    solnState = Final()

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

    # if typeof(alg) <: Trapezoidal{Continuous}
    #     algType = "Trapezoidal_Continuous_input"
    #  elseif typeof(alg) <: Trapezoidal{Discontinuous}
    #     algType = "Trapezoidal_Discontinuous_input"
    #  elseif typeof(alg) <: HermiteSimpson{Continuous}
    #     algType = "Hermite-Simpson_Continuous_input"
    # elseif typeof(alg) <: HermiteSimpson{Discontinuous}
    #     algType = "Hermite-Simpson_Discontinuous_input"
    # end

    if typeof(alg) <: Trapezoidal{Continuous}
        algType = "TRPuC"
     elseif typeof(alg) <: Trapezoidal{Discontinuous}
        algType = "TRPuD"
     elseif typeof(alg) <: HermiteSimpson{Continuous}
        algType = "HSSuC"
    elseif typeof(alg) <: HermiteSimpson{Discontinuous}
        algType = "HSSuD"
    end

    pltResN = plot!(pltResN, soln.solnError.numInterv, soln.solnError.errorTNIR, label = "$(algType)", dpi = 1000)

    algSolnType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\pltResN\\$algSolnType Residual")
    display(pltResN)

    return pltResN;
end
