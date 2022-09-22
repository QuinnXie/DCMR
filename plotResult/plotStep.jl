function plotStep(soln::DOSolution, alg::Collocation, meshRefi::MeshRefi, solnState::SolnState=Final())
    step = soln.solnTime
    # τ = (1:size(step,1))./size(step,1);
    # τ = (1:size(soln.solnTime,1));
    # τ = [0.0; τ]; 
    τ₁ = cumsum(soln.solnTime); τ₂ = [0.0; τ₁[1:end-1]];
    τ = (τ₁ + τ₂)./2;
 
    step = [step[1]; step[1:end-1]; step[end-1]];
    # pltStep = plot(τ, step, label = "$(typeof(solnState))", linetype=:steppre)
    pltStep = plot(τ, step, label = "$(typeof(solnState))", seriestype=:scatter, dpi =1000)
   
    if typeof(alg) <: Trapezoidal{Continuous}
        algType = "TRPuC"
        title!(pltStep, "(TRPuC) Time step size")  
    elseif typeof(alg) <: Trapezoidal{Discontinuous}
        algType = "TRPuD"
        title!(pltStep, "(TRPuD) Time step size")  
    elseif typeof(alg) <: HermiteSimpson{Continuous}
        algType = "HSSuC"
        title!(pltStep, "(HSSuC) Time step size")  
    elseif typeof(alg) <: HermiteSimpson{Discontinuous}
        algType = "HSSuD"
        title!(pltStep, "(HSSuD) Time step size")  
    end  

    xlabel!(pltStep,"Whole time scale (s)")
    xlims!(pltStep, (-0.05,2.05))
    xticks!(pltStep,[0.0:0.25:2.0;], 
    ["0", "0.25", "0.5", "0.75", "1.0", "1.25", "1.50", "1.75", "2.0"])

    ylabel!(pltStep,"Length of time in each interval (s)")
    # ylims!(pltStep, (1e-4,1))
    # yaxis!(pltStep, yaxis=:log)
    # yticks!(pltStep,[1, 1e-1, 1e-2, 1e-3, 1e-4], 
    # ["1", "0.1", "0.01", "1e-3", "1e-4"])
    
    algType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\$algType StepSize")
    display(pltStep)

    return pltStep;
end

function plotStep(pltStep, soln::DOSolution, alg::Collocation, meshRefi::MeshRefi)
    step = soln.solnTime
    # τ = (1:size(step,1))./size(step,1);
    # τ = (1:size(soln.solnTime,1));
    # τ = [0.0; τ]; 
    τ₁ = cumsum(soln.solnTime); τ₂ = [0.0; τ₁[1:end-1]];
    τ = (τ₁ + τ₂)./2;
 
    step = [step[1]; step[1:end-1]; step[end-1]];
    # pltStep = plot!(pltStep, τ, step, label = "Final", legend = :bottomright, linetype=:steppre) seriestype=:scatter
    pltStep = plot!(pltStep, τ, step, label = "Final", legend = :bottomright, seriestype=:scatter, dpi =1000) 

    ylabel!(pltStep,"Length of time in each interval (s)")
    ylims!(pltStep, (1e-4,1.1))
    yaxis!(pltStep, yaxis=:log)
    yticks!(pltStep,[1, 1e-1, 1e-2, 1e-3, 1e-4], 
    ["1", "0.1", "0.01", "1e-3", "1e-4"])
    
    algType = "$(typeof(alg))"
    savefig("Cart_Pole_v4_plot\\$algType StepSizeComp")

    display(pltStep)

    return pltStep
end


function plotStep(prob_alg, solnState::SolnState=Final())
    if solnState == Final()
        solnInit, algInit, soln, alg, meshRefi= prob_alg;
    elseif solnState == Initial()
        soln, alg, solnFina, algFina, meshRefi= prob_alg;
    end
    
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

    step = soln.solnTime
    # τ = (1:size(step,1))./size(step,1);
    # τ = (1:size(soln.solnTime,1));
    # τ = [0.0; τ]; 
    τ₁ = cumsum(soln.solnTime); τ₂ = [0.0; τ₁[1:end-1]];
    τ = (τ₁ + τ₂)./2;

    step = [step[1]; step[1:end-1]; step[end-1]];
    # pltStep = plot(τ, step, label = "Initial", linetype=:steppre)
    pltStep = plot(τ, step, label = "Initial", seriestype=:scatter, dpi =1000)

    xlabel!(pltStep,"τ")
    xlims!(pltStep, (-0.05,2.05))
    xticks!(pltStep,[0.0:0.25:2.0;], 
    ["0", "0.25", "0.5", "0.75", "1.0", "1.25", "1.50", "1.75", "2.0"])

    ylabel!(pltStep,"Length of time in each interval (s)")
    # ylims!(pltStep, (1e-4,1))
    # yaxis!(pltStep, yaxis=:log)
    # yticks!(pltStep,[1, 1e-1, 1e-2, 1e-3, 1e-4], 
    # ["1", "0.1", "0.01", "1e-3", "1e-4"])
    
    if typeof(alg) <: Trapezoidal{Continuous}
        algType = "TRPuC"
    elseif typeof(alg) <: Trapezoidal{Discontinuous}
        algType = "TRPuD"
    elseif typeof(alg) <: HermiteSimpson{Continuous}
        algType = "HSSuC"
    elseif typeof(alg) <: HermiteSimpson{Discontinuous}
        algType = "HSSuD"
    end

    title!(pltStep,"($(algType)) Time step size")   
    algType = "$(typeof(alg))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\pltStep\\$algType StepSize")
    display(pltStep)

    return pltStep;
end

function plotStep(prob_alg, pltStep::Plots.Plot{Plots.GRBackend}, solnState::SolnState=Final())
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

    step = soln.solnTime
    # τ = (1:size(step,1))./size(step,1);
    # τ = (1:size(soln.solnTime,1));
    # τ = [0.0; τ]; 
    τ₁ = cumsum(soln.solnTime); τ₂ = [0.0; τ₁[1:end-1]];
    τ = (τ₁ + τ₂)./2;

    step = [step[1]; step[1:end-1]; step[end-1]];
    if typeof(alg) <: Trapezoidal
        # pltStep = plot!(pltStep, τ, step, label = "$MRType", legend = :topright, linetype=:steppre)
        pltStep = plot!(pltStep, τ, step, label = "$MRType", legend = :topright, seriestype=:scatter, dpi =1000)
    elseif typeof(alg) <: HermiteSimpson
        # pltStep = plot!(pltStep, τ, step, label = "$MRType", legend = :bottomright, linetype=:steppre)
        pltStep = plot!(pltStep, τ, step, label = "$MRType", legend = :bottomright, seriestype=:scatter, dpi =1000)
    end
    # ply = plot(τ, step, label = "$MRType")

    ylabel!(pltStep,"Length of time in each interval (s)")
    ylims!(pltStep, (1e-4,1.1))
    yaxis!(pltStep, yaxis=:log)
    yticks!(pltStep,[1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5], 
    ["1", "0.1", "0.01", "1e-3", "1e-4", "1e-5"])
    
    algType = "$(typeof(alg))"
    savefig("Cart_Pole_v4_plot\\$algType StepSizeComp")

    display(pltStep)

    return pltStep
end