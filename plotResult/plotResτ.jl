function plotResτ(soln::DOSolution, alg::Collocation, meshRefi::MeshRefi, solnState::SolnState=Final())
    Resi = soln.solnError.errorSqua# errorMean
    # Resi = soln.solnError.errorMean
    # τ = (1:size(Resi,1))./size(Resi,1);
    # τ = [0.0; τ]; Resi = [Resi[1]; Resi];

    # τ = (1:size(soln.solnTime,1));
    τ₁ = cumsum(soln.solnTime[1:end-1]); 
    τ₂ = [0.0; τ₁[1:end-1]];
    τ = (τ₁ + τ₂)./2;
    Resi = [Resi[1]; Resi[1:end-1];];

    pltResτ = plot(τ, Resi, label = "$(typeof(solnState))", seriestype=:scatter, dpi = 1000)
    # linetype=:steppre
    if typeof(alg) <: Trapezoidal{Continuous}
        algType = "TRPuC"
        title!(pltResτ, "(TRPuC) Local Integrated Residuals")  
    elseif typeof(alg) <: Trapezoidal{Discontinuous}
        algType = "TRPuD"
        title!(pltResτ, "(TRPuD) Local Integrated Residuals")  
    elseif typeof(alg) <: HermiteSimpson{Continuous}
        algType = "HSSuC"
        title!(pltResτ, "(HSSuC) Local Integrated Residuals")  
    elseif typeof(alg) <: HermiteSimpson{Discontinuous}
        algType = "HSSuD"
        title!(pltResτ, "(HSSuD) Local Integrated Residuals")  
    end  

    # title!(pltResτ,"Local Integrated Residuals")

    xlabel!(pltResτ,"Whole time scale (s)")
    xlims!(pltResτ, (-0.05,2.05))
    xticks!(pltResτ,[0.0:0.25:2.0;], 
    ["0", "0.25", "0.5", "0.75", "1.0", "1.25", "1.50", "1.75", "2.0"])
    ylabel!(pltResτ,"Residuals")
    yaxis!(pltResτ, yaxis=:log)
    # ylims!(pltResτ, (1e-13,1e-8))
    # yticks!(pltResτ,[1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13], 
    # ["1e-8", "1e-9", "1e-10", "1e-11", "1e-12", "1e-13"])

    algType = "$(typeof(algType))$(typeof(solnState))"
    savefig("Cart_Pole_v4_plot\\pltResτ\\$algType ResiDistri")
    display(pltResτ)

    return pltResτ;
end

function plotResτ(pltResτ, soln::DOSolution, alg::Collocation, meshRefi::MeshRefi)
    Resi = soln.solnError.errorSqua # errorSqua
    # τ = (1:size(Resi,1))./size(Resi,1);
    # τ = [0.0; τ]; Resi = [Resi[1]; Resi];

    # τ = (1:size(soln.solnTime,1));
    τ₁ = cumsum(soln.solnTime[1:end-1]); 
    τ₂ = [0.0; τ₁[1:end-1]];
    τ = (τ₁ + τ₂)./2;
    Resi = [Resi[1]; Resi[1:end-1];];

    pltResτ = plot!(pltResτ, τ, Resi, label = "Final", seriestype=:scatter, dpi = 1000)


    if typeof(alg) <: Trapezoidal{Continuous}
        algType = "TRPuC"
        title!(pltResτ, "(TRPuC) Local Integrated Residuals")  
    elseif typeof(alg) <: Trapezoidal{Discontinuous}
        algType = "TRPuD"
        title!(pltResτ, "(TRPuD) Local Integrated Residuals")  
    elseif typeof(alg) <: HermiteSimpson{Continuous}
        algType = "HSSuC"
        title!(pltResτ, "(HSSuC) Local Integrated Residuals")  
    elseif typeof(alg) <: HermiteSimpson{Discontinuous}
        algType = "HSSuD"
        title!(pltResτ, "(HSSuD) Local Integrated Residuals")  
    end  

    yaxis!(pltResτ, yaxis=:log)
    ylims!(pltResτ, (1e-12,1e4))
    yticks!(pltResτ,[1e2, 1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12], 
    ["100", "1", "1e-2", "1e-4", "1e-6", "1e-8", "1e-10", "1e-12"])
    # ["10", "0.1", "1e-3", "1e-5", "1e-7", "1e-9", "1e-11", "1e-13"]
    algType = "$(typeof(algType))"
    savefig("Cart_Pole_v4_plot\\pltResτ\\$algType ResiDistComp")
    display(pltResτ)
    
    return pltResτ;
end