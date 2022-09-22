include("plotSoln.jl")
include("plotErro.jl")
include("plotResN.jl")
include("plotResτ.jl")
include("plotStep.jl")

function plotResults(soln, alg, meshRefi::MeshRefi, solnState::SolnState = Final())
    pltSoln = plotSoln(soln, alg, meshRefi, solnState); # plot state trajectory
    pltErro = plotErro(soln, alg, meshRefi, solnState); # plot absolute local error
    pltResN = plotResN(soln, alg, meshRefi, solnState); # plot integrated residuals with number of intervals
    pltResτ = plotResτ(soln, alg, meshRefi, solnState); # plot integrated residuals distribution
    pltStep = plotStep(soln, alg, meshRefi, solnState); # plot step size distribution
end

function plotResults(solnInit, algInit, solnFina, algFina, meshRefi::MeshRefi)
    pltResτ = plotResτ(solnInit, algInit, meshRefi, Initial());
    pltResτ = plotResτ(pltResτ, solnFina, algFina, meshRefi);

    pltStep = plotStep(solnInit, algInit, meshRefi, Initial());
    pltStep = plotStep(pltStep, solnFina, algFina, meshRefi);

    pltSolnState, pltSolnInput = plotSoln(solnInit, algInit, meshRefi, Initial())
    pltSoln = plotSoln(pltSolnState, pltSolnInput, solnFina, algFina, meshRefi, Final());
    # eps(pltSoln, "pltSoln")
    return pltResτ, pltStep, pltSoln;
end

function plotResults(prob_alg)
    solnInit, algInit,solnFina, algFina, meshRefi = prob_alg;
    plotResults(solnFina, algFina, meshRefi, Final()); plotResults(solnInit, algInit, meshRefi, Initial());
    pltResτ, pltStep, pltSoln = plotResults(solnInit, algInit, solnFina, algFina, meshRefi);
    return pltResτ, pltStep, pltSoln;
end
