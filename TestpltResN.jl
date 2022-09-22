include("CollocationMeshRefine.jl")
include("cartPoleDynamics.jl")
include("plotResult/plotResN_MeshRefine.jl")

## Global variables

## Physical parameters for the cart-pole example (Standard units
const m₁ = 1.0      # mass of cart
const m₂ = 0.3      # mass of pole
const 𝓁 = 0.5       # pole length
const g = 9.81      # grabity acceleration

## Initial conditions
const u_max = 20.0  # maximum actuator force
const d_max = 2.0   # extent of the rail that cart travels on
const d = 1.0       # distance traveled during swing-up
const αₒ = π         # angle traveled during swing-up
const T = 2.0       # duration of swing-up

const t_s = 0.0;                        # initial time
const t_t = T;                          # terminal time
const q_s = [0.0, 0.0, 0.0, 0.0];       # initial state
const q_t = [d, αₒ, 0.0, 0.0];          # termial state
const q_lb = -[d_max, Inf, Inf, Inf];   # lower bound of state
const q_ub = [d_max, Inf, Inf, Inf];    # upper bound of state
const u_lb = -u_max;                    # lower bound of input
const u_ub = u_max;                     # upper bound of input

## plot all results in contrast
prob_alg_Group = []
for collocationMethod in ["Trapezoidal"; "Hermite-Simpson";]
    for input_continuity in [Continuous(); Discontinuous()]
        MRMethodRun = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionOne"]
        # MRMethodRun = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionAll"; "MRBisectionOne"]
        for MRMethod in MRMethodRun
            for (start,termination) in [(warmStart(), earlyTerm());]# [(coldStart(), lateTerm()); (warmStart(), earlyTerm());]
                optimizer = Ipopt.Optimizer();
                numInterval = 4;
                resiThreshold = 1e-6;
                alg = Collocation(optimizer, numInterval, input_continuity, collocationMethod)
                meshRefi = MeshRefi(resiThreshold, start, termination, MRMethod)
                prob_alg = solveMeshRefi(alg, meshRefi);
                # plotResults(prob_alg);
                prob_alg_Group = [prob_alg_Group; prob_alg]
                println(prob_alg_Group[end][2]," ", prob_alg_Group[end][5])
            end
        end
    end
end

pltResN = Array{Plots.Plot{Plots.GRBackend},1}(undef, 4)

pltResN[1] = plotResN(prob_alg_Group[1], Final())
for i in 2:4
    pltResN[1] = plotResN(prob_alg_Group[i], pltResN[1])
end
savefig(pltResN[1], "Cart_Pole_v4_plot\\pltResN\\pltResN_TRPuC.png")

pltResN[2] = plotResN(prob_alg_Group[5], Final())
for i in 6:8
    pltResN[2] = plotResN(prob_alg_Group[i], pltResN[2])
end
savefig(pltResN[2], "Cart_Pole_v4_plot\\pltResN\\pltResN_TRPuD.png")

N = 20
pltResN[3] = plotResN(prob_alg_Group[9], Final())
for i in 10:12
    pltResN[3] = plotResN(prob_alg_Group[i], pltResN[3])
end
savefig(pltResN[3], "Cart_Pole_v4_plot\\pltResN\\pltResN_HSSuC.png")

pltResN[4] = plotResN(prob_alg_Group[13], Final())
for i in 14:16
    pltResN[4] = plotResN(prob_alg_Group[i], pltResN[4])
end
savefig(pltResN[4], "Cart_Pole_v4_plot\\pltResN\\pltResN_HSSuD.png")

pltResNall = plot(
    pltResN[1], 
    pltResN[2],
    pltResN[3],
    pltResN[4],
    layout=grid(2,2),
    linewidth=2,
    size=(900, 900), 
    dpi = 1000
)

savefig(pltResNall, "Cart_Pole_v4_plot\\pltResN\\pltResNall.png")
savefig(pltResNall, "Cart_Pole_v4_plot\\pltResN\\pltResNall.svg")

pltResNmeshRefine = Array{Plots.Plot{Plots.GRBackend},1}(undef, 5)
for m in 1:5
    N = 10; j = 2*m;
    pltResNmeshRefine[m] = plotResN_MeshRefine(prob_alg_Group[j], Final())
    for i in N+j:10:40
        pltResNmeshRefine[m] = plotResN_MeshRefine(prob_alg_Group[i], pltResNmeshRefine[m])
    end
end

pltResNmeshRefineall = plot(
    pltResNmeshRefine[1], 
    pltResNmeshRefine[2],
    pltResNmeshRefine[3],
    pltResNmeshRefine[4],
    pltResNmeshRefine[5],
    layout=grid(3,2),
    linewidth=2,
    size=(1000, 880),
    dpi = 1000
)
savefig("Cart_Pole_v4_plot\\pltResNmeshRefineall_Residual")




