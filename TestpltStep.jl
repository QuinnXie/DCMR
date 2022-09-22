include("CollocationMeshRefine.jl")
include("cartPoleDynamics.jl")

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
global prob_alg_Group = []
for collocationMethod in ["Trapezoidal", "Hermite-Simpson";]
    for input_continuity in [Continuous(); Discontinuous()]
        for MRMethod in ["MREquidistributed"; "MRMinimax"; "MRBisectionMulti"; "MRBisectionAll"; "MRBisectionOne"]
            for (start,termination) in [(coldStart(), lateTerm()); (warmStart(), earlyTerm()); ]
                optimizer = Ipopt.Optimizer();
                numInterval = 5;
                resiThreshold = 1e-6;
                alg = Collocation(optimizer, numInterval, input_continuity, collocationMethod)
                meshRefi = MeshRefi(resiThreshold, start, termination, MRMethod)
                prob_alg = solveMeshRefi(alg, meshRefi);
                prob_alg_Group = [prob_alg_Group; prob_alg]
                println(prob_alg_Group[end][2]," ", prob_alg_Group[end][5])
            end
        end
    end
end

# plotResults(prob_alg_Group[1]);

pltStepContrast = Array{Plots.Plot{Plots.GRBackend},1}(undef, 4)

pltStepContrast[1] = plotStep(prob_alg_Group[2], Initial(), dpi = 1000)
for i in 2:2:8
    println(typeof(pltStepContrast[1]))
    pltStepContrast[1] = plotStep(prob_alg_Group[i], pltStepContrast[1], Final())
end
savefig(pltStepContrast[1], "Cart_Pole_v4_plot\\pltStep\\pltStepContrast_TRPuC.png")

pltStepContrast[2] = plotStep(prob_alg_Group[10], Initial(), dpi = 1000)
for i in 10:2:16
    pltStepContrast[2] = plotStep(prob_alg_Group[i], pltStepContrast[2], Final())
end
savefig(pltStepContrast[2], "Cart_Pole_v4_plot\\pltStep\\pltStepContrast_TRPuD.png")

N = 16
pltStepContrast[3] = plotStep(prob_alg_Group[N+2], Initial(), dpi = 1000)
for i in N+2:2:N+10
    pltStepContrast[3] = plotStep(prob_alg_Group[i], pltStepContrast[3], Final())
end
savefig(pltStepContrast[3], "Cart_Pole_v4_plot\\pltStep\\pltStepContrast_HSCuC.png")

pltStepContrast[4] = plotStep(prob_alg_Group[N+12], Initial(), dpi = 1000)
for i in N+12:2:N+20
    pltStepContrast[4] = plotStep(prob_alg_Group[i], pltStepContrast[4], Final())
end
savefig(pltStepContrast[4], "Cart_Pole_v4_plot\\pltStep\\pltStepContrast_HSCuD.png")

pltStepContrastall = plot(
    pltStepContrast[1], 
    pltStepContrast[2],
    pltStepContrast[3],
    pltStepContrast[4],
    layout=grid(2,2),
    linewidth=2,
    size=(800, 700),
    dpi = 1000
)

display(pltStepContrastall)
savefig(pltStepContrastall, "Cart_Pole_v4_plot\\pltStep\\pltStepContrastall.png")





