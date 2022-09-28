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


## Test running time
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10;
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60;
io = IOContext(stdout)

medianTime = []; meanTime = [];
for collocationMethod in ["Trapezoidal";"Hermite-Simpson"; ]
    for input_continuity in [Continuous(); Discontinuous()]
        MRMethodRun = ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionOne"]
        for MRMethod in MRMethodRun
            MRstart = [coldStart(); onceWarmStart(); quasiWarmStart(); warmStart()]
            for start in MRstart
                # MRTermRun = [lateTerm(); earlyTerm()]
                MRTermRun = [lateTerm();]
                for termination in MRTermRun
                    optimizer = Ipopt.Optimizer();
                    numInterval = 4
                    resiThreshold = 1e-6;
                    @show alg = Collocation(optimizer, numInterval, input_continuity, collocationMethod)
                    @show meshRefi = MeshRefi(resiThreshold, start, termination, MRMethod)

                    b = @benchmark (prob_alg = solveMeshRefi(alg, meshRefi)); show(io, MIME("text/plain"), b); println("\n");

                    medianTime = [medianTime; median(b.times).*1e-9]; # Change the time unit from ns to s
                    meanTime = [meanTime; mean(b.times).*1e-9]; # Change the time unit from ns to s
                    # prob_alg_Group = [prob_alg_Group; prob_alg]
                    # println("$(typeof(prob_alg[2]))\n")
                end
            end
        end
    end
end


## record running time to txt files

using DelimitedFiles

touch("plotResult\\medianTimeStart.txt")
open("plotResult\\medianTimeStart.txt", "w") do io
    writedlm(io, medianTime) 
end
touch("plotResult\\meanTimeStart.txt")
open("plotResult\\meanTimeStart.txt", "w") do io
    writedlm(io, meanTime)
end

# touch("plotResult\\RunningTimeWithType.txt")
# open("plotResult\\RunningTimeWithType.txt", "w") do io
#     i = 1
#     for collocationMethod in ["Hermite-Simpson"; "Trapezoidal"]
#         for input_continuity in [Continuous(); Discontinuous()]
#             for MRMethod in ["MRMinimax"; "MREquidistributed"; "MRBisectionMulti"; "MRBisectionAll"; "MRBisectionOne"]
#                 for start in [coldStart(); onceWarmStart(); quasiWarmStart(); warmStart()]
#                     for termination in [lateTerm()] # [lateTerm(); earlyTerm()]
#                         optimizer = Ipopt.Optimizer();
#                         numInterval = 5;
#                         resiThreshold = 1e-6;
#                         alg = Collocation(optimizer, numInterval, input_continuity, collocationMethod)
#                         meshRefi = MeshRefi(resiThreshold, start, termination, MRMethod)

#                         # write(io, "$(typeof(alg)) $(typeof(meshRefi))\n")
#                         # write(io, "meanTime\n")
#                         # writedlm(io, meanTime[i])

#                         write(io, "medianTime\n")
#                         writedlm(io, medianTime[i])
#                         write(io, "\n")

#                         i = i + 1
#                     end
#                 end
#             end
#         end
#     end
# end;

# for collocationMethod in ["Trapezoidal";"Hermite-Simpson"; ]
#     for input_continuity in [Continuous(); Discontinuous()]
#         for start in [coldStart(); onceWarmStart(); quasiWarmStart(); warmStart()]
#             for termination in [lateTerm();]# [lateTerm(); earlyTerm()]
#                 optimizer = Ipopt.Optimizer();
#                 MRMethod = "MRMinimax"
#                 numInterval = 107
#                 resiThreshold = 1e-6;
#                 @show alg = Collocation(optimizer, numInterval, input_continuity, collocationMethod)
#                 # @show alg
#                 @show meshRefi = MeshRefi(resiThreshold, start, termination, MRMethod)
#                 # @show meshRefi

#                 ## Dynamic optimization problem defined
#                 time  = JuDOTime(t_s, t_t);               # DOVariable(t_s, t_t)
#                 state = DOVariables(q_s, q_s, q_t, q_t, q_lb, q_ub);
#                 input = DOVariables(u_lb, u_ub);
#                 guess = initial_guess(time, state, input, alg);
#                 funcDynamics = CartPole!;                  # system dynamics

#                 prob = DOProblem(time, state, input, funcDynamics, guess);
                
#                 b = @benchmark (prob_alg = solve(prob, alg, meshRefi)); show(io, MIME("text/plain"), b); println("\n");
#                 # @time prob_alg = solveMeshRefi(alg, meshRefi); # plotResults(prob_alg);
#                 medianTime = [medianTime; median(b.times)];
#                 meanTime = [meanTime; mean(b.times)]
#                 prob_alg_Group = [prob_alg_Group; prob_alg]
#                 # println("$(typeof(prob_alg[2]))\n")
#             end
#         end
#     end
# end

