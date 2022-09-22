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

## define property of collocation method
optimizer = Ipopt.Optimizer(); # NLP solver

# collocation methods
collocationMethod = "Trapezoidal";
collocationMethod = "Hermite-Simpson";

# continuity of inputs
input_continuity = Continuous();
input_continuity = Discontinuous();

# Minimum N is selected for probelm feasibility
# TRP: Minimum numInterval = 4
# HSS: Minimum numInterval = 3
numInterval = 4

## define property of meshRefinement method
# target Mean Integrated Residual Norm Squared 
resiThreshold = 1e-6;

# initial start setting
start = coldStart(); 
start = warmStart(); 
start = quasiWarmStart();
start = onceWarmStart();

# termination setting
termination = lateTerm();
termination = earlyTerm();

# mesh refinement methods
MRMethod = "MRBisectionOne";
MRMethod = "MRMinimax";
MRMethod = "MRBisectionMulti";
MRMethod = "MREquidistributed";
# MRMethod = "MRBisectionAll";

## define collocation algorithm applied
alg = Collocation(optimizer, numInterval, input_continuity, collocationMethod)
meshRefi = MeshRefi(resiThreshold, start, termination, MRMethod)

## solve
@benchmark solveMeshRefi(alg, meshRefi)
@time prob_alg = solveMeshRefi(alg, meshRefi);
plotResults(prob_alg);
# @time prob_alg = solveMeshRefi(prob, alg, meshRefi);

BenchmarkTools.DEFAULT_PARAMETERS.samples = 10;
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60;
io = IOContext(stdout)
medianTimeEarly2 = []
b = @benchmark prob_alg = solveMeshRefi(alg, meshRefi); show(io, MIME("text/plain"), b); println("\n");
medianTimeEarly2 = [medianTimeEarly2; median(b.times)];


solnInit, algInit, soln, alg, meshRefi = prob_alg;
plotSoln(solnInit, algInit, meshRefi)
plotSoln(soln, alg, meshRefi)
plotResults(prob_alg);

solnInit, algInit, soln, alg, meshRefi = prob_alg;
# soln, alg = prob_alg; plotResults(soln, alg);
# plotResults(prob_alg);





