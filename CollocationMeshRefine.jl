import Interpolations
import Ipopt
import MathOptInterface

using JuMP
using BenchmarkTools
using DifferentialEquations
using Polynomials
using Plots

## Mesh refinemnet definition
include("meshRefine\\typeMeshRefinement.jl")

## algorithm definition
include("collocation\\typeCollocation.jl")

## problem definition
include("collocation\\typeDOProblem.jl")

## create initial guess for problem
include("collocation\\add_initialGuess.jl")

## build each part of NLP problem for dynamic optimization probelm
include("collocation\\add_earlyTerm.jl")
include("collocation\\add_variables.jl")
include("collocation\\add_guess.jl")
include("collocation\\add_constraints.jl")
include("collocation\\add_objective.jl")

include("collocation\\interpolation.jl")
include("collocation\\simpsonQuadrature.jl")
include("collocation\\soln_trajectory.jl")

## solve the probelm
include("collocation\\solveMeshRefi.jl")

## mesh refinement function
include("meshRefine\\MRMinimaxError.jl")
include("meshRefine\\MREquidistributedError.jl")
include("meshRefine\\MRBisection.jl")

include("meshRefine\\MRQuasiWarmStart.jl")
include("meshRefine\\MRWarmStart.jl")
include("meshRefine\\MRColdStart.jl")

## plot the results
include("plotResult\\plotResults.jl")