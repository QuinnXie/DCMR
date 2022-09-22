const MOI = MathOptInterface

# define type of collocation method
abstract type Collocation end

# define type of continuity
abstract type Continuity end
struct Continuous <: Continuity end
struct Discontinuous <: Continuity end

# define Trapezoidal Collocation parameters
# mutable struct Trapezoidal{C<:Continuity, M<:MeshRefi} <: Collocation
mutable struct Trapezoidal{C<:Continuity} <: Collocation
    optimizer::MOI.AbstractOptimizer
    numInterval::Integer
    numPoint::Integer
    input_continuity::C
end

# build function for Trapezoidal Collocation method
function Trapezoidal(
        optimizer_factory = Ipopt.Optimizer(), 
        N::Integer = 5, 
        input_continuity::Continuity = Discontinuous()
        )
    return Trapezoidal(optimizer_factory, N, N+1, input_continuity)
end

# define Hermite-Simpson Collocation parameters
mutable struct HermiteSimpson{C<:Continuity} <: Collocation  # O<:MOI.AbstractOptimizer
    optimizer::MOI.AbstractOptimizer
    numInterval::Integer
    numPoint::Integer
    input_continuity::C
end
# build function for Hermite-Simpson Collocation method
function HermiteSimpson(
        optimizer_factory = Ipopt.Optimizer(), 
        N::Integer = 25, 
        input_continuity::Continuity = Discontinuous(), 
        )
    return HermiteSimpson(optimizer_factory, N, 2*N+1, input_continuity)
end

# build function for changing the number of interval in collocation method
function Collocation(alg::Trapezoidal, numInterval::Integer)
    return Trapezoidal(
        alg.optimizer, 
        numInterval, 
        numInterval+1, 
        alg.input_continuity)
end

function Collocation(alg::HermiteSimpson, numInterval::Integer)
    return HermiteSimpson(
        alg.optimizer, 
        numInterval, 
        2*numInterval+1, 
        alg.input_continuity) 
end

# both Trapezoidal Collocation and HermiteSimpson Collocation included in Collocation function
function Collocation(
    optimizer_factory = Ipopt.Optimizer(), 
    N::Integer = 5, 
    input_continuity::Continuity = Discontinuous(), 
    collocMethod::String = "HermiteSimpson")

    if collocMethod == "Trapezoidal"
        return Trapezoidal(optimizer_factory, N, N+1, input_continuity)
    elseif collocMethod == "Hermite-Simpson"
        return HermiteSimpson(optimizer_factory, N, 2*N+1, input_continuity)   
    else
        @error "Unexpected collocation method '$(collocMethod)'"
    end
end

## type of solution state to distinguish initial solution and final solution
abstract type SolnState end
struct Initial <: SolnState end
struct Final <: SolnState end