abstract type AbstractProblem end
abstract type JuDOTime end

function JuDOTime(initial::Real,terminal::Real)
    return fixedTime(initial,terminal)
end

function JuDOTime(initial::Real,terminal::Tuple{Real,Real})
    return variableTime(initial,terminal)
end

mutable struct fixedTime{} <: JuDOTime
    initial::Real
    terminal::Real
end

mutable struct variableTime{} <: JuDOTime
    initial::Real
    terminal::Tuple{Real,Real}
end

mutable struct DOVariables{F<:Union{AbstractFloat, AbstractArray}}
    initialLow::F
    initialUpp::F

    terminalLow::F
    terminalUpp::F

    boundsLow::F
    boundsUpp::F
end

function DOVariables(
    initialLow::Union{Real,AbstractArray},
    initialUpp::Union{Real,AbstractArray},
    terminalLow::Union{Real,AbstractArray},
    terminalUpp::Union{Real,AbstractArray},
    boundsLow::Union{Real,AbstractArray},
    boundsUpp::Union{Real,AbstractArray}
    )
    return DOVariables(
        Float64.(initialLow),
        Float64.(initialUpp),
        Float64.(terminalLow),
        Float64.(terminalUpp),
        Float64.(boundsLow),
        Float64.(boundsUpp)
    )
end

function DOVariables(boundsLow::T, boundsUpp::T) where T<:Union{Real,AbstractArray}
    return DOVariables(boundsLow, boundsUpp, boundsLow, boundsUpp, boundsLow, boundsUpp)
end

function DOVariables(initial::T, terminal::T) where T<:Union{Real,AbstractArray}
    return DOVariables(initial, initial, terminal, terminal, -Inf, Inf)
end

mutable struct DOSolnFunc{}
    FuncTime::Array{Float64,1}
    FuncState::Array{Polynomial{Float64},2}
    FuncInput::Array{Polynomial{Float64},2}
    FuncDynamic::Array{Polynomial{Float64},2}
end

mutable struct DOSolnError{}
    errorTime::Array{Float64,1}
    errorQuad::AbstractArray
    errorSqua::AbstractArray
    errorMean::AbstractArray
    errorTNIR::Vector{Float64}
    numInterv::Vector{Integer}
    errorFunc::Array{Function}
end

mutable struct DOSolution{}
    solnTime::AbstractArray
    solnState::AbstractArray
    solnInput::AbstractArray
    solnDynamic::AbstractArray
    solnFunc::DOSolnFunc
    solnError::DOSolnError
end

mutable struct DOProblem{T<:JuDOTime} <: AbstractProblem
    time::T;
    state::DOVariables;
    input::DOVariables;

    # ## function
    funcDynamics::Function;
    # MayerCost::Expr
    # runningCost::Expr # Force-squared cost function

    # guess::DOGuess
    solution::DOSolution; # initial guess given here
end