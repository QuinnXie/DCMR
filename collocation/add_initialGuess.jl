function initial_guess(time::T, state::DOVariables, input::DOVariables, alg::Collocation) where T<:fixedTime
    N = alg.numInterval
    n = alg.numPoint
    nstate = size(state.initialLow,1)
    ninput = size(input.initialLow,1)

    duartion = time.terminal - time.initial;
    guessTime = ones(N+1).*duartion/N

    interp_linearState = Interpolations.LinearInterpolation([1, n], [state.initialLow, state.terminalLow])
    guessState = mapreduce(transpose, vcat, interp_linearState.(1:n))

    guessInput = initial_guess(alg) # if both are Integer, there would be an error!
    guessDynamic = initial_guess(n, 0.0, 0.0)
    guessFunc = DOSolnFunc(zeros(N+1), 
                        Array{Polynomial{Float64},2}(undef,N,nstate),
                        Array{Polynomial{Float64},2}(undef,N,ninput),
                        Array{Polynomial{Float64},2}(undef,N,nstate)
                        )
    guessError = DOSolnError(zeros(N+1), 
                        Array{Float64,2}(undef,N,nstate),
                        Array{Float64,2}(undef,N,1),
                        Array{Float64,2}(undef,N,1),
                        [],
                        [],
                        Array{Function}(undef, N)
                        )
    return DOSolution(guessTime, guessState, guessInput, guessDynamic, guessFunc, guessError)
end

## now the setting is the same with that of T<:fixedTime
function initial_guess(time::T, state::DOVariables, input::DOVariables, alg::Collocation) where T<:variableTime
    N = alg.numInterval
    n = alg.numPoint
    nstate = size(state.initialLow,1)
    ninput = size(input.initialLow,1)
    
    guessTime  = ones(N+1).*((time.terminal - time.initial)/N);
    guessState = initial_guess(n, Float64.(state.initialLow), Float64.(state.terminalLow))
    guessInput = initial_guess(n, Float64.(input.initialLow), Float64.(input.initialLow)) # if both are Integer, there would be an error!
    guessDynamic = initial_guess(n, 0.0, 0.0)
    guessFunc = DOSolnFunc(zeros(N+1), 
                        Array{Polynomial{Float64},2}(undef,N,nstate),
                        Array{Polynomial{Float64},2}(undef,N,ninput),
                        Array{Polynomial{Float64},2}(undef,N,nstate)
                        )
    guessError = DOSolnError(zeros(N+1), 
                        Array{Float64,2}(undef,N,nstate),
                        Array{Float64,2}(undef,N,1),
                        Array{Float64,2}(undef,N,1),
                        [0.0],
                        [0],
                        Array{Function}(undef, N)
                        )
    return DOSolution(guessTime, guessState, guessInput,guessDynamic,guessFunc,guessError)
end

function initial_guess(n::Integer, low::Union{Real,AbstractArray}, upp::Union{Real,AbstractArray}) # where T<:Union{Real,AbstractArray}
    interp_linear = Interpolations.LinearInterpolation([1, n], [low, upp])
    return mapreduce(transpose, vcat, interp_linear.(1:n))
end

function initial_guess(alg::Trapezoidal{Continuous})
    n = alg.numPoint
    return zeros(n,1)
end 

function initial_guess(alg::Trapezoidal{Discontinuous})
    N = alg.numInterval
    return zeros(2*N,1)
end 

function initial_guess(alg::HermiteSimpson{Continuous})
    n = alg.numPoint
    return zeros(n,1)
end 

function initial_guess(alg::HermiteSimpson{Discontinuous})
    N = alg.numInterval
    return zeros(3*N,1)
end
