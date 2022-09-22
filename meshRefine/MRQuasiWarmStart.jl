## mesh refinement method without warm start 
function meshRefinement(numInsert::AbstractArray, 
                        prob::DOProblem, 
                        alg::Collocation, 
                        meshRefi::T)   where T <: Union{MRBisectionAll{quasiWarmStart},
                                                        MRBisectionMulti{quasiWarmStart}, 
                                                        MRBisectionOne{quasiWarmStart},
                                                        MRMinimax{quasiWarmStart},
                                                        MREquidistributed{quasiWarmStart}}
    N = alg.numInterval
    n = alg.numPoint
    nstate = size(prob.state.initialLow,1)
    ninput = size(prob.input.initialLow,1)

    indexMesh = []
    for i in eachindex(numInsert)
        indexMesh = numInsert[i] >= 1 ? [indexMesh; i] : indexMesh
    end

    # new setting of time
    guessTime = prob.solution.solnTime;
    indexMesh = reverse(indexMesh)
    for i in indexMesh
        # time multidivide
        guessTime[i] = guessTime[i]/(numInsert[i]+1)
        [insert!(guessTime, i, guessTime[i]) for k in 1:numInsert[i]]
    end

    # new number of intervals
    N = alg.numInterval + Int64(sum(numInsert))
    alg = Collocation(alg, N)
    n = alg.numPoint
    
    interp_linearState = Interpolations.LinearInterpolation([1, n], [prob.state.initialLow, prob.state.terminalLow])

    prob.solution.solnTime = guessTime    
    prob.solution.solnState = mapreduce(transpose, vcat, interp_linearState.(1:n))
    prob.solution.solnInput = initial_guess(alg) # if both are Integer, there would be an error!
    prob.solution.solnDynamic = initial_guess(n, 0.0, 0.0)
    prob.solution.solnFunc = DOSolnFunc(zeros(N+1), 
                        Array{Polynomial{Float64},2}(undef,N,nstate),
                        Array{Polynomial{Float64},2}(undef,N,ninput),
                        Array{Polynomial{Float64},2}(undef,N,nstate))

    # println("quasiWarmStart")
    return prob.solution, alg
end



