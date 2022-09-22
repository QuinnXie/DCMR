function meshRefinement(numInsert::AbstractArray, 
    prob::DOProblem, 
    alg::Collocation, 
    meshRefi::T)   where T <: Union{MRBisectionAll{onceWarmStart}, MRBisectionAll{coldStart},
                                    MRBisectionMulti{onceWarmStart}, MRBisectionMulti{coldStart}, 
                                    MRBisectionOne{onceWarmStart}, MRBisectionOne{coldStart},
                                    MRMinimax{onceWarmStart}, MRMinimax{coldStart},
                                    MREquidistributed{onceWarmStart}, MREquidistributed{coldStart}}
    N = alg.numInterval

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

    return prob.solution, alg
end