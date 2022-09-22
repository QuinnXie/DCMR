## mesh refinement for Trapezoidal{Continuous} with warm start 
function meshRefinement(numInsert::AbstractArray, 
                        prob::DOProblem, 
                        alg::Trapezoidal{Continuous},
                        meshRefi::T)   where T <: Union{MRBisectionAll{warmStart},
                                                        MRBisectionMulti{warmStart}, 
                                                        MRBisectionOne{warmStart},
                                                        MRMinimax{warmStart},
                                                        MREquidistributed{warmStart}}
    ts = prob.solution.solnTime
    nstate = size(prob.solution.solnState,2)
    ninput = size(prob.solution.solnInput,2)

    indexMesh = []
    for i in eachindex(numInsert)
        indexMesh = numInsert[i] >= 1 ? [indexMesh; i] : indexMesh;
    end

    vecSolnState = reshape(prob.solution.solnState',:,1)
    vecSolnState = mapreduce(transpose,vcat,vecSolnState)
    vecSolnInput = reshape(prob.solution.solnInput',:,1)
    vecSolnInput = mapreduce(transpose,vcat,vecSolnInput)

    indexMesh = reverse(indexMesh)
    for i in indexMesh
        addState = [];
        addInput = [];
        M₁ = numInsert[i]
        N₁ = M₁+1
        # time subdivide
        ts[i] = ts[i]/(M₁+1)
        for k in 1:M₁
            insert!(ts, i, ts[i])                    
        end

        for k in 1:M₁
            # state interpolation
            addState = [addState; (prob.solution.solnState[i,:]*(N₁-k)/N₁+prob.solution.solnState[i+1,:]*k/N₁)]
            # input interpolation
            addInput = [addInput; (prob.solution.solnInput[i,:]*(N₁-k)/N₁+prob.solution.solnInput[i+1,:]*k/N₁)]   
        end
        vecSolnState = [vecSolnState[1:nstate*i]; addState; vecSolnState[nstate*i+1:end]]
        vecSolnInput = [vecSolnInput[1:ninput*i]; addInput; vecSolnInput[ninput*i+1:end]]
    end

    numInterval = alg.numInterval + Int64(sum(numInsert))
    alg = Collocation(alg, numInterval)
    
    prob.solution.solnTime = ts
    prob.solution.solnState = reshape(vecSolnState,nstate,:)'
    prob.solution.solnInput = reshape(vecSolnInput,ninput,:)'

    return prob.solution, alg
end

## mesh refinement for Trapezoidal{Discontinuous} with warm start 
function meshRefinement(numInsert::AbstractArray,
                        prob::DOProblem, 
                        alg::Trapezoidal{Discontinuous},
                        meshRefi::T)   where T <: Union{MRBisectionAll{warmStart},
                                                        MRBisectionMulti{warmStart}, 
                                                        MRBisectionOne{warmStart},
                                                        MRMinimax{warmStart},
                                                        MREquidistributed{warmStart}}
    ts = prob.solution.solnTime
    nstate = size(prob.solution.solnState,2)
    ninput = size(prob.solution.solnInput,2)

    indexMesh = []
    for i in eachindex(numInsert)
        indexMesh = numInsert[i] >= 1 ? [indexMesh; i] : indexMesh;
    end

    vecSolnState = reshape(prob.solution.solnState',:,1)
    vecSolnState = mapreduce(transpose,vcat,vecSolnState)
    vecSolnInput = reshape(prob.solution.solnInput',:,1)
    vecSolnInput = mapreduce(transpose,vcat,vecSolnInput)

    indexMesh = reverse(indexMesh)
    for i in indexMesh
        addState = [];
        addInput = [];
        M₁ = numInsert[i]
        N₁ = M₁+1
        # time subdivide
        ts[i] = ts[i]/(M₁+1)
        for k in 1:M₁
            insert!(ts, i, ts[i])                    
        end

        for k in 1:M₁
            # state interpolation
            addState = [addState; (prob.solution.solnState[i,:]*(N₁-k)/N₁+prob.solution.solnState[i+1,:]*k/N₁)]
            # input interpolation
            addInput = [addInput; (prob.solution.solnInput[i,1]*(N₁-k)/N₁+prob.solution.solnInput[i,2]*k/N₁);
                                  (prob.solution.solnInput[i,1]*(N₁-k)/N₁+prob.solution.solnInput[i,2]*k/N₁)]
        end
        vecSolnState = [vecSolnState[1:nstate*i]; addState; vecSolnState[nstate*i+1:end]]
        vecSolnInput = [vecSolnInput[1:ninput*i]; addInput; vecSolnInput[ninput*i+1:end]]
    end

    numInterval = alg.numInterval + Int64(sum(numInsert))
    alg = Collocation(alg, numInterval)
    
    prob.solution.solnTime = ts
    prob.solution.solnState = reshape(vecSolnState,nstate,:)'
    prob.solution.solnInput = reshape(vecSolnInput,ninput,:)'

    return prob.solution, alg
end

## mesh refinement for HermiteSimpson{Continuous} with warm start 
function meshRefinement(numInsert::AbstractArray, 
                        prob::DOProblem, 
                        alg::HermiteSimpson{Continuous},
                        meshRefi::T)   where T <: Union{MRBisectionAll{warmStart},
                                                        MRBisectionMulti{warmStart}, 
                                                        MRBisectionOne{warmStart},
                                                        MRMinimax{warmStart},
                                                        MREquidistributed{warmStart}}
    ts = prob.solution.solnTime
    nstate = size(prob.solution.solnState,2)
    ninput = size(prob.solution.solnInput,2)

    indexMesh = []
    for i in eachindex(numInsert)
        indexMesh = numInsert[i] >= 1 ? [indexMesh; i] : indexMesh;
    end

    vecSolnState = reshape(prob.solution.solnState',:,1)
    vecSolnState = mapreduce(transpose,vcat,vecSolnState)
    vecSolnInput = reshape(prob.solution.solnInput',:,1)
    vecSolnInput = mapreduce(transpose,vcat,vecSolnInput)

    indexMesh = reverse(indexMesh)
    for i in indexMesh
        addStateL = []; addStateR = [];
        addInputL = []; addInputR = [];
        M₁ = numInsert[i]
        N₁ = M₁+1
        # time subdivide
        ts[i] = ts[i]/(M₁+1)
        for k in 1:M₁
            insert!(ts, i, ts[i])                    
        end

        for k in 1:M₁
            # state interpolation
            addStateL = [addStateL; (prob.solution.solnState[2*i-1,:]*(N₁-k)/N₁+prob.solution.solnState[2*i,:]*k/N₁)]
            addStateR = [addStateR; (prob.solution.solnState[2*i,:]*(N₁-k)/N₁+prob.solution.solnState[2*i+1,:]*k/N₁)]
            # input interpolation
            addInputL = [addInputL; (prob.solution.solnInput[2*i-1,:]*(N₁-k)/N₁+prob.solution.solnInput[2*i,:]*k/N₁)]
            addInputR = [addInputR; (prob.solution.solnInput[2*i,:]*(N₁-k)/N₁+prob.solution.solnInput[2*i+1,:]*k/N₁)] 
        end
        vecSolnState = [vecSolnState[1:nstate*(2*i-1)]; 
                        addStateL; vecSolnState[(nstate*(2*i-1)+1):nstate*(2*i)]; addStateR; 
                        vecSolnState[nstate*(2*i)+1:end]]
        vecSolnInput = [vecSolnInput[1:ninput*(2*i-1)]; 
                        addInputL; vecSolnInput[ninput*(2*i-1)+1:ninput*(2*i)]; addInputR;
                        vecSolnInput[ninput*(2*i)+1:end]]
    end

    numInterval = alg.numInterval + Int64(sum(numInsert))
    alg = Collocation(alg, numInterval)
    
    prob.solution.solnTime = ts
    prob.solution.solnState = reshape(vecSolnState,nstate,:)'
    prob.solution.solnInput = reshape(vecSolnInput,ninput,:)'

    return prob.solution, alg
end

## mesh refinement for HermiteSimpson{Discontinuous} with warm start 
function meshRefinement(numInsert::AbstractArray,
                        prob::DOProblem, 
                        alg::HermiteSimpson{Discontinuous},
                        meshRefi::T)   where T <: Union{MRBisectionAll{warmStart},
                                                        MRBisectionMulti{warmStart}, 
                                                        MRBisectionOne{warmStart},
                                                        MRMinimax{warmStart},
                                                        MREquidistributed{warmStart}}
    ts = prob.solution.solnTime
    nstate = size(prob.solution.solnState,2)
    ninput = size(prob.solution.solnInput,2)

    indexMesh = []
    for i in eachindex(numInsert)
        indexMesh = numInsert[i] >= 1 ? [indexMesh; i] : indexMesh;
    end

    vecSolnState = reshape(prob.solution.solnState',:,1)
    vecSolnState = mapreduce(transpose,vcat,vecSolnState)
    vecSolnInput = reshape(prob.solution.solnInput',:,1)
    vecSolnInput = mapreduce(transpose,vcat,vecSolnInput)

    indexMesh = reverse(indexMesh)
    for i in indexMesh
        addStateL = []; addStateR = [];
        addInputL = []; addInputM = []; addInputR = [];
        M₁ = numInsert[i]
        N₁ = M₁+1
        # time subdivide
        ts[i] = ts[i]/(M₁+1)
        for k in 1:M₁
            insert!(ts, i, ts[i])                    
        end
        for k in 1:M₁
            # state interpolation
            addStateL = [addStateL; (prob.solution.solnState[2*i-1,:]*(N₁-k)/N₁+prob.solution.solnState[2*i,:]*k/N₁)]
            addStateR = [addStateR; (prob.solution.solnState[2*i,:]*(N₁-k)/N₁+prob.solution.solnState[2*i+1,:]*k/N₁)]
        end
        # input interpolation
        if mod(M₁,2) == 1
            for k in 1:M₁
                if mod(k,2) == 1
                    addInputL = [addInputL; (prob.solution.solnInput[i,1]*(N₁-k)/N₁+prob.solution.solnInput[i,2]*k/N₁)]
                    addInputR = [addInputR; (prob.solution.solnInput[i,2]*(N₁-k)/N₁+prob.solution.solnInput[i,3]*k/N₁)]     
                else
                    addInputL = [addInputL; (prob.solution.solnInput[i,1]*(N₁-k)/N₁+prob.solution.solnInput[i,2]*k/N₁);
                                            (prob.solution.solnInput[i,1]*(N₁-k)/N₁+prob.solution.solnInput[i,2]*k/N₁)]
                    addInputR = [addInputR; (prob.solution.solnInput[i,2]*(N₁-k)/N₁+prob.solution.solnInput[i,3]*k/N₁)
                                            (prob.solution.solnInput[i,2]*(N₁-k)/N₁+prob.solution.solnInput[i,3]*k/N₁)]                              
                end                     
            end
            addInputM = [prob.solution.solnInput[i,2]; prob.solution.solnInput[i,2]]
        else
            for k in 1:M₁
                if mod(k,2) == 1
                    addInputL = [addInputL; (prob.solution.solnInput[i,1]*(N₁-k)/N₁+prob.solution.solnInput[i,2]*k/N₁)]
                    addInputR = [addInputR; (prob.solution.solnInput[i,2]*(N₁-k)/N₁+prob.solution.solnInput[i,3]*k/N₁)]     
                else
                    addInputL = [addInputL; (prob.solution.solnInput[i,1]*(N₁-k)/N₁+prob.solution.solnInput[i,2]*k/N₁);
                                            (prob.solution.solnInput[i,1]*(N₁-k)/N₁+prob.solution.solnInput[i,2]*k/N₁)]
                    addInputR = [addInputR; (prob.solution.solnInput[i,2]*(N₁-k)/N₁+prob.solution.solnInput[i,3]*k/N₁)
                                            (prob.solution.solnInput[i,2]*(N₁-k)/N₁+prob.solution.solnInput[i,3]*k/N₁)]                              
                end                      
            end
            addInputM = prob.solution.solnInput[i,2]
       end

        vecSolnState = [vecSolnState[1:nstate*(2*i-1)]; 
                        addStateL; vecSolnState[(nstate*(2*i-1)+1):nstate*(2*i)]; addStateR; 
                        vecSolnState[nstate*(2*i)+1:end]]
        vecSolnInput = [vecSolnInput[1:ninput*(i-1)+1]; 
                        addInputL; addInputM; addInputR;
                        vecSolnInput[ninput*i:end]]
    end

    numInterval = alg.numInterval + Int64(sum(numInsert))
    alg = Collocation(alg, numInterval)
    
    prob.solution.solnTime = ts
    prob.solution.solnState = reshape(vecSolnState,nstate,:)'
    prob.solution.solnInput = reshape(vecSolnInput,ninput,:)'

    return prob.solution, alg
end