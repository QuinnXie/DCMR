function DOSolnError(solnFunc::DOSolnFunc, prob::DOProblem, alg::Collocation)
    N = alg.numInterval
    nstate = size(prob.state.initialLow, 1)

    ## Error - Simpson quadrature
    # Interpolation for checking collocation constraint along trajectory:
    # collocation constraint = (dynamics) - (derivative of state trajectory)
    ts = solnFunc.FuncTime
    solnFuncState = solnFunc.FuncState
    solnFuncDynamic = solnFunc.FuncDynamic
    solnFuncInput = solnFunc.FuncInput

    errorQuad = zeros(N, nstate)
    errorSqua = zeros(N)
    errorMean = zeros(N)
    errorFunc = Array{Function}(undef, N)

    ## Collocation constraint error Careful!!! Not Dynamic!
    for i in 1:N
        cartPoleDynamics(t) = map(t -> CartPole!([solnFuncState[i,1](t),solnFuncState[i,2](t),solnFuncState[i,3](t),solnFuncState[i,4](t)],solnFuncInput[i](t)), t)
        cartPoleDynamicsSoln(t) = map(t -> [solnFuncDynamic[i,1](t); solnFuncDynamic[i,2](t); solnFuncDynamic[i,3](t); solnFuncDynamic[i,4](t)], t)
        cartPoleDynamicsError(t) = cartPoleDynamics(t) - cartPoleDynamicsSoln(t)
        errorFunc[i] = cartPoleDynamicsError
    end
    
    # simpsonQuadrature
    nQuadSegment = 9;
    [errorQuad[i, :] = simpsonQuadrature(errorFunc[i], ts[i], ts[i+1], nQuadSegment) for i in 1:N]
    # [errorSqua[i] = sum(errorQuad[i,:].^2) for i in 1:N]
    errorSqua = sum(errorQuad.^2, dims = 2)./(nstate*(ts[end] - ts[1])) # sqrt.(sum(errorQuad.^2, dims = 2))
    [errorMean[i] = errorSqua[i]/(ts[i+1]-ts[i])/2 for i in 1:N]
    
    errorTNIR = sum(errorSqua)

    # append!(prob.solution.solnError.errorTNIR, sum(errorMean)/nstate)
    append!(prob.solution.solnError.errorTNIR, errorTNIR)
    append!(prob.solution.solnError.numInterv, N)

    # @show prob.solution.solnError.errorTNIR[end]
    # @show errorQuad
    # @show errorMean
    # @show errorSqua
    return DOSolnError(ts, errorQuad, errorSqua, errorMean, 
    prob.solution.solnError.errorTNIR, prob.solution.solnError.numInterv, errorFunc)
end

function soln_trajectory(model::Model, prob::DOProblem, value, alg::Trapezoidal{Continuous})
    N = alg.numInterval
    nstate = size(prob.state.initialLow, 1)
    ninput = size(prob.input.initialLow, 1)

    # optimization solution
    ??t = value.(model.obj_dict[:??t])
    q  = value.(model.obj_dict[:q])
    u  = value.(model.obj_dict[:u])
    ??q =   [value.(model.obj_dict[:??q???])'; 
            value.(model.obj_dict[:??q???])'; 
            value.(model.obj_dict[:??q???])'; 
            value.(model.obj_dict[:??q???])']'

    ## initial tspan
    ts = cumsum([0; ??t])[1:end-1]

    ## Interpolation Polynomial matrix
    solnFuncState   = Array{Polynomial{Float64},2}(undef, N, nstate)
    solnFuncInput   = Array{Polynomial{Float64},2}(undef, N, ninput)
    solnFuncDynamic = Array{Polynomial{Float64},2}(undef, N, nstate)

    for i in 1:N
        j = i + 1
        [solnFuncState[i, k]   = interpx(q[i, k], ??q[i, k], ??q[j, k], ts[i], ts[j]) for k in 1:nstate]        
        [solnFuncInput[i, k]   = interpu((u[i, k]), (u[j, k]), ts[i], ts[j]) for k in 1:ninput]
        [solnFuncDynamic[i, k] = interpu(??q[i, k], ??q[j, k], ts[i], ts[j]) for k in 1:nstate] 
    end
    solnFunc = DOSolnFunc(ts, solnFuncState, solnFuncInput, solnFuncDynamic);
    
    ## Error - Simpson quadrature
    solnError = DOSolnError(solnFunc, prob, alg);

    # save and return results
    prob.solution = DOSolution(??t, q, u, ??q, solnFunc, solnError);

    return prob, alg
end

function soln_trajectory(model::Model, prob::DOProblem, value, alg::Trapezoidal{Discontinuous})
    N = alg.numInterval
    nstate = size(prob.state.initialLow, 1)
    ninput = size(prob.input.initialLow, 1)

    # optimization solution
    ??t  = value.(model.obj_dict[:??t])
    q   = value.(model.obj_dict[:q])
    u   = value.(model.obj_dict[:u])
    ??qL = [value.(model.obj_dict[:??q???L])'; value.(model.obj_dict[:??q???L])'; value.(model.obj_dict[:??q???L])'; value.(model.obj_dict[:??q???L])']'
    ??qR = [value.(model.obj_dict[:??q???R])'; value.(model.obj_dict[:??q???R])'; value.(model.obj_dict[:??q???R])'; value.(model.obj_dict[:??q???R])']'
    ??q  = [??qL, ??qR]
    ## initial tspan
    ts  = cumsum([0; ??t])[1:end-1]

    ## Interpolation Polynomial matrix
    solnFuncState   = Array{Polynomial{Float64},2}(undef, N, nstate)
    solnFuncInput   = Array{Polynomial{Float64},2}(undef, N, ninput)
    solnFuncDynamic = Array{Polynomial{Float64},2}(undef, N, nstate)

    for i in 1:N
        j = i + 1
        [solnFuncState[i, k] = interpx(q[i, k], ??q[1][i, k], ??q[2][i, k], ts[i], ts[j]) for k in 1:nstate]        
        solnFuncInput[i, 1] = interpu((u[i, 1]), (u[i, 2]), ts[i], ts[j])
        [solnFuncDynamic[i, k] = interpu(??q[1][i, k], ??q[2][i, k], ts[i], ts[j]) for k in 1:nstate] 
    end
    solnFunc = DOSolnFunc(ts, solnFuncState, solnFuncInput, solnFuncDynamic)
    
    ## Error - Simpson quadrature
    solnError = DOSolnError(solnFunc, prob, alg);

    # save and return results
    prob.solution = DOSolution(??t, q, u, ??q, solnFunc, solnError)

    return prob, alg
end

function soln_trajectory(model::Model, prob::DOProblem, value, alg::HermiteSimpson{Continuous})
    N = alg.numInterval
    nstate = size(prob.state.initialLow, 1)
    ninput = size(prob.input.initialLow, 1)

    # optimization solution
    ??t = value.(model.obj_dict[:??t])
    q  = value.(model.obj_dict[:q])
    u  = value.(model.obj_dict[:u])
    ??q =   [value.(model.obj_dict[:??q???])'; 
            value.(model.obj_dict[:??q???])'; 
            value.(model.obj_dict[:??q???])'; 
            value.(model.obj_dict[:??q???])']'

    ## initial tspan
    ts = cumsum([0; ??t])[1:end-1]

    ## Interpolation Polynomial matrix
    solnFuncState   = Array{Polynomial{Float64},2}(undef, N, nstate)
    solnFuncInput   = Array{Polynomial{Float64},2}(undef, N, ninput)
    solnFuncDynamic = Array{Polynomial{Float64},2}(undef, N, nstate)

    for m in 1:N
        j = 2 * m; i = j - 1; k = j + 1;
        [solnFuncState[m, s]   = interpx(q[i, s], ??q[i, s], ??q[j, s], ??q[k, s], ts[m], ts[m+1]) for s in 1:nstate]
        [solnFuncInput[m, s]   = interpu((u[i, s]), (u[j, s]), (u[k, s]), ts[m], ts[m+1]) for s in 1:ninput]
        [solnFuncDynamic[m, s] = interpu(??q[i, s], ??q[j, s], ??q[k, s], ts[m], ts[m+1]) for s in 1:nstate] 
    end
    solnFunc = DOSolnFunc(ts, solnFuncState, solnFuncInput, solnFuncDynamic)
    
    ## Error - Simpson quadrature
    solnError = DOSolnError(solnFunc, prob, alg);

    # save and return results
    prob.solution = DOSolution(??t, q, u, ??q, solnFunc, solnError)

    return prob, alg
end

function soln_trajectory(model::Model, prob::DOProblem, value, alg::HermiteSimpson{Discontinuous})
    N = alg.numInterval
    nstate = size(prob.state.initialLow, 1)
    ninput = size(prob.input.initialLow, 1)

    # optimization solution
    ??t  = value.(model.obj_dict[:??t])
    q   = value.(model.obj_dict[:q])
    u   = value.(model.obj_dict[:u])
    ??qL = [value.(model.obj_dict[:??q???L])'; value.(model.obj_dict[:??q???L])'; value.(model.obj_dict[:??q???L])'; value.(model.obj_dict[:??q???L])']'
    ??qM = [value.(model.obj_dict[:??q???M])'; value.(model.obj_dict[:??q???M])'; value.(model.obj_dict[:??q???M])'; value.(model.obj_dict[:??q???M])']'
    ??qR = [value.(model.obj_dict[:??q???R])'; value.(model.obj_dict[:??q???R])'; value.(model.obj_dict[:??q???R])'; value.(model.obj_dict[:??q???R])']'
    ??q = [??qL, ??qM, ??qR]

    ## initial tspan
    ts = cumsum([0; ??t])[1:end-1]

    ## Interpolation Polynomial matrix
    solnFuncState   = Array{Polynomial{Float64},2}(undef, N, nstate)
    solnFuncInput   = Array{Polynomial{Float64},2}(undef, N, ninput)
    solnFuncDynamic = Array{Polynomial{Float64},2}(undef, N, nstate)

    for m in 1:N
        i = 2 * m - 1;
        [solnFuncState[m, s] = interpx(q[i, s], ??q[1][m, s], ??q[2][m, s], ??q[3][m, s], ts[m], ts[m+1]) for s in 1:nstate]
        solnFuncInput[m] = interpu(u[m, 1], u[m, 2], u[m, 3], ts[m], ts[m+1])
        [solnFuncDynamic[m, k] = interpu(??q[1][m, k], ??q[2][m, k], ??q[3][m, k], ts[m], ts[m+1]) for k in 1:nstate] 
    end
    solnFunc = DOSolnFunc(ts, solnFuncState, solnFuncInput, solnFuncDynamic)
    
    ## Error - Simpson quadrature
    solnError = DOSolnError(solnFunc, prob, alg);

    # save and return results
    prob.solution = DOSolution(??t, q, u, ??q, solnFunc, solnError);

    return prob, alg
end