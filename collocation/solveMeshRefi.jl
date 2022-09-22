function solve(prob::DOProblem, alg::Collocation, meshRefi::MeshRefi)
    model = Model()
    
    # set_optimizer for model such as Ipopt.Optimizer()
    opt = deepcopy(alg.optimizer)
    set_optimizer(model, ()->opt)

    add_variables(model, prob, alg);
    add_guess(model, prob, meshRefi.start);
    add_constraints(model, alg);
    add_objective(model, prob, alg);
    add_earlyTerm(model, alg, meshRefi.termination) # early termination

    set_silent(model);  # Hide solver's verbose output
    optimize!(model);   # Solve for the control and state 
    @assert termination_status(model) == LOCALLY_SOLVED;

    # trajectory and errorfunction
    prob, alg = soln_trajectory(model, prob, value, alg);

    # @show alg.numInterval
    return prob.solution, alg # objective_value(model)
end

function solve(prob::DOProblem, alg::Collocation)
    model = Model()
    
    # set_optimizer for model such as Ipopt.Optimizer()
    opt = deepcopy(alg.optimizer)
    set_optimizer(model, ()->opt)

    add_variables(model, prob, alg);
    add_guess(model, prob, alg);
    add_constraints(model, alg);
    add_objective(model, prob, alg);
    # add_earlyTerm(model, alg, meshRefi.termination) # early termination

    set_silent(model);  # Hide solver's verbose output
    optimize!(model);   # Solve for the control and state 
    @assert termination_status(model) == LOCALLY_SOLVED;
    
    # trajectory and errorfunction
    prob, alg = soln_trajectory(model, prob, value, alg);

    return prob.solution, alg # objective_value(model)
end

function solveMeshRefi(
    # prob::DOProblem,
    alg::Collocation = Collocation(
        optimizer = Ipopt.Optimizer(), 
        5, 
        "Discontinuous",
        "Hermite-Simpson"),
    meshRefi::MeshRefi = MRBisectionAll(
        1e-6, 
        warmStart(), 
        earlyTerm())
    )
    """
    optimizer = Ipopt.Optimizer()           # optimizer selection
    numInterval = 5                         # number of intervals
    input_continuity = "Continuous"         # continuity of input
    collocationmethod = "Hermite-Simpson"   # "Trapezoidal" or "HermiteSimpson"

    resiThreshold = 1e-6
    start = warmStart()
    termination = earlyTerm()
    MRMethod = "MRMinimax

    alg = Collocation(optimizer, numInterval, input_continuity, collocationMethod)
    meshRefi = MeshRefi(resiThreshold, start, termination, MRMethod)
    """

    ## Dynamic optimization problem defined
    time  = JuDOTime(t_s, t_t);               # DOVariable(t_s, t_t)
    state = DOVariables(q_s, q_s, q_t, q_t, q_lb, q_ub);
    input = DOVariables(u_lb, u_ub);
    guess = initial_guess(time, state, input, alg);
    funcDynamics = CartPole!;                  # system dynamics

    prob = DOProblem(time, state, input, funcDynamics, guess);

    ## solve 
    soln, alg = solve(prob, alg, meshRefi);
    # save initial results
    solnInit = deepcopy(soln); algInit = deepcopy(alg);
    ## mesh refinement
    solnFina, algFina = meshRefinement(prob, alg, meshRefi);

    return solnInit, algInit,solnFina, algFina, meshRefi;
    # return nothing
end

function solveMeshRefi(
    prob::DOProblem,
    alg::Collocation = Collocation(
        optimizer = Ipopt.Optimizer(), 
        5, 
        "Discontinuous",
        "Hermite-Simpson"),
    meshRefi::MeshRefi = MRBisectionAll(
        1e-6, 
        warmStart(), 
        earlyTerm())
    )
    """
    optimizer = Ipopt.Optimizer()           # optimizer selection
    numInterval = 5                         # number of intervals
    input_continuity = "Continuous"         # continuity of input
    collocationmethod = "Hermite-Simpson"   # "Trapezoidal" or "HermiteSimpson"

    resiThreshold = 1e-6
    start = warmStart()
    termination = earlyTerm()
    MRMethod = "MRMinimax

    alg = Collocation(optimizer, numInterval, input_continuity, collocationMethod)
    meshRefi = MeshRefi(resiThreshold, start, termination, MRMethod)
    """

    ## solve 
    soln, alg = solve(prob, alg, meshRefi);
    # save initial results
    solnInit = deepcopy(soln); algInit = deepcopy(alg);
    ## mesh refinement
    solnFina, algFina = meshRefinement(prob, alg, meshRefi);

    return solnInit, algInit,solnFina, algFina, meshRefi;
    # return nothing
end