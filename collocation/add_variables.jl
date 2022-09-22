function add_time(model::Model, prob::DOProblem, alg::Collocation)
    tSample = prob.solution.solnTime
    nTime = size(tSample,1)
    @variable(model, Δt[i = 1:nTime] == tSample[i])     # time step (sec)
end

function add_states(model::Model, state::DOVariables, alg::Collocation)
    n = alg.numPoint;

    nstate = size(state.initialLow,1);
    @variable(model, q[1:n, 1:nstate]);
    [set_lower_bound.(q[i, 1], state.boundsLow[1]) for i in 1:n]
    [set_upper_bound.(q[i, 1], state.boundsUpp[1]) for i in 1:n] 
    ## set flat to show if there is low or upper bound!!!

    ## Fix initial conditions
    [fix(q[1, k], state.initialLow[k]; force=true) for k in 1:nstate]

    ## Fix final conditions
    [fix(q[n, k], state.terminalLow[k]; force=true) for k in 1:nstate]
    # return model
end

function add_inputs(model::Model, input::DOVariables, alg::Trapezoidal{Continuous})
    n = alg.numPoint;
    @variable(model, input.boundsLow ≤ u[1:n] ≤ input.boundsUpp)                # motor force
end

function add_inputs(model::Model, input::DOVariables, alg::Trapezoidal{Discontinuous})
    N = alg.numInterval;
    @variable(model, input.boundsLow ≤ u[1:N, 1:2] ≤ input.boundsUpp)           # motor force
end

function add_inputs(model::Model, input::DOVariables, alg::HermiteSimpson{Continuous})
    n = alg.numPoint;
    @variable(model, input.boundsLow ≤ u[1:n] ≤ input.boundsUpp)                # motor force
end

function add_inputs(model::Model, input::DOVariables, alg::HermiteSimpson{Discontinuous})
    N = alg.numInterval;
    @variable(model, input.boundsLow ≤ u[1:N, 1:3] ≤ input.boundsUpp)           # motor force
end

function add_variables(model::Model, prob::DOProblem{T}, alg::Collocation) where T<:JuDOTime # where C<:Continuity
    add_time(model, prob, alg)
    add_states(model, prob.state, alg)
    add_inputs(model, prob.input, alg)
    return all_variables(model)
end
