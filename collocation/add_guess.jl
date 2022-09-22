function add_guess(model::Model, prob::DOProblem, start::Union{quasiWarmStart, warmStart})
    # Initial guess: time
    vec_guessTime = prob.solution.solnTime
    # Initial guess: state
    initialState = prob.solution.solnState
    vec_guessState = reshape(initialState,:,1)
    # Initial guess: input
    initialInput = prob.solution.solnInput
    vec_guessInput = reshape(initialInput,:,1)

    vec_initial_guess = append!(vec_guessTime, vec_guessState, vec_guessInput)
    set_start_value.(all_variables(model), vec_initial_guess)
end

function add_guess(model::Model, prob::DOProblem, start::onceWarmStart)
    if prob.solution.solnError.errorTNIR == []
        add_guess(model, prob, warmStart())
    end
    # @show "onceWarmStart"
end

function add_guess(model::Model, prob::DOProblem, start::coldStart)
    # @show "no initial guess"
    return "no initial guess"
end