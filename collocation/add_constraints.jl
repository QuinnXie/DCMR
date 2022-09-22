function add_constraints(model::Model, alg::Trapezoidal{Continuous})
    N = alg.numInterval
    n = alg.numPoint

    ## Motion of the vehicle as a differential-algebraic system of equations (DAEs)
    Δt = model.obj_dict[:Δt]
    q = model.obj_dict[:q]
    u = model.obj_dict[:u]
    @NLexpression(model, δq₁[j=1:n], q[j,3])
    @NLexpression(model, δq₂[j=1:n], q[j,4])
    @NLexpression(
        model, 
        δq₃[j=1:n], 
        (𝓁 * m₂ * sin(q[j,2]) * q[j,4]^2 + u[j] + m₂ * g * cos(q[j,2]) * sin(q[j,2])) / (m₁ + m₂ * (1 - cos(q[j,2])^2)))
    @NLexpression(
        model, 
        δq₄[j=1:n], 
        -(𝓁 * m₂ * cos(q[j,2]) * sin(q[j,2]) * q[j,4]^2 + u[j] * cos(q[j,2]) + (m₁ + m₂) * g * sin(q[j,2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(q[j,2])^2)))

    ## System dynamics
    for i in 1:N
        j = i + 1  # index of latter knot
        ## Trapezoidal integration
        @NLconstraint(model, q[j,1] == q[i,1] + 0.5 * Δt[i] * (δq₁[j] + δq₁[i]))
        @NLconstraint(model, q[j,2] == q[i,2] + 0.5 * Δt[i] * (δq₂[j] + δq₂[i]))
        @NLconstraint(model, q[j,3] == q[i,3] + 0.5 * Δt[i] * (δq₃[j] + δq₃[i]))
        @NLconstraint(model, q[j,4] == q[i,4] + 0.5 * Δt[i] * (δq₄[j] + δq₄[i]))
    end
end

function add_constraints(model::Model, alg::Trapezoidal{Discontinuous})
    N = alg.numInterval
    # n = alg.numPoint

    ## Motion of the vehicle as a differential-algebraic system of equations (DAEs)
    Δt = model.obj_dict[:Δt]
    q = model.obj_dict[:q]
    u = model.obj_dict[:u]

    ## Motion of the vehicle as a differential-algebraic system of equations (DAEs)
    @NLexpression(model, δq₁L[j=1:N], q[j,3])
    @NLexpression(model, δq₂L[j=1:N], q[j,4])
    @NLexpression(model, δq₃L[j=1:N], (𝓁 * m₂ * sin(q[j,2]) * q[j,4]^2 + u[j,1] + m₂ * g * cos(q[j,2]) * sin(q[j,2])) / (m₁ + m₂ * (1 - cos(q[j,2])^2)))
    @NLexpression(model, δq₄L[j=1:N], -(𝓁 * m₂ * cos(q[j,2]) * sin(q[j,2]) * q[j,4]^2 + u[j,1] * cos(q[j,2]) + (m₁ + m₂) * g * sin(q[j,2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(q[j,2])^2)))
    @NLexpression(model, δq₁R[j=1:N], q[j+1,3])
    @NLexpression(model, δq₂R[j=1:N], q[j+1,4])
    @NLexpression(model, δq₃R[j=1:N], (𝓁 * m₂ * sin(q[j+1,2]) * q[j+1,4]^2 + u[j,2] + m₂ * g * cos(q[j+1,2]) * sin(q[j+1,2])) / (m₁ + m₂ * (1 - cos(q[j+1,2])^2)))
    @NLexpression(model, δq₄R[j=1:N], -(𝓁 * m₂ * cos(q[j+1,2]) * sin(q[j+1,2]) * q[j+1,4]^2 + u[j,2] * cos(q[j+1,2]) + (m₁ + m₂) * g * sin(q[j+1,2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(q[j+1,2])^2)))

    ## System dynamics
    for i in 1:N
        j = i + 1  # index of latter knot
        ## Trapezoidal integration
        @NLconstraint(model, q[j,1] == q[i,1] + 0.5 * Δt[i] * (δq₁L[i] + δq₁R[i]))
        @NLconstraint(model, q[j,2] == q[i,2] + 0.5 * Δt[i] * (δq₂L[i] + δq₂R[i]))
        @NLconstraint(model, q[j,3] == q[i,3] + 0.5 * Δt[i] * (δq₃L[i] + δq₃R[i]))
        @NLconstraint(model, q[j,4] == q[i,4] + 0.5 * Δt[i] * (δq₄L[i] + δq₄R[i]))
    end
end

function add_constraints(model::Model, alg::HermiteSimpson{Continuous})
    N = alg.numInterval
    n = alg.numPoint
    
    ## Motion of the vehicle as a differential-algebraic system of equations (DAEs)
    Δt = model.obj_dict[:Δt]
    q = model.obj_dict[:q]
    u = model.obj_dict[:u]

    ## Motion of the vehicle as a differential-algebraic system of equations (DAEs)
    @NLexpression(model, δq₁[j=1:n], q[j,3])
    @NLexpression(model, δq₂[j=1:n], q[j,4])
    @NLexpression(model, δq₃[j=1:n], (𝓁 * m₂ * sin(q[j,2]) * q[j,4]^2 + u[j] + m₂ * g * cos(q[j,2]) * sin(q[j,2])) / (m₁ + m₂ * (1 - cos(q[j,2])^2)))
    @NLexpression(model, δq₄[j=1:n], -(𝓁 * m₂ * cos(q[j,2]) * sin(q[j,2]) * q[j,4]^2 + u[j] * cos(q[j,2]) + (m₁ + m₂) * g * sin(q[j,2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(q[j,2])^2)))

    ## System dynamics
    for m in 1:N
        j = 2 * m; i = j - 1; k = j + 1;
        ## interpolation constraints
        @NLconstraint(model, q[j,1] == 0.5 * (q[i,1] + q[k,1]) + 0.125 * Δt[m] * (δq₁[i] - δq₁[k]))
        @NLconstraint(model, q[j,2] == 0.5 * (q[i,2] + q[k,2]) + 0.125 * Δt[m] * (δq₂[i] - δq₂[k]))
        @NLconstraint(model, q[j,3] == 0.5 * (q[i,3] + q[k,3]) + 0.125 * Δt[m] * (δq₃[i] - δq₃[k]))
        @NLconstraint(model, q[j,4] == 0.5 * (q[i,4] + q[k,4]) + 0.125 * Δt[m] * (δq₄[i] - δq₄[k]))
        ## collocation constraints
        @NLconstraint(model, q[k,1] - q[i,1] == 1 / 6 * Δt[m] * (δq₁[i] + 4 * δq₁[j] + δq₁[k]))
        @NLconstraint(model, q[k,2] - q[i,2] == 1 / 6 * Δt[m] * (δq₂[i] + 4 * δq₂[j] + δq₂[k]))
        @NLconstraint(model, q[k,3] - q[i,3] == 1 / 6 * Δt[m] * (δq₃[i] + 4 * δq₃[j] + δq₃[k]))
        @NLconstraint(model, q[k,4] - q[i,4] == 1 / 6 * Δt[m] * (δq₄[i] + 4 * δq₄[j] + δq₄[k]))
    end
end

function add_constraints(model::Model, alg::HermiteSimpson{Discontinuous})
    N = alg.numInterval
    # n = alg.numPoint
    
    ## Motion of the vehicle as a differential-algebraic system of equations (DAEs)
    Δt = model.obj_dict[:Δt]
    q = model.obj_dict[:q]
    u = model.obj_dict[:u]

    ## Motion of the vehicle as a differential-algebraic system of equations (DAEs)
    @NLexpression(model, δq₁L[j=1:N], q[2*j-1,3])
    @NLexpression(model, δq₂L[j=1:N], q[2*j-1,4])
    @NLexpression(model, δq₃L[j=1:N], (𝓁 * m₂ * sin(q[2*j-1,2]) * q[2*j-1,4]^2 + u[j,1] + m₂ * g * cos(q[2*j-1,2]) * sin(q[2*j-1,2])) / (m₁ + m₂ * (1 - cos(q[2*j-1,2])^2)))
    @NLexpression(model, δq₄L[j=1:N], -(𝓁 * m₂ * cos(q[2*j-1,2]) * sin(q[2*j-1,2]) * q[2*j-1,4]^2 + u[j,1] * cos(q[2*j-1,2]) + (m₁ + m₂) * g * sin(q[2*j-1,2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(q[2*j-1,2])^2)))

    @NLexpression(model, δq₁M[j=1:N], q[2*j,3])
    @NLexpression(model, δq₂M[j=1:N], q[2*j,4])
    @NLexpression(model, δq₃M[j=1:N], (𝓁 * m₂ * sin(q[2*j,2]) * q[2*j,4]^2 + u[j,2] + m₂ * g * cos(q[2*j,2]) * sin(q[2*j,2])) / (m₁ + m₂ * (1 - cos(q[2*j,2])^2)))
    @NLexpression(model, δq₄M[j=1:N], -(𝓁 * m₂ * cos(q[2*j,2]) * sin(q[2*j,2]) * q[2*j,4]^2 + u[j,2] * cos(q[2*j,2]) + (m₁ + m₂) * g * sin(q[2*j,2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(q[2*j,2])^2)))

    @NLexpression(model, δq₁R[j=1:N], q[2*j+1,3])
    @NLexpression(model, δq₂R[j=1:N], q[2*j+1,4])
    @NLexpression(model, δq₃R[j=1:N], (𝓁 * m₂ * sin(q[2*j+1,2]) * q[2*j+1,4]^2 + u[j,3] + m₂ * g * cos(q[2*j+1,2]) * sin(q[2*j+1,2])) / (m₁ + m₂ * (1 - cos(q[2*j+1,2])^2)))
    @NLexpression(model, δq₄R[j=1:N], -(𝓁 * m₂ * cos(q[2*j+1,2]) * sin(q[2*j+1,2]) * q[2*j+1,4]^2 + u[j,3] * cos(q[2*j+1,2]) + (m₁ + m₂) * g * sin(q[2*j+1,2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(q[2*j+1,2])^2)))
    
    ## System dynamics
    for m in 1:N
        j = 2 * m; i = j - 1; k = j + 1;
        ## interpolation constraints
        @NLconstraint(model, q[j,1] == 0.5 * (q[i,1] + q[k,1]) + 0.125 * Δt[m] * (δq₁L[m] - δq₁R[m]))
        @NLconstraint(model, q[j,2] == 0.5 * (q[i,2] + q[k,2]) + 0.125 * Δt[m] * (δq₂L[m] - δq₂R[m]))
        @NLconstraint(model, q[j,3] == 0.5 * (q[i,3] + q[k,3]) + 0.125 * Δt[m] * (δq₃L[m] - δq₃R[m]))
        @NLconstraint(model, q[j,4] == 0.5 * (q[i,4] + q[k,4]) + 0.125 * Δt[m] * (δq₄L[m] - δq₄R[m]))
        ## collocation constraints
        @NLconstraint(model, q[k,1] - q[i,1] == 1 / 6 * Δt[m] * (δq₁L[m] + 4 * δq₁M[m] + δq₁R[m]))
        @NLconstraint(model, q[k,2] - q[i,2] == 1 / 6 * Δt[m] * (δq₂L[m] + 4 * δq₂M[m] + δq₂R[m]))
        @NLconstraint(model, q[k,3] - q[i,3] == 1 / 6 * Δt[m] * (δq₃L[m] + 4 * δq₃M[m] + δq₃R[m]))
        @NLconstraint(model, q[k,4] - q[i,4] == 1 / 6 * Δt[m] * (δq₄L[m] + 4 * δq₄M[m] + δq₄R[m]))
    end
end