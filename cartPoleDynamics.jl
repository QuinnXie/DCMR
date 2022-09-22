# dynamic functions
function CartPole!(dx,x,u,t)
    m₁ = 1.0      # mass of cart
    m₂ = 0.3      # mass of pole
    𝓁 = 0.5       # pole length
    g = 9.81      # grabity acceleration

    # dx = zeros(4);
    # ̇dynamics for δq₁
    dx[1] = x[3]
    # ̇dynamics for δq₂
    dx[2] = x[4]
    # ̇dynamics for δ²q₁
    dx[3] = (𝓁 * m₂ * sin(x[2]) * x[4]^2 + u(t) + m₂ * g * cos(x[2]) * sin(x[2])) / (m₁ + m₂ * (1 - cos(x[2])^2))
    # ̇dynamics for δ²q₂
    dx[4] = -(𝓁 * m₂ * cos(x[2]) * sin(x[2]) * x[4]^2 + u(t) * cos(x[2]) + (m₁ + m₂) * g * sin(x[2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(x[2])^2))

    return dx
end

# dynamic functions
function CartPole!(x,u)
    m₁ = 1.0      # mass of cart
    m₂ = 0.3      # mass of pole
    𝓁 = 0.5       # pole length
    g = 9.81      # grabity acceleration

    dx = zeros(size(x,1));
    # dx = x
    # ̇dynamics for δq₁
    dx[1] = x[3]
    # ̇dynamics for δq₂
    dx[2] = x[4]
    # ̇dynamics for δ²q₁
    dx[3] = (𝓁 * m₂ * sin(x[2]) * x[4]^2 + u + m₂ * g * cos(x[2]) * sin(x[2])) / (m₁ + m₂ * (1 - cos(x[2])^2))
    # ̇dynamics for δ²q₂
    dx[4] = - (𝓁 * m₂ * cos(x[2]) * sin(x[2]) * x[4]^2 + u * cos(x[2]) + (m₁ + m₂) * g * sin(x[2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(x[2])^2))

    return dx
end

# dynamic functions
function CartPole!(t,x,u)
    m₁ = 1.0      # mass of cart
    m₂ = 0.3      # mass of pole
    𝓁 = 0.5       # pole length
    g = 9.81      # grabity acceleration

    dx = zeros(size(x,1));
    # ̇dynamics for δq₁
    dx[1] = x[3]
    # ̇dynamics for δq₂
    dx[2] = x[4]
    # ̇dynamics for δ²q₁
    dx[3] = (𝓁 * m₂ * sin(x[2]) * x[4]^2 + u + m₂ * g * cos(x[2]) * sin(x[2])) / (m₁ + m₂ * (1 - cos(x[2])^2))
    # ̇dynamics for δ²q₂
    dx[4] = -(𝓁 * m₂ * cos(x[2]) * sin(x[2]) * x[4]^2 + u * cos(x[2]) + (m₁ + m₂) * g * sin(x[2])) / (𝓁 * m₁ + 𝓁 * m₂ * (1 - cos(x[2])^2))

    return dx
end
