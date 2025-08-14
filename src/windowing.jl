# Utility windowing functions, defined on the interval [0, 1]

"""
    ball_autocorr(x::Float64)::Float64

Spherical ball autocorrelation function

Computes the autocorrelation function of the indicator functions of two 
three-dimensional spherical balls of diameter 1. The result is normalized 
to 1 at the origin.

# Arguments
- `x`: The distance from the center of the balls, normalized to the radius 
(which is 0.5).

# Returns
- `Float64`: The value of the autocorrelation function at distance `x`.
"""
function ball_autocorr(x::Float64)::Float64
    if x < 0.0
        return 0.0
    elseif x < 1.0
        return 0.5*(2+x)*(1-x)^2    
    else
        return 0.0
    end
end


