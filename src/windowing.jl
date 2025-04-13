# Utility windowing functions, defined on the interval [0, 1]

# Spherical ball autocorrelation function. An autocorrelation function of the indicator 
# functions of two sherical balls of diameter 1, normalized to 1 at the origin.
function ball_autocorr(x::Float64)
    if x < 0.0
        return 0.0
    elseif x < 1.0
        return 0.5*(2+x)*(1-x)^2    
    else
        return 0.0
    end
end


