

gauss_2d(x, y, σ; x0=0.0, y0=0.0) = exp(-((x-x0)^2 + (y-y0)^2)/σ^2)
sphere_2d(x, y, r0; x0=0.0, y0=0.0) = (x-x0)^2 + (y-y0)^2 - r0^2

function two_gaussians_2d(x, y)
    w = 3.0
    gauss_2d(x, y, w, x0=2.0, y0=2.0) + 2*gauss_2d(x, y, w, x0=8.0, y0=7.0)
end

function two_gaussians_plus_negative_2d(x, y)
    w = 3.0
    two_gaussians_2d(x, y) - gauss_2d(x, y, w, x0=0.0, y0=9.0)
end