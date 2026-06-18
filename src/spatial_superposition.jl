using LinearAlgebra
using DSP: conv

"""
    successive flux(g)
    successive_flux(t, H, rb, D, ks, Cs, xy)

Iteratively solve spatial superposition for a borefield using the successive flux approach of
Nguyen and Pasquier (2021) to obtain the g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all 
boreholes (Type II). The g-function generated is for an impulse of 1 W/m.
# Arguments
    - `g`: A 3D g-function matrix for all radius of the borefield (nt x nb x nb) [°Cm/W]
        - Each time step (1 x nb x nb) has the borefield response for each radius between boreholes.
    - `t`: Time vector (nt x 1) [s]
    - `H`: Borehole depth (1x1) [m]
    - `rb`: Borehole radius (1x1) [m]
    - `D`: Borehole burried depth (1x1) [m]
    - `ks`: Soil thermal conductivity (1x1) [W/mK]
    - `Cs`: Soil volumetric specific heat (1x1) [J/m³K]
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
        - Can be computed with `borefield_xy()` from utils.jl.
# Output
    - `g`: g-function of the borefield spatial superposition [-]
# Reference
    - Nguyen, A., & Pasquier, P. (2021). A successive flux estimation method for rapid g-function 
        construction of small to large-scale ground heat exchanger. Renewable Energy, 165, 359–368. 
        https://doi.org/10.1016/j.renene.2020.10.074
"""
function successive_flux(g::AbstractArray{<:Real,3})
    # Basic parameters
    nt, nb1, nb2 = size(g)
    @assert nb1 == nb2 "g must be nt × nb × nb"
    @assert all(isfinite, g) "g must contain only finite values"
    nb = nb1

    # First estimation of g-function using block matrix (Eq. 20)
    GG = zeros(eltype(g), nt, nb + 1, nb + 1)
    @views GG[:, 1:nb, 1:nb] .= g
    @views GG[:, nb + 1, 1:nb] .= 1
    @views GG[:, 1:nb, nb + 1] .= 1

    b = zeros(eltype(g), nb + 1)
    b[end] = 1

    x = zeros(eltype(g), nt, nb)
    for it in 1:nt
        M = Matrix(@view GG[it, :, :])
        sol = try
            M \ b
        catch err
            if err isa SingularException
                @warn "Singular initial flux system at time index $it; using an SVD-based minimum-norm solve instead. This usually indicates a rank-deficient g slice caused by a degenerate geometry, duplicate boreholes, or a time step too close to zero."
                svd(M) \ b
            else
                rethrow()
            end
        end
        x[it, :] .= sol[1:nb]
    end

    # Successive flux estimation
    e1 = 10.0
    e2 = Inf
    e3 = Inf
    k  = 0
    kmax = 100
    gi = zeros(eltype(g), nt)

    while e3 > 0.15 && e1 > 1e-3 && e1 < e2 && k < kmax
        k += 1                                  # Iteration counter
        f = vcat(x[1:1, :], diff(x, dims=1))    # Step fluxes (Eq. 4)
        # Temperature response via pairwise convolutions (Eq. 6, 8)
        hh = zeros(eltype(g), nt, nb)
        for j in 1:nb, i in 1:nb
            hh[:, i] .+= conv(f[:, j], g[:, i, j])[1:nt]
        end
        gi = vec(sum(x .* hh, dims=2))          # Eq, 11 ĥ
        c = hh ./ gi .- 1                       # Eq. 12
        x .*= (1 .- c)                          # Eq. 16 (or 15?)
        # Convergence check
        err = maximum(abs, c)
        e3 = abs((e1 - err) / e1)
        e2 = e1
        e1 = err
    end
    return gi
end
function successive_flux(t, H, rb, D, ks, Cs, xy)
    @assert all(>(0), t) "t must be strictly positive to avoid a singular zero-time response"

    # Compute the radius matrix of the borefield
    r, _, _, _, _, _ = borefield_radius(xy, rb)

    # Compute the ground model for all different radius of the borefield
    g = fls(t, H, r, D, ks, Cs) # TODO: Add more models, or choices for the model used

    # Call the bloc matrix function
    return successive_flux(g)
end

"""
    bloc_matrix(g)
    bloc_matrix(t, H, rb, D, ks, Cs, xy)

Function that computes the spatial superposition of a borefield using the bloc matrix approach
of Dusseault et al. (2018) to obtain g-functions of a borefield. This approach assumes that heat
flux is uniform along all the borehole, and that the mean temperature is the same for all 
boreholes (Type II). The g-function generated is for an impulse of 1 W/m.
# Arguments
    - `g`: A 3D g-function matrix for all radius of the borefield (nt x nb x nb) [°Cm/W]
        - Each time step (1 x nb x nb) has the borefield response for each radius between boreholes.
    - `t`: Time vector (nt x 1) [s]
    - `H`: Borehole depth (1x1) [m]
    - `rb`: Borehole radius (1x1) [m]
    - `D`: Borehole burried depth (1x1) [m]
    - `ks`: Soil thermal conductivity (1x1) [W/mK]
    - `Cs`: Soil volumetric specific heat (1x1) [J/m³K]
    - `xy`: Matrix of borehole coordinates where the line source is at (0,0) (nr x 2) [m]
        - Can be computed with `borefield_xy()` from utils.jl.
# Output
    - `g`: g-function of the borefield spatial superposition [-]
# Reference
    - Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient 
        g-function construction. Renewable Energy, 121, 249–260. 
        https://doi.org/10.1016/j.renene.2017.12.092
"""
function bloc_matrix(gm::AbstractArray{<:Real})
    # Basic parameters
    nt, nb = size(gm)

    # Building the convolution matrix
    G = zeros(nt * nb, nt * nb)
    for i in 1:nt
        for ii in 1:nb
            for jj in 1:nb
                G[(ii-1)*nt+i:ii*nt, (jj-1)*nt+i] = gm[i:end, ii, jj]
            end
        end
    end

    # Create inputs to solve the linear system
    Gₕ = [[G; repeat(I(nt), 1, nb)] [repeat(I(nt), nb, 1); zeros(nt, nt)]]
    b = zeros((nb + 1) * nt)
    b[nb*nt+1] = 1

    # Solve the linear system
    sol = Gₕ \ b
    # Output the transfer function
    return -sol[nb*nt+1:end]
end
function bloc_matrix(t, H, rb, D, ks, Cs, xy)
    # Compute the radius matrix of the borefield
    r, _, _, _, _, _ = borefield_radius(xy, rb)

    # Compute the ground model for all different radius of the borefield
    g = fls(t, H, r, D, ks, Cs) # TODO: at some point add more models.

    # Call the bloc matrix function
    return bloc_matrix(g)
end