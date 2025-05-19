using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets
using Plots, Interpolations, LaTeXStrings, Dierckx, BenchmarkTools, LinearAlgebra
using CSV, DataFrames, DelimitedFiles, SavitzkyGolay, LsqFit, Printf, Peaks
using HEOM
const H = HEOM
using OpenQuantumSystems
const OQS = OpenQuantumSystems

function gaussian_mom(x, x0, p0, sigma, h_bar)
    """
    Standard Gaussian but with some momentum
    """
    return @. ((2.0 * pi * sigma^2)^-0.25) *
              exp(-(x - x0)^2.0 / (4 * sigma^2) + (1.0im * p0 * x / h_bar))
end

function v_harm(x, x0, mass, omega)
    return @. 0.5 * mass * omega^2 * (x - x0)^2
end

function v_morse(x, x0, a, v0)
    return @. v0 * (exp(-2.0 * a * (x - x0)) - 2.0 * exp(a * (x - x0)))
end


function prep_SE(x_vec, v_vec, mass, h_bar; apx=2)
    # Calculate the step size
    dx = x_vec[2] - x_vec[1]
    n = length(x_vec)

    # Define the derivative
    d2_dx2 = H.prepare_fd_dq(2, apx, dx, n)

    # Pre allocate the terms
    d2P_dx2 = Array{ComplexF64}(zeros(n))

    # Calculate the terms
    term_kin = @. im * h_bar / (2.0 * mass)
    term_pot = Array{ComplexF64}(@. -im / h_bar * v_vec)

    p = (term_kin=term_kin, term_pot=term_pot, d2_dx2=d2_dx2, d2P_dx2=d2P_dx2)
    return p
end

function SE_fd!(du, u, p, t)
    # Calculate d2P_dx2
    mul!(p.d2P_dx2, p.d2_dx2, u)
    # Main equation
    @. du = p.term_kin * p.d2P_dx2 + p.term_pot * u
end

function prep_VN(x_vec, v_vec, mass, h_bar; apx=2)
    # Calculate the step size
    dx = x_vec[2] - x_vec[1]
    n = length(x_vec)

    # Define the derivative
    d_dx = H.prepare_fd_dq(1, apx, dx, n)
    d2_dx2 = H.prepare_fd_dq(2, apx, dx, n)

    # Pre allocate the terms
    dP_dx = Array{ComplexF64}(zeros(n, n))
    dP_dy = Array{ComplexF64}(zeros(n, n))
    d2P_dx2 = Array{ComplexF64}(zeros(n, n))
    d2P_dy2 = Array{ComplexF64}(zeros(n, n))

    # Calculate the terms
    term_kin = @. im * h_bar / (2.0 * mass)
    term_pot = Array{ComplexF64}(@. -im / h_bar * v_vec)

    p = (
        term_kin=term_kin,
        term_pot=term_pot,
        d_dx=d_dx,
        d2_dx2=d2_dx2,
        dP_dx=dP_dx,
        dP_dy=dP_dy,
        d2P_dx2=d2P_dx2,
        d2P_dy2=d2P_dy2)
    return p
end

function VN_fd!(du, u, p, t)
    # Calculate d2P_dx2
    mul!(p.d2P_dx2, p.d2_dx2, u)
    # Calculate d2P_dy2
    mul!(p.d2P_dy2, u, p.d2_dx2)
    # Main equation
    @. du = p.term_kin * (p.d2P_dx2 - p.d2P_dy2) + p.term_pot * u
end

function prep_CL(x_vec, v_vec, mass, h_bar, gamma, k_b, temperature; apx=2)
    # Calculate the step size
    dx = x_vec[2] - x_vec[1]
    n = length(x_vec)

    # Define the derivative
    d_dx = H.prepare_fd_dq(1, apx, dx, n)
    d2_dx2 = H.prepare_fd_dq(2, apx, dx, n)

    # Pre allocate the terms
    dP_dx = Array{ComplexF64}(zeros(n, n))
    dP_dy = Array{ComplexF64}(zeros(n, n))
    d2P_dx2 = Array{ComplexF64}(zeros(n, n))
    d2P_dy2 = Array{ComplexF64}(zeros(n, n))

    # Calculate the terms
    term_kin = @. im * h_bar / (2.0 * mass)
    term_pot = Array{ComplexF64}(@. -im / h_bar * v_vec)
    term_diss = Array{ComplexF64}(@. -gamma * (x_vec - y_vec))
    term_deco = Array{ComplexF64}(@. -2.0 * mass * gamma * k_b * temperature / h_bar^2 * (x_vec - y_vec)^2)


    p = (
        term_kin=term_kin,
        term_pot=term_pot,
        term_diss=term_diss,
        term_deco=term_deco,
        d_dx=d_dx,
        d2_dx2=d2_dx2,
        dP_dx=dP_dx,
        dP_dy=dP_dy,
        d2P_dx2=d2P_dx2,
        d2P_dy2=d2P_dy2)
    return p
end

function CL_fd!(du, u, p, t)
    # Calculate dP_dx
    mul!(p.dP_dx, p.d_dx, u)
    # mul!(p.dP_dx, u, p.d_dx)
    # Calculate dP_dy
    mul!(p.dP_dy, u, p.d_dx)
    # mul!(p.dP_dy, p.d_dx, u)

    # Calculate d2P_dx2
    mul!(p.d2P_dx2, p.d2_dx2, u)
    # Calculate d2P_dy2
    mul!(p.d2P_dy2, u, p.d2_dx2)

    # Main equation
    @. du = p.term_kin * (p.d2P_dx2 - p.d2P_dy2) + p.term_pot * u + p.term_diss * (p.dP_dx - p.dP_dy) + p.term_deco * u
end

function plot_rho(x_vec, rho; title="", f_name="rho.pdf")
    fig = heatmap(
        x_vec,
        x_vec,
        real(rho),
        seriescolor=:inferno,
        aspect_ratio=:equal,
        title=title,
        xlims=(minimum(x_vec), maximum(x_vec)),
        ylims=(minimum(x_vec), maximum(x_vec)),
        right_margin=10.0Plots.mm)
    # save the figure
    savefig(fig, joinpath(plot_dump, f_name))
    H.os_display(fig)
end

function animate_rho(x_vec, sol; title="", f_name="rho.gif")
    u = real.(sol.u)
    time = sol.t

    # Get the min max range to fix the plot
    basis_low = minimum(x_vec)
    basis_high = maximum(x_vec)
    z_max = maximum([maximum(i[:,:]) for i in u])
    z_min = minimum([minimum(i[:,:]) for i in u])
    @time begin
        anim = @animate for i = 1:length(time)
            fig = heatmap(
                x_vec,
                x_vec,
                u[i],
                seriescolor=:inferno,
                aspect_ratio=:equal,
                title=title,
                xlims=(basis_low, basis_high),
                ylims=(basis_low, basis_high),
                zlims=(z_min, z_max),
                right_margin=10.0Plots.mm)
        end
    end
    p = gif(anim, joinpath(plot_dump, f_name), fps=15)
    H.os_display(p)
end

function calc_entropy(rho)
    # Calculate the von Neumann entropy 
    return -tr(rho * log(rho))
end

function calc_coherence(rho)
    # Calculate the coherence by summing over the off diagonal elements
    return sum(rho) - sum(diag(rho))
end

function plot_coherence(sol; f_name="coherence.pdf")
    # Calculate the coherence
    u = real.(sol.u)
    y = real.([calc_coherence(i) for i in u])
    y = real.([calc_entropy(i) for i in u])

    time_sim = sol.t
    time_sim, prefix, _ = H.best_time_units(time_sim)
    tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    ylab = latexstring("\\mathrm{Coherence}, C(t)")

    # Plot the coherence
    fig = H.plot_general(time_sim, y, tlab, ylab)
    savefig(fig, joinpath(plot_dump, f_name))
    H.os_display(fig)
end

function plot_trace(sol; f_name="trace.pdf")
    # Calculate the coherence
    u = real.(sol.u)
    y = real.([tr(i) for i in u])

    time_sim = sol.t
    time_sim, prefix, _ = H.best_time_units(time_sim)
    tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    ylab = latexstring("\\mathrm{Tr}, ρ(t)")

    # Plot the coherence
    fig = H.plot_general(time_sim, y, tlab, ylab)
    savefig(fig, joinpath(plot_dump, f_name))
    H.os_display(fig)
end

function plot_purity(sol; f_name="purity.pdf")
    # Calculate the coherence
    u = real.(sol.u)
    y = real.([tr(i^2) for i in u])

    time_sim = sol.t
    time_sim, prefix, _ = H.best_time_units(time_sim)
    tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    ylab = latexstring("\\mathrm{Tr}, ρ\\^{2}(t)")

    # Plot the coherence
    fig = H.plot_general(time_sim, y, tlab, ylab)
    savefig(fig, joinpath(plot_dump, f_name))
    H.os_display(fig)
end


if Sys.iswindows()
    global plot_dump =
        joinpath(homedir(), "OneDrive - University of Surrey/DUMP/")
else
    # Set to headless plotting
    ENV["GKSwstype"] = 100
    global plot_dump = pwd()
end

H.printout("Start!")
H.printout("plot_dump = $(plot_dump)")
f_plot_eigen = false

# Units
# https://en.wikipedia.org/wiki/Isotopes_of_sodium
mass_Na23 = H.au_m_dalt * 22.9897692820 # most common
mass_Na24 = H.au_m_dalt * 23.990963012
# https://en.wikipedia.org/wiki/Isotopes_of_potassium
mass_K39 = H.au_m_dalt * 38.963706487 # most common
mass_K40 = H.au_m_dalt * 39.96399817
# https://en.wikipedia.org/wiki/Isotopes_of_lithium
mass_Li6 = H.au_m_dalt * 6.0151228874
mass_Li7 = H.au_m_dalt * 7.0160034344 # most common


# Define the terms
h_bar = H.h_bar
k_b = H.k_b
temperature = 310.15 # 298.15 K
beta = 1.0 / (k_b * temperature)
gamma = (3500 / H.hart_2_inv_cm)
mass = mass_Na23

# Solver settings
tol = 1e-18
apx = 2
algo = VCABM5() #ROCK4() #Tsit5() SSPRK83 VCABM5 Feagin12 Feagin14 Vern9

# Time grid
t0 = 0.0
t1 = 0.001 / H.aut_2_ps # 1 ps
h_step = 5.0 # 1.0
nt = 50
saveat = H.linspace(t0, t1, nt)

# Create a grid
xl = -0.5 # -1.2
xh = 0.5 # 3.5
n = Int(2^8)
x_vec = range(start=xl, stop=xh, length=n) |> Array
y_vec = x_vec'

# Define the potential
x0 = 0.0
omega = gamma
v_vec = v_harm(x_vec, x0, mass, omega)

# set the minimum to zero
v_vec = v_vec .- minimum(v_vec)

# Free particle
#v_vec = zeros(n)

# plot the potential
# fig = plot(x_vec, real(v_vec), label="Potential", xlabel="x", ylabel="V(x)", legend=:topleft)
# H.os_display(fig)

# Solve for the eigenstates
if f_plot_eigen
    e_vals, e_vecs = OQS.eigen_solve(
        x_vec,
        v_vec;
        n_eig=5,
        h_bar=h_bar,
        mass=mass
    )
    e_vecs = OQS.eigen_parity_fixer(x_vec, v_vec, e_vals, e_vecs)

    OQS.plot_eigen(
        x_vec,
        v_vec,
        e_vals,
        e_vecs,
        (-0.0001, 0.4); # -0.0001, 0.004   minimum(v_raw), maximum(v_raw)*0.2
        pmult=0.05, # 0.003
        name=string("pot_eigen_.pdf"),
        f_app="W",
        legend=:outerright,# :topright :bottomright, :outerright
        f_pyplt=false
    )
    @printf "Lowest eigenstate= %.2E Eh \n" e_vals[1]
end

# Calculate the step size
dx = x_vec[2] - x_vec[1]
dy = y_vec[2] - y_vec[1]

P0 = gaussian_mom(x_vec, 0.2, 0.0, 0.05, h_bar) #+ gaussian_mom(x_vec, -0.2, 0.0, 0.05, h_bar)

rho = P0 * P0'
# Normalise
rho = rho ./ tr(rho)

# plot_rho(x_vec, rho)

prep = prep_CL(x_vec, v_vec, mass, h_bar, gamma, k_b, temperature)
prob = ODEProblem(CL_fd!, rho, (t0, t1), prep)
H.printout("Running the simulation...")
@time sol = solve(prob, algo, progress=true, saveat=saveat, abstol=tol, reltol=tol)
sol = sol(saveat)
# plot_rho(x_vec, sol[end])
animate_rho(x_vec, sol)
plot_coherence(sol)
plot_trace(sol)
plot_purity(sol)

H.printout("End!")