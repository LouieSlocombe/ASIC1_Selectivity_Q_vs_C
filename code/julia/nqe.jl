using OrdinaryDiffEq, ModelingToolkit, MethodOfLines, DomainSets
using Plots, Interpolations, LaTeXStrings, Dierckx, BenchmarkTools, LinearAlgebra
using CSV, DataFrames, DelimitedFiles, SavitzkyGolay, LsqFit, Printf, Peaks
using Unitful, UnitfulAtomic
using HEOM
const H = HEOM

function morse_fit(x, p)
    return @. p[1] * (1.0 - exp(-p[2] * (x - p[3])))^2
end

function load_ion_data(file)
    # Read the data
    data = CSV.read(file, DataFrame; delim="	", comment="#", skipto=19, footerskip=1) |> Tables.matrix
    # Extract the position and potential data
    q_raw = data[:, 1]
    v_raw = data[:, 2]
    # Move the minimum to zero
    v_raw = v_raw .- minimum(v_raw)
    # Convert to atomic units
    q_raw = q_raw .* H.si_ns ./ H.au_length
    v_raw = v_raw .* H.kcalmol_2_hart

    return q_raw, v_raw
end

function v_ion_asic(q_vec, file; f_plot=true, idx_l=25, idx_r=15)
    q_raw, v_raw = load_ion_data(file)

    @printf "Lowest raw q= %.2E bohr \n" q_raw[1]
    @printf "Highest raw q= %.2E bohr \n" q_raw[end]

    # # Add 10 points to the right
    # n_r = 2 
    # dq = q_raw[end] - q_raw[end-1]
    # q_raw = vcat(q_raw, q_raw[end] .+ ones(n_r) .*dq)
    # H.printout("q_raw = $(q_raw)")
    # v_raw = vcat(v_raw, v_raw[end] .* ones(n_r))
    # H.printout("v_raw = $(v_raw)")

    # # Filtering
    sg = savitzky_golay(v_raw, 21, 10) # 21, 10
    v_filt = sg.y
    #v_filt = v_raw

    # flatten out the last section of the potential
    v_filt[end-idx_r:end] .= v_filt[end-idx_r]

    q_filt = q_raw

    # Filtering
    # sg = savitzky_golay(v_filt, 9, 3)
    # v_filt = sg.y

    # Interpolate the potential but do not include some initial points
    spl = Spline1D(q_filt[idx_l:end], v_filt[idx_l:end]; k=5, bc="extrapolate")

    # Create a grid
    N = Int(2^9)
    ql = minimum(q_filt .- 0.1)
    qh = maximum(q_filt .+ 0.4)
    pl = -30.0
    ph = 30.0

    q_int, _, _, _, _, _ = H.create_basis(N, ql, qh, pl, ph)
    v_int = spl(q_int)

    chop = findfirst(q_int .>= maximum(q_raw))
    # flatten out the last section of the potential
    v_int[chop:end] .= v_int[chop]

    if f_plot
        y_min = minimum(v_raw .* 1.1)
        y_max = maximum(v_raw .* 1.1)
        q_lab = latexstring("Q \\: \\left[\\alpha_0 \\right]")
        pot_lab = latexstring(
            "\\mathrm{Potential}, \\, V, \\: \\left[E_{\\mathrm{h}} \\right]",
        )
        p = scatter(q_raw, v_raw, lw=2, label="Raw")
        p = plot!(p, q_filt, v_filt, lw=2, label="Filt.")
        p = plot!(p, q_int, v_int, lw=2, label="Int.")
        p = plot!(p,
            xlabel=q_lab,
            ylabel=pot_lab,
            ylims=(y_min, y_max),
            legend=:bottomright,
        )
        savefig(p, joinpath(plot_dump, "compare.pdf"))
        H.os_display(p)
    end
    # Interpolate the potential
    spl = Spline1D(q_int, v_int; k=5, bc="extrapolate")
    v_vec = spl(q_vec)
    v_vec = v_vec .- minimum(v_vec)
    return v_vec, q_raw, v_raw
end

function v_fit()
    path = joinpath(homedir(), "OneDrive - University of Surrey\\Papers\\paper_ion_channel\\data\\latest")
    file = "cedric_asicsf_sod_hw.xvg"
    file = "cedric_asicsf_lit_hw.xvg"

    q_raw, v_raw = load_ion_data(path, file)
    n_l = 20
    n_h = 5
    q_raw = q_raw[n_l:end-n_h]
    v_raw = v_raw[n_l:end-n_h]


    p0 = [v_raw[end], 0.5, 0.0]
    lb = [v_raw[end] * 0.9, -Inf, -Inf]
    ub = [v_raw[end] * 1.1, Inf, Inf]
    H.printout("p0 = $(p0)")
    fit = curve_fit(morse_fit, q_raw, v_raw, p0, lower=lb, upper=ub)
    param = fit.param
    H.printout("param = $(param)")
    v_fit = morse_fit(q_raw, param)

    y_min = minimum(v_raw .* 1.1)
    y_max = maximum(v_raw .* 1.1)
    q_lab = latexstring("Q \\: \\left[\\alpha_0 \\right]")
    pot_lab = latexstring(
        "\\mathrm{Potential}, \\, V, \\: \\left[E_{\\mathrm{h}} \\right]",
    )
    p = scatter(q_raw, v_raw, lw=2, label="Raw")
    p = plot!(p, q_raw, v_fit, lw=2, label="fit")
    p = plot!(p,
        xlabel=q_lab,
        ylabel=pot_lab,
        ylims=(y_min, y_max),
        legend=:bottomright,
    )
    savefig(p, joinpath(plot_dump, "compare.pdf"))
    H.os_display(p)
end


function reduced_mass(m1, m2)
    return (m1 * m2) / (m1 + m2)
end

function water_mass(f_deuteration=false)
    mo = H.au_m_dalt * 15.999
    if f_deuteration
        mh = H.au_m_dalt * 2.014101777844
    else
        mh = H.au_m_dalt * 1.007825
    end
    return mo + 2 * mh
end

function find_local_maxima(vector::Vector, min_prominence)
    n = length(vector)
    local_maxima = []
    
    for i in 2:n-1
        if vector[i] > vector[i-1] && vector[i] > vector[i+1]
            left_min = minimum(vector[1:i-1])
            right_min = minimum(vector[i+1:end])
            prominence = vector[i] - min(left_min, right_min)
            
            if prominence >= min_prominence
                push!(local_maxima, i)
            end
        end
    end
    
    return local_maxima
end

function gaussian_mom(x, x0, p0, sigma, h_bar)
    """
    Standard Gaussian but with some momentum
    """
    # return @. ((2.0 * pi * sigma^2)^-0.25) *
    #           exp(-(x - x0)^2.0 / (4 * sigma^2) + (1.0im * p0 * x / h_bar))
    return @. ((2.0 * pi * sigma^2)^-0.25) * exp(-1.0*(x - x0)^2.0 / (4 * sigma^2)) *exp(1.0im * p0 * x / h_bar)
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

    p = (
        term_kin=term_kin, 
        term_pot=term_pot, 
        d2_dx2=d2_dx2, 
        d2P_dx2=d2P_dx2)
    return p
end

function prob_dens_SE(P)
    return real.(P.* conj.(P))
end

function norm_SE(q_vec, P; k=5)
    return H.int_1d(q_vec, prob_dens_SE(P); k=k)^0.5
end

function normalise_SE(q_vec, P; k=5)
    return P = P ./ norm_SE(q_vec, P;k=k)
end

function SE_fd!(du, u, p, t)
    # Calculate d2P_dx2
    mul!(p.d2P_dx2, p.d2_dx2, u)
    # Main equation
    @. du = p.term_kin * p.d2P_dx2 + p.term_pot * u
end

function animate_SE(q_vec, solu, time; name="Pq.gif", dir=nothing)
    # Get the time loop
    nt = length(time)
    y_vals = [prob_dens_SE(solu[i]) for i = 1:nt]

    # Get the maximum value for plotting
    y_min = minimum([minimum(i) for i in y_vals[:][:]])
    y_max = maximum([maximum(i) for i in y_vals[:][:]])
    @time begin
        # Loop over time
        anim = @animate for i = 1:nt
            plot(
                q_vec,
                y_vals[i],
                ylims=(y_min, y_max),
                lw=3,
                xlabel=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
                ylabel="P(q)",
                legend=false,
                linecolor=:black,
            )
        end
        if dir === nothing
            fig = gif(anim, joinpath(plot_dump, name), fps=60)
        else
            fig = gif(anim, joinpath(dir, name), fps=60)
        end
    end
    H.os_display(fig)
    return fig
end

function animate_SE_log(
    q_vec,
    solu,
    time;
    name="Pq_log.gif",
    max_y_min=-20,
    dir=nothing
)
    # Get the time loop
    nt = length(time)
    y_vals = [prob_dens_SE(solu[i]) for i = 1:nt]
    y_vals = [log10.(abs.(y_vals[i])) for i = 1:nt]
    replace!(y_vals, Inf => NaN)
    # Get the maximum value for plotting
    y_min = minimum([minimum(i) for i in y_vals[:][:]])
    y_max = maximum([maximum(i) for i in y_vals[:][:]])
    # Prevent excessivly small values plotted
    if y_min < max_y_min
        y_min = max_y_min
    end
    @time begin
        # Loop over time
        anim = @animate for i = 1:nt
            plot(
                q_vec,
                y_vals[i],
                ylims=(y_min, y_max),
                lw=3,
                xlabel=latexstring("Q \\: \\left[\\alpha_0 \\right]"),
                ylabel="log P(q)",
                legend=false,
                linecolor=:black,
            )
        end
        if dir === nothing
            fig = gif(anim, joinpath(plot_dump, name), fps=60)
        else
            fig = gif(anim, joinpath(dir, name), fps=60)
        end
    end
    H.os_display(fig)
    return fig
end

function plot_SE_prob_dens(q_vec, P; name="Pq.pdf", dir=nothing)
    P = prob_dens_SE(P)
    fig = H.plot_general(
        q_vec,
        P,
        latexstring("Q \\: \\left[\\alpha_0 \\right]"),
        "P(q)",
    )
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    H.os_display(fig)
    return fig
end

function calc_SE_probability(q_vec, q0, P; k=5, f_renorm=false)
    """
    Calculates the probability of finding the particle in a region of space
    Prob. = ∫ P(q) h_q0(q) dq
    """
    if f_renorm
        P = normalise_SE(q_vec, P; k=k)
    end
    # get probability density
    pd = prob_dens_SE(P)
    # Calculate the step function
    h_q0 = H.step_function(q_vec, pd, q0)
    # Calculate the probability
    return H.int_1d(q_vec, h_q0; k=k)
end

function plot_SE_probability(
    q_vec,
    q0,
    solu,
    time;
    name="SE_probability.pdf",
    f_units="SI",
    dir=nothing)

    prob = [calc_SE_probability(q_vec, q0, solu[i]) for i = 1:length(time)]

    if f_units == "SI"
        time, prefix, _ = H.best_time_units(time)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end
    
    fig = H.plot_general(
        time,
        prob,
        tlab,
        "Occupation probability",
    )
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    H.os_display(fig)
    return fig
end

function plot_QSE_normalisation_log(
    q_vec,
    solu,
    time;
    name="SE_deviation_norm_log.pdf",
    f_units="SI",
    dir=nothing
)
    # Track the deviation from normalisation conditon over time
    norm = [norm_SE(q_vec, solu[i]) for i = 1:length(time)]

    if f_units == "SI"
        time, prefix, _ = H.best_time_units(time)
        tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
    else
        tlab = latexstring("\\mathrm{Time}, t, [AUT]")
    end
    ylab = latexstring("\\log \\: L_{2} \\; \\mathrm{Norm. \\: error}")

    # Loop over the time series
    y_vals = log10.([abs.(i - 1.0) for i in norm])
    replace!(y_vals, Inf => NaN)
    fig = H.plot_general(time, y_vals, tlab, ylab)
    if dir === nothing
        savefig(fig, joinpath(plot_dump, name))
    else
        savefig(fig, joinpath(dir, name))
    end
    H.os_display(fig)
    return fig
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
if Sys.iswindows()
    using OpenQuantumSystems
    const OQS = OpenQuantumSystems
    import PyPlot
    const plt = PyPlot
    rcParams = plt.PyPlot.PyDict(plt.PyPlot.matplotlib."rcParams")
    rcParams["axes.linewidth"] = 2.0

    function pyplot_nice_plotter_helper(
        xlab,
        ylab;
        xs=14,
        ys=14,
        f_legend=true,
        f_leg_loc="best",
        f_leg_size=14
    )
        if f_legend
            plt.legend(loc=f_leg_loc, fontsize=f_leg_size)
        end
        plt.minorticks_on()
        plt.tick_params(
            axis="both",
            which="major",
            labelsize=ys - 2,
            direction="in",
            length=6,
            width=2,
        )
        plt.tick_params(
            axis="both",
            which="minor",
            labelsize=ys - 2,
            direction="in",
            length=4,
            width=2,
        )
        plt.tick_params(axis="both", which="both", top=true, right=true)
        plt.xlabel(xlab, fontsize=xs)
        plt.ylabel(ylab, fontsize=ys)
        plt.axis("tight")
        plt.tight_layout()
        return nothing
    end

    function plot_wigner_s2_entropy_py(
        q,
        p,
        sol,
        plot_dump;
        name="s2_entropy.pdf",
        f_units="SI"
    )
        # Get the time
        time_sim = sol.t
        if f_units == "SI"
            time_sim, prefix, _ = H.best_time_units(time_sim)
            tlab = latexstring("\\mathrm{Time}, t, [$(prefix)s]")
        else
            tlab = latexstring("\\mathrm{Time}, t, [AUT]")
        end

        s2_entropy = [H.calc_wigner_s2_entropy(q, p, sol[i]) for i = 1:length(time_sim)]
        ylab = latexstring("\\mathrm{Entropy} \\, S_{2}")
        plt.plot(time_sim, s2_entropy, linewidth=2, color="black")
        pyplot_nice_plotter_helper(tlab, ylab; f_legend=false)
        plt.savefig(joinpath(plot_dump, name))
        plt.close()
    end
end



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

h_bar = H.h_bar
k_b = H.k_b
temperature = 310.15 # 298.15 K
beta = 1.0 / (k_b * temperature)

# Solver settings
tol = 1e-18
apx = 2
algo = VCABM5() #ROCK4() #Tsit5() SSPRK83 VCABM5 Feagin12 Feagin14 Vern9

# Time grid
t0 = 0.0
t_deco = 5.0 / H.aut_2_fs # 10 fs
# qrate
t_rate = 0.5 / H.aut_2_ps # 1 ps # 0.2 ps
H.printout("t_rate = $(round(t_rate;sigdigits=3)) AUT")
h_step = 5.0 # 1.0
nt = 200

f_plot_hopping = false
f_plot_pot_convert = false
f_eigen = false
f_get_deco = false
f_get_qrate = false
f_qse = false
f_scatter= true
f_scatter_loop = false
f_oqs = false
f_loc = "escape"# "escape" 1

# Na+(H2O)4, for K it's K+(H2O)3 and Li it's Li+(H2O)3
# Single
# 1 = sod_pmf_single
# 2 = sod_hw_pmf_single
# 3 = pot_pmf_single
# 4 = pot_hw_pmf_single
# 5 = lit_pmf_single
# 6 = lit_hw_pmf_single
# Multi
# 7 = sod_pmf_multi
# 8 = sod_hw_pmf_multi
# 9 = pot_pmf_multi
# 10 = pot_hw_pmf_multi
# 11 = lit_pmf_multi
# 12 = lit_hw_pmf_multi
f_choice = 1
N = Int(2^9)
N_qse = Int(2^9)

if f_choice == 1

    file = "sod_pmf_single.xvg"
    idx_l = 10
    idx_r = 1

    # Create a grid
    ql = -10.0 # -1.2
    qh = 20.0 # 3.5
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_Na23, 4 * water_mass(false))
    mass_t = mass_Na23 + 4 * water_mass(false)
    gamma = (3500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 2
    file = "sod_hw_pmf_single.xvg"
    idx_l = 3
    idx_r = 8
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_Na23, 4 * water_mass(true))
    mass_t = mass_Na23 + 4 * water_mass(true)
    gamma = (2500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 3
    file = "pot_pmf_single.xvg"
    idx_l = 3
    idx_r = 2
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_K39, 3 * water_mass(false))
    mass_t = mass_K39 + 3 * water_mass(false)
    gamma = (3500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 4
    file = "pot_hw_pmf_single.xvg"
    idx_l = 6
    idx_r = 1
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_K39, 3 * water_mass(true))
    mass_t = mass_K39 + 3 * water_mass(true)
    gamma = (2500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 5
    file = "lit_pmf_single.xvg"
    idx_l = 5
    idx_r = 6
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_Li7, 3 * water_mass(false))
    mass_t = mass_Li7 + 3 * water_mass(false)
    gamma = (3500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 6
    file = "lit_hw_pmf_single.xvg"
    idx_l = 3
    idx_r = 6
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_Li7, 3 * water_mass(true))
    mass_t = mass_Li7 + 3 * water_mass(true)
    gamma = (2500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 7
    file = "sod_pmf_multi.xvg"
    idx_l = 3
    idx_r = 8

    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_Na23, 4 * water_mass(false))
    mass_t = mass_Na23 + 4 * water_mass(false)
    gamma = (3500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 8
    file = "sod_hw_pmf_multi.xvg"
    idx_l = 1
    idx_r = 14
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_Na23, 4 * water_mass(true))
    mass_t = mass_Na23 + 4 * water_mass(true)
    gamma = (2500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 9
    file = "pot_pmf_multi.xvg"
    idx_l = 4 # 10
    idx_r = 2
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_K39, 3 * water_mass(false))
    mass_t = mass_K39 + 3 * water_mass(false)
    gamma = (3500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 10
    file = "pot_hw_pmf_multi.xvg"
    idx_l = 16
    idx_r = 1
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_K39, 3 * water_mass(true))
    mass_t = mass_K39 + 3 * water_mass(true)
    gamma = (2500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 11
    file = "lit_pmf_multi.xvg"
    idx_l = 3
    idx_r = 1
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_Li7, 3 * water_mass(false))
    mass_t = mass_Li7 + 3 * water_mass(false)
    gamma = (3500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
elseif f_choice == 12
    file = "lit_hw_pmf_multi.xvg"
    idx_l = 2
    idx_r = 4
    # Create a grid
    ql = -10.0
    qh = 20.0
    pl = -50.0
    ph = 50.0
    mass = reduced_mass(mass_Li7, 3 * water_mass(true))
    mass_t = mass_Li7 + 3 * water_mass(true)
    gamma = (2500 / H.hart_2_inv_cm) # H = 3500 cm-1, D = 2500 cm-1
end
H.printout("file = $(file)")
H.printout("gamma=$(gamma)")
if Sys.iswindows()
    path = joinpath(homedir(), "OneDrive - University of Surrey\\Papers\\paper_ion_channel\\data\\v02")
else
    path = joinpath(homedir(), "/mnt/lustre/a2fs-work3/work/e89/e89/louie/paper_ion_channel_data/v02")
end

q_vec, p_vec, Q_vec, P_vec, dq, dp = H.create_basis(N, ql, qh, pl, ph)

v_vec, q_raw, v_raw = v_ion_asic(q_vec, joinpath(path, file); f_plot=f_plot_pot_convert, idx_l=idx_l, idx_r=idx_r)
v_vec = H.trunc_pot(v_vec)

# H.os_display(plot(q_vec,v_vec, xlabel="Q", ylabel="V"))

# Set all values after the end of the raw data to the last value to flatten
idx_r = findfirst(q_vec .>= maximum(q_raw)) - 4
v_vec[idx_r:end] .= v_vec[idx_r]

H.printout("mass = $(round(mass/H.au_m_dalt;sigdigits=3)) Dh")
H.printout("mass_t = $(round(mass_t/H.au_m_dalt;sigdigits=3)) Dh")
# @printf "mass = %.2E Dh \n" mass/H.au_m_dalt

# Get the properties of the potential
loc_min = findmin(v_vec)[2]
dv1 = maximum(v_vec[loc_min:end]) - minimum(v_vec)
H.printout("dv1 = $(dv1)")
@printf "Forward barrier= %.2E Eh \n" dv1
@printf "Forward barrier= %.2E eV \n" dv1 * H.ev_2_hart

k_clas = H.calc_classical_rate(dv1, temperature)

@printf "k_clas = %.2E 1/s \n" k_clas
@printf "I_clas = %.2E q/s \n" k_clas .* 1.602e-19

loc, vals = findmaxima(v_vec, 5, strict=true)
# remove any past idx_r
loc = loc[loc .<= findfirst(q_vec .>= maximum(q_raw.-5.0))]

# get the prominence
loc, proms = peakproms(loc, v_vec)
# filter by prominence
min_prom = 1e-4#1.04E-04#0.0002
min_prom = H.determine_best_local_spring(mass, v_vec, q_vec) ./2
loc = loc[proms .>= min_prom]
# H.printout("proms = $(proms)")
# H.printout("loc = $(loc)")
# H.printout("q = $(q_vec[loc])")
# H.printout("There are $(length(loc)+1) likely states")

if f_plot_hopping
    fig = plot(
        q_vec,
        v_vec,
        lw=2,
        legend=false,
        linecolor=:black,
        xlabel="Q",
        ylabel="V",
        left_margin=2Plots.mm
    )
    # plot the minima if there is one
    fig = scatter!(fig, q_vec[loc], v_vec[loc], markersize=5, markercolor=:red)
    fig = plot!(fig, ylims=(minimum(v_raw), maximum(v_raw) * 1.1))
    # savefig(fig, joinpath(plot_dump, "pot" * string(f_choice)))
    H.os_display(fig)
end

# Pick the hopping rate location 
if f_loc == "escape"
    h_step = 10.0
else
    h_step = q_vec[loc[Int(f_loc)]]
end

if f_eigen
    if Sys.iswindows()
        # Solve for the eigenstates
        e_vals, e_vecs = OQS.eigen_solve(
            q_vec,
            v_vec;
            n_eig=15,
            h_bar=h_bar,
            mass=mass
        )
        e_vecs = OQS.eigen_parity_fixer(q_vec, v_vec, e_vals, e_vecs)
        # select only the values before the end of q_raw
        loc_h = findfirst(q_vec .>= 12.0)
        loc_l = findfirst(q_vec .>= -7.0)
        #loc = length(q_vec)
        q_plot = q_vec[loc_l:loc_h]
        v_plot = v_vec[loc_l:loc_h]
        e_vecs = e_vecs[loc_l:loc_h, :]

        OQS.plot_eigen(
            q_plot,
            v_plot,
            e_vals,
            e_vecs,
            (-0.0001, 0.008); # -0.0001, 0.004   minimum(v_raw), maximum(v_raw)*0.2
            pmult=0.0005, # 0.003
            name=string("pot_eigen_", f_choice, ".pdf"),
            f_app="W",
            legend=:outerright,# :topright :bottomright, :outerright
            f_pyplt=true
        )
        @printf "Lowest eigenstate= %.2E Eh \n" e_vals[1]
    end
end

# Calculate the thermal state
W_thermal = H.w0_thermal(q_vec, p_vec, v_vec, mass, beta)

# Calculate the properties of the thermal state
W_ther_pure = H.calc_wigner_purity(q_vec, p_vec, W_thermal)
W_ther_s2 = H.calc_wigner_s2_entropy(q_vec, p_vec, W_thermal)
W_ther_ener = H.calc_wigner_energy_expect(q_vec, p_vec, v_vec, mass, W_thermal)
W_ther_disp, _ = H.calc_wigner_wqp2_expect(q_vec, p_vec, W_thermal)

# H.printout("Thermal state properties:")
# H.printout("Purity = $(round(W_ther_pure;sigdigits=3))")
# H.printout("Entropy = $(round(W_ther_s2;sigdigits=3))")
# H.printout("Energy = $(round(W_ther_ener;sigdigits=3))")
# H.printout("Dispersion = $(round(W_ther_disp;sigdigits=3)) bohr^2")

# H.printout("Dispersion = $(round(W_ther_disp .* H.au_length^2 ./ H.si_angstrom^2;sigdigits=3)) A^2")

# H.plot_wigner_heatmap(q_vec, p_vec, W_thermal; dir=plot_dump, name="w_thermal.pdf")
# H.plot_wigner_wq(q_vec, p_vec, W_thermal; yscale=:log10, dir=plot_dump, name="wq_thermal.pdf")
# H.plot_wigner_wp(q_vec, p_vec, W_thermal; yscale=:log10, dir=plot_dump, name="wp_thermal.pdf")

if f_get_deco
    # Get the peak of the wavepacket
    W_thermal_wq, _ = H.calc_wigner_wqp(q_vec, p_vec, W_thermal)
    loc_pk = findmax(W_thermal_wq)[2]
    q0 = q_vec[loc_pk]
    sigma = 10.0#5.5 #H.calculate_fwhm(q_vec, psi0)*10.0
    # H.printout("q0 = $(round(q0;sigdigits=3))")
    # H.printout("sigma = $(round(sigma;sigdigits=3))")
    p0 = 0.0
    # Make intial gaussian
    W0 = H.w0_gaussian(q_vec, p_vec, q0, p0, sigma)
    # H.plot_wigner_heatmap(q_vec, p_vec, W0; dir=plot_dump, name="w0.pdf")
    # H.plot_wigner_wq(q_vec, p_vec, W0; yscale=:log10, dir=plot_dump, name="wq0.pdf")
    # H.plot_wigner_wp(q_vec, p_vec, W0; yscale=:log10, dir=plot_dump, name="wp0.pdf")


    # Prepare the time evolution
    saveat = H.linspace(t0, t_deco, nt)
    prep = H.prep_LL_HT_M_fd(q_vec, p_vec, v_vec, mass, h_bar, gamma, beta)
    prob = ODEProblem(H.LL_HT_M_fd_trunc!, W0, (t0, t_deco), prep)
    H.printout("Running the simulation...")
    @time sol = solve(prob, algo, progress=true, saveat=saveat, abstol=tol, reltol=tol)

    H.plot_wigner_normalisation(q_vec, p_vec, sol; yscale=:log10, dir=plot_dump)
    H.plot_wigner_purity(q_vec, p_vec, sol; dir=plot_dump)
    # H.plot_wigner_energy_expect(q_vec, p_vec, mass, v_vec, sol; dir=plot_dump)
    # H.plot_wigner_uncertainty_principle(q_vec, p_vec, sol; dir=plot_dump)

    if Sys.iswindows()
        plot_wigner_s2_entropy_py(q_vec, p_vec, sol, plot_dump)
    else
        H.plot_wigner_s2_entropy(q_vec, p_vec, sol; dir=plot_dump)
    end

    # Calculate the s2 entropy
    s2 = [H.calc_wigner_s2_entropy(q_vec, p_vec, sol[i]) for i = 1:nt]

    # Find the time to reach half the entropy of the thermal state
    idx = findfirst(s2 .>= W_ther_s2 / exp(1.0))
    idx = findfirst(s2 .>= maximum(s2) / exp(1.0))
    t_s2_half = saveat[idx]
    H.printout("Half life = $(round(t_s2_half;sigdigits=3)) AUT")
    H.printout("Half life = $(round(t_s2_half .* H.au_time ./ H.si_fs;sigdigits=3)) fs")

    # Calculate the properties of the final state
    H.printout("Final state properties:")
    H.printout("Purity = $(round(H.calc_wigner_purity(q_vec, p_vec, sol[end]);sigdigits=3))")
    H.printout("Entropy = $(round(H.calc_wigner_s2_entropy(q_vec, p_vec, sol[end]);sigdigits=3))")
    H.printout("Energy = $(round(H.calc_wigner_energy_expect(q_vec, p_vec, v_vec, mass, sol[end]);sigdigits=3))")
    # H.animate_wigner_heatmap(q_vec, p_vec, sol; dir=plot_dump)
    # H.animate_wigner_wq(q_vec, p_vec, sol; yscale=:log10, dir=plot_dump)
    # H.animate_wigner_wp(q_vec, p_vec, sol; yscale=:log10, dir=plot_dump)
end

if f_get_qrate
    W0 = H.w0_thermal_step(q_vec, p_vec, v_vec, mass, beta, h_step)
    H.plot_wigner_heatmap(q_vec, p_vec, W0; dir=plot_dump, name="w0.pdf")
    # H.plot_wigner_wq(q_vec, p_vec, W0; yscale=:log10, dir=plot_dump, name="w0_q.pdf")
    # H.plot_wigner_wp(q_vec, p_vec, W0; yscale=:log10, dir=plot_dump, name="w0_p.pdf")

    # Prepare the time evolution
    saveat = H.linspace(t0, t_rate, nt)
    prep = H.prep_LL_HT_M_fd(q_vec, p_vec, v_vec, mass, h_bar, gamma, beta)
    prob = ODEProblem(H.LL_HT_M_fd_trunc!, W0, (t0, t_rate), prep)
    H.printout("Running the simulation...")
    @time sol = solve(prob, algo, progress=true, saveat=saveat, abstol=tol, reltol=tol)

    H.plot_wigner_normalisation(q_vec, p_vec, sol; yscale=:log10, dir=plot_dump)
    # Calculate the properties of the final state
    H.plot_wigner_step_occupation(q_vec, p_vec, h_step, sol; dir=plot_dump)
    k_qm = H.calc_wigner_k_qm(q_vec, p_vec, h_step, sol)
    @printf "k_qm = %.2E 1/s \n" k_qm[end] ./ H.au_time
    H.plot_QSE_k_qm(sol.t, k_qm, dir=plot_dump)

    # H.animate_wigner_heatmap(q_vec, p_vec, sol; dir=plot_dump)
    # H.animate_wigner_wq(q_vec, p_vec, sol; yscale=:log10, dir=plot_dump)
    # H.animate_wigner_wp(q_vec, p_vec, sol; yscale=:log10, dir=plot_dump)
end


if f_qse
    # Time grid
    t0 = 0.0
    nt = 200
    saveat = H.linspace(t0, t_rate, nt)
    algo = VCABM5() #ROCK4() #Tsit5() SSPRK83 VCABM5 Feagin12 Feagin14 Vern9
    spl = Spline1D(q_vec, v_vec; k=5, bc="extrapolate")
    q_vec, p_vec, Q_vec, P_vec, dq, dp = H.create_basis(N_qse, ql, qh, pl, ph)
    v_vec = spl(q_vec)
    v_vec = v_vec .- minimum(v_vec)

    #W0 = H.w0_gaussian(q_vec, p_vec, q0, p0, sigma)
    W0 = H.w0_thermal_step(q_vec, p_vec, v_vec, mass, beta, h_step)

    P0 = H.convert_W_to_W_q(q_vec, p_vec, W0)
    P0 = H.QSE_normalise(q_vec, P0)

    # H.plot_QSE_P_log(q_vec, P0; dir=plot_dump)
    prep = H.prep_QSE(apx, q_vec, v_vec, mass, temperature, gamma)
    prob = ODEProblem(H.QSE_eq!, P0, (t0, t_rate), prep)

    @time sol = solve(
        prob,
        algo,
        saveat=saveat,
        progress=true,
        reltol=tol,
        abstol=tol,
        maxiters=1e7,
    )
    # Get the solution
    solu = sol.u
    discrete_t = sol.t

    H.plot_QSE_normalisation_log(q_vec, solu, discrete_t)
    #H.plot_QSE_error_compare(q_vec, v_vec, sol, discrete_t; q0=h_step)
    # H.plot_QSE_P_log(q_vec, solu[end])

    k_qm = H.calc_qse_k_qm(q_vec, h_step, sol; f_renorm=true)
    H.printout("k_qm = $(round(k_qm[end] ./ H.au_time;sigdigits=3)) 1/s")
    I_qm = k_qm[end] ./ H.au_time .* 1.602e-19
    H.printout("I_qm = $(round(I_qm;sigdigits=3)) q/s")

    H.plot_QSE_k_qm(discrete_t, k_qm)
    # H.animate_QSE_P(q_vec, solu, discrete_t)
    # H.animate_QSE_P_log(q_vec, solu, discrete_t)
end


if f_scatter
    H.printout("Running scatter process")
    W_thermal_wq, _ = H.calc_wigner_wqp(q_vec, p_vec, W_thermal)
    loc_pk = findmax(W_thermal_wq)[2]
    q0 = q_vec[loc_pk]
    sigma = 10.0
    p0 = 30.0

    # momentum to energy Conversions
    e_incident = 0.5 * p0^2 / mass
    H.printout("e_incident = $(round(e_incident;sigdigits=3)) Eh")
    e_thermal = 0.5 * H.k_b * temperature
    H.printout("e_thermal = $(round(e_thermal;sigdigits=3)) Eh")
    boltzman_weight = exp(-(e_incident) * beta)
    H.printout("boltzman_weight = $(round(boltzman_weight;sigdigits=3))")
    t_scatter = mass *  (h_step-q0) / p0 * 2.0
    H.printout("t_scatter = $(round(t_scatter;sigdigits=3)) AUT")
    H.printout("t_scatter = $(round(t_scatter * H.aut_2_ps;sigdigits=3)) fs")
    t_scatter = minimum([t_scatter, 1.0 / H.aut_2_ps]) # Do not let it go past 1 ps

    # Prepare the time evolution
    saveat = H.linspace(t0, t_scatter, nt)
    if f_oqs
        # Make intial gaussian
        W0 = H.w0_gaussian(q_vec, p_vec, q0, p0, sigma)

        H.plot_wigner_heatmap(q_vec, p_vec, W0; dir=plot_dump, name="w0.pdf")
        # H.plot_wigner_wq(q_vec, p_vec, W0; yscale=:log10, dir=plot_dump, name="w0_q.pdf")
        # H.plot_wigner_wp(q_vec, p_vec, W0; yscale=:log10, dir=plot_dump, name="w0_p.pdf")

        prep = H.prep_LL_HT_M_fd(q_vec, p_vec, v_vec, mass, h_bar, gamma, beta)
        prob = ODEProblem(H.LL_HT_M_fd_trunc!, W0, (t0, t_scatter), prep)
        # prep = H.prep = H.prep_wm_fd(q_vec, p_vec, v_vec, mass, h_bar)
        # prob = ODEProblem(H.wm_fd_trunc!, W0, (t0, t_scatter), prep)
        H.printout("Running the simulation...")
        @time sol = solve(prob, algo, progress=true, saveat=saveat, abstol=tol, reltol=tol)
    
        H.plot_wigner_normalisation(q_vec, p_vec, sol; yscale=:log10, dir=plot_dump)
        # Calculate the properties of the final state
        H.plot_wigner_step_occupation(q_vec, p_vec, h_step, sol; dir=plot_dump)
        # get the transmission probability
        prob = [H.calc_wigner_probability(q_vec, p_vec, h_step, sol[i]) for i = 1:length(saveat)]
        H.printout("Maximum probability = $(round(maximum(prob);sigdigits=3))")
        H.animate_wigner_heatmap(q_vec, p_vec, sol; dir=plot_dump)

    else
        P0 = gaussian_mom(q_vec, q0, p0, sqrt(sigma)/4.0, h_bar)
        # Normalise
        P0 = normalise_SE(q_vec, P0)
        # H.plot_QSE_P(q_vec, real.(P0))
        # H.plot_QSE_P(q_vec, imag.(P0))
        plot_SE_prob_dens(q_vec, P0)

        prep = prep_SE(q_vec, v_vec, mass, h_bar; apx=apx)
        prob = ODEProblem(SE_fd!, P0, (t0, t_scatter), prep)
        H.printout("Running the simulation...")
        @time sol = solve(prob, algo, progress=true, saveat=saveat, abstol=tol, reltol=tol)
    

        prob = [calc_SE_probability(q_vec, h_step, sol[i]) for i = 1:length(saveat)]
        H.printout("Maximum probability = $(round(maximum(prob);sigdigits=3))")
        plot_QSE_normalisation_log(q_vec, sol.u, sol.t)
        plot_SE_probability(q_vec, h_step, sol.u, sol.t)
        # animate_SE(q_vec, sol.u, sol.t)
        # animate_SE_log(q_vec, sol.u, sol.t)

    end

end

H.printout("End")