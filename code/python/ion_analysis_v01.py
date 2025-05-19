import os
import numpy as np
from scipy.fftpack import fft, fftfreq
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve, windows
from scipy import integrate
from ase import units

plt.rcParams['axes.linewidth'] = 2.0


def get_all_files(directory):
    file_list = []
    for root, directories, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            file_list.append(file_path)
    return file_list


# List only the top level folders in a directory
def folder_list(mypath=os.getcwd()):
    """
    List only the top level folders in a directory given by mypath
    :param mypath: specified directory, defaults to current directory
    :return: returns a list of folders
    """
    onlyfolders = [f for f in os.listdir(mypath) if os.path.isdir(os.path.join(mypath, f))]
    return onlyfolders


def file_list(mypath=os.getcwd()):
    """
    List only the files in a directory given by mypath
    :param mypath: specified directory, defaults to current directory
    :return: returns a list of files
    """
    onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
    return onlyfiles


# List only files which contain a substring
def sub_file_list(mypath, sub_str):
    """
    List only files which contain a given substring
    :param mypath: specified directory
    :param sub_str: string to filter by
    :return: list of files which have been filtered
    """
    return [i for i in get_all_files(mypath) if sub_str in i]


def os_display():
    if os.name == 'nt':
        plt.show()
    plt.close()


def n_plot(xlab, ylab, xs=14, ys=14):
    plt.minorticks_on()
    plt.tick_params(axis='both', which='major', labelsize=ys - 2, direction='in', length=6, width=2)
    plt.tick_params(axis='both', which='minor', labelsize=ys - 2, direction='in', length=4, width=2)
    plt.tick_params(axis='both', which='both', top=True, right=True)
    if xlab is not None:
        plt.xlabel(xlab, fontsize=xs)
    if ylab is not None:
        plt.ylabel(ylab, fontsize=ys)
    plt.tight_layout()
    return None

def ax_plot(fig, ax, xlab, ylab, xs=14, ys=14):
    ax.minorticks_on()
    ax.tick_params(axis='both', which='major', labelsize=ys - 2, direction='in', length=6, width=2)
    ax.tick_params(axis='both', which='minor', labelsize=ys - 2, direction='in', length=4, width=2)
    ax.tick_params(axis='both', which='both', top=True, right=True)
    ax.set_xlabel(xlab, fontsize=xs)
    ax.set_ylabel(ylab, fontsize=ys)
    fig.tight_layout()
    return None



def conver_string_2_data(string):
    string = string.split("&")
    string = [float(m) for m in string]
    string = np.array(string)
    # # rearrange the order so that the last two are the first two
    # string = string[[4, 5, 0, 1, 2, 3]]
    return string


def load_vel(file):
    data = np.loadtxt(file, delimiter="\t", dtype=float, comments=["#", "@"])
    time = data[:, 0] * ps_2_fs
    V = np.array([data[:, 1], data[:, 2], data[:, 3]])
    # convert velocity from nm/ps to nm/fs
    V = V / ps_2_fs
    return time, V


def load_force(file):
    data = np.loadtxt(file, delimiter="\t", dtype=float, comments=["#", "@"])
    time = data[:, 0] * ps_2_fs
    F = np.array([data[:, 1], data[:, 2], data[:, 3]])
    # Convert force from kJ/mol/nm to eV/nm
    F = F * kJ_2_eV
    return time, F


def pad_zeros(arr, axis=0, where='end', nadd=None, upto=None, tonext=None,
              tonext_min=None):
    if tonext == False:
        tonext = None
    lst = [nadd, upto, tonext]
    assert lst.count(None) in [2, 3], "`nadd`, `upto` and `tonext` must be " + \
                                      "all None or only one of them not None"
    if nadd is None:
        if upto is None:
            if (tonext is None) or (not tonext):
                # default
                nadd = arr.shape[axis]
            else:
                tonext_min = arr.shape[axis] if (tonext_min is None) \
                    else tonext_min
                # beware of int overflows starting w/ 2**arange(64), but we
                # will never have such long arrays anyway
                two_powers = 2 ** np.arange(30)
                assert tonext_min <= two_powers[-1], ("tonext_min exceeds "
                                                      "max power of 2")
                power = two_powers[np.searchsorted(two_powers,
                                                   tonext_min)]
                nadd = power - arr.shape[axis]
        else:
            nadd = upto - arr.shape[axis]
    if nadd == 0:
        return arr
    add_shape = list(arr.shape)
    add_shape[axis] = nadd
    add_shape = tuple(add_shape)
    if where == 'end':
        return np.concatenate((arr, np.zeros(add_shape, dtype=arr.dtype)), axis=axis)
    elif where == 'start':
        return np.concatenate((np.zeros(add_shape, dtype=arr.dtype), arr), axis=axis)
    else:
        raise Exception("illegal `where` arg: %s" % where)


def slicetake(a, sl, axis=None, copy=False):
    if axis is None:
        slices = sl
    else:
        slices = [slice(None)] * a.ndim
        slices[axis] = sl
    slices = tuple(slices)
    if copy:
        return a[slices].copy()
    else:
        return a[slices]


def smooth(data, kern, axis=0, edge='m', norm=True):
    # https://elcorto.github.io/pwtools/generated/api/pwtools.signal.smooth.html#pwtools.signal.smooth
    N = data.shape[axis]
    M = kern.shape[axis]
    if edge == 'm':
        npad = min(M, N)
        sleft = slice(npad, 0, -1)
        sright = slice(-2, -(npad + 2), -1)
        dleft = slicetake(data, sl=sleft, axis=axis)
        dright = slicetake(data, sl=sright, axis=axis)
        assert dleft.shape == dright.shape
        K = dleft.shape[axis]
        if K < M:
            dleft = pad_zeros(dleft, axis=axis, where='start', nadd=M - K)
            dright = pad_zeros(dright, axis=axis, where='end', nadd=M - K)
    elif edge == 'c':
        sl = [slice(None)] * data.ndim
        sl[axis] = None
        tsl = tuple(sl)
        dleft = np.repeat(slicetake(data, sl=0, axis=axis)[tsl], M, axis=axis)
        dright = np.repeat(slicetake(data, sl=-1, axis=axis)[tsl], M, axis=axis)
        assert dleft.shape == dright.shape
        # 1d special case: (M,1) -> (M,)
        if data.ndim == 1 and dleft.ndim == 2 and dleft.shape[1] == 1:
            dleft = dleft[:, 0]
            dright = dright[:, 0]
    else:
        raise Exception("unknown value for edge")
    sig = np.concatenate((dleft, data, dright), axis=axis)
    kk = kern / float(kern.sum()) if norm else kern
    ret = fftconvolve(sig, kk, 'valid')
    assert ret.shape[axis] == N + M + 1, "unexpected convolve result shape"
    del sig
    if M % 2 == 0:
        ##sl = slice(M//2+1,-(M//2)) # even kernel, shift result to left
        sl = slice(M // 2, -(M // 2) - 1)  # even kernel, shift result to right
    else:
        sl = slice(M // 2 + 1, -(M // 2) - 1)
    ret = slicetake(ret, sl=sl, axis=axis)
    assert ret.shape == data.shape, ("ups, ret.shape (%s)!= data.shape (%s)"
                                     % (ret.shape, data.shape))
    return ret


def calculate_autocorrelation(X, f_pos=False, f_norm=True):
    # Auto-correlation for all degrees of freedom
    auto = [np.correlate(x, x, 'full') for x in X]
    # Normalise
    if f_norm:
        auto /= np.linalg.norm(auto, axis=1)[:, None]
        # auto /= auto[0]
    # Average over all degrees of freedom
    auto = np.mean(auto, axis=0)
    # Remove the background based on the average
    # auto -= np.mean(auto)
    # Only take the positive section of the auto-correlation
    if f_pos:
        auto = auto[int(auto.shape[0] / 2):]
    return auto


def calculate_pdos(vel, dt=1.0, filter=True, avg=False, avg_win=7, filt_win=11):
    # Number of time steps
    n = vel.shape[1]
    # Velocity auto-correlation for all degrees of freedom
    vac = calculate_autocorrelation(vel)
    # Calculate power spectrum (phonon density of states)
    pdos = np.square(np.abs(fft(vac)))
    # Normalise
    pdos /= np.linalg.norm(pdos) / 2  # spectrum is symmetric
    # Frequency axis
    freq = fftfreq(2 * n - 1, dt) * fs_2_cm  # Frequency in cm^-1
    rtn_freq = freq[:n // 2]
    rtn_pdos = pdos[:n // 2]
    # Smoothing of the spectrum removes numerical artifacts due to finite time truncation of the FFT
    if filter:
        # rtn_pdos = gaussian(rtn_pdos, sigma=50)
        # rtn_pdos = signal.savgol_filter(rtn_pdos, 7, 2, mode='nearest')
        # rtn_pdos = signal.savgol_filter(rtn_pdos, 7, 2, mode='nearest')
        k = windows.hann(filt_win)
        rtn_pdos = smooth(rtn_pdos, k)
        pdos /= np.linalg.norm(pdos) / 2  # spectrum is symmetric

    if avg:
        rtn_pdos = np.convolve(rtn_pdos, np.ones(avg_win) / avg_win, mode='same')
    return rtn_freq, rtn_pdos


def plot_pdos(freq, pdos, xlim="default", yscale="log", title="pdos"):
    if xlim == "default":
        xlim = [0.0, 4000.0]
    xlab = 'Frequency [cm$^{-1}$]'
    ylab = 'Density of states [a.u.]'
    plt.plot(freq, pdos)
    plt.xlim(xlim)
    plt.yscale(yscale)
    plt.title('Velocity auto-correlation function')
    n_plot(xlab, ylab)
    plt.savefig(title + '.pdf')
    plt.savefig(title + '.png')
    plt.show()
    return


def plot_pdos_multi(fig, ax, freq, pdos, xlim="default", yscale="log", label=None, f_color=False):
    if xlim == "default":
        xlim = [0.0, 4000.0]
    if f_color:
        ax.plot(freq, pdos, alpha=0.6, label=label)
    else:
        ax.plot(freq, pdos, alpha=0.5, color='k', label=label)
    ax.set_xlim(xlim)
    ax.set_yscale(yscale)
    ax_plot(fig, ax, "Frequency [cm$^{-1}$]", "Density of states [a.u.]")
    return None


def calc_friction(F, dt, temperature=300.0, f_norm=False):

    # Velocity auto-correlation for all degrees of freedom
    fac = calculate_autocorrelation(F, f_pos=True, f_norm=f_norm)
    # Calculate the friction coefficient
    pre_factor = 1.0 / (3.0 * temperature * units.kB)
    lam = fac * pre_factor
    lam = integrate.cumulative_trapezoid(lam, dx=dt, initial=0.0)
    #lam = integrate.cumulative_trapezoid(lam, dx=dt, initial=0.0)
    lam = np.cumsum(lam) * dt
    return lam

def plot_friction(fig, ax, time, fric, label=None, f_color=False):
    if f_color:
        ax.plot(time, fric, alpha=0.6, label=label)
    else:
        ax.plot(time, fric, alpha=0.5, color='k', label=label)
    ax_plot(fig, ax, "Time [fs]", r"Friction coefficient, $\lambda$ [eV/fs]")
    return None


f_rates = False
f_md = True
if f_rates:
    mass = "17.4       & 17.9      & 22.6      & 23.6      & 6.21       & 6.28 "
    mass = conver_string_2_data(mass)

    mass_t = "95.0       & 103.0     & 93.0      & 99.0      & 61.1       & 67.1"
    mass_t = conver_string_2_data(mass_t)

    epsilon = "6.05E-04   & 9.78E-04  & 9.75E-04  & 7.11E-04  & 8.18E-04   & 9.37E-04"
    # epsilon = "5.97E-04   & 1.07E-03  & 1.02E-03  & 7.68E-04  & 1.34E-03   & 9.78E-04"
    epsilon = conver_string_2_data(epsilon)

    dg = "9.88E-03   & 1.43E-02  & 1.25E-02  & 1.39E-02  & 6.03E-03   & 7.55E-03"
    # dg = "1.05E-02   & 1.03E-02  & 1.09E-02  & 1.18E-02  & 1.24E-02   & 9.08E-03"
    dg = conver_string_2_data(dg)

    k_c = "1.78E+08   & 1.71E+06  & 1.14E+07  & 2.63E+06  & 1.05E+10   & 2.10E+09"
    # k_c = "9.40E+07   & 1.08E+08  & 5.86E+07  & 2.39E+07  & 1.25E+07   & 4.12E+08"
    k_c = conver_string_2_data(k_c)

    k_q = "1.57E+07   & 5.20E+05  & 6.79E+06  & 7.14E+05  & 1.71E+09   & 4.60E+08"
    # k_q = "8.48E+06   & 1.99E+07  & 8.72E+06  & 2.88E+06  & 2.30E+06   & 6.99E+07"
    k_q = conver_string_2_data(k_q)

    t_d = "0.126      & 0.151     & 0.101     & 0.126     & 0.251      & 0.327"
    # t_d = "0.126      & 0.151     & 0.101     & 0.126     & 0.226      & 0.302"
    t_d = conver_string_2_data(t_d)

    text_labs = ["Na(H)", "Na(D)", "K(H)", "K(D)", "Li(H)", "Li(D)"]

    mass_plot = mass_t
    idx = np.argsort(mass_plot)

    text_labs = [text_labs[m] for m in idx]
    mass_plot = mass_plot[idx]
    epsilon = epsilon[idx]
    dg = dg[idx]
    k_c = k_c[idx]
    k_q = k_q[idx]
    t_d = t_d[idx]

    plt.scatter(mass_plot, epsilon)
    plt.plot(mass_plot, epsilon)
    n_plot("Mass [Da]", "Lowest eigenstate energy, $\epsilon_0$ [Ha]")
    plt.show()

    plt.scatter(mass_plot, dg)
    plt.plot(mass_plot, dg)
    n_plot("Mass [Da]", "Trapping potential barrier, $\Delta G$ [Ha]")
    plt.show()

    plt.scatter(mass_plot, k_c)
    plt.plot(mass_plot, k_c)
    plt.yscale("log")
    n_plot("Mass [Da]", r"Escape rate, $k_{\mathrm{class}}$  [$s^{-1}$]")
    plt.show()

    plt.scatter(mass_plot, k_q)
    plt.plot(mass_plot, k_q)
    plt.yscale("log")
    n_plot("Mass [Da]", r"$Escape rate, k_{\mathrm{QM}}$  [$s^{-1}$]")
    plt.show()

    plt.scatter(mass_plot, k_c, label="Classical")
    plt.plot(mass_plot, k_c)
    plt.scatter(mass_plot, k_q, label="Quantum")
    plt.plot(mass_plot, k_q)
    # Annotate the plot
    for i, txt in enumerate(text_labs):
        plt.annotate(txt, (mass_plot[i], k_c[i]), fontsize=12)

    plt.yscale("log")
    plt.legend()
    n_plot("Mass [Da]", r"Escape rate, $k$ [$s^{-1}$]")
    plt.savefig("escape_rates.pdf")
    plt.show()

    plt.scatter(mass_plot, t_d)
    plt.plot(mass_plot, t_d)
    n_plot("Mass [Da]", r"Decoherence timescale, $\tau_{\mathrm{deco}}$ [fs]")
    plt.show()

if f_md:
    # convert ps to fs
    ps_2_fs = 1e3
    # convert kJ/mol to eV
    kJ_2_eV = 0.010364
    # convert fs to cm^-1
    fs_2_cm = 33356.4095198152
    super_dir = r"C:\Users\ls0067\OneDrive - University of Surrey\Papers\paper_ion_channel\md_data"
    dirs = folder_list(super_dir)
    dirs = ['Li', 'Li_HW', 'Na', 'Na_HW', 'K', 'K_HW']
    # dirs = ['K', 'K_HW']
    # dirs = ['Na', 'Na_HW']
    # dirs = ['Li', 'Li_HW']
    # dirs = ['K', 'K_HW', 'Na', 'Na_HW']
    print(dirs)
    fig_1, ax_1 = plt.subplots()
    fig_2, ax_2 = plt.subplots()
    fig_3, ax_3 = plt.subplots()
    for j, directory in enumerate(dirs):
        # directory = r"C:\Users\ls0067\OneDrive - University of Surrey\Papers\paper_ion_channel\md_data\K"
        files = sub_file_list(os.path.join(super_dir, directory), "vel")
        array_pdos = []

        for i, file in enumerate(files):
            time, V = load_vel(file)
            freq, pdos = calculate_pdos(V, dt=1.0, filter=True, avg=True, avg_win=7, filt_win=11)
            array_pdos.append(pdos)

        array_fric = []
        files = sub_file_list(os.path.join(super_dir, directory), "for")
        for i, file in enumerate(files):
            time, F = load_force(file)
            fric = calc_friction(F, dt=1.0, temperature=300.0, f_norm=True)
            array_fric.append(fric)
        avg_fric = np.mean(np.array(array_fric), axis=0)
        plot_friction(fig_3, ax_3, time, avg_fric, label=directory, f_color=True)

        avg_pdos = np.mean(np.array(array_pdos), axis=0)
        print("Integral of the spectrum: {}".format(integrate.simps(avg_pdos, freq)))
        plot_pdos_multi(fig_1, ax_1, freq, avg_pdos, xlim=[0.0, 3000.0], yscale="log", label=directory, f_color=True)
        plot_pdos_multi(fig_2, ax_2, freq, avg_pdos, xlim=[0.0, 400.0], yscale="linear", label=directory, f_color=True)
    ax_1.legend(loc="lower left")
    fig_1.savefig("pdos_log_all.pdf")
    fig_1.savefig("pdos_log_all.png")
    ax_2.legend(loc="upper right")
    fig_2.savefig("pdos_all.pdf")
    fig_2.savefig("pdos_all.png")

    ax_3.legend(loc="upper left")
    fig_3.savefig("friction_all.pdf")
    fig_3.savefig("friction_all.png")

    plt.show()
    plt.close("all")