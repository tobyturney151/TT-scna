"""
Takes as input a text file containing
electrophysiological data in two conditions, and
computes the difference in single-channel noise between the two conditions.

Must contain one column, with time, expressed in seconds,
and a second column, with current, expressed in picoamps.

Inputs
----------
file_name: str
    Path and File_name to text data file.

on_start: float
    Start time of 'on' condition.

off_start: float
    Start time of 'off' condition.

sweep_len: float
    Length of time, in seconds, of each condition to analyze.

sample_rate: float
    Sampling rate of the time series, in Hz.

skip_rows: int
    Rows to skip in the header of the data file.

sep: str
    Delimiter used in each column for the data file.

fmin: float
    Minimum frequency to fit to, in Hz.

fmax: float
    Maximum frequency to fit to, in Hz.

Output
----------
Figure: .png
    A 4-panel plot containing:
    100 microsecond traces of the On and Off conditions.

    Overlaid power spectra from both conditions.

    The power spectra of the on condition minus the off condition,
    overlaid with a Lorentzian model.

Table: .csv
    A .csv containing the unitary conductance and mean dwell time
    predicted from the data.
"""
# Import modules
# import os
import numpy as np
import pandas as pd

# import scipy.interpolate
import scipy.signal
import scipy.optimize
import matplotlib.pyplot as plt
import argparse

# Define your functions!


def loren(par, x):
    """
    Generates a Lorentzian fit from the given parameters

    Inputs
    ----------
    par: list of floats
        List of parameters to the Lorentzian.
        par[0] = DC Component
        par[1] = Shoulder Frequency

    x: list of floats
        input x-data for the model, in Hz

    Output
    ----------
    lor_model: float
        The Lorentzian model resulting from the inputted parameters.
    """

    lor_model = par[0] / (1 + (x / par[1]))
    return lor_model


def resid(par, x, y, fmin, fmax):
    """
    Calculate the squared error
    from a Lorentzian model with the given parameters.

    Inputs
    ----------
    par: list of floats
        List of parameters to the Lorentzian.
        par[0] = DC Component
        par[1] = Shoulder Frequency

    x: List of floats
        Input x-data for the model

    y: List of floats
        Input y-data for the model

    fmin: float
        Minimum frequency to fit to, in Hz.

    fmax: float
        Maximum frequency to fit to, in Hz.

    Output
    ----------
    SE: float
        The squared errors returned from lor_model with the given parameters.
    """

    # Get indices of min and max frequencies to analyze.
    min_idx = list(x).index(fmin)
    max_idx = list(x).index(fmax)

    # Truncate the x and y-data to the min and max frequencies to analyze.
    x_trunc = x[min_idx:max_idx]
    y_trunc = y[min_idx:max_idx]

    SE = (loren(par, x_trunc) - y_trunc) ** 2
    return SE


def plot_my_data(
    on_time,
    subtracted_on,
    off_time,
    subtracted_off,
    f_on,
    f_off,
    spec_on,
    spec_off,
    delta_spec,
    model_PSD,
    file_name,
):
    """
    Makes a 4-panel plot containing:
        100 microsecond traces of the On and Off conditions.

        Overlaid power spectra from both conditions.

        The power spectra of the on condition minus the off condition,
        overlaid with a Lorentzian model.

    Inputs
    ----------
    on_time: List of floats
        Time (in seconds) of the 'on' condition.

    subtracted_on: List of floats
        Current (in pA) of the 'On' condition
        after subtracting out the mean current.

    off_time: List of floats
        Time (in seconds) of the 'on' condition.

    subtracted_off: List of floats
        Current (in pA) of the 'On' condition
        after subtracting out the mean current.

    f_on: List of floats
        Frequencies for the Power Spectra of
        the 'On' Condition.

    f_off: List of floats
        Frequencies for the Power Spectra of
        the 'Off' Condition.

    spec_on: List of floats
        Power Spectral Density of the 'On' Condition.

    spec_off: List of floats
        Power Spectral Density of the 'Off' Condition.

    delta_spec: List of floats
        Result of subtracting spec_off from spec_on.

    model_PSD: List of floats
        Lorentzian Model fitted to delta_spec.

    file_name:
        Path and File_name to text data file.

    Output
    ----------
    Figure
    """

    plt.rcParams.update({"font.size": 36})
    fig = plt.figure(figsize=(60, 20))

    # For the On and Off traces, only showing a subset
    # of the full trace that was analyzed.

    # On Trace
    on = fig.add_subplot(2, 2, 1)
    on.spines["top"].set_visible(False)
    on.spines["right"].set_visible(False)
    on.spines["bottom"].set_color("k")
    on.spines["bottom"].set_linewidth(5)
    on.spines["left"].set_color("k")
    on.spines["left"].set_linewidth(5)
    on.plot(on_time[0:2000], subtracted_on[0:2000], c="k", linewidth=5)
    on.tick_params(direction="out", length=15, width=5, colors="k")
    on.set(
        ylabel="Current (pA)",
        xlabel="Time (\u03bcs)",
        title="On",
        ylim=[
            min(min(subtracted_on), min(subtracted_off)),
            max(max(subtracted_on), max(subtracted_off)),
        ],
        xticks=[
            on_time[0],
            on_time[400],
            on_time[800],
            on_time[1200],
            on_time[1600],
            on_time[2000],
        ],
        xticklabels=["0", "20", "40", "60", "80", "100"],
    )

    # Off trace
    off = fig.add_subplot(2, 2, 3)
    off.spines["top"].set_visible(False)
    off.spines["right"].set_visible(False)
    off.spines["bottom"].set_color("k")
    off.spines["bottom"].set_linewidth(5)
    off.spines["left"].set_color("k")
    off.spines["left"].set_linewidth(5)
    off.plot(off_time[0:2000], subtracted_off[0:2000], c="k", linewidth=5)
    off.tick_params(direction="out", length=15, width=5, colors="k")
    off.set(
        ylabel="Current (pA)",
        xlabel="Time (\u03bcs)",
        title="Off",
        ylim=[
            min(min(subtracted_on), min(subtracted_off)),
            max(max(subtracted_on), max(subtracted_off)),
        ],
        xticks=[
            off_time[0],
            off_time[400],
            off_time[800],
            off_time[1200],
            off_time[1600],
            off_time[2000],
        ],
        xticklabels=["0", "20", "40", "60", "80", "100"],
    )

    # Overlaid Power Spectra of both conditions
    ps = fig.add_subplot(2, 2, 2)
    ps.spines["top"].set_visible(False)
    ps.spines["right"].set_visible(False)
    ps.spines["bottom"].set_color("k")
    ps.spines["bottom"].set_linewidth(5)
    ps.spines["left"].set_color("k")
    ps.spines["left"].set_linewidth(5)
    ps.loglog(f_on, spec_on, c="#ff33fc", linewidth=5, label="On")
    ps.loglog(f_off, spec_off, c="#33caff", linewidth=5, label="Off")
    ps.legend()
    ps.tick_params(direction="out", length=15, width=5, colors="k")
    ps.set(
        ylabel="PSD, (pA$^2$/Hz)",
        xlabel="Frequency (Hz)",
        title="Power Spectra",
        ylim=[
            min(min(spec_on[-1000:-1]), min(spec_off[-1000:-1])),
            max(max(spec_off), max(spec_on)),
        ],
    )

    # Delta PSD and Lorentzian Fit
    psd = fig.add_subplot(2, 2, 4)
    psd.spines["top"].set_visible(False)
    psd.spines["right"].set_visible(False)
    psd.spines["bottom"].set_color("k")
    psd.spines["bottom"].set_linewidth(5)
    psd.spines["left"].set_color("k")
    psd.spines["left"].set_linewidth(5)
    psd.loglog(f_on, delta_spec, c="k", linewidth=5, label="Experiment")
    psd.loglog(f_on, model_PSD, c="r", linewidth=7, label="Fitted Model")
    psd.legend()
    psd.tick_params(direction="out", length=15, width=5, colors="k")
    psd.set(
        ylabel="PSD, (pA$^2$/Hz)",
        xlabel="Frequency (Hz)",
        title="Power Spectra - Delta",
        ylim=[model_PSD[-1] / 10, model_PSD[0] * 10],
    )
    fig.tight_layout()

    fig.savefig(file_name[:-4] + "_results.png")


def ttscna(
    file_name,
    on_start,
    off_start,
    sweep_len=5.0,
    sample_rate=20000,
    skip_rows=9,
    sep="\t",
    fmin=20.0,
    fmax=10000.0,
):
    """
    Takes as input a text file containing
    electrophysiological data in two conditions, and
    computes the difference in single-channel noise between the two conditions.

    Must contain one column, with time, expressed in seconds,
    and a second column, with current, expressed in picoamps.

    Inputs
    ----------
    file_name: str
        Path and File_name to text data file.

    on_start: float
        Start time of 'on' condition.

    off_start: float
        Start time of 'off' condition.

    sweep_len: float
        Length of time, in seconds, of each condition to analyze.

    sample_rate: float
        Sampling rate of the time series, in Hz.

    skip_rows: int
        Rows to skip in the header of the data file.

    sep: str
        Delimiter used in each column for the data file.

    fmin: float
        Minimum frequency to fit to, in Hz.

    fmax: float
        Maximum frequency to fit to, in Hz.

    Output
    ----------
    Figure: .png
        A 4-panel plot containing:
        100 microsecond traces of the On and Off conditions.

        Overlaid power spectra from both conditions.

        The power spectra of the on condition minus the off condition,
        overlaid with a Lorentzian model.

    Table: .csv
        A .csv containing the unitary conductance and mean dwell time
        predicted from the data.
    """

    # Open file.
    # Note to self: pd.read_table is approx 30X faster than np.load_txt
    sweeps = pd.read_table(file_name, delimiter=sep, skiprows=skip_rows)

    # Convert data to lists.
    time = list(sweeps["Time (s)"])
    current = list(sweeps["Trace #1 (pA)"])

    # Grab list indices for beginning and end of on and off conditions.
    on_start_idx = time.index(on_start)
    on_stop_idx = time.index(on_start + sweep_len)
    off_start_idx = time.index(off_start)
    off_stop_idx = time.index(off_start + sweep_len)

    # Use list indices to grab the 'on' and 'off' traces.
    on_time = time[on_start_idx:on_stop_idx]
    on_curr = current[on_start_idx:on_stop_idx]
    off_time = time[off_start_idx:off_stop_idx]
    off_curr = current[off_start_idx:off_stop_idx]

    # Subtract mean current to remove DC component.
    subtracted_on = on_curr - np.mean(on_curr)
    subtracted_off = off_curr - np.mean(off_curr)

    # Get periodograms.
    f_on, spec_on = scipy.signal.periodogram(subtracted_on, fs=sample_rate)
    f_off, spec_off = scipy.signal.periodogram(subtracted_off, fs=sample_rate)

    # Subtract the 'Off' PSD from the 'On' PSD.
    delta_spec = spec_on - spec_off

    # Now fit the periodogram to a lorentzian!
    rez = scipy.optimize.leastsq(resid,
                                 [1, 1],
                                 args=(f_off, delta_spec, fmin, fmax))[0]
    model_PSD = loren(rez, f_on)

    # Convert the parameters of the Lorentzian
    # to unitary conductances and dwell times.
    # Single-channel Unitary Conductance
    # 1000X because you're converting pS to fS
    u_cond = (
        rez[0]
        * np.pi
        * rez[1]
        / (-0.06 * 2 * (np.mean(on_curr) - np.mean(off_curr)))
        * 1000
    )

    # Mean Dwell Time
    # 1000X because you're converting s to ms.
    dwell_time = 1 / (2 * rez[1] / np.pi) * 1000

    # Average bulk current, subtracting out any leak.
    sub_curr = np.mean(on_curr) - np.mean(off_curr)

    # Plot everything!
    plot_my_data(
        on_time,
        subtracted_on,
        off_time,
        subtracted_off,
        f_on,
        f_off,
        spec_on,
        spec_off,
        delta_spec,
        model_PSD,
        file_name,
    )

    # Export relevant data to a data frame.
    data = pd.DataFrame(
        {
            "Unitary Conductance (fS)": u_cond,
            "Dwell Time (ms)": dwell_time,
            "Mean Bulk Current (pA)": sub_curr,
        },
        index=[0],
    )

    # Now export the metadata from this analysis into a .csv.
    # The '[:-4]' removes the file extension.
    data.to_csv(file_name[:-4] + "_results.csv")


# Parse command line arguments (when using from command line)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="TTscna - Single Channel Noise Analysis"
    )
    parser.add_argument(
        "--file",
        action="store",
        dest="file_name",
        help="Name of the file to read, including path and file extension.",
    )
    parser.add_argument(
        "--on",
        action="store",
        dest="on_start",
        help="Start time of the on condition."
    )
    parser.add_argument(
        "--off",
        action="store",
        dest="off_start",
        help="Start time of the off condition.",
    )
    parser.add_argument(
        "--len",
        action="store",
        dest="sweep_len",
        default=5.0,
        help="Length of time, in seconds, of each condition to analyze.",
    )
    parser.add_argument(
        "--rate",
        action="store",
        dest="sample_rate",
        default=20000.0,
        help="Sampling rate of the time series, in Hz.",
    )
    parser.add_argument(
        "--skip",
        action="store",
        dest="skip_rows",
        default=9,
        help="Rows to skip in the header of the data file.",
    )
    parser.add_argument(
        "--sep",
        action="store",
        dest="sep",
        default="\t",
        help="Delimiter used in each column for the data file.",
    )
    parser.add_argument(
        "--fmin",
        action="store",
        dest="fmin",
        default=20.0,
        help="Minimum frequency (in Hz) to fit to.",
    )
    parser.add_argument(
        "--fmax",
        action="store",
        dest="fmax",
        default=10000.0,
        help="Maximum frequency (in Hz) to fit to.",
    )
    answer = parser.parse_args()

    # Execute Function!
    print(
        ttscna(
            answer.input_string,
            answer.on_start,
            answer.off_start,
            answer.sweep_len,
            answer.sample_rate,
            answer.skip_rows,
            answer.sep,
            answer.fmin,
            answer.fmax,
        )
    )
