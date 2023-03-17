import os

import numpy as np

from ..matlab_engine.matlab_engine import create_engine_for_testing
from ..timestepping import (
    gaussian_pulse_parameters,
    max_allowable_timestep,
    min_timesteps_fdtd,
)

LOCATION_OF_THIS_FILE = os.path.dirname(os.path.abspath(__file__))


def test_gaussian_pulse_parameters() -> None:
    """gaussian_pulse_parameters should be equivalent to it's MATLAB counterpart, fdtdduration"""
    # fdtdduration needs to fetch from an input file, so provide one here.
    MFILE = os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_08.m"
    )
    # Replicate the constants in the input file
    dx = 1300.0e-9 / 8.0
    dt = 2.0 * dx * 0.95 / (np.sqrt(3.0) * np.pi * 3.0e8 * 1.35)
    wavelength_width = 120.0e-9
    f_an = np.arcsin(2.0 * np.pi * 2.997924580105029e8 * dt / (2 * 1300e-9)) / (
        np.pi * dt
    )
    epsr = 1.35**2

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    mat_t_0, mat_HWHM = engine.fdtdduration(MFILE, nargout=2)
    mat_t_0 = float(mat_t_0)
    mat_HWHM = float(mat_HWHM)
    engine.quit()

    py_t_0, py_HWHM = gaussian_pulse_parameters(wavelength_width, f_an, epsr)

    assert np.isclose(py_t_0, mat_t_0)
    assert np.isclose(py_HWHM, mat_HWHM)

    return


def test_max_allowable_timestep() -> None:
    """max_allowable_timestep should be equivalent to it's MATLAB counterpart, fdtdts"""
    # fdtdts needs to fetch from an input file, so provide one here.
    # arc_12 chosen because the extents are different
    MFILE = os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_12.m"
    )
    delta = {
        "x": 1300.0e-9 / 16.0,
        "y": 1300.0e-9 / 20.0,
        "z": 1300.0e-9 / 15.0,
    }

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    mat_dt_upper = engine.fdtdts(MFILE, nargout=1)
    mat_dt_upper = float(mat_dt_upper)
    engine.quit()

    py_dt_upper = max_allowable_timestep(delta)

    assert np.isclose(py_dt_upper, mat_dt_upper)

    return


def test_min_timesteps_fdtd() -> None:
    """min_timesteps with dt not provided should replicate the behaviour of it's MATLAB counterpart, minsteps_fdtd."""
    # minsteps_fdtd needs to fetch from an input file, so provide one here.
    # arc_12 chosen because it is an FDTD simulation
    MFILE = os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_12.m"
    )
    # Replicate file information
    delta = {
        "x": 1300.0e-9 / 16.0,
        "y": 1300.0e-9 / 20.0,
        "z": 1300.0e-9 / 15.0,
    }
    total_cells_in_K = 64
    K0_cell_index = 10
    wavelength_width = 120.0e-9
    eps_relative = 1.35**2
    dt = 0.95 / (
        np.sqrt(1.0 / delta["x"] ** 2 + 1.0 / delta["y"] ** 2 + 1.0 / delta["z"] ** 2)
        * (3e8 / 1.35)
    )
    f_an = np.arcsin(2.0 * np.pi * 2.997924580105029e8 * dt / (2.0 * 1300.0e-9)) / (
        np.pi * dt
    )

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    mat_min_timesteps = engine.minsteps_fdtd(MFILE)
    mat_min_timesteps = int(mat_min_timesteps)
    engine.quit()

    py_min_timesteps = min_timesteps_fdtd(
        wavelength_width,
        f_an,
        eps_relative,
        delta,
        "z",
        total_cells_in_K,
        K0_cell_index,
    )

    assert py_min_timesteps == mat_min_timesteps

    return


def test_min_timesteps_pstd() -> None:
    """min_timesteps with a value of dt provided should replicate the behaviour of it's MATLAB counterpart, minsteps_pstd."""
    # minsteps_pstd needs to fetch from an input file, so provide one here.
    # arc_01 chosen because it is a PSTD simulation
    MFILE = os.path.abspath(
        LOCATION_OF_THIS_FILE
        + "../../../data/input_generation/input_files/input_file_01.m"
    )
    # Replicate file information
    delta = {
        "x": 1300.0e-9 / 4.0,
        "y": 1300.0e-9 / 4.0,
        "z": 1300.0e-9 / 4.0,
    }
    total_cells_in_K = 256
    K0_cell_index = 10
    wavelength_width = 120.0e-9
    eps_relative = 1.35**2
    dt = 2.0 * 0.95 * 1.35 * delta["x"] / (np.pi * np.sqrt(2) * 3e8)
    f_an = np.arcsin(2.0 * np.pi * 2.997924580105029e8 * dt / (2.0 * 1300.0e-9)) / (
        np.pi * dt
    )

    engine = create_engine_for_testing(LOCATION_OF_THIS_FILE)
    mat_min_timesteps = engine.minsteps_fdtd(MFILE)
    mat_min_timesteps = int(mat_min_timesteps)
    engine.quit()

    py_min_timesteps = min_timesteps_fdtd(
        wavelength_width,
        f_an,
        eps_relative,
        delta,
        "z",
        total_cells_in_K,
        K0_cell_index,
    )

    assert py_min_timesteps == mat_min_timesteps

    return
