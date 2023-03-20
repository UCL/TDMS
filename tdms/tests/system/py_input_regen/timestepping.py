from typing import Literal, Tuple, Union

import numpy as np

from .physical_constants import C


def gaussian_pulse_parameters(
    wavelength_width: float, f_an: float, eps_relative: float
) -> Tuple[float, float]:
    """
    Calculates the parameters HWHM and t_0 for the Gaussian pulse employed in the FDTD code according to:
        G(t) = exp( - pi * ((t - t_0) / (HWHM))**2 ).

    HWHM is the Half-Width at Half-Maximum of the pulse.
    t_0 is the time delay: this delay ensures that when the pulse starts to be introduced, the envelope of the pulse has a magnitude of 10^{-8}. In turn, this ensures that there aren't any erroneous spectral components introduced as a result of having a sharp discontinuity in time.

    :param wavelength_width: FIXME
    :param f_an: FIXME
    :param eps_relative: Relative electric permittivity of the medium
    :returns: Tuple, (t_0, HWHM).
    """
    refractive_index = np.sqrt(eps_relative)
    lambda_an = C / (f_an * refractive_index)

    HWHM = (
        2.0
        * np.sqrt(np.log(2.0) / np.pi)
        * refractive_index
        * lambda_an**2
        / (C * wavelength_width)
    )
    t_0 = HWHM * np.sqrt(np.log(1e8) / np.pi)

    return t_0, HWHM


def max_allowable_timestep(delta: dict[str, float]) -> float:
    """
    Calculates the maximum allowable timestep (absolute time) subject to the stability criterion. This is
        dt_upper = ( c * sqrt(dx^{-2} + dy^{-2} + dz^{-2}) )^{-1}.

    :param delta: Yee cell extents in each of the coordinate dimensions
    :returns: Maximum allowable timestep
    """
    d = 1.0 / (np.array([delta["x"], delta["y"], delta["z"]]) ** 2)
    dt_upper = 1.0 / (C * np.sqrt(np.sum(d)))

    return dt_upper


def min_timesteps_fdtd(
    wavelength_width: float,
    f_an: float,
    eps_relative: float,
    delta: dict[str, float],
    travel_direction: Literal["x", "y", "z"],
    total_cells: int,
    interface_cell: int,
    dt: Union[None, float] = None,
) -> int:
    """
    Estimate the minimum number of (time)steps required in a simulation with timestep-interval dt. If not provided, dt is assumed to be .95 times the maximum allowable timestep.

    Assumes that the trailing edge of the gaussian pulse travels at the speed of light (dispersion free).

    Suppose that we introduce a pulsed source at the interface P (where P is one of I0, J0, or K0) which is perpendicular to the p-axis (where p is x, y, or z as corresponds to I0, J0, or K0). This pulse must travel from the interface to the opposite side of the grid, then reflect and travel back (through the interface) to the opposite side of the grid.

    This means that the pulse must travel through (2*total_cells - interface_cell_index) Yee cells before the simulation ends.

    :param wavelength_width: FIXME
    :param f_an: FIXME
    :param eps_relative: Relative electric permittivity of the medium.
    :param delta: Yee cell extents in each of the coordinate dimensions.
    :param travel_direction: The coordinate direction "p" that the pulse will be travelling in - this should also be a key of delta.
    :param total_cells: The total number of cells in the "P"-dimension.
    :param interface_cell: The P-index of the Yee-cell at which the interface sits.
    """
    pulse_time_delay, _ = gaussian_pulse_parameters(
        wavelength_width, f_an, eps_relative
    )
    max_timestep = max_allowable_timestep(delta)
    if dt is None:
        timestep_being_used = 0.95 * max_timestep
    else:
        timestep_being_used = dt

    # Compute timestep estimate, by computing the time for the pulse to propagate from the interface, to the other side of the grid, and back again.
    # We must travel through 2*total_cells - interface_cell Yee cells,
    # each with an extent of delta[travel_direction] in the p-direction,
    # at the speed of light in the medium
    distance_covered = (2 * total_cells - interface_cell) * delta[travel_direction]
    travel_speed = C / np.sqrt(eps_relative)

    min_time = distance_covered / travel_speed
    # Add an additional 2 * t_0 to the ellapsed time, to account for the possibiltity that the "front" edge of the pulse crosses the interface at time 0, and the "back" edge of the reflected pulse crosses reaches the edge of the grid at the final time.
    min_time += 2 * pulse_time_delay

    # Compute the minimum number of timesteps from this information, given the timestep
    min_timesteps = np.ceil(min_time / timestep_being_used)

    return min_timesteps
