import numpy as np


def zero_field(field_size: tuple[int]) -> list[np.ndarray]:
    """Return a zero-field of the requested size.

    :param field_size: The size of each of the component arrays
    :returns: A (complex) zero-field of the requested size
    """
    field = [np.zeros(field_size, dtype=complex)] * 3
    return field
