"""Selectable spatial length normalizations for Iseult.

Historically Iseult plotted every spatial axis in electron skin depths,
c/omega_pe, by dividing cell indices by ``c_omp`` -- Tristan's electron skin
depth measured in cells.  This module generalizes that single hard-coded
factor into a selectable normalization.

Each normalization is expressed as a **length scale measured in simulation
cells**, so call sites keep their original arithmetic and simply divide by
``scale`` where they used to divide by ``c_omp``::

    grid arrays:        np.arange(n) / scale * istep
    particle positions: x / scale
    and the inverse, norm -> cells:   x_norm * scale

What may be used here
---------------------
Only quantities that Tristan writes into its **output** files, and that
``data_loading`` maps for both Tristan v1 and v2.  That restricts this module
to ``c_omp``, ``mi`` and ``me``, which all load cleanly from either version.

Input-only parameters are deliberately not used.

Conventions
-----------
``c_omp`` is a Tristan *input* parameter that is echoed into the output param
file: the cold, non-relativistic electron skin depth of the initial upstream
plasma, in cells.  Both scales below derive from it and therefore inherit
that reference state.  They are fixed reference lengths, not local ones.

Note that only *spatial* uses of ``c_omp`` belong here.  ``c_omp`` is also
used to convert a timestep into cells (the ``omega_pe^-1`` time axis, the ExB
drift speed and the inductive gauge integral in ``streamlines.py``); those
must keep using the true electron skin depth regardless of what the user
selects for the spatial axes.
"""

import warnings

import numpy as np

# Normalization ids.  These are persisted in the .yml view files, so once a
# value has shipped its meaning must not change -- append new options rather
# than renumbering these.
ELECTRON_SKIN_DEPTH = 0
ION_SKIN_DEPTH = 1
CELLS = 2

DEFAULT_NORM = ELECTRON_SKIN_DEPTH

# Short names for the GUI radio buttons, indexed by normalization id.
NORM_NAMES = [
    'Electron skin depth',
    'Ion skin depth',
    'Cells',
]

# The unit string that goes inside the axis label's square brackets, as LaTeX
# (without the enclosing $...$).
_UNIT_LATEX = {
    ELECTRON_SKIN_DEPTH: r'c/\omega_{\rm pe}',
    ION_SKIN_DEPTH: r'c/\omega_{\rm pi}',
    CELLS: r'{\rm cells}',
}

# The reciprocal unit, for Fourier axes.
_INVERSE_UNIT_LATEX = {
    ELECTRON_SKIN_DEPTH: r'\omega_{\rm pe}/c',
    ION_SKIN_DEPTH: r'\omega_{\rm pi}/c',
    CELLS: r'{\rm cells}^{-1}',
}

# A short plain-text unit, for Tk labels that cannot render LaTeX.
_UNIT_PLAIN = {
    ELECTRON_SKIN_DEPTH: 'c/w_pe',
    ION_SKIN_DEPTH: 'c/w_pi',
    CELLS: 'cells',
}


def _scalar(val, default=np.nan):
    """Coerce a Tristan parameter to a plain float.

    Parameters are read out of the hdf5 files as length-1 arrays about as often
    as they are read out as scalars, and some of them can be missing entirely.
    """
    if val is None:
        return default
    try:
        arr = np.asarray(val, dtype=float).ravel()
    except (TypeError, ValueError):
        return default
    if arr.size == 0:
        return default
    out = float(arr[0])
    if not np.isfinite(out):
        return default
    return out


def is_valid(norm):
    return norm in _UNIT_LATEX


def sanitize(norm):
    """Map anything unexpected (bad .yml, old config) onto the default."""
    try:
        norm = int(norm)
    except (TypeError, ValueError):
        return DEFAULT_NORM
    return norm if is_valid(norm) else DEFAULT_NORM


def norm_from_name(name):
    """The id whose NORM_NAMES entry is ``name``, else the default.

    The Tk option menu works in display strings, so it hands back a name
    rather than an id.
    """
    try:
        return NORM_NAMES.index(name)
    except ValueError:
        return DEFAULT_NORM


def norm_name(norm):
    """The display string for a normalization id."""
    return NORM_NAMES[sanitize(norm)]


def mass_ratio(mi, me):
    """``mi/me`` from the output param file, or NaN if unusable."""
    mi_f, me_f = _scalar(mi), _scalar(me)
    if not np.isfinite(mi_f) or not np.isfinite(me_f):
        return np.nan
    if mi_f <= 0 or me_f <= 0:
        return np.nan
    return mi_f / me_f


def length_scale(norm, c_omp, mi=None, me=None):
    """The selected length scale, measured in simulation cells.

    Returns ``(scale, used_norm, message)``.  ``used_norm`` is the
    normalization actually applied: if the requested one cannot be computed
    it falls back to the electron skin depth and ``message`` explains why.
    ``scale`` is always a finite positive float.
    """
    norm = sanitize(norm)
    c_omp = _scalar(c_omp)
    if not np.isfinite(c_omp) or c_omp <= 0:
        # Without a usable c_omp the only honest axis is raw cells.
        return 1.0, CELLS, 'c_omp is missing or non-positive; plotting in cells.'

    if norm == ELECTRON_SKIN_DEPTH:
        return c_omp, norm, None

    if norm == CELLS:
        return 1.0, norm, None

    if norm == ION_SKIN_DEPTH:
        # d_i = c/omega_pi = d_e * sqrt(mi/me), since omega_p goes as
        # 1/sqrt(m) at fixed density and charge.
        ratio = mass_ratio(mi, me)
        if not np.isfinite(ratio):
            return (c_omp, ELECTRON_SKIN_DEPTH,
                    'Cannot read mi/me; falling back to electron skin depths.')
        return c_omp * np.sqrt(ratio), norm, None

    return c_omp, ELECTRON_SKIN_DEPTH, None


def unit_latex(norm):
    """The unit on its own, e.g. ``c/\\omega_{\\rm pe}`` (no ``$``)."""
    return _UNIT_LATEX[sanitize(norm)]


def inverse_unit_latex(norm):
    """The reciprocal unit, for Fourier-space axes."""
    return _INVERSE_UNIT_LATEX[sanitize(norm)]


def unit_plain(norm):
    """A short plain-text unit, for Tk labels that cannot render LaTeX."""
    return _UNIT_PLAIN[sanitize(norm)]


def axis_label(norm, axis='x', prime=False):
    """A full matplotlib axis label, e.g. ``$x\\ [c/\\omega_{\\rm pe}]$``.

    ``axis`` is the coordinate name ('x', 'y', 'z'); ``prime`` adds the prime
    used by the boosted-frame phase plots.
    """
    name = axis + (r'\prime' if prime else '')
    return r'$%s\ [%s]$' % (name, unit_latex(norm))


def k_axis_label(norm):
    """The Fourier-axis label matching this normalization."""
    return r'$k \ [%s]$' % inverse_unit_latex(norm)


def warn_about(message):
    """Emit a one-shot warning about how a normalization was resolved."""
    if message is not None:
        warnings.warn(message, RuntimeWarning)
