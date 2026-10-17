"""Selectable temporal normalizations for Iseult.

Iseult reports time in inverse electron plasma frequencies throughout -- the
figure title and the total-energy time axis both hard-code
``omega_pe^{-1}``.  This module generalizes that into a selectable unit, the
same way :mod:`length_norms` does for the spatial axes.

Both Tristan versions already hand Iseult a time in ``omega_pe^{-1}``, so one
conversion layer covers both: v2 derives it in ``data_loading`` as
``timestep * c / c_omp``, and v1 stores the same quantity directly.

Each normalization is expressed as **how many omega_pe^{-1} make up one unit
of it**, so a stored time converts with a single division::

    t_display = t_stored / scale

Unlike the spatial scale, this one must be applied at the *display* boundary
only.  ``time`` also feeds physics: ``shock_speed`` is derived from it as
cells per ``omega_pe^{-1}`` and multiplied back out in ``shock_loc``.
Rescaling the stored value would silently move the shock marker.

What may be used here
---------------------
Only quantities that Tristan writes into its **output** files, and that
``data_loading`` maps for both v1 and v2.  That restricts this module to
``c``, ``c_omp``, ``mi``, ``me`` and ``sigma``, which all load cleanly from
either version.  Input-only parameters are deliberately not used.

The gyration times assume Tristan's standard species convention: electrons
as the species with ``|ch/m| = 1`` and an upstream field of unit internal
amplitude, which is what makes ``sigma`` the electron magnetization.  A setup
that renumbered species or scaled the initial field differently would need
``ch`` and ``m`` read per species rather than assumed; ``particles:ch1..ch4``
are in the v2 output if that is ever wanted.
"""

import warnings

import numpy as np

import length_norms

# The parameter coercion helper is shared with length_norms rather than
# duplicated: the two modules face exactly the same hdf5 quirks (scalars that
# arrive as length-1 arrays, parameters that are absent entirely).
_scalar = length_norms._scalar

# Normalization ids.  The menu is ordered by id, so these also fix the order
# the options appear in: the physical units first, then time steps last.
#
# They are persisted in the .yml view files, so once a value has shipped its
# meaning must not change -- append new options rather than renumbering these.
# This numbering was settled before the feature merged and while no saved view
# carried a 'TimeNorm' key, which is the only reason reordering was free.
INVERSE_OMEGA_PE = 0
INVERSE_OMEGA_PI = 1
INVERSE_OMEGA_CE = 2
INVERSE_OMEGA_CI = 3
TIMESTEPS = 4

DEFAULT_NORM = INVERSE_OMEGA_PE

# Short names for the GUI option menu, indexed by normalization id.
NORM_NAMES = [
    'Electron plasma time',
    'Ion plasma time',
    'Electron gyration time',
    'Ion gyration time',
    'Time steps',
]

# The unit string that goes inside the axis label's square brackets, as LaTeX
# (without the enclosing $...$).
_UNIT_LATEX = {
    INVERSE_OMEGA_PE: r'\omega_{\rm pe}^{-1}',
    INVERSE_OMEGA_PI: r'\omega_{\rm pi}^{-1}',
    TIMESTEPS: r'{\rm steps}',
    INVERSE_OMEGA_CE: r'\Omega_{\rm ce}^{-1}',
    INVERSE_OMEGA_CI: r'\Omega_{\rm ci}^{-1}',
}

# A short plain-text unit, for Tk labels and titles that cannot render LaTeX.
_UNIT_PLAIN = {
    INVERSE_OMEGA_PE: '1/w_pe',
    INVERSE_OMEGA_PI: '1/w_pi',
    TIMESTEPS: 'steps',
    INVERSE_OMEGA_CE: '1/W_ce',
    INVERSE_OMEGA_CI: '1/W_ci',
}


# How many decimals the figure title shows, by normalization id; anything
# absent gets none.
#
# The rule is the unit's own coarseness: a unit whose ``scale`` exceeds one
# is coarser than ``omega_pe^{-1}``, so the same run spans fewer of them and
# rounding to an integer throws away real precision.  That is every derived
# unit except time steps, which are *finer* than ``omega_pe^{-1}`` (one step
# is ``c/c_omp`` of them, well under one) and so always print large.
#
#   omega_pi^{-1}  t/sqrt(mi/me)        /43 at the true mass ratio
#   Omega_ce^{-1}  t*sqrt(sigma)        /100 at sigma = 1e-4
#   Omega_ci^{-1}  t*sqrt(sigma)*me/mi  tens, or less early in a run
#
# This is per-unit rather than adaptive on the value: a rule like "one decimal
# below 100" would change the title's format partway through a run, which
# looks like flicker when the frames are played as a movie.
_UNIT_DECIMALS = {
    INVERSE_OMEGA_PI: 1,
    INVERSE_OMEGA_CE: 1,
    INVERSE_OMEGA_CI: 1,
}


def decimals(norm):
    """How many decimals the title shows for this unit."""
    return _UNIT_DECIMALS.get(sanitize(norm), 0)


# Which output parameters each unit's scale actually needs.  The call sites
# read these one at a time, out of an hdf5 file, every time an output is
# loaded -- so reading only what the selected unit uses keeps the default
# unit off the disk entirely, and stops a unit that does not care about a
# parameter from warning when it happens to be absent.
#
# Kept beside time_scale rather than at the call sites so the two cannot
# drift apart; test_required_params_matches_time_scale pins them together.
_UNIT_PARAMS = {
    INVERSE_OMEGA_PE: (),
    INVERSE_OMEGA_PI: ('mi', 'me'),
    INVERSE_OMEGA_CE: ('sigma',),
    INVERSE_OMEGA_CI: ('sigma', 'mi', 'me'),
    TIMESTEPS: ('c', 'c_omp'),
}


def required_params(norm):
    """The output parameters ``time_scale`` needs to compute this unit."""
    return _UNIT_PARAMS[sanitize(norm)]


def electron_magnetization(sigma):
    """Tristan's ``sigma``, which *is* the electron magnetization.

    Both versions use this convention.  For v2 it follows from the code: the
    Boris rotation angle per step is ``(ch/m) * sqrt(sigma) * B`` (combining
    ``boris_push.F08`` with ``B_norm = c**2 sqrt(sigma)/c_omp`` from
    ``global_namespace.F90``), and since masses and charges are given in units
    of ``m_e`` and ``q_e``, electrons have ``|ch/m| = 1``, leaving
    ``Omega_ce/omega_pe = sqrt(sigma)`` for a unit upstream field.

    Returns NaN when it cannot be used, which includes an unmagnetized run.
    A negative sigma is taken as a magnitude; Tristan writes one for some
    setups.
    """
    sig = abs(_scalar(sigma))
    if not np.isfinite(sig) or sig <= 0.0:
        return np.nan
    return sig


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


def time_scale(norm, c=None, c_omp=None, mi=None, me=None, sigma=None):
    """The selected time unit, measured in ``omega_pe^{-1}``.

    Returns ``(scale, used_norm, message)``.  ``used_norm`` is the
    normalization actually applied: if the requested one cannot be computed it
    falls back to ``omega_pe^{-1}`` and ``message`` explains why.  ``scale`` is
    always a finite positive float.
    """
    norm = sanitize(norm)

    if norm == INVERSE_OMEGA_PE:
        return 1.0, norm, None

    if norm == INVERSE_OMEGA_PI:
        # omega_p goes as 1/sqrt(m) at fixed density and charge, so one ion
        # plasma time is sqrt(mi/me) electron plasma times.
        ratio = length_norms.mass_ratio(mi, me)
        if not np.isfinite(ratio):
            return (1.0, INVERSE_OMEGA_PE,
                    'Cannot read mi/me; falling back to electron plasma times.')
        return np.sqrt(ratio), norm, None

    if norm == TIMESTEPS:
        # Both versions store t = step * c / c_omp, so one step is c/c_omp
        # electron plasma times.
        c_f, c_omp_f = _scalar(c), _scalar(c_omp)
        if (not np.isfinite(c_f) or not np.isfinite(c_omp_f)
                or c_f <= 0 or c_omp_f <= 0):
            return (1.0, INVERSE_OMEGA_PE,
                    'Cannot read c and c_omp; falling back to electron '
                    'plasma times.')
        return c_f / c_omp_f, norm, None

    if norm in (INVERSE_OMEGA_CE, INVERSE_OMEGA_CI):
        # Omega_ce/omega_pe = sqrt(sigma), so one electron gyration time is
        # 1/sqrt(sigma) electron plasma times.
        sig = electron_magnetization(sigma)
        if not np.isfinite(sig):
            return (1.0, INVERSE_OMEGA_PE,
                    'sigma is zero or unavailable, so the gyration time is '
                    'undefined; falling back to electron plasma times.')
        scale = 1.0 / np.sqrt(sig)

        if norm == INVERSE_OMEGA_CE:
            return scale, norm, None

        # Omega goes as q/m, so the ion gyration time is mi/me times longer.
        # Equivalently sigma_i = sigma * me/mi and one works through omega_pi;
        # that route is left out because it cancels to the same factor and two
        # expressions for one quantity can drift apart.
        ratio = length_norms.mass_ratio(mi, me)
        if not np.isfinite(ratio):
            return (1.0, INVERSE_OMEGA_PE,
                    'Cannot read mi/me; falling back to electron plasma times.')
        return scale * ratio, norm, None

    return 1.0, INVERSE_OMEGA_PE, None


def unit_latex(norm):
    """The unit on its own, e.g. ``\\omega_{\\rm pe}^{-1}`` (no ``$``)."""
    return _UNIT_LATEX[sanitize(norm)]


def unit_plain(norm):
    """A short plain-text unit, for labels that cannot render LaTeX."""
    return _UNIT_PLAIN[sanitize(norm)]


def axis_label(norm):
    """A full matplotlib axis label, e.g. ``$t\\ [\\omega_{\\rm pe}^{-1}]$``."""
    return r'$t\ \  [%s]$' % unit_latex(norm)


def title_time(norm, t_display):
    """The ``t = ... [unit]`` fragment used in the figure title."""
    return r't = %.*f $%s$' % (decimals(norm), float(t_display),
                               unit_latex(norm))


def warn_about(message):
    """Emit a one-shot warning about how a normalization was resolved."""
    if message is not None:
        warnings.warn(message, RuntimeWarning)
