import numpy as np
import pytest

import length_norms as ln


C_OMP = 10.0
MI, ME = 100.0, 1.0


def scale(norm, **kw):
    params = dict(c_omp=C_OMP, mi=MI, me=ME)
    params.update(kw)
    s, used, _ = ln.length_scale(norm, **params)
    return s, used


# --- the scales themselves ------------------------------------------------

def test_electron_skin_depth_is_c_omp():
    s, used = scale(ln.ELECTRON_SKIN_DEPTH)
    assert s == pytest.approx(C_OMP)
    assert used == ln.ELECTRON_SKIN_DEPTH


def test_ion_skin_depth_scales_as_sqrt_mass_ratio():
    s, used = scale(ln.ION_SKIN_DEPTH)
    assert s == pytest.approx(C_OMP * np.sqrt(MI / ME))
    assert used == ln.ION_SKIN_DEPTH


def test_cells_is_unity():
    s, used = scale(ln.CELLS)
    assert s == 1.0
    assert used == ln.CELLS


def test_pair_plasma_ion_and_electron_skin_depths_agree():
    s_e, _ = scale(ln.ELECTRON_SKIN_DEPTH, mi=1.0, me=1.0)
    s_i, _ = scale(ln.ION_SKIN_DEPTH, mi=1.0, me=1.0)
    assert s_e == pytest.approx(s_i)


# --- only output quantities may be used ----------------------------------

def test_length_scale_takes_only_output_parameters():
    """The signature is the guard: anything beyond c_omp/mi/me would mean
    reading a parameter Tristan does not reliably write to its output."""
    import inspect
    assert set(inspect.signature(ln.length_scale).parameters) == \
        {'norm', 'c_omp', 'mi', 'me'}


def test_exactly_three_normalizations():
    valid = [n for n in range(16) if ln.is_valid(n)]
    assert valid == [ln.ELECTRON_SKIN_DEPTH, ln.ION_SKIN_DEPTH, ln.CELLS]
    assert len(ln.NORM_NAMES) == len(valid)


# --- fallbacks: a scale must always come back finite and positive ---------

def test_ion_skin_depth_falls_back_when_mass_ratio_missing():
    s, used, msg = ln.length_scale(ln.ION_SKIN_DEPTH, C_OMP, mi=None, me=None)
    assert s == pytest.approx(C_OMP)
    assert used == ln.ELECTRON_SKIN_DEPTH
    assert msg is not None


@pytest.mark.parametrize('bad', [0.0, -1.0, np.nan, None, np.array([])])
def test_degenerate_masses_fall_back(bad):
    s, used, msg = ln.length_scale(ln.ION_SKIN_DEPTH, C_OMP, mi=bad, me=ME)
    assert s == pytest.approx(C_OMP)
    assert used == ln.ELECTRON_SKIN_DEPTH
    assert msg is not None


def test_missing_c_omp_degrades_to_cells():
    for bad in (None, 0.0, -1.0, np.nan):
        s, used, msg = ln.length_scale(ln.ELECTRON_SKIN_DEPTH, bad)
        assert s == 1.0
        assert used == ln.CELLS
        assert msg is not None


@pytest.mark.parametrize('norm', [0, 1, 2])
def test_scale_is_always_finite_and_positive(norm):
    s, _, _ = ln.length_scale(norm, C_OMP, mi=MI, me=ME)
    assert np.isfinite(s) and s > 0


# --- parameters arrive as length-1 arrays as often as scalars -------------

def test_array_valued_parameters():
    s, used, _ = ln.length_scale(ln.ION_SKIN_DEPTH,
                                 np.array([C_OMP]),
                                 mi=np.array([MI]),
                                 me=np.array([ME]))
    assert s == pytest.approx(C_OMP * np.sqrt(MI / ME))
    assert used == ln.ION_SKIN_DEPTH


def test_mass_ratio_helper():
    assert ln.mass_ratio(MI, ME) == pytest.approx(MI / ME)
    assert ln.mass_ratio(np.array([MI]), np.array([ME])) == pytest.approx(MI / ME)
    for bad in (None, 0.0, -1.0, np.nan):
        assert not np.isfinite(ln.mass_ratio(bad, ME))
        assert not np.isfinite(ln.mass_ratio(MI, bad))


# --- sanitizing what comes out of a .yml view file -----------------------

@pytest.mark.parametrize('junk', [None, 'nonsense', 99, -1, 3, 4, object()])
def test_sanitize_rejects_junk(junk):
    assert ln.is_valid(ln.sanitize(junk))


def test_sanitize_passes_through_valid_ids():
    for norm in (ln.ELECTRON_SKIN_DEPTH, ln.ION_SKIN_DEPTH, ln.CELLS):
        assert ln.sanitize(norm) == norm


def test_missing_key_defaults_to_electron_skin_depth():
    # Views saved before this feature existed have no LengthNorm key at all,
    # so the default must preserve Iseult's historical behaviour.
    assert ln.DEFAULT_NORM == ln.ELECTRON_SKIN_DEPTH


# --- labels ---------------------------------------------------------------

def test_axis_label_tracks_the_normalization():
    assert ln.axis_label(ln.ELECTRON_SKIN_DEPTH, 'x') == r'$x\ [c/\omega_{\rm pe}]$'
    assert ln.axis_label(ln.ION_SKIN_DEPTH, 'y') == r'$y\ [c/\omega_{\rm pi}]$'
    assert ln.axis_label(ln.CELLS, 'z') == r'$z\ [{\rm cells}]$'


def test_axis_label_prime():
    assert ln.axis_label(ln.ELECTRON_SKIN_DEPTH, 'x', prime=True) == \
        r'$x\prime\ [c/\omega_{\rm pe}]$'


def test_k_label_is_the_reciprocal_unit():
    assert ln.k_axis_label(ln.ELECTRON_SKIN_DEPTH) == r'$k \ [\omega_{\rm pe}/c]$'
    assert ln.k_axis_label(ln.ION_SKIN_DEPTH) == r'$k \ [\omega_{\rm pi}/c]$'


@pytest.mark.parametrize('norm', [0, 1, 2])
def test_every_norm_has_a_name_and_labels(norm):
    assert ln.NORM_NAMES[norm]
    for fn in (ln.unit_latex, ln.inverse_unit_latex, ln.unit_plain):
        assert fn(norm)
    assert ln.axis_label(norm, 'x').startswith('$')


# --- the conversion the call sites actually perform ----------------------

@pytest.mark.parametrize('norm', [0, 1, 2])
def test_axis_round_trips_through_cells(norm):
    """Panels go cells -> normalized and the viewport goes back again; the
    two must be exact inverses or zooming drifts."""
    s, _, _ = ln.length_scale(norm, C_OMP, mi=MI, me=ME)
    istep = 4.0
    cells = np.arange(20) * istep
    assert np.allclose((cells / s) * s, cells)


def test_changing_norm_rescales_a_viewport_consistently():
    """The GUI converts stored limits by old_scale/new_scale; check that this
    keeps the same physical region."""
    s_old, _, _ = ln.length_scale(ln.ELECTRON_SKIN_DEPTH, C_OMP, mi=MI, me=ME)
    s_new, _, _ = ln.length_scale(ln.ION_SKIN_DEPTH, C_OMP, mi=MI, me=ME)
    x_left_old, x_right_old = 20.0, 80.0
    ratio = s_old / s_new
    x_left_new, x_right_new = x_left_old * ratio, x_right_old * ratio
    # Same region, measured in cells, either way.
    assert x_left_old * s_old == pytest.approx(x_left_new * s_new)
    assert x_right_old * s_old == pytest.approx(x_right_new * s_new)
    # mi/me = 100 -> d_i = 10 d_e, so the numbers shrink tenfold.
    assert x_right_new == pytest.approx(8.0)


# --- the Tk option menu round-trips through display names ----------------

@pytest.mark.parametrize('norm', [0, 1, 2])
def test_name_round_trips(norm):
    assert ln.norm_from_name(ln.norm_name(norm)) == norm


@pytest.mark.parametrize('junk', ['', 'Gyroradius', None, 'cells'])
def test_unknown_name_falls_back_to_the_default(junk):
    assert ln.norm_from_name(junk) == ln.DEFAULT_NORM


def test_norm_name_sanitizes():
    assert ln.norm_name(99) == ln.NORM_NAMES[ln.DEFAULT_NORM]
