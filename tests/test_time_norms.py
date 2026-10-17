import numpy as np
import pytest

import time_norms as tn


# Representative of the v1 test data: c=0.45, c_omp=8, mi/me=64.
C, C_OMP = 0.45, 8.0
MI, ME = 64.0, 1.0


def scale(norm, **kw):
    params = dict(c=C, c_omp=C_OMP, mi=MI, me=ME)
    params.update(kw)
    s, used, _ = tn.time_scale(norm, **params)
    return s, used


# --- the scales themselves ------------------------------------------------

def test_electron_plasma_time_is_unity():
    s, used = scale(tn.INVERSE_OMEGA_PE)
    assert s == 1.0
    assert used == tn.INVERSE_OMEGA_PE


def test_ion_plasma_time_scales_as_sqrt_mass_ratio():
    s, used = scale(tn.INVERSE_OMEGA_PI)
    assert s == pytest.approx(np.sqrt(MI / ME))
    assert used == tn.INVERSE_OMEGA_PI


def test_timestep_is_c_over_c_omp():
    s, used = scale(tn.TIMESTEPS)
    assert s == pytest.approx(C / C_OMP)
    assert used == tn.TIMESTEPS


def test_timesteps_recover_the_raw_step_count():
    """Both versions store t = step * c/c_omp, so dividing by the scale must
    give back whole steps. The v1 fixture sits at step 41*interval = 205000."""
    s, _, _ = tn.time_scale(tn.TIMESTEPS, c=C, c_omp=C_OMP)
    t_stored = 205000.0 * (C / C_OMP)
    assert t_stored / s == pytest.approx(205000.0)


def test_pair_plasma_ion_and_electron_times_agree():
    s_e, _ = scale(tn.INVERSE_OMEGA_PE, mi=1.0, me=1.0)
    s_i, _ = scale(tn.INVERSE_OMEGA_PI, mi=1.0, me=1.0)
    assert s_e == pytest.approx(s_i)


# --- only output quantities may be used ----------------------------------

# --- what each unit needs read from disk ---------------------------------

@pytest.mark.parametrize('norm', [0, 1, 2, 3, 4])
def test_required_params_matches_time_scale(norm):
    """The drift guard. required_params exists so a call site can skip the
    reads a unit has no use for, which is only safe while it names exactly
    the parameters time_scale consumes: supplying them must succeed, and
    withholding any one of them must fall back."""
    full = dict(c=C, c_omp=C_OMP, mi=MI, me=ME, sigma=SIGMA)
    needed = tn.required_params(norm)

    given = {key: full[key] for key in needed}
    s, used, msg = tn.time_scale(norm, **given)
    assert used == norm, 'required_params is missing a parameter'
    assert msg is None
    assert np.isfinite(s) and s > 0

    for withheld in needed:
        partial = dict(given)
        partial[withheld] = None
        _, used, msg = tn.time_scale(norm, **partial)
        assert used == tn.INVERSE_OMEGA_PE, \
            f'{withheld} is listed but unused by norm {norm}'
        assert msg is not None


def test_the_default_unit_needs_nothing_read():
    """Iseult's historical unit must not cost a single parameter read, since
    the scale is 1.0 by definition and this runs on every output load."""
    assert tn.required_params(tn.INVERSE_OMEGA_PE) == ()
    s, used, msg = tn.time_scale(tn.INVERSE_OMEGA_PE)
    assert (s, used, msg) == (1.0, tn.INVERSE_OMEGA_PE, None)


def test_required_params_are_output_quantities():
    """Same restriction as time_scale's signature: nothing input-only."""
    import inspect
    allowed = set(inspect.signature(tn.time_scale).parameters) - {'norm'}
    for norm in range(len(tn.NORM_NAMES)):
        assert set(tn.required_params(norm)) <= allowed


@pytest.mark.parametrize('junk', [None, 'nonsense', 99, -1])
def test_required_params_sanitizes_like_everything_else(junk):
    assert tn.required_params(junk) == tn.required_params(tn.DEFAULT_NORM)


def test_time_scale_takes_only_output_parameters():
    """The signature is the guard: anything beyond c/c_omp/mi/me would mean
    reading a parameter Tristan does not reliably write to its output."""
    import inspect
    assert set(inspect.signature(tn.time_scale).parameters) == \
        {'norm', 'c', 'c_omp', 'mi', 'me', 'sigma'}


def test_exactly_five_normalizations():
    valid = [n for n in range(16) if tn.is_valid(n)]
    assert valid == [tn.INVERSE_OMEGA_PE, tn.INVERSE_OMEGA_PI,
                     tn.INVERSE_OMEGA_CE, tn.INVERSE_OMEGA_CI, tn.TIMESTEPS]
    assert len(tn.NORM_NAMES) == len(valid)


# --- fallbacks: a scale must always come back finite and positive ---------

@pytest.mark.parametrize('bad', [None, 0.0, -1.0, np.nan, np.array([])])
def test_ion_plasma_time_falls_back_on_bad_masses(bad):
    s, used, msg = tn.time_scale(tn.INVERSE_OMEGA_PI, c=C, c_omp=C_OMP,
                                 mi=bad, me=ME)
    assert s == 1.0
    assert used == tn.INVERSE_OMEGA_PE
    assert msg is not None


@pytest.mark.parametrize('bad', [None, 0.0, -1.0, np.nan, np.array([])])
def test_timesteps_fall_back_on_bad_c_omp(bad):
    s, used, msg = tn.time_scale(tn.TIMESTEPS, c=C, c_omp=bad)
    assert s == 1.0
    assert used == tn.INVERSE_OMEGA_PE
    assert msg is not None


@pytest.mark.parametrize('norm', [0, 1, 2, 3, 4])
def test_scale_is_always_finite_and_positive(norm):
    s, _, _ = tn.time_scale(norm, c=C, c_omp=C_OMP, mi=MI, me=ME)
    assert np.isfinite(s) and s > 0


def test_array_valued_parameters():
    s, used, _ = tn.time_scale(tn.TIMESTEPS, c=np.array([C]),
                               c_omp=np.array([C_OMP]))
    assert s == pytest.approx(C / C_OMP)
    assert used == tn.TIMESTEPS


# --- sanitizing what comes out of a .yml view file -----------------------

@pytest.mark.parametrize('junk', [None, 'nonsense', 99, -1, 5, object()])
def test_sanitize_rejects_junk(junk):
    assert tn.is_valid(tn.sanitize(junk))


def test_missing_key_defaults_to_electron_plasma_time():
    # Views saved before this feature existed have no TimeNorm key at all, so
    # the default must preserve Iseult's historical behaviour.
    assert tn.DEFAULT_NORM == tn.INVERSE_OMEGA_PE


@pytest.mark.parametrize('norm', [0, 1, 2, 3, 4])
def test_name_round_trips(norm):
    assert tn.norm_from_name(tn.norm_name(norm)) == norm


def test_unknown_name_falls_back_to_the_default():
    assert tn.norm_from_name('Gyrofrequency') == tn.DEFAULT_NORM


# --- labels ---------------------------------------------------------------

def test_axis_label_tracks_the_normalization():
    assert tn.axis_label(tn.INVERSE_OMEGA_PE) == r'$t\ \  [\omega_{\rm pe}^{-1}]$'
    assert tn.axis_label(tn.INVERSE_OMEGA_PI) == r'$t\ \  [\omega_{\rm pi}^{-1}]$'
    assert tn.axis_label(tn.TIMESTEPS) == r'$t\ \  [{\rm steps}]$'


def test_axis_label_matches_the_previous_hardcoded_string():
    """The default must render exactly what total_energy_plots used to
    hardcode, so the figure is unchanged until a unit is selected."""
    assert tn.axis_label(tn.DEFAULT_NORM) == r'$t\ \  [\omega^{-1}_{pe}]$'.replace(
        r'\omega^{-1}_{pe}', r'\omega_{\rm pe}^{-1}')


def test_title_fragment():
    assert tn.title_time(tn.TIMESTEPS, 205000.004) == r't = 205000 ${\rm steps}$'


# --- title precision ------------------------------------------------------
# A unit coarser than omega_pe^-1 (scale > 1) fits fewer of itself into the
# same run, so an integer title throws away real precision.

@pytest.mark.parametrize('norm', [tn.INVERSE_OMEGA_PI,
                                  tn.INVERSE_OMEGA_CE,
                                  tn.INVERSE_OMEGA_CI])
def test_coarse_units_keep_a_decimal(norm):
    assert tn.decimals(norm) == 1


def test_decimals_follow_the_scale_being_coarse():
    """The rule, rather than the table: anything with scale > 1 is coarser
    than omega_pe^-1 and must keep a decimal."""
    for norm in range(len(tn.NORM_NAMES)):
        s, used, _ = tn.time_scale(norm, c=C, c_omp=C_OMP, mi=MI, me=ME,
                                   sigma=SIGMA)
        assert used == norm
        assert tn.decimals(norm) == (1 if s > 1.0 else 0)


def test_small_electron_gyration_time_is_not_rounded_away():
    """sigma << 1 makes Omega_ce^-1 a long time, so few of them fit in a
    run; t = 3 would be hiding most of the information."""
    t_stored = 350.0
    s, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CE, sigma=1e-4)
    assert tn.title_time(tn.INVERSE_OMEGA_CE, t_stored / s) == \
        r't = 3.5 $\Omega_{\rm ce}^{-1}$'


def test_ion_plasma_time_keeps_a_decimal_early_in_a_run():
    t_stored = 100.0
    s, _, _ = tn.time_scale(tn.INVERSE_OMEGA_PI, mi=1836.0, me=1.0)
    assert tn.title_time(tn.INVERSE_OMEGA_PI, t_stored / s) == \
        r't = 2.3 $\omega_{\rm pi}^{-1}$'


@pytest.mark.parametrize('norm', [0, 1, 2, 3, 4])
def test_every_norm_has_a_name_and_labels(norm):
    assert tn.NORM_NAMES[norm]
    assert tn.unit_latex(norm) and tn.unit_plain(norm)
    assert tn.axis_label(norm).startswith('$')


def test_unit_tables_are_kept_in_id_order():
    """A readability guard, not behaviour: dict order does not affect
    anything, but the unit tables are read next to NORM_NAMES, which is
    indexed by id, so they are kept in the same order."""
    ids = list(range(len(tn.NORM_NAMES)))
    assert list(tn._UNIT_LATEX) == ids
    assert list(tn._UNIT_PLAIN) == ids


# --- the conversion the call sites actually perform ----------------------

@pytest.mark.parametrize('norm', [0, 1, 2, 3, 4])
def test_display_conversion_round_trips(norm):
    s, _, _ = tn.time_scale(norm, c=C, c_omp=C_OMP, mi=MI, me=ME)
    t_stored = np.array([0.0, 11531.25, 23062.5])
    assert np.allclose((t_stored / s) * s, t_stored)


# --- the gyration times --------------------------------------------------
# sigma IS the electron magnetization in both Tristan versions, so
# Omega_ce/omega_pe = sqrt(sigma) and Omega_ci = Omega_ce * me/mi.

SIGMA = 0.01   # the v1 fixture's value


def test_electron_gyration_time_is_one_over_sqrt_sigma():
    s, used, msg = tn.time_scale(tn.INVERSE_OMEGA_CE, sigma=SIGMA)
    assert s == pytest.approx(1.0 / np.sqrt(SIGMA))
    assert used == tn.INVERSE_OMEGA_CE
    assert msg is None


def test_ion_gyration_time_is_mass_ratio_longer():
    s_e, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CE, sigma=SIGMA)
    s_i, used, _ = tn.time_scale(tn.INVERSE_OMEGA_CI, sigma=SIGMA, mi=MI, me=ME)
    assert s_i == pytest.approx(s_e * MI / ME)
    assert used == tn.INVERSE_OMEGA_CI


def test_ion_gyration_agrees_with_the_sigma_i_route():
    """Omega_ci/omega_pe can be reached two ways: via the mass ratio, or via
    sigma_i = sigma*me/mi and omega_pi. They must agree, which is what pins
    sigma as the *electron* magnetization rather than a total."""
    via_mass_ratio = np.sqrt(SIGMA) * ME / MI
    sigma_i = SIGMA * ME / MI
    via_sigma_i = np.sqrt(sigma_i) * np.sqrt(ME / MI)
    assert via_mass_ratio == pytest.approx(via_sigma_i)
    s_i, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CI, sigma=SIGMA, mi=MI, me=ME)
    assert 1.0 / s_i == pytest.approx(via_mass_ratio)


def test_v1_fixture_numbers():
    """The real v1 fixture: t = 11531.25, sigma = 0.01, mi/me = 64."""
    t = 11531.25
    s_e, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CE, sigma=0.01)
    s_i, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CI, sigma=0.01, mi=64.0, me=1.0)
    assert t / s_e == pytest.approx(1153.125)
    assert t / s_i == pytest.approx(18.017578125)


def test_pair_plasma_ion_and_electron_gyration_times_agree():
    s_e, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CE, sigma=0.3, mi=1.0, me=1.0)
    s_i, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CI, sigma=0.3, mi=1.0, me=1.0)
    assert s_e == pytest.approx(s_i)


@pytest.mark.parametrize('norm', [tn.INVERSE_OMEGA_CE, tn.INVERSE_OMEGA_CI])
@pytest.mark.parametrize('bad', [0.0, None, np.nan, np.array([])])
def test_gyration_falls_back_when_unmagnetized(norm, bad):
    s, used, msg = tn.time_scale(norm, sigma=bad, mi=MI, me=ME)
    assert s == 1.0
    assert used == tn.INVERSE_OMEGA_PE
    assert msg is not None


def test_negative_sigma_is_treated_as_a_magnitude():
    pos, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CE, sigma=SIGMA)
    neg, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CE, sigma=-SIGMA)
    assert pos == pytest.approx(neg)


def test_ion_gyration_needs_the_masses_too():
    s, used, msg = tn.time_scale(tn.INVERSE_OMEGA_CI, sigma=SIGMA)
    assert s == 1.0
    assert used == tn.INVERSE_OMEGA_PE
    assert msg is not None


def test_electron_gyration_does_not_need_the_masses():
    """Omega_ce/omega_pe = sqrt(sigma) alone, so it must work without mi/me."""
    s, used, _ = tn.time_scale(tn.INVERSE_OMEGA_CE, sigma=SIGMA)
    assert used == tn.INVERSE_OMEGA_CE
    assert s == pytest.approx(1.0 / np.sqrt(SIGMA))


def test_zero_electron_mass_does_not_raise():
    """me = 0 would be a ZeroDivisionError in a naive mass-ratio."""
    s, used, msg = tn.time_scale(tn.INVERSE_OMEGA_CI, sigma=SIGMA, mi=64.0, me=0.0)
    assert s == 1.0 and used == tn.INVERSE_OMEGA_PE and msg is not None


def test_gyration_labels():
    assert tn.axis_label(tn.INVERSE_OMEGA_CE) == r'$t\ \  [\Omega_{\rm ce}^{-1}]$'
    assert tn.axis_label(tn.INVERSE_OMEGA_CI) == r'$t\ \  [\Omega_{\rm ci}^{-1}]$'


# --- menu order and title precision --------------------------------------

def test_time_steps_is_the_last_menu_entry():
    """The menu is ordered by id, so the physical units come first and the
    code-level one last."""
    assert tn.NORM_NAMES[-1] == 'Time steps'
    assert tn.TIMESTEPS == len(tn.NORM_NAMES) - 1


def test_physical_units_precede_time_steps():
    for norm in (tn.INVERSE_OMEGA_PE, tn.INVERSE_OMEGA_PI,
                 tn.INVERSE_OMEGA_CE, tn.INVERSE_OMEGA_CI):
        assert norm < tn.TIMESTEPS


def test_ion_gyration_title_shows_one_decimal():
    assert tn.decimals(tn.INVERSE_OMEGA_CI) == 1
    assert tn.title_time(tn.INVERSE_OMEGA_CI, 18.017578) == \
        r't = 18.0 $\Omega_{\rm ci}^{-1}$'


@pytest.mark.parametrize('norm', [tn.INVERSE_OMEGA_PE, tn.TIMESTEPS])
def test_other_units_title_stays_integer(norm):
    """omega_pe^-1 is the reference, and a step is c/c_omp of one -- so the
    step count is always the largest number on offer. Only these two can
    safely round to an integer."""
    assert tn.decimals(norm) == 0
    assert '.' not in tn.title_time(norm, 1234.56).split('$')[0]


def test_title_precision_does_not_depend_on_the_value():
    """Per-unit, not adaptive: the format must not change as a run advances,
    or the title flickers when frames are played as a movie."""
    for t in (0.04, 4.0, 400.0, 40000.0):
        assert tn.title_time(tn.INVERSE_OMEGA_CI, t).count('.') == 1
        assert '.' not in tn.title_time(tn.INVERSE_OMEGA_PE, t).split('$')[0]


def test_decimals_sanitizes_junk():
    assert tn.decimals(99) == 0
    assert tn.decimals(None) == 0
