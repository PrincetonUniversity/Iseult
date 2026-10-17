"""MainApp and Oengus time-normalization plumbing.

Stored limits that live on a time axis must follow the selected unit.

The total-energy panels are the only ones that plot time on x, and they keep
x_min/x_max in the *displayed* unit. Switching the unit without converting
them leaves the plot framed on an interval that no longer exists -- the
ratio between omega_pe^-1 and time steps is c/c_omp, so the numbers are off
by more than an order of magnitude for a typical run.

These drive MainApp.RescaleTimeLimits and SettingsFrame.TimeNormChanged as
unbound functions against stand-ins: both only touch a handful of attributes,
and a real MainApp needs a display.
"""

import numpy as np
import pytest

import oengus
import time_norms as tn
from main_app import MainApp, SettingsFrame


class FakeTkVar:
    """Enough of a Tk StringVar for the entry-refreshing path."""

    def __init__(self, value=''):
        self.value = value

    def get(self):
        return self.value

    def set(self, value):
        self.value = value


class FakeSettingsWindow:
    def __init__(self, x_min, x_max):
        self.Xmin = FakeTkVar(str(x_min))
        self.Xmax = FakeTkVar(str(x_max))


class FakeGraph:
    def __init__(self, settings_window=None):
        self.settings_window = settings_window


class FakeSubPlot:
    """Stands in for SubPlotWrapper: a per-chart-type parameter store."""

    def __init__(self, chart_type='TotalEnergyPlot', x_min=0.0, x_max=10.0,
                 settings_window=None):
        self.chartType = chart_type
        self.graph = FakeGraph(settings_window)
        self.PlotParamsDict = {
            'TotalEnergyPlot': {'x_min': x_min, 'x_max': x_max,
                                'y_min': 1e-5, 'y_max': 1.0},
            'PhasePlot': {'x_min': 3.0, 'x_max': 7.0},
        }


def app(*subplots):
    """A stand-in carrying only what RescaleTimeLimits reads."""
    obj = type('Stub', (), {})()
    obj.SubPlotList = [list(subplots)]
    return obj


# --- the conversion itself ------------------------------------------------

def test_limits_are_converted_into_the_new_unit():
    panel = FakeSubPlot(x_min=100.0, x_max=2000.0)
    MainApp.RescaleTimeLimits(app(panel), 0.5)
    params = panel.PlotParamsDict['TotalEnergyPlot']
    assert params['x_min'] == pytest.approx(50.0)
    assert params['x_max'] == pytest.approx(1000.0)


def test_the_same_physical_interval_is_kept():
    """The point of the conversion: the limits must still bracket the same
    stored time, which is what the user is actually looking at."""
    t_stored = 11531.25
    old_scale, _, _ = tn.time_scale(tn.INVERSE_OMEGA_PE)
    new_scale, _, _ = tn.time_scale(tn.TIMESTEPS, c=0.45, c_omp=8.0)

    panel = FakeSubPlot(x_min=0.0, x_max=t_stored / old_scale)
    MainApp.RescaleTimeLimits(app(panel), old_scale / new_scale)

    assert panel.PlotParamsDict['TotalEnergyPlot']['x_max'] == \
        pytest.approx(t_stored / new_scale)


def test_switching_there_and_back_restores_the_limits():
    panel = FakeSubPlot(x_min=12.5, x_max=987.5)
    pe, _, _ = tn.time_scale(tn.INVERSE_OMEGA_PE)
    ci, _, _ = tn.time_scale(tn.INVERSE_OMEGA_CI, sigma=0.01, mi=64.0, me=1.0)

    MainApp.RescaleTimeLimits(app(panel), pe / ci)
    MainApp.RescaleTimeLimits(app(panel), ci / pe)

    params = panel.PlotParamsDict['TotalEnergyPlot']
    assert params['x_min'] == pytest.approx(12.5)
    assert params['x_max'] == pytest.approx(987.5)


# --- what must be left alone ---------------------------------------------

def test_y_limits_and_other_chart_types_are_untouched():
    panel = FakeSubPlot()
    MainApp.RescaleTimeLimits(app(panel), 4.0)
    assert panel.PlotParamsDict['TotalEnergyPlot']['y_min'] == pytest.approx(1e-5)
    assert panel.PlotParamsDict['TotalEnergyPlot']['y_max'] == pytest.approx(1.0)
    # A phase plot's x is a length, so the time unit must not touch it.
    assert panel.PlotParamsDict['PhasePlot']['x_min'] == pytest.approx(3.0)
    assert panel.PlotParamsDict['PhasePlot']['x_max'] == pytest.approx(7.0)


@pytest.mark.parametrize('ratio', [0.0, -2.0, np.nan, np.inf])
def test_an_unusable_ratio_changes_nothing(ratio):
    """time_scale can only fail by falling back, but UpdateTimeScale is
    called before the first output is loaded too; a bad ratio must be inert
    rather than zeroing the limits."""
    panel = FakeSubPlot(x_min=100.0, x_max=2000.0)
    MainApp.RescaleTimeLimits(app(panel), ratio)
    params = panel.PlotParamsDict['TotalEnergyPlot']
    assert params['x_min'] == pytest.approx(100.0)
    assert params['x_max'] == pytest.approx(2000.0)


def test_panels_not_showing_total_energy_are_still_converted():
    """Limits are stored per chart type, so a panel switched to total energy
    after the unit change must not come up framed in the old unit."""
    panel = FakeSubPlot(chart_type='PhasePlot', x_min=100.0, x_max=2000.0)
    MainApp.RescaleTimeLimits(app(panel), 0.5)
    assert panel.PlotParamsDict['TotalEnergyPlot']['x_min'] == pytest.approx(50.0)


# --- the open settings window --------------------------------------------

def test_an_open_settings_window_is_refreshed():
    window = FakeSettingsWindow(100.0, 2000.0)
    panel = FakeSubPlot(x_min=100.0, x_max=2000.0, settings_window=window)
    MainApp.RescaleTimeLimits(app(panel), 0.5)
    assert float(window.Xmin.get()) == pytest.approx(50.0)
    assert float(window.Xmax.get()) == pytest.approx(1000.0)


def test_a_closed_settings_window_is_not_an_error():
    panel = FakeSubPlot(settings_window=None)
    MainApp.RescaleTimeLimits(app(panel), 2.0)
    assert panel.PlotParamsDict['TotalEnergyPlot']['x_max'] == pytest.approx(20.0)


# --- the handler that drives it ------------------------------------------

class FakeNavStack:
    def __init__(self):
        self.cleared = False

    def clear(self):
        self.cleared = True


class FakeToolbar:
    def __init__(self):
        self._nav_stack = FakeNavStack()


class FakeParent:
    """Stands in for MainApp as TimeNormChanged sees it."""

    def __init__(self, norm=tn.INVERSE_OMEGA_PE, scale=1.0, new_scale=1.0,
                 resolved_norm=None):
        self.MainParamDict = {'TimeNorm': norm}
        self.time_scale = scale
        self.time_norm = norm
        self._new_scale = new_scale
        self._resolved = resolved_norm
        self.toolbar = FakeToolbar()
        self.ratios = []
        self.renew_calls = []

    def UpdateTimeScale(self):
        requested = self.MainParamDict['TimeNorm']
        self.time_norm = requested if self._resolved is None else self._resolved
        self.time_scale = self._new_scale
        return self.time_scale

    def RescaleTimeLimits(self, ratio):
        self.ratios.append(ratio)

    def RenewCanvas(self, keep_view=True, ForceRedraw=False):
        self.renew_calls.append((keep_view, ForceRedraw))


def frame(parent, chosen):
    obj = type('Stub', (), {})()
    obj.parent = parent
    obj.timeNormVar = FakeTkVar(tn.norm_name(chosen))
    return obj


def test_handler_passes_the_scale_ratio_and_drops_the_stale_view():
    parent = FakeParent(scale=1.0, new_scale=0.45 / 8.0)
    SettingsFrame.TimeNormChanged(frame(parent, tn.TIMESTEPS))

    assert parent.ratios == [pytest.approx(1.0 / (0.45 / 8.0))]
    assert parent.toolbar._nav_stack.cleared
    # keep_view would restore the old time limits verbatim.
    assert parent.renew_calls == [(False, True)]


def test_handler_does_nothing_when_the_unit_is_unchanged():
    parent = FakeParent(norm=tn.TIMESTEPS)
    SettingsFrame.TimeNormChanged(frame(parent, tn.TIMESTEPS))
    assert parent.ratios == []
    assert parent.renew_calls == []
    assert not parent.toolbar._nav_stack.cleared


def test_handler_uses_the_unit_actually_applied_after_a_fallback():
    """An unmagnetized run cannot show gyration times, so time_norms falls
    back; the ratio must come from the scale really in use (1.0), not from
    the one that was asked for."""
    parent = FakeParent(scale=1.0, new_scale=1.0,
                        resolved_norm=tn.INVERSE_OMEGA_PE)
    settings = frame(parent, tn.INVERSE_OMEGA_CE)
    SettingsFrame.TimeNormChanged(settings)

    assert parent.MainParamDict['TimeNorm'] == tn.INVERSE_OMEGA_PE
    assert settings.timeNormVar.get() == tn.norm_name(tn.INVERSE_OMEGA_PE)
    assert parent.ratios == [pytest.approx(1.0)]


# --- only the parameters the selected unit needs are read ----------------

class CountingOutput:
    """Stands in for an OutputPoint: every attribute is an hdf5 read."""

    def __init__(self):
        self.reads = []

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        self.reads.append(name)
        return 1.0


def oengus_stub(norm):
    obj = type('Stub', (), {})()
    obj.MainParamDict = {'TimeNorm': norm}
    return obj


@pytest.mark.parametrize('norm', [0, 1, 2, 3, 4])
def test_update_time_scale_reads_only_what_the_unit_needs(norm):
    """UpdateTimeScale runs on every output load, so a unit must not pull
    parameters it has no use for -- and must not warn about them going
    missing either."""
    output = CountingOutput()
    oengus.Oengus.UpdateTimeScale(oengus_stub(norm), output)
    assert set(output.reads) == set(tn.required_params(norm))


def test_the_default_unit_touches_no_parameter():
    output = CountingOutput()
    stub = oengus_stub(tn.INVERSE_OMEGA_PE)
    assert oengus.Oengus.UpdateTimeScale(stub, output) == 1.0
    assert output.reads == []
    assert stub.time_norm == tn.INVERSE_OMEGA_PE
