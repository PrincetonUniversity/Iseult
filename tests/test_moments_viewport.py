import numpy as np
import pytest
import sys
import pathlib

# Add src/ to python path
sys.path.append(str(pathlib.Path(__file__).resolve().parents[1] / 'src'))

from moments_panel import MomentsPanel
import length_norms

class MockParent:
    def __init__(self):
        self.MainParamDict = {
            'ColorMap': 'viridis',
            'NumOfRows': 1,
            'NumOfCols': 1,
            '2DSlicePlane': 0, # x-y plane
            'SetxLim': False,
            'SetyLim': False,
            'LinkSpatial': 1,
            'DoLorentzBoost': False,
            'GammaBoost': 0.0,
            'Average1D': False,
            'xSlice': 0,
            'ySlice': 0,
            'zSlice': 0,
            'shock_loc': 0.0,
            'shock_color': 'r',
            'LengthNorm': length_norms.DEFAULT_NORM,
        }
        self.ion_color = 'r'
        self.electron_color = 'b'
        self.figure = None
        self.SubPlotList = None
        # The panels divide by parent.length_scale rather than by c_omp. With
        # the default normalization this is just MockOutput.c_omp, so every
        # expectation below is unchanged.
        self.length_norm = length_norms.DEFAULT_NORM
        self.length_scale = 2.0

    def LengthLabel(self, axis = 'x', prime = False):
        return length_norms.axis_label(self.length_norm, axis, prime = prime)

    def KLabel(self):
        return length_norms.k_axis_label(self.length_norm)

class MockOutput:
    def __init__(self):
        self.c_omp = 2.0
        self.istep = 1.0
        self.me = 1.0
        self.mi = 16.0
        self.bx = np.zeros((10, 10, 10)) # shape[-1] is 10
        self.xe = np.array([2.0, 4.0, 6.0, 8.0, 10.0]) # in units of c_omp: 1, 2, 3, 4, 5
        self.ye = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
        self.ze = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        self.ue = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        self.ve = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        self.we = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        self.che = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
        
        self.xi = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
        self.yi = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
        self.zi = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        self.ui = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        self.vi = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        self.wi = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
        self.chi = np.array([1.0, 1.0, 1.0, 1.0, 1.0])

def test_moments_panel_no_viewport():
    parent = MockParent()
    panel = MomentsPanel(parent, (0,0), {'filter_by_viewport': False})
    output = MockOutput()
    panel.update_data(output)
    
    # Without viewport, all particles should be kept
    assert len(panel.xe) == 5
    assert len(panel.xi) == 5

def test_moments_panel_with_viewport():
    parent = MockParent()
    panel = MomentsPanel(parent, (0,0), {'filter_by_viewport': True})
    output = MockOutput()
    
    # We mock the active viewport on parent.SubPlotList
    class MockSubPlot:
        def __init__(self):
            self.chartType = 'DensityPlot'
            self.param_dict = {'twoD': True}
            
            class MockAxes:
                def get_xlim(self):
                    return (1.5, 3.5) # only xe/c_omp in [1.5, 3.5], i.e., index 2 and 3 (xe = 4.0, 6.0 -> 2.0, 3.0)
                def get_ylim(self):
                    return (1.5, 3.5) # only ye/c_omp in [1.5, 3.5]
            self.axes = MockAxes()

    parent.SubPlotList = [[MockSubPlot()]]
    panel.update_data(output)
    
    # viewport limits should filter to particles at indexes 1 and 2 (xe/c_omp = 2.0 and 3.0)
    assert len(panel.xe) == 2
    assert len(panel.xi) == 2
    assert np.allclose(panel.xe / panel.c_omp, [2.0, 3.0])
    
    # Moments should be correctly binned and not all zeros
    assert np.any(panel.ex != 0)
    assert np.any(panel.ey != 0)
    assert np.any(panel.ez != 0)
    assert np.any(panel.ix != 0)
    
    # Verify that particles are binned at the expected bin indices
    # xmin * c_omp = 3.0, bin_width = 0.04
    # x = 4.0: (4.0 - 3.0)/0.04 = 25
    # x = 6.0: (6.0 - 3.0)/0.04 = 75
    assert panel.ecounts[25] == 1.0
    assert panel.ecounts[75] == 1.0
    assert panel.ecounts[0] == 0.0

def test_moments_panel_with_viewport_high_offset():
    # A test case where xlim_min is large enough that the unit mismatch bug
    # would bin all particles out of bounds, resulting in all-zero moments.
    parent = MockParent()
    panel = MomentsPanel(parent, (0,0), {'filter_by_viewport': True})
    output = MockOutput()
    
    class MockSubPlot:
        def __init__(self):
            self.chartType = 'DensityPlot'
            self.param_dict = {'twoD': True}
            
            class MockAxes:
                def get_xlim(self):
                    return (3.5, 4.5) # only xe/c_omp in [3.5, 4.5], i.e., index 3 (xe = 8.0 -> 4.0)
                def get_ylim(self):
                    return (3.5, 4.5) # ye/c_omp = 4.0
            self.axes = MockAxes()

    parent.SubPlotList = [[MockSubPlot()]]
    panel.update_data(output)
    
    # Only index 3 particle is kept
    assert len(panel.xe) == 1
    assert len(panel.xi) == 1
    
    # Without the fix, this particle would bin to:
    # l = int((8.0 - 3.5)/((4.5-3.5)/100*2.0)) = int(4.5 / 0.02) = 225 >= 100 (out of bounds)
    # With the fix, it bins to:
    # l = int((8.0 - 3.5*2.0)/0.02) = int(1.0 / 0.02) = 50 (in bounds)
    assert panel.ecounts[50] == 1.0
    assert np.any(panel.ex != 0)

def test_moments_panel_with_unzoomed_viewport():
    parent = MockParent()
    panel = MomentsPanel(parent, (0,0), {'filter_by_viewport': True})
    output = MockOutput()
    
    # We mock the active viewport on parent.SubPlotList to span the entire box size
    class MockSubPlot:
        def __init__(self):
            self.chartType = 'DensityPlot'
            self.param_dict = {'twoD': True}
            
            class MockAxes:
                def get_xlim(self):
                    return (0.0, 5.0) # spans full x range of bx (bx.shape[-1] = 10 -> 10/2*1 = 5.0)
                def get_ylim(self):
                    return (0.0, 5.0) # spans full y range
            self.axes = MockAxes()

    parent.SubPlotList = [[MockSubPlot()]]
    panel.update_data(output)
    
    # Since it's unzoomed, filtering should be bypassed and keep all 5 particles
    assert len(panel.xe) == 5
    assert len(panel.xi) == 5


def test_moments_panel_vs_y_no_viewport():
    parent = MockParent()
    panel = MomentsPanel(parent, (0,0), {'oneD_axis': 'y', 'filter_by_viewport': False})
    output = MockOutput()
    output.xe = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    output.ye = np.array([2.0, 4.0, 6.0, 8.0, 10.0]) # in units of c_omp: 1, 2, 3, 4, 5
    output.xi = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    output.yi = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    panel.update_data(output)
    
    assert panel.GetPlotParam('oneD_axis') == 'y'
    assert len(panel.ye) == 5
    assert len(panel.yi) == 5
    assert panel.ymin == 1.0
    assert panel.ymax == 5.0
    assert np.isclose(panel.x_bins[0], 1.0)
    assert np.isclose(panel.x_bins[-1], 5.0)
    assert np.any(panel.ex != 0)
    assert np.any(panel.ey != 0)


def test_moments_panel_vs_y_with_viewport():
    parent = MockParent()
    panel = MomentsPanel(parent, (0,0), {'oneD_axis': 'y', 'filter_by_viewport': True})
    output = MockOutput()
    output.xe = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    output.ye = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    output.xi = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    output.yi = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    
    class MockSubPlot:
        def __init__(self):
            self.chartType = 'DensityPlot'
            self.param_dict = {'twoD': True}
            
            class MockAxes:
                def get_xlim(self):
                    return (1.5, 3.5)
                def get_ylim(self):
                    return (1.5, 3.5)
            self.axes = MockAxes()

    parent.SubPlotList = [[MockSubPlot()]]
    panel.update_data(output)
    
    assert len(panel.ye) == 2
    assert len(panel.yi) == 2
    assert np.allclose(panel.ye / panel.c_omp, [2.0, 3.0])
    assert panel.ymin == 1.5
    assert panel.ymax == 3.5
    assert np.isclose(panel.x_bins[0], 1.5)
    assert np.isclose(panel.x_bins[-1], 3.5)
    assert np.any(panel.ex != 0)


def test_moments_panel_set_v_min_max():
    import matplotlib.pyplot as plt
    parent = MockParent()
    parent.figure = plt.figure()
    parent.gs0 = [parent.figure.add_gridspec(1, 1)[0, 0]]
    parent.axes_extent = [0, 100, 0, 100]
    parent.MainParamDict.update({
        'xLabelPad': 5,
        'yLabelPad': 5,
        'AxLabelSize': 10,
        'NumFontSize': 10,
        'legendLabelSize': 10,
        'shock_loc': 0.0,
        'shock_color': 'r',
        'LinkSpatial': 0,
        'xLeft': 0,
        'xRight': 10,
        'xLimsRelative': False,
    })
    panel = MomentsPanel(parent, 0, {'set_v_min': True, 'v_min': -0.5, 'set_v_max': True, 'v_max': 0.8})
    output = MockOutput()
    panel.update_data(output)
    panel.draw()
    assert np.isclose(panel.axes.get_ylim()[0], -0.5)
    assert np.isclose(panel.axes.get_ylim()[1], 0.8)
    plt.close(parent.figure)


from fields_panel import FieldsPanel
from density_panel import DensPanel
from mag_panel import BPanel

def test_fields_density_mag_panel_vs_y():
    class FieldParent:
        def __init__(self):
            self.MainParamDict = {
                'ColorMap': 'viridis',
                'NumOfRows': 1,
                'NumOfCols': 1,
                '2DSlicePlane': 0,
                'SetxLim': False,
                'SetyLim': False,
                'LinkSpatial': 0,
                'Average1D': False,
                'xSlice': 3,
                'ySlice': 2,
                'zSlice': 0,
                'xLabelPad': 5,
                'yLabelPad': 5,
                'AxLabelSize': 10,
                'NumFontSize': 10,
                'HorizontalCbars': False,
                'ImageAspect': False,
                'CbarTickNum': 5,
                'shock_loc': 0.0,
                'shock_color': 'r',
                'FFT_color': 'g',
                'FFTLeft': 0,
                'FFTRight': 0,
                'FFTRelative': 0,
                'annotateTextSize': 10,
                'DivColorMap': 'PuOr'
            }
            self.btheta = 0.0
            self.bx0 = 1.0
            self.by0 = 0.0
            self.bz0 = 0.0
            self.b0 = 1.0
            self.SubPlotList = []
            self.figure = None
            # The panels divide by parent.length_scale rather than by c_omp;
            # with the default normalization it is just FieldOutput.c_omp.
            self.MainParamDict['LengthNorm'] = length_norms.DEFAULT_NORM
            self.length_norm = length_norms.DEFAULT_NORM
            self.length_scale = 2.0

        def LengthLabel(self, axis = 'x', prime = False):
            return length_norms.axis_label(self.length_norm, axis, prime = prime)

        def KLabel(self):
            return length_norms.k_axis_label(self.length_norm)

    class FieldOutput:
        def __init__(self):
            self.c_omp = 2.0
            self.istep = 1.0
            self.bx = np.arange(48, dtype=float).reshape((2, 4, 6))
            self.by = np.zeros((2, 4, 6), dtype=float)
            self.bz = np.zeros((2, 4, 6), dtype=float)
            self.ex = np.zeros((2, 4, 6), dtype=float)
            self.ey = np.zeros((2, 4, 6), dtype=float)
            self.ez = np.zeros((2, 4, 6), dtype=float)
            self.jx = np.zeros((2, 4, 6), dtype=float)
            self.jy = np.zeros((2, 4, 6), dtype=float)
            self.jz = np.zeros((2, 4, 6), dtype=float)
            self.dens = np.ones((2, 4, 6), dtype=float)
            self.densi = np.ones((2, 4, 6), dtype=float)
            self.ppc0 = 1.0

    parent = FieldParent()
    output = FieldOutput()

    # FieldsPanel vs y
    fpanel = FieldsPanel(parent, (0,0), {'twoD': 0, 'oneD_axis': 'y', 'fld_plot_type': 0})
    fpanel.update_data(output)
    assert fpanel.GetPlotParam('oneD_axis') == 'y'
    assert len(fpanel.yaxis_values) == 4
    assert np.allclose(fpanel.yaxis_values, np.arange(4) / 2.0)
    assert fpanel.GetPlotParam('spatial_y') == False
    assert fpanel.GetPlotParam('spatial_x') == False

    # DensPanel vs y
    dpanel = DensPanel(parent, (0,0), {'twoD': 0, 'oneD_axis': 'y'})
    dpanel.update_data(output)
    assert dpanel.GetPlotParam('oneD_axis') == 'y'
    assert len(dpanel.yaxis_values) == 4
    assert np.allclose(dpanel.yaxis_values, np.arange(4) / 2.0)
    assert dpanel.GetPlotParam('spatial_y') == False
    assert dpanel.GetPlotParam('spatial_x') == False

    # BPanel vs y
    bpanel = BPanel(parent, (0,0), {'twoD': 0, 'oneD_axis': 'y'})
    bpanel.update_data(output)
    assert bpanel.GetPlotParam('oneD_axis') == 'y'
    assert len(bpanel.yaxis_values) == 4
    assert np.allclose(bpanel.yaxis_values, np.arange(4) / 2.0)
    assert bpanel.GetPlotParam('spatial_y') == False
    assert bpanel.GetPlotParam('spatial_x') == False


from phase_panel import PhasePanel

def test_phase_panel_vs_y_no_viewport():
    parent = MockParent()
    panel = PhasePanel(parent, (0, 0), {'oneD_axis': 'y', 'filter_by_viewport': False, 'prtl_type': 0, 'mom_dim': 0})
    output = MockOutput()
    panel.update_data(output)

    assert panel.GetPlotParam('oneD_axis') == 'y'
    assert panel.GetPlotParam('spatial_x') == False
    assert panel.GetPlotParam('spatial_y') == False
    assert len(panel.x_values) == 5
    assert np.allclose(panel.x_values, output.yi / output.c_omp)
    assert panel.hist2d is not None
    assert panel.xmin == 0.0
    assert np.isclose(panel.xmax, output.bx.shape[0] / output.c_omp * output.istep)


def test_phase_panel_vs_y_with_viewport():
    parent = MockParent()
    panel = PhasePanel(parent, (0, 0), {'oneD_axis': 'y', 'filter_by_viewport': True, 'prtl_type': 0, 'mom_dim': 0})
    output = MockOutput()
    output.xe = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    output.ye = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    output.xi = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    output.yi = np.array([2.0, 4.0, 6.0, 8.0, 10.0])

    class MockSubPlot:
        def __init__(self):
            self.chartType = 'DensityPlot'
            self.param_dict = {'twoD': True}

            class MockAxes:
                def get_xlim(self):
                    return (1.5, 3.5)
                def get_ylim(self):
                    return (2.5, 4.5)
            self.axes = MockAxes()

    parent.SubPlotList = [[MockSubPlot()]]
    panel.update_data(output)

    # Particle index 2 (xi=6.0 -> 3.0, yi=6.0 -> 3.0) satisfies both x in [1.5, 3.5] and y in [2.5, 4.5]
    assert panel.xmin == 2.5
    assert panel.xmax == 4.5
    assert np.isclose(panel.hist2d[2][0], 2.5)
    assert np.isclose(panel.hist2d[2][1], 4.5)


def test_main_app_spatial_y_axis_sync():
    import matplotlib.pyplot as plt
    from main_app import MainApp

    class DummyApp:
        pass
    app = DummyApp()
    app.MainParamDict = {
        'LinkSpatial': 1,
        '2DSlicePlane': 0,
    }
    app.first_y = (0, 0)
    app.first_y_1d = (0, 1)
    app._syncing_y_axes = False
    app._spatial_y_cids = []

    fig, (ax_2d, ax_1d) = plt.subplots(1, 2)
    ax_2d.set_xlim(0, 10)
    ax_2d.set_ylim(2, 8)
    ax_1d.set_xlim(0, 20)
    ax_1d.set_ylim(-5, 5)

    class DummySubplot:
        def __init__(self, ax):
            class DummyGraph:
                def __init__(self, axes):
                    self.axes = axes
            self.graph = DummyGraph(ax)

    app.SubPlotList = [[DummySubplot(ax_2d), DummySubplot(ax_1d)]]
    app.canvas = None

    # Bind the methods from MainApp
    app._sync_2d_y_to_1d_x = MainApp._sync_2d_y_to_1d_x.__get__(app)
    app._sync_1d_x_to_2d_y = MainApp._sync_1d_x_to_2d_y.__get__(app)
    app.LinkSpatialAxes = MainApp.LinkSpatialAxes.__get__(app)

    # Initial linking should sync 1D x-lim to 2D y-lim
    app.LinkSpatialAxes()
    assert np.allclose(ax_1d.get_xlim(), (2, 8))

    # Zooming in 2D plot in y should update 1D plot in x
    ax_2d.set_ylim(3, 6)
    assert np.allclose(ax_1d.get_xlim(), (3, 6))

    # Zooming in 1D plot in x should update 2D plot in y
    ax_1d.set_xlim(4, 5)
    assert np.allclose(ax_2d.get_ylim(), (4, 5))

    # Test slice plane 2 (y is horizontal axis of 2D plot)
    app.MainParamDict['2DSlicePlane'] = 2
    ax_2d.set_xlim(1, 9)
    app.LinkSpatialAxes()
    assert np.allclose(ax_1d.get_xlim(), (1, 9))

    ax_2d.set_xlim(2.5, 7.5)
    assert np.allclose(ax_1d.get_xlim(), (2.5, 7.5))

    ax_1d.set_xlim(3.5, 6.5)
    assert np.allclose(ax_2d.get_xlim(), (3.5, 6.5))

    plt.close(fig)


import streamlines


def test_add_streamline_params():
    param_dict = {}
    streamlines.add_streamline_params(param_dict)
    assert "az_contours_fiducial_axis" in param_dict
    assert param_dict["az_contours_fiducial_axis"] == "X"
    assert "az_contours_fiducial_x" in param_dict
    assert param_dict["az_contours_fiducial_x"] == 200.0
    assert "show_az_contours" in param_dict
    assert "az_contours_count" in param_dict


def test_compute_vector_potential_2d_derivatives_x():
    """Verify that Az computed with fiducial_axis='X' satisfies curl(Az z) = B."""
    ny, nx = 40, 50
    y, x = np.mgrid[0:ny, 0:nx].astype(np.float64)
    kx = 2.0 * np.pi / nx
    ky = 2.0 * np.pi / ny
    bx = (-ky * np.sin(kx * x) * np.sin(ky * y)).astype(np.float32)
    by = (-kx * np.cos(kx * x) * np.cos(ky * y)).astype(np.float32)

    Az = streamlines.compute_vector_potential_2d(bx, by, stride=1, dx=1.0, fiducial_axis="X")

    dAz_dy = 0.5 * (Az[2:, 1:-1] - Az[:-2, 1:-1])
    dAz_dx = 0.5 * (Az[1:-1, 2:] - Az[1:-1, :-2])
    bx_int = bx[1:-1, 1:-1]
    by_int = by[1:-1, 1:-1]

    np.testing.assert_allclose(dAz_dy, bx_int, atol=1e-2, rtol=1e-2)
    np.testing.assert_allclose(-dAz_dx, by_int, atol=1e-2, rtol=1e-2)


def test_compute_vector_potential_2d_derivatives_y():
    """Verify that Az computed with fiducial_axis='Y' satisfies curl(Az z) = B."""
    ny, nx = 40, 50
    y, x = np.mgrid[0:ny, 0:nx].astype(np.float64)
    kx = 2.0 * np.pi / nx
    ky = 2.0 * np.pi / ny
    bx = (-ky * np.sin(kx * x) * np.sin(ky * y)).astype(np.float32)
    by = (-kx * np.cos(kx * x) * np.cos(ky * y)).astype(np.float32)

    Az = streamlines.compute_vector_potential_2d(bx, by, stride=1, dx=1.0, fiducial_axis="Y", x_fiducial=200.0)

    dAz_dy = 0.5 * (Az[2:, 1:-1] - Az[:-2, 1:-1])
    dAz_dx = 0.5 * (Az[1:-1, 2:] - Az[1:-1, :-2])
    bx_int = bx[1:-1, 1:-1]
    by_int = by[1:-1, 1:-1]

    np.testing.assert_allclose(dAz_dy, bx_int, atol=1e-2, rtol=1e-2)
    np.testing.assert_allclose(-dAz_dx, by_int, atol=1e-2, rtol=1e-2)


def test_lagrangian_tracker_fiducial_axes():
    # Test X axis
    tracker_x = streamlines.LagrangianFieldLineTracker(n_contours=10, fiducial_axis="X")
    xs_x, ys_x = tracker_x.init_step(1, nx=100, ny=80)
    assert len(xs_x) == 10
    assert len(ys_x) == 10
    np.testing.assert_allclose(ys_x, 40.0)
    assert xs_x[0] < xs_x[-1]

    # Test Y axis: with x_fiducial=200 c/omp, c_omp=1, istep=1, stride=1 in domain nx=300
    tracker_y = streamlines.LagrangianFieldLineTracker(n_contours=10, fiducial_axis="Y", x_fiducial=200.0)
    xs_y, ys_y = tracker_y.init_step(1, nx=300, ny=80, c_omp=1.0, istep=1.0, stride=1)
    assert len(xs_y) == 10
    assert len(ys_y) == 10
    np.testing.assert_allclose(xs_y, 200.0)
    assert ys_y[0] < ys_y[-1]

    # Test Y axis with downsampling/stride: c_omp=10, istep=2, stride=2 -> eff_dx=4 -> 200*10/4 = 500
    xs_y_scaled, _ = tracker_y.init_step(1, nx=1000, ny=80, c_omp=10.0, istep=2.0, stride=2)
    np.testing.assert_allclose(xs_y_scaled, 500.0)

    # Test Y axis clamping when domain nx is smaller than fiducial point
    xs_y_clamped, _ = tracker_y.init_step(1, nx=100, ny=80, c_omp=1.0, istep=1.0, stride=1)
    np.testing.assert_allclose(xs_y_clamped, 99.0)


def test_fiducial_point_fixed_in_expanding_domain():
    """Verify that as the active domain nx increases with time, the fiducial index remains at x=200 c/omp."""
    c_omp = 10.0
    istep = 1.0
    stride = 1
    x_fiducial = 200.0

    # Active domain expanding across timesteps from nx=2500 to nx=5000
    idx_t1 = streamlines._get_fiducial_x_idx(2500, c_omp=c_omp, istep=istep, stride=stride, x_fiducial=x_fiducial)
    idx_t2 = streamlines._get_fiducial_x_idx(3500, c_omp=c_omp, istep=istep, stride=stride, x_fiducial=x_fiducial)
    idx_t3 = streamlines._get_fiducial_x_idx(5000, c_omp=c_omp, istep=istep, stride=stride, x_fiducial=x_fiducial)

    # In middle of domain, indices would have been 1250, 1750, 2500 (drifting with expansion)
    # With fixed fiducial point at x=200 c/omp, index must stay exactly 2000!
    assert idx_t1 == 2000
    assert idx_t2 == 2000
    assert idx_t3 == 2000


def test_gauge_integral_cache_invalidation():
    class DummyParent:
        def __init__(self):
            self.PathDict = {"Flds": ["f1", "f2", "f3"], "Param": []}
            self.MainParamDict = {"2DSlicePlane": 0}
            self._az_gauge_cache = {}

    parent = DummyParent()
    # Cache initially empty, reading real files fails -> returns 0.0
    g_x = streamlines._get_gauge_integral(parent, 1, 50, 60, 0, stride=1, fiducial_axis="X")
    assert g_x == 0.0

    parent._az_gauge_cache = {
        "G": np.array([0.0, 1.0, 2.0]),
        "xref": 50,
        "ymid": 25,
        "fiducial_axis": "X",
        "x_fiducial": 200.0,
    }
    # When asking for X, cached value for step 2 should be 1.0
    val_x = streamlines._get_gauge_integral(parent, 2, 50, 60, 0, stride=1, fiducial_axis="X", x_fiducial=200.0)
    assert val_x == 1.0

    # When asking for Y with dummy fields, cache should invalidate and recompute
    val_y = streamlines._get_gauge_integral(parent, 2, 50, 60, 0, stride=1, fiducial_axis="Y", x_fiducial=200.0)
    assert val_y == 0.0
    assert parent._az_gauge_cache.get("fiducial_axis") == "Y"

    # When x_fiducial changes, cache should invalidate and recompute
    parent._az_gauge_cache = {
        "G": np.array([0.0, 1.0, 2.0]),
        "xref": 50,
        "ymid": 25,
        "fiducial_axis": "Y",
        "x_fiducial": 200.0,
    }
    val_y_new_fid = streamlines._get_gauge_integral(parent, 2, 50, 60, 0, stride=1, fiducial_axis="Y", x_fiducial=150.0)
    assert val_y_new_fid == 0.0
    assert parent._az_gauge_cache.get("x_fiducial") == 150.0


import vector_arrows


def test_vector_arrows_execute_update_handles_plots_without_show_vectors():
    class DummyMainApp:
        def __init__(self):
            self.MainParamDict = {'NumOfRows': 1, 'NumOfCols': 2}
            self.canvas = None
            self._zoom_timer_id = None
        def after(self, ms, func):
            func()
            return 1
        def after_cancel(self, timer_id):
            pass

    main_app = DummyMainApp()

    class DummySubplotPhase:
        chartType = 'PhasePlot'
        class DummyGraph:
            def GetPlotParam(self, pname):
                raise KeyError(pname)
        graph = DummyGraph()

    class DummySubplotDensity:
        chartType = 'DensityPlot'
        class DummyGraph:
            c_omp = 1.0
            istep = 1.0
            def GetPlotParam(self, pname):
                if pname == 'show_vectors':
                    return False
                raise KeyError(pname)
        graph = DummyGraph()

    main_app.SubPlotList = [[DummySubplotPhase(), DummySubplotDensity()]]

    class DummyPanel:
        parent = main_app
        _in_refresh = False

    class DummyCanvas:
        def draw_idle(self):
            pass

    main_app.canvas = DummyCanvas()

    # Calling on_limits_changed should not raise KeyError
    vector_arrows.on_limits_changed(DummyPanel())


def test_phase_plot_param_dict_has_set_p_max():
    import phase_plots
    import phase_panel

    assert 'set_p_max' in phase_plots.PhasePanel.plot_param_dict
    assert phase_plots.PhasePanel.plot_param_dict['set_p_max'] is False
    assert 'set_p_max' in phase_panel.PhasePanel.plot_param_dict
    assert phase_panel.PhasePanel.plot_param_dict['set_p_max'] is False


def test_subplot_wrapper_get_plot_param_fallback():
    from main_app import SubPlotWrapper

    class DummyParent:
        pass

    class DummyClass:
        plot_param_dict = {'key_in_class': 'val_from_class', 'set_something': False}

    # Instantiate SubPlotWrapper without full Tk initialization
    wrapper = SubPlotWrapper.__new__(SubPlotWrapper)
    wrapper.chartType = 'TestPlot'
    wrapper.PlotTypeDict = {'TestPlot': DummyClass}
    wrapper.PlotParamsDict = {'TestPlot': {'existing_key': 42}}

    # 1. Existing key returns directly
    assert wrapper.GetPlotParam('existing_key') == 42

    # 2. Missing key present in class plot_param_dict is fetched and cached
    assert wrapper.GetPlotParam('key_in_class') == 'val_from_class'
    assert wrapper.PlotParamsDict['TestPlot']['key_in_class'] == 'val_from_class'

    # 3. Missing boolean flag starting with set_ / show_ defaults to False
    assert wrapper.GetPlotParam('set_unknown_flag') is False
    assert wrapper.GetPlotParam('show_unknown_flag') is False

    # 4. Explicit default parameter
    assert wrapper.GetPlotParam('nonexistent_key', default='my_default') == 'my_default'




# --- the spatial normalization actually reaching a panel -------------------
# MockOutput has c_omp = 2 and mi/me = 16, so the ion skin depth is
# c_omp*sqrt(16) = 8 cells, four times the electron skin depth. Switching the
# normalization must therefore divide every plotted length by four.

def _panel_at(norm, scale):
    parent = MockParent()
    parent.MainParamDict['LengthNorm'] = norm
    parent.length_norm = norm
    parent.length_scale = scale
    panel = MomentsPanel(parent, (0, 0), {'filter_by_viewport': False})
    panel.update_data(MockOutput())
    return panel


def test_panel_picks_up_the_parents_scale():
    panel = _panel_at(length_norms.ION_SKIN_DEPTH, 8.0)
    assert panel.length_scale == 8.0


def test_ion_skin_depth_shrinks_the_x_extent_fourfold():
    d_e = _panel_at(length_norms.ELECTRON_SKIN_DEPTH, 2.0)
    d_i = _panel_at(length_norms.ION_SKIN_DEPTH, 8.0)
    # xe runs 2..10 cells, i.e. 1..5 d_e but only 0.25..1.25 d_i
    assert d_e.xmax == pytest.approx(5.0)
    assert d_i.xmax == pytest.approx(1.25)
    assert d_e.xmax / d_i.xmax == pytest.approx(4.0)


def test_cells_normalization_plots_raw_cell_coordinates():
    panel = _panel_at(length_norms.CELLS, 1.0)
    assert panel.xmax == pytest.approx(10.0)
    assert panel.xmin == pytest.approx(2.0)


def test_histogram_bins_stay_in_cells_whatever_the_unit():
    """The numba kernels take raw particle coordinates, so bin_width and xmin
    are handed to them in cells. The binned result must be identical however
    the axis is labelled."""
    d_e = _panel_at(length_norms.ELECTRON_SKIN_DEPTH, 2.0)
    d_i = _panel_at(length_norms.ION_SKIN_DEPTH, 8.0)
    assert np.allclose(d_e.ecounts, d_i.ecounts)
    assert np.allclose(d_e.ex, d_i.ex)
    assert np.allclose(d_e.icounts, d_i.icounts)


def test_x_bins_follow_the_unit():
    d_e = _panel_at(length_norms.ELECTRON_SKIN_DEPTH, 2.0)
    d_i = _panel_at(length_norms.ION_SKIN_DEPTH, 8.0)
    # same physical edges, expressed in different units
    assert np.allclose(d_e.x_bins * 2.0, d_i.x_bins * 8.0)


def test_axis_label_follows_the_unit():
    parent = MockParent()
    parent.length_norm = length_norms.ION_SKIN_DEPTH
    assert parent.LengthLabel('x') == r'$x\ [c/\omega_{\rm pi}]$'
    parent.length_norm = length_norms.CELLS
    assert parent.LengthLabel('x') == r'$x\ [{\rm cells}]$'
