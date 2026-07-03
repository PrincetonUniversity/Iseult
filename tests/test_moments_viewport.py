import numpy as np
import pytest
import sys
import pathlib

# Add src/ to python path
sys.path.append(str(pathlib.Path(__file__).resolve().parents[1] / 'src'))

from moments_panel import MomentsPanel

class MockParent:
    def __init__(self):
        self.MainParamDict = {
            'ColorMap': 'viridis',
            'NumOfRows': 1,
            'NumOfCols': 1,
            '2DSlicePlane': 0, # x-y plane
            'SetxLim': False,
            'SetyLim': False,
            'LinkSpatial': 1
        }
        self.figure = None
        self.SubPlotList = None

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
