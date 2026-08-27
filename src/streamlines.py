#!/usr/bin/env python
"""This file contains the functions to control and generate magnetic
streamline and vector potential (Az) contour overlays for Iseult.
"""

import tkinter as Tk
import numpy as np
import matplotlib
import matplotlib.patches


def add_streamline_plot_keys(panel):
    """Add magnetic and electric field components needed for streamlines or Az contours."""
    slice_plane = panel.parent.MainParamDict["2DSlicePlane"]
    if slice_plane == 0:  # x-y plane
        panel.arrs_needed.extend(["bx", "by"])
        if panel.GetPlotParam("show_az_contours") and panel.GetPlotParam("az_contours_lagrangian"):
            panel.arrs_needed.append("ez")
    elif slice_plane == 1:  # x-z plane
        panel.arrs_needed.extend(["bx", "bz"])
        if panel.GetPlotParam("show_az_contours") and panel.GetPlotParam("az_contours_lagrangian"):
            panel.arrs_needed.append("ey")
    elif slice_plane == 2:  # y-z plane
        panel.arrs_needed.extend(["by", "bz"])
        if panel.GetPlotParam("show_az_contours") and panel.GetPlotParam("az_contours_lagrangian"):
            panel.arrs_needed.append("ex")


def add_streamline_params(param_dictionary):
    """Add data to the parameter dictionary for controlling streamlines and Az contours.

    Parameters
    ----------
    param_dictionary : dict
        The dictionary to add elements to.
    """
    # Streamlines parameters
    param_dictionary["show_streamlines"] = False
    param_dictionary["streamlines_stride"] = 10
    param_dictionary["streamlines_density"] = 1
    param_dictionary["streamlines_color"] = "black"

    # Az Contours parameters
    param_dictionary["show_az_contours"] = False
    param_dictionary["az_contours_count"] = 18
    param_dictionary["az_contours_color"] = "black"
    param_dictionary["az_contours_width"] = 1.0
    param_dictionary["az_contours_stride"] = 1
    param_dictionary["az_contours_gauge_tracking"] = False
    param_dictionary["az_contours_lagrangian"] = True


def add_streamline_buttons(settings, panel, starting_row):
    """Add buttons and fields for streamlines and Az contours to the settings pane.

    Parameters
    ----------
    settings :
        The settings object for the panel being drawn in
    panel :
        The Panel object the streamlines/contours are being drawn in
    starting_row :
        The row to start adding settings at
    """

    # --- 1. Streamlines Section ---
    Tk.ttk.Label(settings.frm, text="Streamline Settings: (2D only)").grid(
        row=starting_row, columnspan=2, sticky=Tk.W
    )

    settings.show_streamlines = Tk.BooleanVar()
    settings.show_streamlines.set(settings.parent.GetPlotParam("show_streamlines"))
    Tk.ttk.Checkbutton(
        settings.frm,
        text="Display Streamlines",
        variable=settings.show_streamlines,
        command=lambda: __show_streamline_handler(settings, panel),
    ).grid(row=starting_row + 1, column=0, sticky=Tk.W)

    settings.streamlines_stride = Tk.IntVar(
        value=settings.parent.GetPlotParam("streamlines_stride")
    )
    Tk.ttk.Label(settings.frm, text="Stride:").grid(row=starting_row + 2, column=0, sticky=Tk.W)
    Tk.ttk.Entry(settings.frm, textvariable=settings.streamlines_stride, width=7).grid(
        row=starting_row + 2, column=0, sticky=Tk.E
    )

    settings.streamlines_density = Tk.DoubleVar(
        value=settings.parent.GetPlotParam("streamlines_density")
    )
    Tk.ttk.Label(settings.frm, text="Line Density:").grid(
        row=starting_row + 3, column=0, sticky=Tk.W
    )
    Tk.ttk.Entry(settings.frm, textvariable=settings.streamlines_density, width=7).grid(
        row=starting_row + 3, column=0, sticky=Tk.E
    )

    settings.streamlines_color = Tk.StringVar(
        value=settings.parent.GetPlotParam("streamlines_color")
    )
    Tk.ttk.Label(settings.frm, text="Line Color:").grid(
        row=starting_row + 4, column=0, sticky=Tk.W
    )
    Tk.ttk.Entry(settings.frm, textvariable=settings.streamlines_color, width=7).grid(
        row=starting_row + 4, column=0, sticky=Tk.E
    )

    # --- 2. Az Contours Section ---
    az_row = starting_row + 5
    Tk.ttk.Label(settings.frm, text="Az Contour Settings: (2D only)").grid(
        row=az_row, columnspan=2, sticky=Tk.W
    )

    settings.show_az_contours = Tk.BooleanVar()
    settings.show_az_contours.set(settings.parent.GetPlotParam("show_az_contours"))
    Tk.ttk.Checkbutton(
        settings.frm,
        text="Display Az Contours",
        variable=settings.show_az_contours,
        command=lambda: __show_az_contours_handler(settings, panel),
    ).grid(row=az_row + 1, column=0, sticky=Tk.W)

    settings.az_contours_lagrangian = Tk.BooleanVar()
    settings.az_contours_lagrangian.set(settings.parent.GetPlotParam("az_contours_lagrangian"))
    Tk.ttk.Checkbutton(
        settings.frm,
        text="Lagrangian Tracking",
        variable=settings.az_contours_lagrangian,
        command=lambda: __show_az_contours_handler(settings, panel),
    ).grid(row=az_row + 1, column=1, sticky=Tk.W)

    settings.az_contours_count = Tk.IntVar(
        value=settings.parent.GetPlotParam("az_contours_count")
    )
    Tk.ttk.Label(settings.frm, text="Count:").grid(row=az_row + 2, column=0, sticky=Tk.W)
    Tk.ttk.Entry(settings.frm, textvariable=settings.az_contours_count, width=7).grid(
        row=az_row + 2, column=0, sticky=Tk.E
    )

    settings.az_contours_width = Tk.DoubleVar(
        value=settings.parent.GetPlotParam("az_contours_width")
    )
    Tk.ttk.Label(settings.frm, text="Line Width:").grid(row=az_row + 2, column=1, sticky=Tk.W)
    Tk.ttk.Entry(settings.frm, textvariable=settings.az_contours_width, width=7).grid(
        row=az_row + 2, column=1, sticky=Tk.E
    )

    settings.az_contours_color = Tk.StringVar(
        value=settings.parent.GetPlotParam("az_contours_color")
    )
    Tk.ttk.Label(settings.frm, text="Line Color:").grid(row=az_row + 3, column=0, sticky=Tk.W)
    Tk.ttk.Entry(settings.frm, textvariable=settings.az_contours_color, width=7).grid(
        row=az_row + 3, column=0, sticky=Tk.E
    )

    settings.az_contours_stride = Tk.IntVar(
        value=settings.parent.GetPlotParam("az_contours_stride")
    )
    Tk.ttk.Label(settings.frm, text="Stride:").grid(row=az_row + 3, column=1, sticky=Tk.W)
    Tk.ttk.Entry(settings.frm, textvariable=settings.az_contours_stride, width=7).grid(
        row=az_row + 3, column=1, sticky=Tk.E
    )


def __show_streamline_handler(settings, panel):
    """Handle what happens when the `show_streamlines` button is toggled."""
    settings.parent.SetPlotParam(
        "show_streamlines", settings.show_streamlines.get(), update_plot=False, NeedsRedraw=True
    )
    if not settings.parent.GetPlotParam("show_streamlines"):
        remove_streamlines(panel)
    else:
        settings.parent.parent.LoadAllKeys()

    settings.parent.parent.canvas.draw()
    settings.parent.parent.canvas.get_tk_widget().update_idletasks()


def __show_az_contours_handler(settings, panel):
    """Handle what happens when the `show_az_contours` or `az_contours_lagrangian` button is toggled."""
    settings.parent.SetPlotParam(
        "show_az_contours", settings.show_az_contours.get(), update_plot=False, NeedsRedraw=True
    )
    settings.parent.SetPlotParam(
        "az_contours_lagrangian", settings.az_contours_lagrangian.get(), update_plot=False, NeedsRedraw=True
    )
    if not settings.parent.GetPlotParam("show_az_contours"):
        remove_az_contours(panel)
    else:
        # Invalidate base level cache and tracker so levels recompute cleanly
        if hasattr(settings.parent.parent, "_az_base_level_cache"):
            settings.parent.parent._az_base_level_cache = None
        if hasattr(settings.parent.parent, "_lagrangian_tracker"):
            settings.parent.parent._lagrangian_tracker = None
        settings.parent.parent.LoadAllKeys()

    settings.parent.parent.canvas.draw()
    settings.parent.parent.canvas.get_tk_widget().update_idletasks()


def streamlines_callback(settings, update_plot=True):
    """Update streamline settings from GUI entries."""
    if not settings.parent.GetPlotParam("show_streamlines"):
        update_plot = False

    if settings.streamlines_stride.get() != settings.parent.plot_param_dict["streamlines_stride"]:
        settings.parent.plot_param_dict["streamlines_stride"] = settings.streamlines_stride.get()
        settings.parent.SetPlotParam("streamlines_stride", settings.streamlines_stride.get(), update_plot=update_plot)

    if settings.streamlines_density.get() != settings.parent.plot_param_dict["streamlines_density"]:
        settings.parent.plot_param_dict["streamlines_density"] = settings.streamlines_density.get()
        settings.parent.SetPlotParam("streamlines_density", settings.streamlines_density.get(), update_plot=update_plot)

    if settings.streamlines_color.get() != settings.parent.plot_param_dict["streamlines_color"]:
        settings.parent.plot_param_dict["streamlines_color"] = settings.streamlines_color.get()
        settings.parent.SetPlotParam("streamlines_color", settings.streamlines_color.get(), update_plot=update_plot)


def az_contours_callback(settings, update_plot=True):
    """Update Az contour settings from GUI entries."""
    if not settings.parent.GetPlotParam("show_az_contours"):
        update_plot = False

    # Count
    if hasattr(settings, "az_contours_count") and settings.az_contours_count.get() != settings.parent.plot_param_dict.get("az_contours_count"):
        try:
            val = int(settings.az_contours_count.get())
            settings.parent.plot_param_dict["az_contours_count"] = val
            settings.parent.SetPlotParam("az_contours_count", val, update_plot=update_plot)
            if hasattr(settings.parent.parent, "_az_base_level_cache"):
                settings.parent.parent._az_base_level_cache = None
            if hasattr(settings.parent.parent, "_lagrangian_tracker"):
                settings.parent.parent._lagrangian_tracker = None
        except ValueError:
            pass

    # Color
    if hasattr(settings, "az_contours_color") and settings.az_contours_color.get() != settings.parent.plot_param_dict.get("az_contours_color"):
        val = settings.az_contours_color.get()
        settings.parent.plot_param_dict["az_contours_color"] = val
        settings.parent.SetPlotParam("az_contours_color", val, update_plot=update_plot)

    # Width
    if hasattr(settings, "az_contours_width") and settings.az_contours_width.get() != settings.parent.plot_param_dict.get("az_contours_width"):
        try:
            val = float(settings.az_contours_width.get())
            settings.parent.plot_param_dict["az_contours_width"] = val
            settings.parent.SetPlotParam("az_contours_width", val, update_plot=update_plot)
        except ValueError:
            pass

    # Stride
    if hasattr(settings, "az_contours_stride") and settings.az_contours_stride.get() != settings.parent.plot_param_dict.get("az_contours_stride"):
        try:
            val = int(settings.az_contours_stride.get())
            settings.parent.plot_param_dict["az_contours_stride"] = val
            settings.parent.SetPlotParam("az_contours_stride", val, update_plot=update_plot)
            if hasattr(settings.parent.parent, "_az_base_level_cache"):
                settings.parent.parent._az_base_level_cache = None
            if hasattr(settings.parent.parent, "_lagrangian_tracker"):
                settings.parent.parent._lagrangian_tracker = None
        except ValueError:
            pass

    # Lagrangian Tracking Toggle
    if hasattr(settings, "az_contours_lagrangian") and settings.az_contours_lagrangian.get() != settings.parent.plot_param_dict.get("az_contours_lagrangian"):
        val = bool(settings.az_contours_lagrangian.get())
        settings.parent.plot_param_dict["az_contours_lagrangian"] = val
        settings.parent.SetPlotParam("az_contours_lagrangian", val, update_plot=update_plot)
        if hasattr(settings.parent.parent, "_az_base_level_cache"):
            settings.parent.parent._az_base_level_cache = None
        if hasattr(settings.parent.parent, "_lagrangian_tracker"):
            settings.parent.parent._lagrangian_tracker = None


# ------------------------------------------------------------------------------
# 3. STREAMLINES DRAWING
# ------------------------------------------------------------------------------

def draw_streamlines(panel):
    """Draw streamlines using matplotlib.streamplot."""
    stride = panel.GetPlotParam("streamlines_stride")
    slice_plane = panel.parent.MainParamDict["2DSlicePlane"]
    if slice_plane == 0:  # x-y plane
        bx_name, by_name = "bx", "by"
        slice_tuple = np.s_[panel.parent.zSlice, ::stride, ::stride]
    elif slice_plane == 1:  # x-z plane
        bx_name, by_name = "bx", "bz"
        slice_tuple = np.s_[::stride, panel.parent.ySlice, ::stride]
    else:
        bx_name, by_name = "by", "bz"
        slice_tuple = np.s_[::stride, ::stride, panel.parent.xSlice]

    if bx_name not in panel.parent.DataDict or by_name not in panel.parent.DataDict:
        return

    bx = panel.parent.DataDict[bx_name][slice_tuple]
    by = panel.parent.DataDict[by_name][slice_tuple]

    if bx.ndim != 2 or by.ndim != 2:
        return

    xmin = getattr(panel, "xmin", 0.0)
    xmax = getattr(panel, "xmax", float(bx.shape[1]))
    ymin = getattr(panel, "ymin", 0.0)
    ymax = getattr(panel, "ymax", float(bx.shape[0]))

    coords_x = np.linspace(xmin, xmax, bx.shape[1])
    coords_y = np.linspace(ymin, ymax, bx.shape[0])
    coords_x, coords_y = np.meshgrid(coords_x, coords_y)

    panel.FigWrap.streamlines = panel.FigWrap.graph.axes.streamplot(
        coords_x,
        coords_y,
        bx,
        by,
        density=panel.GetPlotParam("streamlines_density"),
        color=panel.GetPlotParam("streamlines_color"),
    )


def refresh_streamlines(panel):
    """Refresh the streamlines."""
    remove_streamlines(panel)
    draw_streamlines(panel)


def remove_streamlines(panel):
    """Remove streamlines."""
    if hasattr(panel.FigWrap, "streamlines") and panel.FigWrap.streamlines is not None:
        try:
            panel.FigWrap.streamlines.lines.remove()
        except Exception:
            pass
        for artist in panel.FigWrap.graph.axes.get_children():
            if isinstance(artist, matplotlib.patches.FancyArrowPatch):
                try:
                    artist.remove()
                except Exception:
                    pass
        panel.FigWrap.streamlines = None


# ------------------------------------------------------------------------------
# 4. AZ (VECTOR POTENTIAL) CONTOURS
# ------------------------------------------------------------------------------

def compute_vector_potential_2d(bx, by, stride=1):
    """Computes the 2D magnetic vector potential Az(x,y) from Bx, By components.

    By definition:
        Bx =  ∂Az / ∂y
        By = -∂Az / ∂x

    Parameters
    ----------
    bx : np.ndarray (ny, nx)
    by : np.ndarray (ny, nx)
    stride : int or float
        Grid stride factor

    Returns
    -------
    np.ndarray (ny, nx)
        Vector potential Az on the 2D grid
    """
    ny, nx = bx.shape
    ymid = ny // 2
    Az = np.zeros((ny, nx), dtype=np.float32)
    dx = float(stride)

    # 1. Integrate along midplane x: dAz = -By dx -> Az(x, ymid) = - ∫ By(x, ymid) dx
    Az[ymid, 1:] = -np.cumsum(0.5 * (by[ymid, 1:] + by[ymid, :-1]) * dx)

    # 2. Integrate upward in y: dAz = Bx dy -> Az(y, x) = Az(ymid, x) + ∫ Bx dy
    if ymid < ny - 1:
        dAz_up = 0.5 * (bx[ymid+1:, :] + bx[ymid:-1, :]) * dx
        Az[ymid+1:, :] = Az[ymid, :] + np.cumsum(dAz_up, axis=0)

    # 3. Integrate downward in y: dAz = -Bx dy
    if ymid > 0:
        dAz_down = -0.5 * (bx[ymid-1::-1, :] + bx[ymid:0:-1, :]) * dx
        Az[ymid-1::-1, :] = Az[ymid, :] + np.cumsum(dAz_down, axis=0)

    return Az


def find_dynamic_levels(level0, delta, minval, maxval):
    """Finds all levels in { level0 + n * delta } within [minval, maxval]."""
    if delta <= 0 or maxval <= minval:
        return np.array([level0], dtype=np.float32)
    n_start = int(np.floor((minval - level0) / delta))
    n_end = int(np.ceil((maxval - level0) / delta))
    n_vals = np.arange(n_start, n_end + 1)
    lvls = level0 + n_vals * delta
    lvls = lvls[(lvls >= minval) & (lvls <= maxval)]
    if len(lvls) == 0:
        lvls = np.linspace(minval, maxval, 5)
    return np.sort(lvls)


def sample_bilinear_vec(grid, xs, ys):
    """Vectorized 2D bilinear interpolation on a (ny, nx) scalar grid."""
    ny, nx = grid.shape
    xc = np.clip(xs, 0, nx - 1)
    yc = np.clip(ys, 0, ny - 1)
    ix = np.clip(np.floor(xc).astype(int), 0, nx - 2)
    iy = np.clip(np.floor(yc).astype(int), 0, ny - 2)
    fx = xc - ix
    fy = yc - iy
    return ((1.0 - fx) * (1.0 - fy) * grid[iy, ix] +
            fx * (1.0 - fy) * grid[iy, ix + 1] +
            (1.0 - fx) * fy * grid[iy + 1, ix] +
            fx * fy * grid[iy + 1, ix + 1])


class LagrangianFieldLineTracker:
    """Tracks Lagrangian fluid markers via sub-stepped ExB advection across timesteps."""

    def __init__(self, n_contours=18):
        self.n_contours = n_contours
        self.history = {}  # step -> (xs, ys)
        self.last_step = None
        self.last_flds = None

    def reset(self, n_contours=18):
        self.n_contours = n_contours
        self.history = {}
        self.last_step = None
        self.last_flds = None

    def init_step(self, step, nx, ny):
        ymid = ny // 2
        xs = np.linspace(nx * 0.04, nx * 0.48, self.n_contours).astype(np.float32)
        ys = np.full(len(xs), float(ymid), dtype=np.float32)
        self.history[step] = (xs, ys)
        return xs, ys

    def get_markers(self, parent, cur_step, bx, by, ez, c_omp):
        ny, nx = bx.shape
        ymid = ny // 2
        xmid = nx // 2

        if cur_step in self.history:
            self.last_step = cur_step
            self.last_flds = (bx, by, ez)
            return self.history[cur_step]

        # Try to advect from last_step if available
        if self.last_step is not None and self.last_step in self.history and self.last_flds is not None:
            prev_xs, prev_ys = self.history[self.last_step]
            prev_bx, prev_by, prev_ez = self.last_flds

            # Retrieve physical time difference
            param_paths = parent.PathDict.get("Param", [])
            t_cur, t_prev = None, None
            try:
                import h5py
                if cur_step - 1 < len(param_paths):
                    with h5py.File(param_paths[cur_step - 1], "r") as fp:
                        if "time" in fp:
                            t_cur = float(fp["time"][0])
                if self.last_step - 1 < len(param_paths):
                    with h5py.File(param_paths[self.last_step - 1], "r") as fp:
                        if "time" in fp:
                            t_prev = float(fp["time"][0])
            except Exception:
                pass

            if t_cur is not None and t_prev is not None:
                dt_phys = t_cur - t_prev
            else:
                dt_phys = float(cur_step - self.last_step) * 3.5156

            direction = 1.0 if cur_step > self.last_step else -1.0
            dt = abs(dt_phys) * direction
            n_sub = 10
            dt_sub = dt / n_sub

            cur_xs = prev_xs.copy()
            cur_ys = prev_ys.copy()

            for sub in range(n_sub):
                frac = (sub + 0.5) / n_sub
                bx_sub = (1.0 - frac) * sample_bilinear_vec(prev_bx, cur_xs, cur_ys) + frac * sample_bilinear_vec(bx, cur_xs, cur_ys)
                by_sub = (1.0 - frac) * sample_bilinear_vec(prev_by, cur_xs, cur_ys) + frac * sample_bilinear_vec(by, cur_xs, cur_ys)
                ez_sub = (1.0 - frac) * sample_bilinear_vec(prev_ez, cur_xs, cur_ys) + frac * sample_bilinear_vec(ez, cur_xs, cur_ys)
                b2 = np.maximum(bx_sub**2 + by_sub**2, 1e-8)
                vx = - ez_sub * by_sub / b2 * c_omp * dt_sub
                vy =   ez_sub * bx_sub / b2 * c_omp * dt_sub
                cur_xs = np.clip(cur_xs + vx, 0.0, nx - 1.0)
                cur_ys = np.clip(cur_ys + vy, 0.0, ny - 1.0)

            # Replenish left boundary markers if needed
            new_xs, new_ys = list(cur_xs), list(cur_ys)
            left_xs = [x for x in cur_xs if x < xmid]
            if len(left_xs) == 0 or min(left_xs) > nx * 0.10:
                new_xs.append(nx * 0.02)
                new_ys.append(float(ymid))

            xs = np.array(new_xs, dtype=np.float32)
            ys = np.array(new_ys, dtype=np.float32)
            self.history[cur_step] = (xs, ys)
            self.last_step = cur_step
            self.last_flds = (bx, by, ez)
            return xs, ys
        else:
            xs, ys = self.init_step(cur_step, nx, ny)
            self.last_step = cur_step
            self.last_flds = (bx, by, ez)
            return xs, ys


def _compute_frame0_base_levels(parent, n_contours, zSlice, stride):
    """Computes base A0 and delta using the first available snapshot in PathDict['Flds']."""
    num_flds = len(parent.PathDict.get("Flds", []))
    if num_flds == 0:
        return 0.0, 0.05

    first_file = parent.PathDict["Flds"][0]
    try:
        import h5py
        with h5py.File(first_file, "r") as f:
            slice_plane = parent.MainParamDict["2DSlicePlane"]
            if slice_plane == 0:
                bx_name, by_name = "bx", "by"
                sl = np.s_[min(zSlice, f[bx_name].shape[0] - 1), ::stride, ::stride]
            elif slice_plane == 1:
                bx_name, by_name = "bx", "bz"
                sl = np.s_[::stride, min(zSlice, f[bx_name].shape[1] - 1), ::stride]
            else:
                bx_name, by_name = "by", "bz"
                sl = np.s_[::stride, ::stride, min(zSlice, f[bx_name].shape[2] - 1)]

            bx0 = f[bx_name][sl]
            by0 = f[by_name][sl]

        Az0 = compute_vector_potential_2d(bx0, by0, stride=stride)
        az_min0, az_max0 = float(Az0.min()), float(Az0.max())
        width0 = max(az_max0 - az_min0, 1e-4)
        left_level = az_min0 + 0.04 * width0
        right_level = az_max0 - 0.04 * width0
        delta = (right_level - left_level) / max(1, n_contours - 1)
        A0 = left_level
        return A0, delta
    except Exception as err:
        return 0.0, 0.05


def _get_az_contour_levels(panel, Az, ny, nx, zSlice, n_contours, lagrangian, stride=1):
    """Computes contour levels for Az (Lagrangian fluid tracking or Eulerian fixed flux)."""
    parent = panel.parent
    cur_step = parent.TimeStep.value

    if not lagrangian:
        # Eulerian Fixed Flux Surfaces (Option 1)
        cache = getattr(parent, "_az_base_level_cache", None)
        if cache is None or cache.get("n_contours") != n_contours:
            A0, delta = _compute_frame0_base_levels(parent, n_contours, zSlice, stride)
            parent._az_base_level_cache = {"A0": A0, "delta": delta, "n_contours": n_contours}
        else:
            A0 = parent._az_base_level_cache["A0"]
            delta = parent._az_base_level_cache["delta"]

        az_min_s, az_max_s = float(Az.min()), float(Az.max())
        levels = find_dynamic_levels(A0, delta, az_min_s, az_max_s)
        return levels, Az

    # Lagrangian Fluid Tracking
    tracker = getattr(parent, "_lagrangian_tracker", None)
    if tracker is None or tracker.n_contours != n_contours:
        tracker = LagrangianFieldLineTracker(n_contours=n_contours)
        parent._lagrangian_tracker = tracker

    c_omp = getattr(parent, "c_omp", 1.0)
    if isinstance(c_omp, np.ndarray) and c_omp.size > 0:
        c_omp = float(c_omp.flat[0])
    elif not isinstance(c_omp, (int, float)):
        c_omp = 1.0

    slice_plane = parent.MainParamDict["2DSlicePlane"]
    if slice_plane == 0:
        bx_name, by_name, ez_name = "bx", "by", "ez"
        sl = np.s_[min(zSlice, parent.DataDict[bx_name].shape[0] - 1), ::stride, ::stride]
    elif slice_plane == 1:
        bx_name, by_name, ez_name = "bx", "bz", "ey"
        sl = np.s_[::stride, min(zSlice, parent.DataDict[bx_name].shape[1] - 1), ::stride]
    else:
        bx_name, by_name, ez_name = "by", "bz", "ex"
        sl = np.s_[::stride, ::stride, min(zSlice, parent.DataDict[bx_name].shape[2] - 1)]

    bx = parent.DataDict[bx_name][sl]
    by = parent.DataDict[by_name][sl]
    ez = parent.DataDict[ez_name][sl] if ez_name in parent.DataDict else np.zeros_like(bx)

    xs, ys = tracker.get_markers(parent, cur_step, bx, by, ez, c_omp)
    sampled_lvls = sample_bilinear_vec(Az, xs, ys)
    unique_lvls = np.unique(np.round(sampled_lvls, 5))
    if len(unique_lvls) == 0:
        unique_lvls = np.linspace(Az.min(), Az.max(), 5)
    return np.sort(unique_lvls), Az


def draw_az_contours(panel):
    """Draws magnetic vector potential Az contours on the 2D panel."""
    remove_az_contours(panel)

    stride = max(1, int(panel.GetPlotParam("az_contours_stride")))
    slice_plane = panel.parent.MainParamDict["2DSlicePlane"]

    if slice_plane == 0:  # x-y plane
        bx_name, by_name = "bx", "by"
        zSlice = panel.parent.zSlice
        slice_tuple = np.s_[zSlice, ::stride, ::stride]
    elif slice_plane == 1:  # x-z plane
        bx_name, by_name = "bx", "bz"
        zSlice = panel.parent.ySlice
        slice_tuple = np.s_[::stride, zSlice, ::stride]
    else:  # y-z plane
        bx_name, by_name = "by", "bz"
        zSlice = panel.parent.xSlice
        slice_tuple = np.s_[::stride, ::stride, zSlice]

    if bx_name not in panel.parent.DataDict or by_name not in panel.parent.DataDict:
        return

    bx = panel.parent.DataDict[bx_name][slice_tuple]
    by = panel.parent.DataDict[by_name][slice_tuple]

    if bx.ndim != 2 or by.ndim != 2:
        return

    Az = compute_vector_potential_2d(bx, by, stride=stride)

    n_contours = max(2, int(panel.GetPlotParam("az_contours_count")))
    lagrangian = bool(panel.GetPlotParam("az_contours_lagrangian"))

    levels, Az_plot = _get_az_contour_levels(panel, Az, bx.shape[0], bx.shape[1], zSlice, n_contours, lagrangian, stride=stride)

    xmin = getattr(panel, "xmin", 0.0)
    xmax = getattr(panel, "xmax", float(bx.shape[1]))
    ymin = getattr(panel, "ymin", 0.0)
    ymax = getattr(panel, "ymax", float(bx.shape[0]))

    coords_x = np.linspace(xmin, xmax, bx.shape[1])
    coords_y = np.linspace(ymin, ymax, bx.shape[0])

    color = panel.GetPlotParam("az_contours_color")
    width = float(panel.GetPlotParam("az_contours_width"))

    panel.FigWrap.az_contours = panel.FigWrap.graph.axes.contour(
        coords_x,
        coords_y,
        Az_plot,
        levels=levels,
        colors=color,
        linewidths=width,
        linestyles="solid",
    )


def refresh_az_contours(panel):
    """Refreshes Az contours by clearing and redrawing."""
    remove_az_contours(panel)
    draw_az_contours(panel)


def remove_az_contours(panel):
    """Removes Az contours from the figure."""
    if hasattr(panel.FigWrap, "az_contours") and panel.FigWrap.az_contours is not None:
        try:
            panel.FigWrap.az_contours.remove()
        except (AttributeError, TypeError):
            for coll in getattr(panel.FigWrap.az_contours, "collections", []):
                try:
                    coll.remove()
                except Exception:
                    pass
        panel.FigWrap.az_contours = None
