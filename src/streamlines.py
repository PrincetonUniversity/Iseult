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
        if panel.GetPlotParam("show_az_contours") and panel.GetPlotParam("az_contours_gauge_tracking"):
            panel.arrs_needed.append("ez")
    elif slice_plane == 1:  # x-z plane
        panel.arrs_needed.extend(["bx", "bz"])
        if panel.GetPlotParam("show_az_contours") and panel.GetPlotParam("az_contours_gauge_tracking"):
            panel.arrs_needed.append("ey")
    elif slice_plane == 2:  # y-z plane
        panel.arrs_needed.extend(["by", "bz"])
        if panel.GetPlotParam("show_az_contours") and panel.GetPlotParam("az_contours_gauge_tracking"):
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

    settings.az_contours_gauge_tracking = Tk.BooleanVar()
    settings.az_contours_gauge_tracking.set(settings.parent.GetPlotParam("az_contours_gauge_tracking"))
    Tk.ttk.Checkbutton(
        settings.frm,
        text="Track Gauge (Faraday)",
        variable=settings.az_contours_gauge_tracking,
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
    """Handle what happens when the `show_az_contours` or `az_contours_gauge_tracking` button is toggled."""
    settings.parent.SetPlotParam(
        "show_az_contours", settings.show_az_contours.get(), update_plot=False, NeedsRedraw=True
    )
    settings.parent.SetPlotParam(
        "az_contours_gauge_tracking", settings.az_contours_gauge_tracking.get(), update_plot=False, NeedsRedraw=True
    )
    if not settings.parent.GetPlotParam("show_az_contours"):
        remove_az_contours(panel)
    else:
        # Invalidate base level cache so levels recompute cleanly
        if hasattr(settings.parent.parent, "_az_base_level_cache"):
            settings.parent.parent._az_base_level_cache = None
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
        except ValueError:
            pass

    # Gauge tracking
    if hasattr(settings, "az_contours_gauge_tracking") and settings.az_contours_gauge_tracking.get() != settings.parent.plot_param_dict.get("az_contours_gauge_tracking"):
        val = bool(settings.az_contours_gauge_tracking.get())
        settings.parent.plot_param_dict["az_contours_gauge_tracking"] = val
        settings.parent.SetPlotParam("az_contours_gauge_tracking", val, update_plot=update_plot)


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


def _get_gauge_integral(panel, cur_step, ny, nx, zSlice):
    """Computes or retrieves cached inductive gauge shift G(s) = ∫ Ez(xref, ymid) * c_omp * dt."""
    parent = panel.parent
    if not hasattr(parent, "_az_gauge_cache") or parent._az_gauge_cache is None:
        parent._az_gauge_cache = {}

    num_flds = len(parent.PathDict.get("Flds", []))
    if num_flds == 0:
        return 0.0

    cache = parent._az_gauge_cache
    if "G" not in cache or len(cache["G"]) != num_flds:
        c_omp = getattr(parent, "c_omp", 1.0)
        if isinstance(c_omp, np.ndarray) and c_omp.size > 0:
            c_omp = float(c_omp.flat[0])
        elif not isinstance(c_omp, (int, float)):
            c_omp = 1.0

        ymid = ny // 2
        xref = max(0, nx - max(2, nx // 20))

        times = []
        ez_vals = []

        for idx, fpath in enumerate(parent.PathDict["Flds"]):
            try:
                import h5py
                with h5py.File(fpath, "r") as f:
                    if "ez" in f:
                        ez_arr = f["ez"]
                        if len(ez_arr.shape) == 3:
                            zs = min(zSlice, ez_arr.shape[0] - 1)
                            ym = min(ymid, ez_arr.shape[1] - 1)
                            xr = min(xref, ez_arr.shape[2] - 1)
                            ez_val = float(ez_arr[zs, ym, xr])
                        elif len(ez_arr.shape) == 2:
                            ym = min(ymid, ez_arr.shape[0] - 1)
                            xr = min(xref, ez_arr.shape[1] - 1)
                            ez_val = float(ez_arr[ym, xr])
                        else:
                            ez_val = 0.0
                    else:
                        ez_val = 0.0
                    ez_vals.append(ez_val)
                    if idx == 0 and "c_omp" in f:
                        c_omp = float(f["c_omp"][0])
            except Exception:
                ez_vals.append(0.0)

            param_paths = parent.PathDict.get("Param", [])
            t_val = None
            if idx < len(param_paths):
                try:
                    import h5py
                    with h5py.File(param_paths[idx], "r") as fp:
                        if "time" in fp:
                            t_val = float(fp["time"][0])
                except Exception:
                    pass
            if t_val is None:
                t_val = float(idx)
            times.append(t_val)

        times = np.array(times, dtype=np.float64)
        ez_vals = np.array(ez_vals, dtype=np.float64)
        dt = np.diff(times)
        G = np.zeros(num_flds, dtype=np.float64)
        for i in range(1, num_flds):
            dt_step = dt[i-1] if i-1 < len(dt) else 1.0
            G[i] = G[i-1] + 0.5 * (ez_vals[i] + ez_vals[i-1]) * c_omp * dt_step

        cache["G"] = G
        cache["xref"] = xref
        cache["ymid"] = ymid

    step_idx = max(0, min(cur_step - 1, len(cache["G"]) - 1))
    return float(cache["G"][step_idx])


def _compute_frame0_base_levels(parent, n_contours, zSlice, stride, gauge_tracking=False):
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
        ny0, nx0 = Az0.shape
        if gauge_tracking:
            ymid0 = ny0 // 2
            xref0 = max(0, nx0 - max(2, nx0 // 20))
            Az0 = Az0 - Az0[ymid0, xref0]

        az_min0, az_max0 = float(Az0.min()), float(Az0.max())
        width0 = max(az_max0 - az_min0, 1e-4)
        left_level = az_min0 + 0.04 * width0
        right_level = az_max0 - 0.04 * width0
        delta = (right_level - left_level) / max(1, n_contours - 1)
        A0 = left_level
        return A0, delta
    except Exception as err:
        return 0.0, 0.05


def _get_az_contour_levels(panel, Az, ny, nx, zSlice, n_contours, gauge_tracking):
    """Computes contour levels for Az."""
    parent = panel.parent
    cur_step = parent.TimeStep.value
    stride = max(1, int(panel.GetPlotParam("az_contours_stride")))

    if gauge_tracking:
        ymid = ny // 2
        xref = max(0, nx - max(2, nx // 20))
        ref_val = float(Az[ymid, xref])
        Az_corr = Az - ref_val
    else:
        Az_corr = Az

    # Ensure base level cache is established from Frame 0
    cache = getattr(parent, "_az_base_level_cache", None)
    if cache is None or cache.get("n_contours") != n_contours or cache.get("gauge_tracking") != gauge_tracking:
        A0, delta = _compute_frame0_base_levels(parent, n_contours, zSlice, stride, gauge_tracking=gauge_tracking)
        parent._az_base_level_cache = {"A0": A0, "delta": delta, "n_contours": n_contours, "gauge_tracking": gauge_tracking}
    else:
        A0 = parent._az_base_level_cache["A0"]
        delta = parent._az_base_level_cache["delta"]

    if gauge_tracking:
        G_val = _get_gauge_integral(panel, cur_step, ny, nx, zSlice)
    else:
        G_val = 0.0

    az_min_s, az_max_s = float(Az_corr.min()), float(Az_corr.max())
    level0 = A0 + G_val
    levels = find_dynamic_levels(level0, delta, az_min_s, az_max_s)
    return levels, Az_corr


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
    gauge_tracking = bool(panel.GetPlotParam("az_contours_gauge_tracking"))

    levels, Az_plot = _get_az_contour_levels(panel, Az, bx.shape[0], bx.shape[1], zSlice, n_contours, gauge_tracking)

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
