#!/usr/bin/env python
"""This file contains the functions to control and generate the arrow vector overlays.
"""

import tkinter as Tk
import numpy as np
import matplotlib


def add_vector_params(param_dictionary):
    """Add data to the parameter dictionary for controlling the vectors.
    """
    param_dictionary["show_vectors"] = False
    param_dictionary["vector_type"] = 0  # 0: B, 1: E, 2: J, 3: Vi, 4: Ve


def add_vector_plot_keys(panel):
    """Add components of the selected vector type to arrs_needed.
    """
    vtype = panel.GetPlotParam('vector_type')
    if vtype == 0:  # B Field
        panel.arrs_needed.extend(['bx', 'by', 'bz'])
    elif vtype == 1:  # E field
        panel.arrs_needed.extend(['ex', 'ey', 'ez'])
    elif vtype == 2:  # J [current]
        panel.arrs_needed.extend(['jx', 'jy', 'jz'])
    elif vtype == 3:  # Vi (ion vel)
        panel.arrs_needed.extend(['v3xi', 'v3yi', 'v3zi'])
    elif vtype == 4:  # Ve (electron vel)
        panel.arrs_needed.extend(['v3x', 'v3y', 'v3z', 'v3xi', 'v3yi', 'v3zi', 'dens', 'densi'])


def add_vector_buttons(settings, panel, starting_row):
    """Add the vectors checkbox and dropdown selection to the settings window next to streamlines.
    """
    settings.VectorList = ['B Field', 'E field', 'J [current]', 'Vi (ion vel)', 'Ve (electron vel)']

    settings.show_vectors = Tk.BooleanVar()
    settings.show_vectors.set(settings.parent.GetPlotParam("show_vectors"))
    Tk.ttk.Checkbutton(
        settings.frm,
        text="Display Vectors",
        variable=settings.show_vectors,
        command=lambda: __show_vector_handler(settings, panel),
    ).grid(row=starting_row + 1, column=1, sticky=Tk.W)

    settings.vector_type = Tk.StringVar()
    vtype_val = settings.parent.GetPlotParam("vector_type")
    # Clamp in case it's out of range of the list
    if vtype_val >= len(settings.VectorList):
        vtype_val = 0
    settings.vector_type.set(settings.VectorList[vtype_val])

    vector_chooser = Tk.ttk.OptionMenu(
        settings.frm,
        settings.vector_type,
        settings.VectorList[vtype_val],
        *tuple(settings.VectorList),
        command=lambda val: __change_vector_type(settings, panel)
    )
    vector_chooser.grid(row=starting_row + 1, column=2, sticky=Tk.W + Tk.E)


def __show_vector_handler(settings, panel):
    """Handle what happens when the `show_vectors` button is toggled.
    """
    settings.parent.SetPlotParam(
        "show_vectors", settings.show_vectors.get(), update_plot=False, NeedsRedraw=True
    )

    if not settings.parent.GetPlotParam("show_vectors"):
        remove_vectors(panel)
    else:
        # Reload keys first to ensure we have the correct vector data loaded
        settings.parent.parent.LoadAllKeys()
        register_zoom_callback(panel)

    settings.parent.parent.canvas.draw_idle()


def __change_vector_type(settings, panel):
    """Handle when the vector selection dropdown is changed.
    """
    selected_val = settings.vector_type.get()
    try:
        vtype_idx = settings.VectorList.index(selected_val)
    except ValueError:
        vtype_idx = 0

    settings.parent.SetPlotParam("vector_type", vtype_idx, update_plot=False, NeedsRedraw=True)

    # Reload keys since we might need new datasets (e.g. if we switched from B to Vi)
    settings.parent.parent.LoadAllKeys()

    if settings.parent.GetPlotParam("show_vectors"):
        refresh_vectors(panel)

    settings.parent.parent.canvas.draw_idle()


def remove_vectors(panel):
    """Remove quiver vector arrows if they exist.
    """
    if hasattr(panel, 'vector_quiver') and panel.vector_quiver is not None:
        try:
            panel.vector_quiver.remove()
        except Exception:
            pass
        panel.vector_quiver = None


def register_zoom_callback(panel):
    """Connect zoom/pan events on axis limits to redrawing vectors.
    """
    if not hasattr(panel, 'vector_cid_x') or panel.vector_cid_x is None:
        panel.vector_cid_x = panel.axes.callbacks.connect('xlim_changed', lambda ax: on_limits_changed(panel))
    if not hasattr(panel, 'vector_cid_y') or panel.vector_cid_y is None:
        panel.vector_cid_y = panel.axes.callbacks.connect('ylim_changed', lambda ax: on_limits_changed(panel))


def on_limits_changed(panel):
    """Triggered on axes zoom/pan. Recalculates and updates the quiver grid for this panel.
    """
    if panel.GetPlotParam("show_vectors"):
        refresh_vectors(panel)
        panel.parent.canvas.draw_idle()


def draw_vectors(panel):
    """Draw vectors on the panel's axes.
    """
    remove_vectors(panel)

    slice_plane = panel.parent.MainParamDict["2DSlicePlane"]
    vtype = panel.GetPlotParam('vector_type')

    try:
        if vtype == 0:  # B Field
            if slice_plane == 0:
                U_full = panel.FigWrap.LoadKey('bx')[panel.parent.zSlice, :, :]
                V_full = panel.FigWrap.LoadKey('by')[panel.parent.zSlice, :, :]
            elif slice_plane == 1:
                U_full = panel.FigWrap.LoadKey('bx')[:, panel.parent.ySlice, :]
                V_full = panel.FigWrap.LoadKey('bz')[:, panel.parent.ySlice, :]
            elif slice_plane == 2:
                U_full = panel.FigWrap.LoadKey('by')[:, :, panel.parent.xSlice]
                V_full = panel.FigWrap.LoadKey('bz')[:, :, panel.parent.xSlice]
        elif vtype == 1:  # E Field
            if slice_plane == 0:
                U_full = panel.FigWrap.LoadKey('ex')[panel.parent.zSlice, :, :]
                V_full = panel.FigWrap.LoadKey('ey')[panel.parent.zSlice, :, :]
            elif slice_plane == 1:
                U_full = panel.FigWrap.LoadKey('ex')[:, panel.parent.ySlice, :]
                V_full = panel.FigWrap.LoadKey('ez')[:, panel.parent.ySlice, :]
            elif slice_plane == 2:
                U_full = panel.FigWrap.LoadKey('ey')[:, :, panel.parent.xSlice]
                V_full = panel.FigWrap.LoadKey('ez')[:, :, panel.parent.xSlice]
        elif vtype == 2:  # J Field
            if slice_plane == 0:
                U_full = panel.FigWrap.LoadKey('jx')[panel.parent.zSlice, :, :]
                V_full = panel.FigWrap.LoadKey('jy')[panel.parent.zSlice, :, :]
            elif slice_plane == 1:
                U_full = panel.FigWrap.LoadKey('jx')[:, panel.parent.ySlice, :]
                V_full = panel.FigWrap.LoadKey('jz')[:, panel.parent.ySlice, :]
            elif slice_plane == 2:
                U_full = panel.FigWrap.LoadKey('jy')[:, :, panel.parent.xSlice]
                V_full = panel.FigWrap.LoadKey('jz')[:, :, panel.parent.xSlice]
        elif vtype == 3:  # Vi Field
            if slice_plane == 0:
                U_full = panel.FigWrap.LoadKey('v3xi')[panel.parent.zSlice, :, :]
                V_full = panel.FigWrap.LoadKey('v3yi')[panel.parent.zSlice, :, :]
            elif slice_plane == 1:
                U_full = panel.FigWrap.LoadKey('v3xi')[:, panel.parent.ySlice, :]
                V_full = panel.FigWrap.LoadKey('v3zi')[:, panel.parent.ySlice, :]
            elif slice_plane == 2:
                U_full = panel.FigWrap.LoadKey('v3yi')[:, :, panel.parent.xSlice]
                V_full = panel.FigWrap.LoadKey('v3zi')[:, :, panel.parent.xSlice]
        elif vtype == 4:  # Ve Field
            dens = panel.FigWrap.LoadKey('dens')
            densi = panel.FigWrap.LoadKey('densi')
            dense = np.maximum(dens - densi, 1e-5)

            if slice_plane == 0:
                v3x = panel.FigWrap.LoadKey('v3x')[panel.parent.zSlice, :, :]
                v3xi = panel.FigWrap.LoadKey('v3xi')[panel.parent.zSlice, :, :]
                v3y = panel.FigWrap.LoadKey('v3y')[panel.parent.zSlice, :, :]
                v3yi = panel.FigWrap.LoadKey('v3yi')[panel.parent.zSlice, :, :]

                U_full = (dens[panel.parent.zSlice, :, :] * v3x - densi[panel.parent.zSlice, :, :] * v3xi) / dense[panel.parent.zSlice, :, :]
                V_full = (dens[panel.parent.zSlice, :, :] * v3y - densi[panel.parent.zSlice, :, :] * v3yi) / dense[panel.parent.zSlice, :, :]
            elif slice_plane == 1:
                v3x = panel.FigWrap.LoadKey('v3x')[:, panel.parent.ySlice, :]
                v3xi = panel.FigWrap.LoadKey('v3xi')[:, panel.parent.ySlice, :]
                v3z = panel.FigWrap.LoadKey('v3z')[:, panel.parent.ySlice, :]
                v3zi = panel.FigWrap.LoadKey('v3zi')[:, panel.parent.ySlice, :]

                U_full = (dens[:, panel.parent.ySlice, :] * v3x - densi[:, panel.parent.ySlice, :] * v3xi) / dense[:, panel.parent.ySlice, :]
                V_full = (dens[:, panel.parent.ySlice, :] * v3z - densi[:, panel.parent.ySlice, :] * v3zi) / dense[:, panel.parent.ySlice, :]
            elif slice_plane == 2:
                v3y = panel.FigWrap.LoadKey('v3y')[:, :, panel.parent.xSlice]
                v3yi = panel.FigWrap.LoadKey('v3yi')[:, :, panel.parent.xSlice]
                v3z = panel.FigWrap.LoadKey('v3z')[:, :, panel.parent.xSlice]
                v3zi = panel.FigWrap.LoadKey('v3zi')[:, :, panel.parent.xSlice]

                U_full = (dens[:, :, panel.parent.xSlice] * v3y - densi[:, :, panel.parent.xSlice] * v3yi) / dense[:, :, panel.parent.xSlice]
                V_full = (dens[:, :, panel.parent.xSlice] * v3z - densi[:, :, panel.parent.xSlice] * v3zi) / dense[:, :, panel.parent.xSlice]
        else:
            return
    except (AttributeError, KeyError, TypeError):
        return

    register_zoom_callback(panel)

    xlim = panel.axes.get_xlim()
    ylim = panel.axes.get_ylim()

    c_omp = panel.c_omp
    istep = panel.istep

    Ny_grid, Nx_grid = U_full.shape

    i_min = int(np.floor(xlim[0] * c_omp / istep))
    i_max = int(np.ceil(xlim[1] * c_omp / istep))
    j_min = int(np.floor(ylim[0] * c_omp / istep))
    j_max = int(np.ceil(ylim[1] * c_omp / istep))

    # Clamp indices
    i_min = max(0, min(i_min, Nx_grid - 1))
    i_max = max(0, min(i_max, Nx_grid - 1))
    j_min = max(0, min(j_min, Ny_grid - 1))
    j_max = max(0, min(j_max, Ny_grid - 1))

    if i_max <= i_min:
        i_max = i_min + 1
    if j_max <= j_min:
        j_max = j_min + 1

    # Target 50 vector points in each direction
    stride_x = max(1, (i_max - i_min) // 50)
    stride_y = max(1, (j_max - j_min) // 50)

    ix = np.arange(i_min, i_max + 1, stride_x)
    iy = np.arange(j_min, j_max + 1, stride_y)

    ix = ix[ix < Nx_grid]
    iy = iy[iy < Ny_grid]

    if len(ix) == 0 or len(iy) == 0:
        return

    IX, IY = np.meshgrid(ix, iy)

    X = IX * (istep / c_omp)
    Y = IY * (istep / c_omp)

    U = U_full[IY, IX]
    V = V_full[IY, IX]

    # Temporarily turn off autoscale so quiver doesn't alter axes limits
    autoscale_on = panel.axes.get_autoscale_on()
    panel.axes.set_autoscale_on(False)
    try:
        panel.vector_quiver = panel.axes.quiver(X, Y, U, V, pivot='middle', color='black')
    finally:
        panel.axes.set_autoscale_on(autoscale_on)


def refresh_vectors(panel):
    """Refresh the vector overlay by removing the old quiver and drawing a new one.
    """
    remove_vectors(panel)
    draw_vectors(panel)
