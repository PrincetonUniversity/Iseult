#!/usr/bin/env python
"""This file contains the functions to control and generate the magnetic
streamline overlays. Due to the highly coupled nature of the classes that
control the figures these functions have many side effects. This isn't ideal
but it's better than having this code duplicated in half a dozen different
places.
"""

import tkinter as Tk
import numpy as np
import matplotlib


def add_streamline_params(param_dictionary):
    """Add data to the parameter dictionary for controlling the streamlines.

    Parameters
    ----------
    param_dictionary : dict
        The dictionary to add elements to.
    """
    param_dictionary["show_streamlines"] = False
    param_dictionary["streamlines_fields"] = "bxby"
    param_dictionary["streamlines_stride"] = 10
    param_dictionary["streamlines_density"] = 1
    param_dictionary["streamlines_color"] = "black"


def add_streamline_buttons(settings, panel, starting_row):
    """Add the various buttons, fields, etc for the streamlines to the settings pane

    Parameters
    ----------
    settings :
        The settings object for the panel being drawn in
    panel :
        The Panel object the streamlines are being drawn in
    starting_row :
        The row to start adding settings at
    """

    # Label
    Tk.ttk.Label(settings.frm, text="Streamline Settings: (2D only)").grid(
        row=starting_row, columnspan=2, sticky=Tk.W
    )

    # Create checkbox to toggle the streamlines
    settings.show_streamlines = Tk.BooleanVar()
    settings.show_streamlines.set(settings.parent.GetPlotParam("show_streamlines"))
    Tk.ttk.Checkbutton(
        settings.frm,
        text="Display Streamlines",
        variable=settings.show_streamlines,
        command=lambda: __show_streamline_handler(settings, panel),
    ).grid(row=starting_row + 1, sticky=Tk.W)

    # Choose which streamlines to display
    settings.streamlines_fields = Tk.StringVar()
    settings.streamlines_fields.set(settings.parent.GetPlotParam("streamlines_fields"))
    Tk.ttk.Radiobutton(
        settings.frm,
        text="Bx/By",
        variable=settings.streamlines_fields,
        value="bxby",
        command=lambda: __update_parameter(
            settings, "streamlines_fields", settings.streamlines_fields.get()
        ),
    ).grid(row=starting_row + 2, sticky=Tk.W)
    Tk.ttk.Radiobutton(
        settings.frm,
        text="By/Bz",
        variable=settings.streamlines_fields,
        value="bybz",
        command=lambda: __update_parameter(
            settings, "streamlines_fields", settings.streamlines_fields.get()
        ),
    ).grid(row=starting_row + 3, sticky=Tk.W)
    Tk.ttk.Radiobutton(
        settings.frm,
        text="Bx/Bz",
        variable=settings.streamlines_fields,
        value="bxbz",
        command=lambda: __update_parameter(
            settings, "streamlines_fields", settings.streamlines_fields.get()
        ),
    ).grid(row=starting_row + 4, sticky=Tk.W)

    # Set stride
    settings.streamlines_stride = Tk.IntVar(
        value=settings.parent.GetPlotParam("streamlines_stride")
    )
    Tk.ttk.Label(settings.frm, text="Stride:").grid(row=starting_row + 5, sticky=Tk.W)
    Tk.ttk.Entry(settings.frm, textvariable=settings.streamlines_stride, width=7).grid(
        row=starting_row + 5, sticky=Tk.E
    )

    # Set line density
    settings.streamlines_density = Tk.DoubleVar(
        value=settings.parent.GetPlotParam("streamlines_density")
    )
    Tk.ttk.Label(settings.frm, text="Line Density:").grid(
        row=starting_row + 6, sticky=Tk.W
    )
    Tk.ttk.Entry(settings.frm, textvariable=settings.streamlines_density, width=7).grid(
        row=starting_row + 6, sticky=Tk.E
    )

    # Set line color
    settings.streamlines_color = Tk.StringVar(
        value=settings.parent.GetPlotParam("streamlines_color")
    )
    Tk.ttk.Label(settings.frm, text="Line Color:").grid(
        row=starting_row + 7, sticky=Tk.W
    )
    Tk.ttk.Entry(settings.frm, textvariable=settings.streamlines_color, width=7).grid(
        row=starting_row + 7, sticky=Tk.E
    )


def __show_streamline_handler(settings, panel):
    """Handle what happens when the `show_streamlines` button is toggledself.

    Parameters
    ----------
    settings :
        The settings object for the panel being drawn in
    panel :
        The Panel object the streamlines are being drawn in
    """
    # Write value to the settings dictionary
    settings.parent.SetPlotParam(
        "show_streamlines", settings.show_streamlines.get(), update_plot=False
    )
    # Either create the streamlines or remove them depending on the state of the checkbox
    if settings.parent.GetPlotParam("show_streamlines"):
        draw_streamlines(panel)
    else:
        remove_streamlines(panel)

    # Update everything
    settings.parent.parent.canvas.draw()
    settings.parent.parent.canvas.get_tk_widget().update_idletasks()


def __update_parameter(settings, param_name, value):
    settings.parent.SetPlotParam(param_name, value, update_plot=True)


def streamlines_callback(settings, update_plot=True):
    # Don't update the plot if there's no streamplot active
    if not settings.parent.GetPlotParam("show_streamlines"):
        update_plot = False

    # Handle streamline stride
    if settings.streamlines_stride.get() != settings.streamlines_stride:
        settings.parent.SetPlotParam(
            "streamlines_stride",
            settings.streamlines_stride.get(),
            update_plot=update_plot,
        )

    # Handle streamline density
    if settings.streamlines_density.get() != settings.streamlines_density:
        settings.parent.SetPlotParam(
            "streamlines_density",
            settings.streamlines_density.get(),
            update_plot=update_plot,
        )

    # Handle streamline color
    if settings.streamlines_color.get() != settings.streamlines_color:
        settings.parent.SetPlotParam(
            "streamlines_color",
            settings.streamlines_color.get(),
            update_plot=update_plot,
        )


def draw_streamlines(panel):
    """Draw streamlines.

    Parameters
    ----------
    panel :
        The Panel object the streamlines are being drawn in
    """
    # Get settings from panel
    fields = panel.GetPlotParam("streamlines_fields")
    stride = panel.GetPlotParam("streamlines_stride")
    line_density = panel.GetPlotParam("streamlines_density")
    line_color = panel.GetPlotParam("streamlines_color")

    # Grab the data, use aliases to make it easier to work with
    if fields == "bxby":
        bx_name, by_name = "bx", "by"
    elif fields == "bybz":
        bx_name, by_name = "by", "bz"
    elif fields == "bxbz":
        bx_name, by_name = "bx", "bz"
    else:
        raise ValueError(
            f"Improper fields value of {fields} passed to draw_streamlines. "
        )
    bx = panel.parent.DataDict[bx_name][0, ::stride, ::stride]
    by = panel.parent.DataDict[by_name][0, ::stride, ::stride]

    # Create meshgrid
    x_min, x_max = panel.FigWrap.graph.axes.get_xlim()
    y_min, y_max = panel.FigWrap.graph.axes.get_ylim()
    coords_x = np.linspace(x_min, x_max, bx.shape[1])
    coords_y = np.linspace(y_min, y_max, bx.shape[0])
    coords_x, coords_y = np.meshgrid(coords_x, coords_y)

    # Draw plots
    panel.FigWrap.streamlines = panel.FigWrap.graph.axes.streamplot(
        coords_x, coords_y, bx, by, density=line_density, color=line_color
    )


def refresh_streamlines(panel):
    """Refresh the streamlines. Since the StreamplotSet object has no set_data() method we will just clear the streamlines and redraw them

    Parameters
    ----------
    panel :
        The panel object the streamlines are being drawn in
    """
    remove_streamlines(panel)
    draw_streamlines(panel)


def remove_streamlines(panel):
    """Remove streamlines.

    Parameters
    ----------
    panel :
        The Panel object that the streamlines are being drawn in
    """
    # Remove lines
    panel.FigWrap.streamlines.lines.remove()

    # Remove arrows
    for artist in panel.FigWrap.graph.axes.get_children():
        if isinstance(artist, matplotlib.patches.FancyArrowPatch):
            artist.remove()
