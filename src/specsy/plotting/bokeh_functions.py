try:
    from bokeh.plotting import figure, show
    from bokeh.models import ColumnDataSource, Whisker, HoverTool
    from bokeh.models import WheelZoomTool, PanTool, HoverTool
    bokeh_check = True
except ImportError:
    bokeh_check = False

def update_bokeh_figure(figure_obj, config_dict):

    if bokeh_check:

        # Set general figure properties
        for key, value in config_dict.items():

            # Dictionary based entries
            if isinstance(value, dict):
                match key:
                    case "xaxis":
                        for axis in figure_obj.xaxis:  # Update all x-axes
                            for attr, val in value.items():
                                setattr(axis, attr, val)

                    case "yaxis":
                        for axis in figure_obj.yaxis:  # Update all y-axes
                            for attr, val in value.items():
                                setattr(axis, attr, val)

                    case "title":
                        for attr, val in value.items():
                            setattr(figure_obj.title, attr, val)

                    case "xgrid":
                        for grid in figure_obj.xgrid:  # Update all x-grids
                            for attr, val in value.items():
                                val = None if val == 'None' else val
                                setattr(grid, attr, val)

                    case "ygrid":
                        for grid in figure_obj.ygrid:  # Update all y-grids
                            for attr, val in value.items():
                                val = None if val == 'None' else val
                                setattr(grid, attr, val)

            # Single value entries
            else:
                if key != 'tools':
                    setattr(figure_obj, key, value)

        # Set zoom and pan as active
        figure_obj.toolbar.active_scroll = figure_obj.select_one(WheelZoomTool)  # Activate zoom wheel
        figure_obj.toolbar.active_drag = figure_obj.select_one(PanTool)  # Activate pan tool

        # Legend
        figure_obj.legend.background_fill_color = config_dict['background_fill_color']
        figure_obj.legend.label_text_color = config_dict['outline_line_color']
        figure_obj.legend.background_fill_alpha = config_dict['legend']['background_fill_alpha']

        return figure_obj

    else:

        return