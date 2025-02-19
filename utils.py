from dash import dcc, html
import dash_bootstrap_components as dbc


def slider_with_label(label_text, slider_id, min_val, max_val, step, value, 
                      marks=None, tooltip_text=None, label_style=None):
    """
    Helper to create a labeled slider with optional tooltip.
    """
    if label_style is None:
        label_style = {'fontSize': '24px'}

    return html.Div([
        html.Label(label_text, style=label_style),
        dcc.Slider(id=slider_id, 
                   min=min_val, 
                   max=max_val, 
                   step=step, 
                   value=value, 
                   marks=marks),
        dbc.Tooltip(tooltip_text, target=slider_id, placement="top") if tooltip_text else None
    ], style={'marginBottom': '15px'})

def slider_input_with_label(label_text, slider_id, input_id, min_val, max_val, step, value, 
                            marks=None, tooltip_text=None, label_style=None):
    """
    Creates a Div that contains a slider and a numeric input, which will be synchronized.
    """
    if label_style is None:
        label_style = {'fontSize': '24px'}
    return html.Div([
        html.Label(label_text, style=label_style),
        dcc.Slider(id=slider_id, 
                   min=min_val, 
                   max=max_val, 
                   step=step, 
                   value=value, 
                   marks=marks,
                   tooltip={"always_visible": True}),
        dcc.Input(id=input_id, type='number', value=value, 
                  style={'marginTop': '10px', 'width': '100%'})
    ], style={'marginBottom': '20px'})
