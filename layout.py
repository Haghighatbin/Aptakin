# layout.py

import dash_bootstrap_components as dbc
from dash import dcc, html, dash_table
from dash_bootstrap_templates import ThemeSwitchAIO
import data

def create_layout(app):
    """
    Constructs the layout for the Dash application including theme toggles,
    data tables, input controls, and graph components.

    Args:
        app: The Dash app instance.

    Returns:
        dbc.Container: A Bootstrap-styled container holding all UI elements.
    """
    return dbc.Container([
        ThemeSwitchAIO(
            aio_id="theme",
            themes=[dbc.themes.DARKLY, dbc.themes.BOOTSTRAP],
            icons={"left": "fa fa-sun", "right": "fa fa-moon"},
        ),

        html.H1("Electrochemical Aptamer-based Biosensors (E-AB) / SWV Simulation"),

        dash_table.DataTable(
            id='swv-table',
            columns=[
                {"name": "Frequency / (Hz)", "id": "Frequency", "type": "numeric"},
                {"name": "Unbound-state Charge / (C)", "id": "Unbound", "type": "numeric"},
                {"name": "Bound-state Charge / (C)", "id": "Bound", "type": "numeric"},
            ],
            data=data.default_data,
            editable=True,
            row_deletable=True,
            style_cell={'textAlign': 'center', 'minWidth': '120px'},
            style_header={'fontWeight': 'bold'},
            style_table={'overflowX': 'auto'},
        ),

        dbc.Button(
            "Add row",
            id="add-row-button",
            n_clicks=0,
            color="primary",
            style={"marginBottom": 20, "marginTop": 20}
        ),

        dcc.Graph(id='swv-graph'),

        dbc.Row([
            dbc.Col(
                dbc.Button("Fit Model", id="fit-button", color="primary"),
                width="auto"
            ),
            dbc.Col(
                dcc.Dropdown(
                    id='swv-model-dropdown',
                    options=[
                        {'label': 'Exponential Model', 'value': 'exponential'},
                        {'label': 'Butler-Volmer Model', 'value': 'butler-volmer'},
                        {'label': 'Marcus Theory Model', 'value': 'marcus'},
                    ],
                    value='exponential',
                    className='custom-dropdown',
                    style={
                        'borderRadius': '8px',
                        'height': '38px',
                        'width': '100%',
                    },
                ),
                width=2
            )
        ]),

        dash_table.DataTable(
            id='results-table',
            columns=[
                {"name": "Fit Parameters", "id": "Parameters", "type": "text"},
                {"name": "Unbound-state", "id": "Unbound", "type": "numeric"},
                {"name": "Bound-state", "id": "Bound", "type": "numeric"},
            ],
            data=data.results_data,
            editable=False,
            row_deletable=False,
            style_cell={'textAlign': 'center', 'minWidth': '120px'},
            style_header={'fontWeight': 'bold'},
            style_table={'overflowX': 'auto'},
        ),
    ], fluid=True)

