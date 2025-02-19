import dash
import dash_bootstrap_components as dbc

THEMES = [dbc.themes.COSMO, dbc.themes.DARKLY]
app = dash.Dash(__name__, external_stylesheets=THEMES)

app.title = "E-AB Sensor Simulation"
