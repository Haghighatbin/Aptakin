# =======================================================================================
# Title: E-AB / SWV Biosensor Simulation Platform
# Description:
#     This Dash-based Python application provides a modular simulation
#     environment for analysing electrochemical aptamer-based (E-AB) sensors
#     using Square Wave Voltammetry (SWV). It includes built-in support for
#     empirical double exponential, Butler–Volmer, and Marcus theory models,
#     and visualises real and fitted electrochemical data in an interactive UI.
#
# Author: Dr Amin Haghighatbin
# Version: 3.0.0
# License: Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 (CC BY-NC-ND 4.0)
# =======================================================================================

# app.py

from dash import Dash
import dash_bootstrap_components as dbc

THEMES = [dbc.themes.COSMO, dbc.themes.DARKLY]
app = Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])

app.title = "E-AB/SWV Bioensor Simulation"
