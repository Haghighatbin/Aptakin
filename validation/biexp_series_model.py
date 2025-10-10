import dash
from dash import Dash, dcc, html, dash_table, Input, Output, Patch, clientside_callback, callback
import plotly.express as px
import plotly.io as pio
from dash_bootstrap_templates import load_figure_template, ThemeSwitchAIO
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
import numpy as np
from scipy.optimize import curve_fit
from config import Config

# ---------------------------------------------------------------------------------
# 1. Define the Dash app
# ---------------------------------------------------------------------------------
app = Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])
server = app.server

# ---------------------------------------------------------------------------------
# 2. Helper function: Exponential SWV model
# ---------------------------------------------------------------------------------
def double_exponential_swv_model_series(f, I0, k_et, f_decay1, f_decay2):
    """
    Double-exponential SWV model (series/multiplicative approach):

    I_net(f) = I0 * exp(-k_et/(2f)) * exp(-f/f_decay1) * exp(-f/f_decay2)

    Here:
      - I0 is a scaling factor.
      - exp(-k_et/(2f)) is the electron-transfer gating factor often used in simple SWV models.
      - exp(-f/f_decay1) might represent, e.g., aptamer conformational gating (or one timescale).
      - exp(-f/f_decay2) might lump together other processes (double-layer, mass transport, adsorption).
    
    All these exponentials multiply together, implying they act in series on a single pathway.
    """
    return (
        I0
        * np.exp(-k_et / (2.0 * f))
        * np.exp(-f / f_decay1)
        * np.exp(-f / f_decay2)
    )

# ---------------------------------------------------------------------------------
# 3. Default Table Data
# ---------------------------------------------------------------------------------
default_data = [
    {"Frequency": 10,  "Unbound": 3.77e-07, "Bound": 5.14e-07},
    {"Frequency": 16,  "Unbound": 4.66e-07, "Bound": 5.84e-07},
    {"Frequency": 20,  "Unbound": 4.66e-07, "Bound": 6.14e-07},
    {"Frequency": 25,  "Unbound": 4.48e-07, "Bound": 6.58e-07},
    {"Frequency": 40,  "Unbound": 3.77e-07, "Bound": 8.35e-07},
    {"Frequency": 60,  "Unbound": 3.03e-07, "Bound": 8.16e-07},
    {"Frequency": 100, "Unbound": 2.13e-07, "Bound": 5.94e-07},
    {"Frequency": 150, "Unbound": 1.58e-07, "Bound": 3.63e-07},
    {"Frequency": 180, "Unbound": 1.44e-07, "Bound": 3.25e-07},
    {"Frequency": 210, "Unbound": 1.09e-07, "Bound": 2.95e-07},
    {"Frequency": 240, "Unbound": 1.05e-07, "Bound": 2.75e-07},
]
results_data = [
    {"Parameters": "I0",         "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "k_et",       "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "f_decay1",   "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "f_decay2",   "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "SSR",        "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "R²",         "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "Reduced-X²", "Unbound": 0.0, "Bound": 0.0},
    ]
# ---------------------------------------------------------------------------------
# 4. Layout: Theme switcher + DataTable + "Fit" button + graph + parameter output
# ---------------------------------------------------------------------------------
app.layout = dbc.Container([
    ThemeSwitchAIO(
        aio_id="theme", 
        themes=[dbc.themes.DARKLY, dbc.themes.BOOTSTRAP],
        icons={"left": "fa fa-sun", "right": "fa fa-moon"},
    ),
    
    html.H1("SWV Exponential Model Fitting", style={"marginTop": 20}),
    
    dash_table.DataTable(
        id='swv-table',
        columns=[
            {"name": "Frequency/(Hz)", "id": "Frequency", "type": "numeric"},
            {"name": "Unbound-state Charge/(C)",     "id": "Unbound",   "type": "numeric"},
            {"name": "Bound-state Charge/(C)",      "id": "Bound",     "type": "numeric"},
        ],
        data=default_data,
        editable=True,
        row_deletable=True,
        style_cell={'textAlign': 'center', 'minWidth': '120px'},
        style_header={'fontWeight': 'bold'},
        style_table={'overflowX': 'auto'},
    ),

    dbc.Button("Add row", id="add-row-button", n_clicks=0, color="primary", style={"marginButtom": 20, "marginTop": 20}),    

    html.Br(),
    
    dcc.Graph(id='swv-graph', figure={}),
    
    dbc.Button("Fit Model", id="fit-button", color="primary", style={"marginBottom": 20, "marginTop": 20}),

    dash_table.DataTable(
        id='results-table',
        columns=[
            {"name": "Fit Parameters/Goodness", "id": "Parameters", "type": "text"},
            {"name": "Unbound-state",     "id": "Unbound",   "type": "numeric"},
            {"name": "Bound-state",      "id": "Bound",     "type": "numeric"},
        ],
        data=results_data,
        editable=True,
        row_deletable=False,
        style_cell={'textAlign': 'center', 'minWidth': '120px'},
        style_header={'fontWeight': 'bold'},
        style_table={'overflowX': 'auto'},
    ),
], fluid=True)

# ---------------------------------------------------------------------------------
# Callback to update the DataTable styling based on the global theme
# ---------------------------------------------------------------------------------
@app.callback(
    Output('swv-table', 'data'),
    Input('add-row-button', 'n_clicks'),
    State('swv-table', 'data'),
    prevent_initial_call=True
)
def add_row(n_clicks, rows):
    if rows is None:
        rows = []
    rows.append({"Frequency": 0, "Unbound": 0, "Bound": 0})
    return rows

@app.callback(
    [Output('results-table', 'style_cell'),
     Output('results-table', 'style_header')],
    Input(ThemeSwitchAIO.ids.switch("theme"), "value")
)
def update_table_styles(current_theme):
    if current_theme:
        style_cell = {
            'textAlign': 'center',
            'minWidth': '120px',
            'backgroundColor': '#303030',
            'color': 'white'
        }
        style_header = {
            'fontWeight': 'bold',
            'backgroundColor': '#1a1a1a',
            'color': 'white'
        }
    else:
        style_cell = {
            'textAlign': 'center',
            'minWidth': '120px',
            'backgroundColor': 'white',
            'color': 'black'
        }
        style_header = {
            'fontWeight': 'bold',
            'backgroundColor': 'lightgrey',
            'color': 'black'
        }
    return style_cell, style_header

@app.callback(
    [Output('swv-table', 'style_cell'),
     Output('swv-table', 'style_header')],
    Input(ThemeSwitchAIO.ids.switch("theme"), "value")
)
def update_table_styles(current_theme):
    if current_theme:
        style_cell = {
            'textAlign': 'center',
            'minWidth': '120px',
            'backgroundColor': '#303030',
            'color': 'white'
        }
        style_header = {
            'fontWeight': 'bold',
            'backgroundColor': '#1a1a1a',
            'color': 'white'
        }
    else:
        style_cell = {
            'textAlign': 'center',
            'minWidth': '120px',
            'backgroundColor': 'white',
            'color': 'black'
        }
        style_header = {
            'fontWeight': 'bold',
            'backgroundColor': 'lightgrey',
            'color': 'black'
        }
    return style_cell, style_header

# ---------------------------------------------------------------------------------
# 5. Callback: "Fit" button and theme switcher -> 
# perform separate fits for unbound & bound -> update graph & text
# ---------------------------------------------------------------------------------
@app.callback(
    [Output("swv-graph", "figure"),
     Output("results-table", "data")],
    [Input("fit-button", "n_clicks"),
     Input(ThemeSwitchAIO.ids.switch("theme"), "value")],
    State("swv-table", "data"),
)
def perform_fit(n_clicks, current_theme, table_data):
    template = "plotly_white" if not current_theme else "plotly_dark"

    if not n_clicks:
        fig = go.Figure()
        fig.update_layout(template=template)

        return fig, []
    
    # ---------------------------------------------------------------------------------
    # A) Extract data from the table
    # ---------------------------------------------------------------------------------
    freqs = []
    unbound_vals = []
    bound_vals = []
    for row in table_data:
        try:
            f = float(row["Frequency"])
            ub = float(row["Unbound"])
            bd = float(row["Bound"])
            # Filter out any incomplete rows
            if f > 0 and ub >= 0 and bd >= 0:
                freqs.append(f)
                unbound_vals.append(ub)
                bound_vals.append(bd)
        except:
            pass
    
    freqs = np.array(freqs)
    unbound_vals = np.array(unbound_vals)
    bound_vals = np.array(bound_vals)
    
    if len(freqs) < 2:
        # Need at least two points to do a basic fit
        fig = go.Figure()
        fig.update_layout(template=template)
        return fig, "Not enough data points to perform fit."
    
    # ---------------------------------------------------------------------------------
    # B) Fit the unbound data
    # ---------------------------------------------------------------------------------
    p0_unbound = [
        np.max(unbound_vals) * Config.CF_UNBOUND_I0_MUL,    # I0
        Config.UNBOUND_K_ET,                                # k_et
        Config.UNBOUND_F_DECAY1,                            # f_decay1
        Config.UNBOUND_F_DECAY2,                            # f_decay2
    ]
    # Double-exponential
    try:
        bounds_u = (Config.CF_UNBOUND_L_BNDS, Config.CF_UNBOUND_H_BNDS)
        sigma_vals = Config.CF_UNBOUND_PCT_SIGMA * unbound_vals
        method = Config.CF_UNBOUND_METHOD
        maxfev = Config.CF_UNBOUND_MAXFEV
        popt_unbound, _ = curve_fit(double_exponential_swv_model_series,
                                                freqs,
                                                unbound_vals,
                                                p0=p0_unbound,
                                                bounds=bounds_u,
                                                sigma=sigma_vals, 
                                                method=method,
                                                maxfev=maxfev
                                                )
        I0_u, k_et_u, fd1_u, fd2_u = popt_unbound

        print("Optimised parameters (UnBound):")
        param_names = ["I0", "k_et", "f_decay1", "f_decay2"]
        for name, val in zip(param_names, popt_unbound):
            print(f"  {name} = {val:.6g}")

        residuals = unbound_vals - double_exponential_swv_model_series(freqs, *popt_unbound)    # 1) Residuals
        unbound_ss_res = np.sum(residuals**2)                                                   # 2) Sum of squared residuals
        unbound_ss_tot = np.sum((unbound_vals - np.mean(unbound_vals))**2)                      # 3) Total sum of squares (proportional to data variance)
        unbound_r_squared = 1 - (unbound_ss_res / unbound_ss_tot)                               # 4) R^2
        sigma_array = Config.CF_UNBOUND_PCT_SIGMA * unbound_vals                                # 5% Relative Error

        # sigma_array = np.full_like(bound_vals, Config.CF_BOUND_RAW_VAL_SIGMA)
        dof = len(freqs) - len(popt_unbound)                                                    # Degrees of freedom
        unbound_chi2 = np.sum((residuals / sigma_array)**2)                                     # 5% for each point instead of 5% overall
        unbound_reduced_chi2 = unbound_chi2 / dof

        print(f"\nSum of squared residuals (SSR) = {unbound_ss_res:.6g}")
        print(f"R² = {unbound_r_squared:.6g}")
        print(f"Reduced chi-square = {unbound_reduced_chi2:.6g}")

    except Exception as e:
        print(e)
        I0_u, k_et_u, fd1_u, fd2_u = (None,) * 4    
    
    # ---------------------------------------------------------------------------------
    # C) Fit the bound data
    # ---------------------------------------------------------------------------------
    p0_bound = [
        np.max(bound_vals) * Config.CF_BOUND_I0_MUL,    # I0
        Config.BOUND_K_ET,                              # k_et
        Config.BOUND_F_DECAY1,                          # f_decay1
        Config.BOUND_F_DECAY2,                          # f_decay2
    ]
    try:
        bounds_b = (Config.CF_BOUND_L_BNDS, Config.CF_BOUND_H_BNDS)
        sigma_vals = Config.CF_BOUND_PCT_SIGMA * bound_vals
        method = Config.CF_BOUND_METHOD
        maxfev = Config.CF_BOUND_MAXFEV
        popt_bound, _ = curve_fit(double_exponential_swv_model_series,
                                                freqs,
                                                bound_vals,
                                                p0=p0_bound,
                                                bounds=bounds_b,
                                                sigma=sigma_vals, 
                                                method=method,
                                                maxfev=maxfev
                                                )
        I0_b, k_et_b, fd1_b, fd2_b = popt_bound

        print("Optimised parameters (Bound):")
        param_names = ["I0", "k_et", "f_decay1", "f_decay2"]
        for name, val in zip(param_names, popt_bound):
            print(f"  {name} = {val:.6g}")

        residuals = bound_vals - double_exponential_swv_model_series(freqs, *popt_bound)  # 1) Residuals
        bound_ss_res = np.sum(residuals**2)                                               # 2) Sum of squared residuals
        bound_ss_tot = np.sum((bound_vals - np.mean(bound_vals))**2)                      # 3) Total sum of squares (proportional to data variance)
        bound_r_squared = 1 - (bound_ss_res / bound_ss_tot)                               # 4) R^2
        sigma_array = Config.CF_BOUND_PCT_SIGMA * bound_vals                              # 5% Relative Error

        # sigma_array = np.full_like(bound_vals, Config.CF_BOUND_RAW_VAL_SIGMA)
        dof = len(freqs) - len(popt_bound)                                                # Degrees of freedom
        bound_chi2 = np.sum((residuals / sigma_array)**2)                                 # 5% for each point instead of 5% overall
        bound_reduced_chi2 = bound_chi2 / dof

        print(f"\nSum of squared residuals (SSR) = {bound_ss_res:.6g}")
        print(f"R² = {bound_r_squared:.6g}")
        print(f"Reduced chi-square = {bound_reduced_chi2:.6g}")

    except Exception as e:
        print(e)
        I0_b, k_et_b, fd1_b, fd2_b = (None,) * 4
    
    # ---------------------------------------------------------------------------------
    # D) Prepare a frequency range for plotting the fitted curves
    # ---------------------------------------------------------------------------------
    f_min, f_max = np.min(freqs), np.max(freqs)
    f_plot = np.linspace(f_min, f_max, Config.PLT_LNSPC_NUM)
    step_size = (f_max - f_min)/(Config.PLT_LNSPC_NUM - 1)
    
    # unbound fit
    if I0_u is not None:
        unbound_fit = double_exponential_swv_model_series(f_plot, I0_u, k_et_u, fd1_u, fd2_u)
    else:
        unbound_fit = None
    
    # bound fit
    if I0_b is not None:
        bound_fit = double_exponential_swv_model_series(f_plot, I0_b, k_et_b, fd1_b, fd2_b)
    else:
        bound_fit = None

    # Extrapolated
    num_new_points = int(round((Config.PLT_EXTP_MAXFREQ - Config.PLT_EXTP_MINFREQ) / step_size)) + 1

    f_extp = np.linspace(Config.PLT_EXTP_MINFREQ, Config.PLT_EXTP_MAXFREQ, num_new_points)
    if I0_u is not None:
        unbound_fit_extp = double_exponential_swv_model_series(f_extp, I0_u, k_et_u, fd1_u, fd2_u)
    else:
        unbound_fit_extp = None
    
    # bound fit
    if I0_b is not None:
        bound_fit_extp = double_exponential_swv_model_series(f_extp, I0_b, k_et_b, fd1_b, fd2_b)
    else:
        bound_fit_extp = None

    # ---------------------------------------------------------------------------------
    # E) Normalise data & fits for plotting: "Normalised charge" vs frequency
    # ---------------------------------------------------------------------------------
    unbound_norm_data = unbound_vals / (np.max(unbound_vals) + Config.PLT_NODEVZERO_VAL)
    bound_norm_data   = bound_vals   / (np.max(bound_vals)   + Config.PLT_NODEVZERO_VAL)
    
    unbound_norm_fit = None
    bound_norm_fit   = None

    if unbound_fit is not None:
        unbound_norm_fit = unbound_fit / (np.max(unbound_fit) + Config.PLT_NODEVZERO_VAL)
    if bound_fit is not None:
        bound_norm_fit = bound_fit / (np.max(bound_fit) + Config.PLT_NODEVZERO_VAL)

    # EXTP
    unbound_norm_fit_extp = None
    bound_norm_fit_extp   = None

    if unbound_fit_extp is not None:
        unbound_norm_fit_extp = unbound_fit_extp / (np.max(unbound_fit_extp) + Config.PLT_NODEVZERO_VAL)
    if bound_fit_extp is not None:
        bound_norm_fit_extp = bound_fit_extp / (np.max(bound_fit_extp) + Config.PLT_NODEVZERO_VAL)

    # ---------------------------------------------------------------------------------
    # F) Build the figure
    # ---------------------------------------------------------------------------------
    fig = go.Figure()
    # Extrapolated fit
    if unbound_norm_fit_extp is not None:
        fig.add_trace(go.Scatter(
            x=f_extp, y=unbound_norm_fit_extp, mode='lines',
            name='EXTP Unbound (Fit)', line=dict(color='cyan', dash='dash')
        ))
    if bound_norm_fit_extp is not None:
        fig.add_trace(go.Scatter(
            x=f_extp, y=bound_norm_fit_extp, mode='lines',
            name='EXTP Bound (Fit)', line=dict(color='orange', dash='dash')
        ))

    fig.add_trace(go.Scatter(
        x=freqs, y=unbound_norm_data, mode='markers',
        name='Unbound (Data)', marker=dict(color='blue', symbol='circle')
    ))
    fig.add_trace(go.Scatter(
        x=freqs, y=bound_norm_data, mode='markers',
        name='Bound (Data)', marker=dict(color='red', symbol='square')
    ))
    
    # Plot fitted curves
    if unbound_norm_fit is not None:
        fig.add_trace(go.Scatter(
            x=f_plot, y=unbound_norm_fit, mode='lines',
            name='Unbound (Fit)', line=dict(color='blue')
        ))
    if bound_norm_fit is not None:
        fig.add_trace(go.Scatter(
            x=f_plot, y=bound_norm_fit, mode='lines',
            name='Bound (Fit)', line=dict(color='red')
        ))

    fig.update_layout(
        title="Double-Exponential Model Fit (Normalised Charge vs Frequency)",
        xaxis_title="Frequency / (Hz)",
        yaxis_title="Normalised Charge / (C)",
        template=template
    )
    
    # ---------------------------------------------------------------------------------
    # G) Prepare output text
    # ---------------------------------------------------------------------------------
    if any(val is None for val in (I0_u, k_et_u, fd1_u, fd2_u)):
        print('Failed to fit the Unbound-State!')
    
    if any(val is None for val in (I0_b, k_et_b, fd1_b, fd2_b)):
        print('Failed to fit the Bound-State!')
    
    results_table_data = [
    {"Parameters": "I0",         "Unbound": f'{I0_u:.10f}',                 "Bound": f'{I0_b:.10f}'},
    {"Parameters": "k_et",       "Unbound": f'{k_et_u:.3f}',                "Bound": f'{k_et_b:.3f}'},
    {"Parameters": "f_decay1",   "Unbound": f'{fd1_u:.3f}',                 "Bound": f'{fd1_b:.3f}'},
    {"Parameters": "f_decay2",   "Unbound": f'{fd2_u:.3f}',                 "Bound": f'{fd2_b:.3f}'},
    {"Parameters": "SSR",        "Unbound": f'{unbound_ss_res:.6g}',        "Bound": f'{bound_ss_res:.6g}'},
    {"Parameters": "R²",         "Unbound": f'{unbound_r_squared:.6g}',     "Bound": f'{bound_r_squared:.6g}'},
    {"Parameters": "Reduced-X²", "Unbound": f'{unbound_reduced_chi2:.6g}',  "Bound": f'{bound_reduced_chi2:.6g}'},
    ]
    return fig, results_table_data

# ---------------------------------------------------------------------------------
# 6. Run the App
# ---------------------------------------------------------------------------------
if __name__ == "__main__":
    app.run_server(debug=True)
