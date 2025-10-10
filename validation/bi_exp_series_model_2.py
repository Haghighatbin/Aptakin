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

# Suppose you store various config parameters in a separate file
from config import Config

app = Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])
server = app.server

################################################################################
# 1. MODEL FUNCTIONS
################################################################################

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

def butler_volmer_swv_model_series(f, I0, j0, alpha, eta0, k_freq, f_decay1, f_decay2):
    """
    Butler-Volmer-based SWV model (series).
    Example approach:
      I_net(f) = I0 * j_BV(eta_eff(f)) * exp(-f/f_decay1) * exp(-f/f_decay2)

    Where:
      j_BV(eta) = j0 * [exp((1-alpha) * F * eta / (R*T)) - exp(-alpha * F * eta / (R*T))]
      eta_eff(f) is a simple phenomenological guess, e.g. eta0 * exp(-k_et/(2f)) 
        or just treat 'eta0' as a constant if you want a fixed overpotential.

    For simplicity here, we'll do:
      eta_eff = eta0 * np.exp(-f/f_decay1 * 0.0)  # or something trivial
    REV1:
        Butler-Volmer-based SWV model (series) with a frequency-dependent overpotential:
      eta_eff(f) = eta0 * [1 - exp(-k_freq / f)]
    so that at low freq, overpotential ~ 0, at higher freq, it approaches eta0.

    I_net(f) = I0 * j_BV( eta_eff(f) ) * exp(-f/f_decay1) * exp(-f/f_decay2)
    """
    # In reality, you might define a frequency-dependent overpotential.
    # For demonstration, let's just use a fixed eta0 for all frequencies:
    # j_BV(eta0) = j0 * [exp((1-alpha)*F*eta0/(R*T)) - exp(-alpha*F*eta0/(R*T))]
    # We'll read physical constants from Config, or treat them as global.
    F  = Config.FARADAY_CONST
    R  = Config.UNI_GAS_CONST
    T  = Config.STRD_TEMP

    eta_eff = eta0 * (1.0 - np.exp(-k_freq / f))

    term_forward = np.exp((1 - alpha) * F * eta_eff / (R * T))
    term_reverse = np.exp(-alpha * F * eta_eff / (R * T))
    j_bv = j0 * (term_forward - term_reverse)

    # Then multiply by the "penalty" exponentials for frequency:
    return (
        I0
        * j_bv
        * np.exp(-f / f_decay1)
        * np.exp(-f / f_decay2)
    )

def marcus_swv_model_series(f, I0, A_marcus, lambda_marcus, deltaG0, k_freq, f_decay1, f_decay2):
    """
    Marcus-based SWV model (series).
    Example approach:
      I_net(f) = I0 * k_et_marcus(...) * exp(-f/f_decay1) * exp(-f/f_decay2)

    Where:
      k_et_marcus = A_marcus * exp( -((lambda + deltaG0)^2 / (4 * lambda * kB_eV * T)) )

    We'll read kB_eV, T from config. This is still simplified, ignoring frequency-dependence of ET.
    """
    # Compute k_et from Marcus
    k_B_eV = Config.BOLTZMANN_EVK_CONST  # eV/K
    T      = Config.STRD_TEMP
    denom  = 4.0 * lambda_marcus * k_B_eV * T
    numerator = (lambda_marcus + deltaG0)**2
    exponent  = -(numerator / denom)
    # k_et = A_marcus * np.exp(exponent)  # a single number

    # REV for exponetial
    k_et_base = A_marcus * np.exp(exponent)  # e.g. at f -> infinity ?

    # Frequency-dependent factor
    # (At low freq => near zero => rate is small, at mid freq => it saturates)
    # factor_freq = 1.0 - np.exp(-k_freq / f)
    factor_freq = 1.0 - np.exp(-f / k_freq)

    k_et = k_et_base * factor_freq
    # Multiply by exponentials for gating/other processes
    return (
        I0
        * k_et
        * np.exp(-f / f_decay1)
        * np.exp(-f / f_decay2)
    )

################################################################################
# 2. Helper: A single function to evaluate the chosen model
################################################################################

def swv_model_eval(f, model_type, params):
    """
    Evaluate the SWV current at frequencies 'f' for the chosen model type.
    'params' is a tuple of parameters in the correct order for each model.

    model_type in ['exponential', 'butler-volmer', 'marcus']
    """
    if model_type == 'exponential':
        # Expect params = (I0, k_et, f_decay1, f_decay2)
        return double_exponential_swv_model_series(f, *params)
    elif model_type == 'butler-volmer':
        # Expect params = (I0, j0, alpha, eta0, f_decay1, f_decay2)
        return butler_volmer_swv_model_series(f, *params)
    elif model_type == 'marcus':
        # Expect params = (I0, A_marcus, lambda_marcus, deltaG0, f_decay1, f_decay2)
        return marcus_swv_model_series(f, *params)
    else:
        # Fallback
        return np.zeros_like(f)

################################################################################
# 3. Layout
################################################################################

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
    {"Parameters": "Param1", "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "Param2", "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "Param3", "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "SSR",    "Unbound": 0.0, "Bound": 0.0},
    {"Parameters": "R²",     "Unbound": 0.0, "Bound": 0.0},
]

app.layout = dbc.Container([
    ThemeSwitchAIO(
        aio_id="theme", 
        themes=[dbc.themes.DARKLY, dbc.themes.BOOTSTRAP],
        icons={"left": "fa fa-sun", "right": "fa fa-moon"},
    ),

    html.H1("SWV Model Fitting with Multiple Options"),

    dash_table.DataTable(
        id='swv-table',
        columns=[
            {"name": "Frequency/(Hz)", "id": "Frequency", "type": "numeric"},
            {"name": "Unbound-state Charge/(C)",    "id": "Unbound",   "type": "numeric"},
            {"name": "Bound-state Charge/(C)",     "id": "Bound",     "type": "numeric"},
        ],
        data=default_data,
        editable=True,
        row_deletable=True,
        style_cell={'textAlign': 'center', 'minWidth': '120px'},
        style_header={'fontWeight': 'bold'},
        style_table={'overflowX': 'auto'},
    ),
    dbc.Button("Add row", id="add-row-button", n_clicks=0, color="primary", style={"marginButtom": 20, "marginTop": 20}),    

    dcc.Graph(id='swv-graph'),
    dbc.Row(
    [
        dbc.Col(dbc.Button("Fit Model", id="fit-button", color="primary"), width="auto"),
        dbc.Col(
            dcc.Dropdown(
                id='swv-model-dropdown',
                options=[
                    {'label': 'Exponential Model', 'value': 'exponential'},
                    {'label': 'Butler-Volmer Model', 'value': 'butler-volmer'},
                    {'label': 'Marcus Theory Model', 'value': 'marcus'},
                ],
                value='exponential',
                style={
                'backgroundColor': 'green',
                'color': 'orange',
                'borderRadius': '8px',
                'height': '38px',
                'width': '100%',
                    }
            ), width=3
        )
    ],
    className="mb-3",  # Optional margin
),
    dash_table.DataTable(
        id='results-table',
        columns=[
            {"name": "Fit Parameters", "id": "Parameters", "type": "text"},
            {"name": "Unbound-state",    "id": "Unbound",    "type": "numeric"},
            {"name": "Bound-state",      "id": "Bound",      "type": "numeric"},
        ],
        data=results_data,
        editable=False,
        row_deletable=False,
        style_cell={'textAlign': 'center', 'minWidth': '120px'},
        style_header={'fontWeight': 'bold'},
        style_table={'overflowX': 'auto'},
    ),
], fluid=True)

################################################################################
# 4. The Fitting Callback
################################################################################
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

@app.callback(
    [Output('swv-graph', 'figure'),
     Output('results-table', 'data')],
    [Input('fit-button', 'n_clicks'),
     Input('swv-model-dropdown', 'value'),
     Input(ThemeSwitchAIO.ids.switch("theme"), "value")],
    [State('swv-table', 'data')]
)
def perform_fit(n_clicks, model_type, current_theme, table_data):
    template = "plotly_white" if not current_theme else "plotly_dark"
    if not n_clicks:
        fig = go.Figure()
        fig.update_layout(template=template)
        return fig, []

    # Extract data
    freqs, unbound_vals, bound_vals = [], [], []
    for row in table_data:
        try:
            f = float(row["Frequency"])
            ub = float(row["Unbound"])
            bd = float(row["Bound"])
            if f > 0 and ub >= 0 and bd >= 0:
                freqs.append(f)
                unbound_vals.append(ub)
                bound_vals.append(bd)
        except:
            pass

    freqs = np.array(freqs)
    unbound_vals = np.array(unbound_vals)
    bound_vals   = np.array(bound_vals)

    if len(freqs) < 2:
        # Need at least two points to do a basic fit
        fig = go.Figure()
        fig.update_layout(template=template)
        return fig, "Not enough data points to perform fit."

    # Depending on model_type, define initial guesses & param order
    if model_type == 'exponential':
        # (I0, k_et, f_decay1, f_decay2)
        p0_unbound = [
        np.max(unbound_vals) * Config.CF_UNBOUND_I0_MUL,    # I0
        Config.UNBOUND_K_ET,                                # k_et
        Config.UNBOUND_F_DECAY1,                            # f_decay1
        Config.UNBOUND_F_DECAY2,                            # f_decay2
        ]
        p0_bound = [
        np.max(bound_vals) * Config.CF_BOUND_I0_MUL,        # I0
        Config.BOUND_K_ET,                                  # k_et
        Config.BOUND_F_DECAY1,                              # f_decay1
        Config.BOUND_F_DECAY2,                              # f_decay2
        ]
    
    elif model_type == 'butler-volmer':
        # (I0, j0, alpha, eta0, k_freq, f_decay1, f_decay2)
        # Just an example
        p0_unbound = [np.max(unbound_vals), 
                    Config.BV_CURR_DENS_J0,
                    Config.BV_ALPHA,
                    Config.BV_ETA_0,
                    Config.BV_UNBOUND_K_EFF, 
                    Config.BV_UNBOUND_F_DECAY1,
                    Config.BV_UNBOUND_F_DECAY2,
                    ]
        p0_bound   = [np.max(bound_vals), 
                    Config.BV_CURR_DENS_J0,
                    Config.BV_ALPHA,
                    Config.BV_ETA_0,
                    Config.BV_BOUND_K_EFF,
                    Config.BV_BOUND_F_DECAY1,
                    Config.BV_BOUND_F_DECAY2,
                    ]
    elif model_type == 'marcus':
        # (I0, A_marcus, lambda_marcus, deltaG0, k_freq, f_decay1, f_decay2)
        p0_unbound = [np.max(unbound_vals), 
                    Config.MARCUS_PREEXP_A,
                    Config.MARCUS_REORG_E_LAMBDA,
                    Config.MARCUS_DELTAG,
                    Config.MARCUS_UNBOUND_K_FREQ,
                    Config.MARCUS_UNBOUND_F_DECAY1,
                    Config.MARCUS_UNBOUND_F_DECAY2,
                    ]
        p0_bound   = [np.max(bound_vals),                    
                    Config.MARCUS_PREEXP_A,
                    Config.MARCUS_REORG_E_LAMBDA,
                    Config.MARCUS_DELTAG,
                    Config.MARCUS_BOUND_K_FREQ,
                    Config.MARCUS_BOUND_F_DECAY1,
                    Config.MARCUS_BOUND_F_DECAY2,
                    ]
    else:
        return go.Figure(), []

    # Define a wrapper for curve_fit
    def model_wrapper(f, *params):
        return swv_model_eval(f, model_type, params)

    # Fit unbound
    try:
        popt_u, _ = curve_fit(model_wrapper, freqs, unbound_vals, p0=p0_unbound, maxfev=10000)
        unbound_fit = swv_model_eval(freqs, model_type, popt_u)
        residuals_u = unbound_vals - unbound_fit
        ss_res_u    = np.sum(residuals_u**2)
        ss_tot_u    = np.sum((unbound_vals - np.mean(unbound_vals))**2)
        r2_u        = 1 - (ss_res_u/ss_tot_u)
        print(f"""
                Model type: {model_type}
                Status    : Unbound
            """)
        if model_type == 'exponential':
            if popt_u is not None:
                I0u, keu, fd1u, fd2u = popt_u
                print(f"""
                        I0u: {I0u}
                        keu: {keu}
                        fd1u: {fd1u}
                        fd2u: {fd2u}
                    """)
            else:
                print(f'exponential fit failed: popt_u: {popt_u}')

        elif model_type == 'butler-volmer':
            if popt_u is not None:
                I0_u, j0_u, alpha_u, eta0_u, kfreq_u, fd1_u, fd2_u = popt_u
                print(f"""
                        I0_u: {I0_u}
                        j0_u: {j0_u}
                        alpha_u: {alpha_u}
                        eta0_u: {eta0_u}
                        kfreq_u: {kfreq_u}
                        fd1_u: {fd1_u}
                        fd2_u: {fd2_u}
                    """)
            else:
                print(f'butler-volmer fit failed: popt_u: {popt_u}')  

        elif model_type == 'marcus':
            if popt_u is not None:
                I0_u, A_u, lambda_u, deltag0_u, kfreq_u, fd1_u, fd2_u = popt_u
                print(f"""
                        I0_u: {I0_u}
                        A_u: {A_u}
                        lambda_u: {alpha_u}
                        deltag0_u: {eta0_u}
                        kfreq_u: {kfreq_u}
                        fd1_u: {fd1_u}
                        fd2_u: {fd2_u}
                    """)
            else:
                print(f'marcus fit failed: popt_u: {popt_b}') 
        else:
            pass

    except Exception as e:
        print(e)
        popt_u = None
        ss_res_u = r2_u = 0

    # Fit bound
    try:
        popt_b, _ = curve_fit(model_wrapper, freqs, bound_vals, p0=p0_bound, maxfev=10000)
        bound_fit = swv_model_eval(freqs, model_type, popt_b)
        residuals_b = bound_vals - bound_fit
        ss_res_b    = np.sum(residuals_b**2)
        ss_tot_b    = np.sum((bound_vals - np.mean(bound_vals))**2)
        r2_b        = 1 - (ss_res_b/ss_tot_b)
        print(f"""
                Model type: {model_type}
                Status    : Bound
            """)
        if model_type == 'exponential':
            if popt_b is not None:
                I0b, keb, fd1b, fd2b = popt_b
                print(f"""
                        I0b: {I0b}
                        keb: {keb}
                        fd1b: {fd1b}
                        fd2u: {fd2b}
                    """)
            else:
                print(f'exponential fit failed: popt_u: {popt_b}')

        elif model_type == 'butler-volmer':
            if popt_b is not None:
                I0_b, j0_b, alpha_b, eta0_b, kfreq_b, fd1_b, fd2_b = popt_b
                print(f"""
                        I0_b: {I0_b}
                        j0_b: {j0_b}
                        alpha_b: {alpha_b}
                        eta0_b: {eta0_b}
                        kfreq_b: {kfreq_b}
                        fd1_b: {fd1_b}
                        fd2_b: {fd2_b}
                    """)
            else:
                print(f'butler-volmer fit failed: popt_u: {popt_b}')  

        elif model_type == 'marcus':
            if popt_b is not None:
                I0_b, A_b, lambda_b, deltag0_b, kfreq_b, fd1_b, fd2_b = popt_b
                print(f"""
                        I0_b: {I0_b}
                        A_b: {A_b}
                        lambda_b: {alpha_b}
                        deltag0_b: {eta0_b}
                        kfreq_b: {kfreq_b}
                        fd1_b: {fd1_b}
                        fd2_b: {fd2_b}
                    """)
            else:
                print(f'marcus fit failed: popt_u: {popt_b}')  
        else:
            pass

    except Exception as e:
        print(e)
        popt_b = None
        ss_res_b = r2_b = 0

    f_min, f_max = np.min(freqs), np.max(freqs)
    f_plot = np.linspace(f_min, f_max, Config.PLT_LNSPC_NUM)
    step_size = (f_max - f_min)/(Config.PLT_LNSPC_NUM - 1)
    num_new_points = int(round((Config.PLT_EXTP_MAXFREQ - Config.PLT_EXTP_MINFREQ) / step_size)) + 1
    f_extp = np.linspace(Config.PLT_EXTP_MINFREQ, Config.PLT_EXTP_MAXFREQ, num_new_points)

    # Build figure
    fig = go.Figure()
    # Normalised data
    unbound_norm = unbound_vals / (np.max(unbound_vals) + Config.PLT_NODEVZERO_VAL)
    bound_norm   = bound_vals   / (np.max(bound_vals) + Config.PLT_NODEVZERO_VAL)
    fig.add_trace(go.Scatter(x=freqs, y=unbound_norm, mode='markers', name='Unbound (Data)', marker=dict(color='blue')))
    fig.add_trace(go.Scatter(x=freqs, y=bound_norm,   mode='markers', name='Bound (Data)',   marker=dict(color='red')))

    # Plot fitted curves
    f_plot = np.linspace(np.min(freqs), np.max(freqs), Config.PLT_LNSPC_NUM)
    if popt_u is not None:
        fit_u_plot = swv_model_eval(f_plot, model_type, popt_u)
        fit_u_norm = fit_u_plot / (np.max(fit_u_plot) + Config.PLT_NODEVZERO_VAL)

        extp_fit_u_plot = swv_model_eval(f_extp, model_type, popt_u)
        extp_fit_u_norm = extp_fit_u_plot / (np.max(extp_fit_u_plot) + Config.PLT_NODEVZERO_VAL)

        fig.add_trace(go.Scatter(x=f_plot, y=fit_u_norm, mode='lines', name='Unbound (Fit)', line=dict(color='blue')))
        fig.add_trace(go.Scatter(x=f_extp, y=extp_fit_u_norm, mode='lines', name='EXTP Unbound (Fit)', line=dict(color='cyan', dash='dash')))

    if popt_b is not None:
        fit_b_plot = swv_model_eval(f_plot, model_type, popt_b)
        fit_b_norm = fit_b_plot / (np.max(fit_b_plot) + Config.PLT_NODEVZERO_VAL)

        extp_fit_b_plot = swv_model_eval(f_extp, model_type, popt_b)
        extp_fit_b_norm = extp_fit_b_plot / (np.max(extp_fit_b_plot) + Config.PLT_NODEVZERO_VAL)

        fig.add_trace(go.Scatter(x=f_plot, y=fit_b_norm, mode='lines', name='Bound (Fit)', line=dict(color='red')))
        fig.add_trace(go.Scatter(x=f_extp, y=extp_fit_b_norm, mode='lines', name='EXTP Bound (Fit)', line=dict(color='orange', dash='dash')))

    fig.update_layout(
        title=f"{model_type.capitalize()} Model Fit",
        xaxis_title="Frequency / (Hz)",
        yaxis_title="Normalised Charge / (C)",
        template=template
    )

    # Prepare results table
    results_table = []
    if model_type == 'exponential':
        # popt_u, popt_b each has 4 params: (I0, k_et, f_decay1, f_decay2)
        if popt_u is not None:
            I0u, keu, fd1u, fd2u = popt_u
        else:
            I0u, keu, fd1u, fd2u = (0,0,0,0)
        if popt_b is not None:
            I0b, keb, fd1b, fd2b = popt_b
        else:
            I0b, keb, fd1b, fd2b = (0,0,0,0)
        results_table = [
            {"Parameters": "I0",       "Unbound": f"{I0u:.4g}", "Bound": f"{I0b:.4g}"},
            {"Parameters": "k_et",     "Unbound": f"{keu:.4g}", "Bound": f"{keb:.4g}"},
            {"Parameters": "f_decay1", "Unbound": f"{fd1u:.4g}","Bound": f"{fd1b:.4g}"},
            {"Parameters": "f_decay2", "Unbound": f"{fd2u:.4g}","Bound": f"{fd2b:.4g}"},
            {"Parameters": "SSR",      "Unbound": f"{ss_res_u:.4g}", "Bound": f"{ss_res_b:.4g}"},
            {"Parameters": "R²",       "Unbound": f"{r2_u:.4g}",     "Bound": f"{r2_b:.4g}"},
        ]
    elif model_type == 'butler-volmer':
        # popt has 7 params: (I0, j0, alpha, eta0, k_freq, f_decay1, f_decay2)
        if popt_u is not None:
            I0_u, j0_u, alpha_u, eta0_u, kfreq_u, fd1_u, fd2_u = popt_u
        else:
            I0_u, j0_u, alpha_u, eta0_u, kfreq_u, fd1_u, fd2_u = (0,0,0,0,0,0,0)

        if popt_b is not None:
            I0_b, j0_b, alpha_b, eta0_b, kfreq_b, fd1_b, fd2_b = popt_b
        else:
            I0_b, j0_b, alpha_b, eta0_b, kfreq_b, fd1_b, fd2_b = (0,0,0,0,0,0,0)

        results_table = [
            {"Parameters": "I0",       "Unbound": f"{I0_u:.4g}", "Bound": f"{I0_b:.4g}"},
            {"Parameters": "J0",       "Unbound": f"{j0_u:.4g}", "Bound": f"{j0_b:.4g}"},
            {"Parameters": "alpha",    "Unbound": f"{alpha_u:.4g}", "Bound": f"{alpha_b:.4g}"},
            {"Parameters": "eta0",     "Unbound": f"{eta0_u:.4g}", "Bound": f"{eta0_b:.4g}"},
            {"Parameters": "k_freq",   "Unbound": f"{kfreq_u:.4g}", "Bound": f"{kfreq_b:.4g}"},
            {"Parameters": "f_decay1", "Unbound": f"{fd1_u:.4g}","Bound": f"{fd1_b:.4g}"},
            {"Parameters": "f_decay2", "Unbound": f"{fd2_u:.4g}","Bound": f"{fd2_b:.4g}"},
            {"Parameters": "SSR",      "Unbound": f"{ss_res_u:.4g}", "Bound": f"{ss_res_b:.4g}"},
            {"Parameters": "R²",       "Unbound": f"{r2_u:.4g}",     "Bound": f"{r2_b:.4g}"},
        ]

    elif model_type == 'marcus':
        # popt has 6 params: (I0, A_marcus, lambda_marcus, deltaG0, k_freq, f_decay1, f_decay2)
        if popt_u is not None:
            I0_u, A_u, lambda_u, deltag0_u, kfreq_u, fd1_u, fd2_u = popt_u
        else:
            I0_u, A_u, lambda_u, deltag0_u, kfreq_u, fd1_u, fd2_u = (0,0,0,0,0,0,0)

        if popt_b is not None:
            I0_b, A_b, lambda_b, deltag0_b, kfreq_b, fd1_b, fd2_b = popt_b
        else:
            I0_b, A_b, lambda_b, deltag0_b, kfreq_b, fd1_b, fd2_b = (0,0,0,0,0,0,0)

        results_table = [
            {"Parameters": "I0",       "Unbound": f"{I0_u:.4g}", "Bound": f"{I0_b:.4g}"},
            {"Parameters": "A",       "Unbound": f"{A_u:.4g}", "Bound": f"{A_b:.4g}"},
            {"Parameters": "lambda",    "Unbound": f"{lambda_u:.4g}", "Bound": f"{lambda_b:.4g}"},
            {"Parameters": "deltag0",     "Unbound": f"{deltag0_u:.4g}", "Bound": f"{deltag0_b:.4g}"},
            {"Parameters": "k_freq",   "Unbound": f"{kfreq_u:.4g}", "Bound": f"{kfreq_b:.4g}"},
            {"Parameters": "f_decay1", "Unbound": f"{fd1_u:.4g}","Bound": f"{fd1_b:.4g}"},
            {"Parameters": "f_decay2", "Unbound": f"{fd2_u:.4g}","Bound": f"{fd2_b:.4g}"},
            {"Parameters": "SSR",      "Unbound": f"{ss_res_u:.4g}", "Bound": f"{ss_res_b:.4g}"},
            {"Parameters": "R²",       "Unbound": f"{r2_u:.4g}",     "Bound": f"{r2_b:.4g}"},
        ]

    return fig, results_table

################################################################################
# Run the app
################################################################################

if __name__ == "__main__":
    app.run_server(debug=True)
