import numpy as np
import dash
from dash import dcc, html, callback_context
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.graph_objects as go
from scipy.integrate import odeint
from dash_bootstrap_templates import ThemeSwitchAIO

###############################################################################
# Global Constants and Helper Functions
###############################################################################

EPSILON = 1e-12

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

###############################################################################
# AptamerBindingModel Class
###############################################################################

class AptamerBindingModel:
    """
    Encapsulates the math/physics/chemistry for an aptamer-based sensor:
    - Binding kinetics
    - Hill–Langmuir model
    - SWV frequency response
    - Butler–Volmer equation
    - Marcus theory (in eV)
    """
    
    F = 96485.3329       # Faraday's constant (C/mol)
    R = 8.314            # Universal gas constant (J/(mol·K))
    kB = 1.380649e-23    # Boltzmann constant (J/K)

    def __init__(self):
        # -------------------------------------------------
        # 1. Binding Kinetics Parameters
        # -------------------------------------------------
        self.A_total = 1.0        # Total aptamer concentration (arbitrary units)
        self.T_conc = 10.0        # Target concentration
        self.k_on = 1.0           # Binding rate constant (1/s)
        self.k_off = 0.5          # Unbinding rate constant (1/s)
        self.I_min = 0.1          # Minimum current
        self.I_max = 1.0          # Maximum current
        self.AT0 = 0.0            # Initial bound aptamer concentration
        self.n = 1.0              # Hill coefficient

        # -------------------------------------------------
        # 2. SWV Simulation Parameters
        # -------------------------------------------------
        self.k_et_unbound = 10
        self.k_et_bound   = 50
        self.f_decay      = 100
        self.I0           = 1.0

        # -------------------------------------------------
        # 3. Butler–Volmer Parameters
        # -------------------------------------------------
        self.j0    = 1e-3   # Exchange current density (A/cm²)
        self.alpha = 0.5    # Charge transfer coefficient
        self.T     = 298.15 # Temperature (K)
        # For integration: baseline overpotential (V) when no target is bound.
        self.eta0  = 0.5

        # -------------------------------------------------
        # 4. Marcus Theory Parameters (in eV)
        # -------------------------------------------------
        self.A_marcus      = 1e12      # Pre-exponential factor (s⁻¹)
        self.lambda_marcus = 0.2       # Reorganisation energy (eV)
        self.deltaG0       = -0.1      # Standard free energy change (eV)
        # For integration: shift in free energy upon binding (eV)
        self.deltaG_shift  = 0.1

        self.eps = EPSILON

        # Experimental data
        self.freq_low = 16       # Frequency below the critical frequency (e.g. 16 Hz)
        self.freq_high = 240     # Frequency above the critical frequency (e.g. 240 Hz)
        # Calculated charge (from integration) at these frequencies (in, say, µC)
        self.charge_low_unbound  = 5.0
        self.charge_low_bound    = 5.0
        self.charge_high_unbound = 2.0
        self.charge_high_bound   = 2.0
        # Electrode geometry (nominal area in mm²)
        self.electrode_area = 15.0
        # Number of electrons transferred during MB redox (default 2)
        self.mb_electron = 2

    ###########################################################################
    # 1. Aptamer Binding Kinetics
    ###########################################################################
    def _binding_ode(self, AT, t):
        A_free = self.A_total - AT
        return self.k_on * A_free * self.T_conc - self.k_off * AT

    def compute_binding_curve(self, t_end=4, num_points=200):
        t = np.linspace(0, t_end, num_points)
        AT = odeint(self._binding_ode, self.AT0, t).flatten()
        fraction_bound = AT / (self.A_total + self.eps)
        I = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return t, fraction_bound, I

    ###########################################################################
    # 2. Hill–Langmuir Isotherm
    ###########################################################################
    def hill_langmuir(self, T, Kd, n):
        return np.power(T, n) / (np.power(Kd, n) + np.power(T, n) + self.eps)

    def compute_hill_langmuir_curve(self, t_min=0.01, t_max=10, num_points=1000):
        Kd = self.k_off / (self.k_on + self.eps)
        T_conc_range = np.logspace(np.log10(t_min), np.log10(t_max), num_points)
        fraction_bound = self.hill_langmuir(T_conc_range, Kd, self.n)
        I_hl = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return T_conc_range, fraction_bound, I_hl

    ###########################################################################
    # 3. Square Wave Voltammetry (SWV) - inverse frequency transform
    ###########################################################################

    def compute_swv_curve(self, num_points=200):
        # Use the experimental parameters for frequency bounds
        frequencies = np.linspace(self.freq_low, self.freq_high, num_points)
        
        # Calculate net SWV current (for unbound vs. bound)
        denom_unbound = (2 * frequencies + self.eps)
        denom_bound   = (2 * frequencies + self.eps)
        I_net_unbound = self.I0 * np.exp(-self.k_et_unbound / denom_unbound) * np.exp(-frequencies / (self.f_decay + self.eps))
        I_net_bound   = self.I0 * np.exp(-self.k_et_bound   / denom_bound)   * np.exp(-frequencies / (self.f_decay + self.eps))
        
        # Transform: (peak current) / frequency
        I_u_F = I_net_unbound / frequencies
        I_b_F = I_net_bound   / frequencies
        
        # x-axis: 1/frequency
        inv_freq = 1 / frequencies
        
        # Normalize for plotting
        I_u_norm = I_u_F / (np.max(I_u_F) + self.eps)
        I_b_norm = I_b_F / (np.max(I_b_F) + self.eps)
        
        # Find crossover point
        diff_array = np.abs(I_u_norm - I_b_norm)
        crossover_index = np.argmin(diff_array)
        crossover_freq = frequencies[crossover_index]
        crossover_inv_freq = inv_freq[crossover_index]
        
        return inv_freq, I_u_F, I_b_F, I_u_norm, I_b_norm, crossover_freq, crossover_inv_freq


    ###########################################################################
    # 4. Butler–Volmer Equation
    ###########################################################################
    def butler_volmer(self, eta):
        term_forward = np.exp((1 - self.alpha) * self.F * eta / (self.R * self.T))
        term_reverse = np.exp(-self.alpha * self.F * eta / (self.R * self.T))
        return self.j0 * (term_forward - term_reverse)

    def compute_butler_volmer_curve(self, eta_min=-0.5, eta_max=0.5, num_points=200):
        eta_range = np.linspace(eta_min, eta_max, num_points)
        j_values = self.butler_volmer(eta_range)
        return eta_range, j_values

    ###########################################################################
    # 5. Marcus Theory (All in eV)
    ###########################################################################
    def marcus_rate(self, deltaG_eV: float) -> float:
        k_B_eV = 8.617333262145e-5
        denom = 4.0 * max(self.lambda_marcus, self.eps) * k_B_eV * self.T + self.eps
        numerator = (self.lambda_marcus + deltaG_eV)**2
        exponent = - (numerator / denom)
        return self.A_marcus * np.exp(exponent)

    def compute_marcus_curve(self, deltaG0_eV: float = None, window_eV=0.5, num_points=200):
        if deltaG0_eV is None:
            deltaG0_eV = self.deltaG0
        deltaG_min = deltaG0_eV - window_eV
        deltaG_max = deltaG0_eV + window_eV
        deltaG_range = np.linspace(deltaG_min, deltaG_max, num_points)
        k_et = np.array([self.marcus_rate(dG) for dG in deltaG_range])
        return deltaG_range, k_et

    ###########################################################################
    # 6. Integrated Model
    ###########################################################################
    def compute_integrated_current_curve(self, method='butler-volmer', 
                                         t_min=0.01, t_max=10, num_points=1000):
        Kd = self.k_off / (self.k_on + self.eps)
        T_conc_range = np.logspace(np.log10(t_min), np.log10(t_max), num_points)
        fraction_bound = self.hill_langmuir(T_conc_range, Kd, self.n)
        integrated_current = np.zeros_like(T_conc_range)
        
        if method == 'butler-volmer':
            # reference
            j_max = self.butler_volmer(self.eta0)
            for i, theta in enumerate(fraction_bound):
                eta_eff = self.eta0 * (1 - theta)
                j = self.butler_volmer(eta_eff)
                j_norm = j / (j_max + self.eps)
                integrated_current[i] = self.I_min + (self.I_max - self.I_min) * j_norm

        elif method == 'marcus':
            # reference
            k_max = self.marcus_rate(self.deltaG0)
            for i, theta in enumerate(fraction_bound):
                effective_deltaG = self.deltaG0 - self.deltaG_shift * theta
                k_et = self.marcus_rate(effective_deltaG)
                k_norm = k_et / (k_max + self.eps)
                integrated_current[i] = self.I_min + (self.I_max - self.I_min) * k_norm

        else:
            # fallback: linear HL response
            integrated_current = self.I_min + (self.I_max - self.I_min) * fraction_bound
        
        return T_conc_range, fraction_bound, integrated_current

###############################################################################
# Dash App Initialisation and Layout
###############################################################################
THEMES = [dbc.themes.COSMO, dbc.themes.DARKLY]
app = dash.Dash(__name__, external_stylesheets=THEMES)
model = AptamerBindingModel()

app.layout = dbc.Container([
    html.H1("Electrochemical Aptamer-Based (E-AB) Sensor Simulation",
            style={'marginTop': '5px', 'fontSize': '40px', 
                   'marginBottom': '5px', 'padding': '5px'}),

    dbc.Row([
        dbc.Col(ThemeSwitchAIO(aio_id="theme", themes=THEMES), width=12)
    ]),

    # -------------------------------
    # 1. Binding Kinetics Layout
    # -------------------------------
    dbc.Row([
        dbc.Col([dcc.Graph(id='kinetics-plot', style={'height': '40vh'})], width=9, 
                style={'flex': '0 0 40%', 'maxWidth': '40%'}),
        dbc.Col([dcc.Graph(id='hill-langmuir-plot', style={'height': '40vh'})], width=9, 
                style={'flex': '0 0 40%', 'maxWidth': '40%'}),
        dbc.Col([
            html.H2("Binding Kinetics Parameters", style={'color': 'orange'}),
            slider_with_label("k_on (Binding Rate Constant, 1/s)", 
                              'kon-slider', 0.1, 5.0, 0.01, model.k_on,
                              marks={1: "1", 5: "5"}, 
                              tooltip_text="Controls how fast aptamer binds the target."),
            slider_with_label("k_off (Unbinding Rate Constant, 1/s)", 
                              'koff-slider', 0.1, 5.0, 0.01, model.k_off,
                              marks={1: "1", 5: "5"}, 
                              tooltip_text="Controls how fast the aptamer dissociates."),
            slider_with_label("Target Concentration (uM)", 
                              'Tconc-slider', 1, 50, 0.1, model.T_conc,
                              marks={1: "1", 50: "50"}),
            slider_with_label("I_min (Minimum Current)", 
                              'Imin-slider', 0.0, 0.5, 0.01, model.I_min,
                              marks={0.0: "0.0", 0.5: "0.5"}),
            slider_with_label("I_max (Maximum Current)", 
                              'Imax-slider', 0.5, 2.0, 0.01, model.I_max,
                              marks={0.5: "0.5", 2.0: "2.0"}),
            slider_with_label("Hill Coefficient (n)", 
                              'hill-n-slider', 0.5, 4.0, 0.01, model.n,
                              marks={0.5: "0.5", 4.0: "4.0"}, 
                              tooltip_text="Controls cooperativity in binding.")
        ], width=3, style={'flex': '0 0 20%', 'maxWidth': '20%'})
    ], style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}),
    
    html.Hr(),

    # -------------------------------
    # 3. Integrated Model Layout
    # -------------------------------
    dbc.Row([
        dbc.Col([dcc.Graph(id='integrated-model-plot')], width=9,
                style={'flex': '0 0 80%', 'maxWidth': '80%'}),
        dbc.Col([
            html.H2("Integrated Model Parameters", style={'color': 'orange'}),
            dcc.Dropdown(
                id='integrated-method-dropdown',
                options=[
                    {'label': 'Butler-Volmer', 'value': 'butler-volmer'},
                    {'label': 'Marcus', 'value': 'marcus'}
                ],
                value='butler-volmer', 
                style={'background': 'cyan', 'height': '40px', 'font-size': '120%'}
            ),
        ], width=3, style={'flex': '0 0 20%', 'maxWidth': '20%'})
    ], style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}),

    html.Hr(),

    # -------------------------------
    # 4. Butler-Volmer Layout
    # -------------------------------
    dbc.Row([
        dbc.Col([dcc.Graph(id='butler-volmer-plot')], width=9,
                style={'flex': '0 0 40%', 'maxWidth': '40%'}),
        dbc.Col([dcc.Graph(id='marcus-plot')], width=9, 
                style={'flex': '0 0 40%', 'maxWidth': '40%'}),
        dbc.Col([
            html.H2("Butler-Volmer Parameters", style={'color': 'orange'}),
            slider_with_label("Exchange Current Density (j0)", 'j0-slider', 
                              1e-4, 1e-2, 1e-4, model.j0,
                              marks={1e-4: "1e-4", 1e-3: "1e-3", 1e-2: "1e-2"}),
            slider_with_label("Charge Transfer Coefficient (alpha)", 'alpha-slider', 
                              0.1, 1.0, 0.01, model.alpha,
                              marks={round(x,1): str(round(x,1)) for x in 
                                     [0.1 * i for i in range(1, 11)]}),
            slider_with_label("Reorganisation Energy (eV)", 'lambda-slider', 
                              0.1, 1.0, 0.01, model.lambda_marcus,
                              marks={round(x,1): str(round(x,1)) 
                                     for x in [0.1 * i for i in range(1, 11)]}),
            slider_with_label("Standard Free Energy (eV)", 'deltaG-slider', 
                              -0.5, 0.5, 0.01, model.deltaG0,
                              marks={round(x,1): str(round(x,1)) 
                                     for x in [0.1 * i for i in range(-5, 6)]}),
        ], width=3, style={'flex': '0 0 20%', 'maxWidth': '20%'})
    ], style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}),
    
    html.Hr(),

    html.H2("SWV Simulation | Experimental Parameters", style={'color': 'orange'}),
    dbc.Row([
        dbc.Col([dcc.Graph(id='swv-plot', style={'height': '100vh'})], width=9, 
                style={'flex': '0 0 80%', 'maxWidth': '80%'}),
        dbc.Col([
            slider_with_label("k_et Unbound (s⁻¹)", 'ket-unbound-slider', 
                              1, 100, 0.1, model.k_et_unbound,
                              marks={1: "1", 100: "100"}),
            slider_with_label("k_et Bound (s⁻¹)", 'ket-bound-slider', 
                              1, 100, 0.1, model.k_et_bound,
                              marks={1: "1", 100: "100"}),
            slider_with_label("f_decay (Hz)", 'fdecay-slider', 
                              50, 300, 1, model.f_decay,
                              marks={50: "50", 300: "300"}),
            slider_with_label("I0 (Initial Current)", 'I0-slider', 
                              0.1, 2.0, 0.01, model.I0,
                              marks={0.1: "0.1", 2.0: "2.0"}),
            slider_input_with_label("Frequency (Low, Hz)", 
                                    "freq-low-slider", "freq-low-input", 
                                    min_val=1, max_val=50, step=1, 
                                    value=model.freq_low,
                                    marks={1:"1", 50:"50"}),
            slider_input_with_label("Frequency (High, Hz)", 
                                    "freq-high-slider", "freq-high-input", 
                                    min_val=100, max_val=500, step=1, 
                                    value=model.freq_high,
                                    marks={100:"100", 500:"500"}),
            slider_input_with_label("Charge at Low Frequency Unbound (µC)", 
                                    "charge-low-unbound-slider", "charge-low-unbound-input",
                                    min_val=0, max_val=20, step=0.1, 
                                    value=model.charge_low_unbound),
            slider_input_with_label("Charge at Low Frequency Bound (µC)", 
                                    "charge-low-bound-slider", "charge-low-bound-input",
                                    min_val=0, max_val=20, step=0.1, 
                                    value=model.charge_low_bound),
            slider_input_with_label("Charge at High Frequency Unbound (µC)", 
                                    "charge-high-unbound-slider", "charge-high-unbound-input",
                                    min_val=0, max_val=20, step=0.1, 
                                    value=model.charge_high_unbound),
            slider_input_with_label("Charge at High Frequency Bound (µC)", 
                                    "charge-high-bound-slider", "charge-high-bound-input",
                                    min_val=0, max_val=20, step=0.1, 
                                    value=model.charge_high_bound),
            slider_input_with_label("Electrode Area (mm²)", 
                                    "electrode-area-slider", "electrode-area-input",
                                    min_val=1, max_val=100, step=0.1, 
                                    value=model.electrode_area),
            slider_input_with_label("MB Electron Transfer (e⁻)", 
                                    "mb-electron-slider", "mb-electron-input",
                                    min_val=1, max_val=5, step=1, 
                                    value=model.mb_electron),
            html.Div(id='crossover-info', style={'fontSize': 20, 'fontWeight': 'bold', 'padding': '10px'})
        ], width=3, style={'flex': '0 0 20%', 'maxWidth': '20%'})
    ], style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'})
], fluid=True)

###############################################################################
# Synchronizing Slider & Input Callbacks (unchanged)
###############################################################################
@app.callback(
    [Output("freq-low-slider", "value"), Output("freq-low-input", "value")],
    [Input("freq-low-slider", "value"), Input("freq-low-input", "value")]
)
def sync_freq_low(slider_val, input_val):
    ctx = callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggered_id == "freq-low-slider":
        return slider_val, slider_val
    else:
        return input_val, input_val

@app.callback(
    [Output("freq-high-slider", "value"), Output("freq-high-input", "value")],
    [Input("freq-high-slider", "value"), Input("freq-high-input", "value")]
)
def sync_freq_high(slider_val, input_val):
    ctx = callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggered_id == "freq-high-slider":
        return slider_val, slider_val
    else:
        return input_val, input_val

@app.callback(
    [Output("charge-low-unbound-slider", "value"), 
     Output("charge-low-unbound-input", "value")],
    [Input("charge-low-unbound-slider", "value"), 
     Input("charge-low-unbound-input", "value")]
)
def sync_charge_low_unbound(slider_val, input_val):
    ctx = callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggered_id == "charge-low-unbound-slider":
        return slider_val, slider_val
    else:
        return input_val, input_val

@app.callback(
    [Output("charge-low-bound-slider", "value"), 
     Output("charge-low-bound-input", "value")],
    [Input("charge-low-bound-slider", "value"), 
     Input("charge-low-bound-input", "value")]
)
def sync_charge_low_bound(slider_val, input_val):
    ctx = callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggered_id == "charge-low-bound-slider":
        return slider_val, slider_val
    else:
        return input_val, input_val

@app.callback(
    [Output("charge-high-unbound-slider", "value"), 
     Output("charge-high-unbound-input", "value")],
    [Input("charge-high-unbound-slider", "value"), 
     Input("charge-high-unbound-input", "value")]
)
def sync_charge_high_unbound(slider_val, input_val):
    ctx = callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggered_id == "charge-high-unbound-slider":
        return slider_val, slider_val
    else:
        return input_val, input_val

@app.callback(
    [Output("charge-high-bound-slider", "value"), 
     Output("charge-high-bound-input", "value")],
    [Input("charge-high-bound-slider", "value"), 
     Input("charge-high-bound-input", "value")]
)
def sync_charge_high_bound(slider_val, input_val):
    ctx = callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggered_id == "charge-high-bound-slider":
        return slider_val, slider_val
    else:
        return input_val, input_val

@app.callback(
    [Output("electrode-area-slider", "value"), 
     Output("electrode-area-input", "value")],
    [Input("electrode-area-slider", "value"), 
     Input("electrode-area-input", "value")]
)
def sync_electrode_area(slider_val, input_val):
    ctx = callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggered_id == "electrode-area-slider":
        return slider_val, slider_val
    else:
        return input_val, input_val

@app.callback(
    [Output("mb-electron-slider", "value"), 
     Output("mb-electron-input", "value")],
    [Input("mb-electron-slider", "value"), 
     Input("mb-electron-input", "value")]
)
def sync_mb_electron(slider_val, input_val):
    ctx = callback_context
    if not ctx.triggered:
        raise dash.exceptions.PreventUpdate
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    if triggered_id == "mb-electron-slider":
        return slider_val, slider_val
    else:
        return input_val, input_val

###############################################################################
# SINGLE Callback Updating All Plots
###############################################################################
@app.callback(
    [
        Output('kinetics-plot', 'figure'),
        Output('hill-langmuir-plot', 'figure'),
        Output('integrated-model-plot', 'figure'),
        Output('butler-volmer-plot', 'figure'),
        Output('marcus-plot', 'figure'),
        Output('swv-plot', 'figure'),
        Output('crossover-info', 'children')
    ],
    [
        Input('kon-slider', 'value'),
        Input('koff-slider', 'value'),
        Input('Tconc-slider', 'value'),
        Input('Imin-slider', 'value'),
        Input('Imax-slider', 'value'),
        Input('hill-n-slider', 'value'),
        Input('j0-slider', 'value'),
        Input('alpha-slider', 'value'),
        Input('lambda-slider', 'value'),
        Input('deltaG-slider', 'value'),
        Input('ket-unbound-slider', 'value'),
        Input('ket-bound-slider', 'value'),
        Input('fdecay-slider', 'value'),
        Input('I0-slider', 'value'),
        # experimental parameters
        Input("freq-low-slider", "value"),
        Input("freq-high-slider", "value"),
        Input("charge-low-unbound-slider", "value"),
        Input("charge-low-bound-slider", "value"),
        Input("charge-high-unbound-slider", "value"),
        Input("charge-high-bound-slider", "value"),
        Input("electrode-area-slider", "value"),
        Input("mb-electron-slider", "value"),
        # integrated model method
        Input('integrated-method-dropdown', 'value'),
        # theme
        Input(ThemeSwitchAIO.ids.switch("theme"), "value")
    ]
)
def update_all_plots(k_on, k_off, T_conc, I_min, I_max, hill_n,
                     j0, alpha, lambda_eV, deltaG_eV,
                     k_et_unbound, k_et_bound, f_decay, I0,
                     freq_low, freq_high, charge_low_unbound, charge_low_bound,
                     charge_high_unbound, charge_high_bound, electrode_area, mb_electron,
                     integrated_method, theme_switch):

    # Set up plotly template based on theme toggle
    plotly_template = "plotly_white" if theme_switch else "plotly_dark"

    # Update model with all slider/input values
    model.k_on = k_on
    model.k_off = k_off
    model.T_conc = T_conc
    model.I_min = I_min
    model.I_max = I_max
    model.n = hill_n
    model.j0 = j0
    model.alpha = alpha
    model.lambda_marcus = lambda_eV
    model.deltaG0 = deltaG_eV
    model.k_et_unbound = k_et_unbound
    model.k_et_bound = k_et_bound
    model.f_decay = f_decay
    model.I0 = I0

    model.freq_low = freq_low
    model.freq_high = freq_high
    model.charge_low_unbound = charge_low_unbound
    model.charge_low_bound = charge_low_bound
    model.charge_high_unbound = charge_high_unbound
    model.charge_high_bound = charge_high_bound
    model.electrode_area = electrode_area
    model.mb_electron = mb_electron

    # -------------------------------------------------
    # 1) Kinetics Plot
    # -------------------------------------------------
    t, frac_bound, I_kin = model.compute_binding_curve()
    fig_kin = go.Figure()
    fig_kin.add_trace(go.Scatter(x=t, y=frac_bound, 
                                 mode='lines', name='Fraction Bound'))
    fig_kin.add_trace(go.Scatter(x=t, y=I_kin, 
                                 mode='lines', name='Current (I)'))
    fig_kin.update_layout(title='Binding Kinetics', 
                          xaxis_title='Time (s)',
                          yaxis_title='Fraction Bound / Current',
                          template=plotly_template)

    # -------------------------------------------------
    # 2) Hill–Langmuir Plot
    # -------------------------------------------------
    T_range, frac_HL, I_HL = model.compute_hill_langmuir_curve()
    fig_hl = go.Figure()
    fig_hl.add_trace(go.Scatter(x=T_range, y=frac_HL, 
                                mode='lines', name='Fraction Bound'))
    fig_hl.add_trace(go.Scatter(x=T_range, y=I_HL, 
                                mode='lines', name='Current (I)', 
                                yaxis='y2'))
    # Double-y axis example
    fig_hl.update_layout(title='Hill-Langmuir Isotherm',
                         xaxis_title='[Target] (µM, log scale)',
                         yaxis=dict(title='Fraction Bound', side='left'),
                         yaxis2=dict(title='Current (I)', overlaying='y', 
                                     side='right', showgrid=False),
                         xaxis_type='log',
                         template=plotly_template)

    # -------------------------------------------------
    # 3) Integrated Model Plot (Butler-Volmer or Marcus)
    # -------------------------------------------------
    T_int, fbound_int, I_int = model.compute_integrated_current_curve(
        method=integrated_method
    )
    fig_int = go.Figure()
    fig_int.add_trace(go.Scatter(x=T_int, y=fbound_int, 
                                 mode='lines', name='Fraction Bound'))
    fig_int.add_trace(go.Scatter(x=T_int, y=I_int, 
                                 mode='lines', name='Integrated Current', 
                                 yaxis='y2'))
    fig_int.update_layout(
        title=f'Integrated Model ({integrated_method})',
        xaxis_title='[Target] (µM, log scale)',
        yaxis=dict(title='Fraction Bound', side='left'),
        yaxis2=dict(title='Current (I)', overlaying='y', 
                    side='right', showgrid=False),
        xaxis_type='log',
        template=plotly_template
    )

    # -------------------------------------------------
    # 4) Butler–Volmer Plot
    # -------------------------------------------------
    eta_range, j_values = model.compute_butler_volmer_curve()
    fig_bv = go.Figure()
    fig_bv.add_trace(go.Scatter(x=eta_range, y=j_values, 
                                mode='lines', name='j(η)'))
    fig_bv.update_layout(title='Butler–Volmer',
                         xaxis_title='Overpotential (V)',
                         yaxis_title='Current Density (A/cm²)',
                         template=plotly_template)

    # -------------------------------------------------
    # 5) Marcus Plot
    # -------------------------------------------------
    dG_range, k_et_vals = model.compute_marcus_curve()
    fig_marcus = go.Figure()
    fig_marcus.add_trace(go.Scatter(x=dG_range, y=k_et_vals, 
                                    mode='lines', name='k_et(ΔG)'))
    fig_marcus.update_layout(title='Marcus Theory',
                             xaxis_title='ΔG (eV)',
                             yaxis_title='ET Rate (s⁻¹)',
                             template=plotly_template)

    # -------------------------------------------------
    # 6) SWV Plot + crossover info
    # -------------------------------------------------
    inv_freqs, I_u_F, I_b_F, I_u_norm, I_b_norm, crossover_freq, crossover_inv_freq = model.compute_swv_curve()
    fig_swv = go.Figure()
    fig_swv.add_trace(go.Scatter(x=inv_freqs, y=I_u_norm, 
                                 mode='lines', name='Unbound (Norm)'))
    fig_swv.add_trace(go.Scatter(x=inv_freqs, y=I_b_norm, 
                                 mode='lines', name='Bound (Norm)'))
    # Add a vertical dashed line at the crossover
    fig_swv.add_trace(go.Scatter(
        x=[crossover_inv_freq, crossover_inv_freq], 
        y=[0, 1], mode='lines',
        line=dict(color='black', dash='dash'),
        name=f'Crossover ~ {crossover_freq:.2f} Hz'
    ))
    fig_swv.update_layout(
        title='SWV Frequency Response',
        xaxis_title='1 / Frequency (s)',
        yaxis_title='(Peak Current)/freq (normalized)',
        template=plotly_template
    )
    crossover_text = f"Crossover Frequency: {crossover_freq:.2f} Hz"

    # Return all 7 outputs
    return fig_kin, fig_hl, fig_int, fig_bv, fig_marcus, fig_swv, crossover_text

###############################################################################
# Run the App
###############################################################################
if __name__ == '__main__':
    app.run_server(debug=True, host='127.0.0.1', port=8050)
