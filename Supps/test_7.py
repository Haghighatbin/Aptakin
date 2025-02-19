import numpy as np
import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
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
    - label_text: Display text for the slider.
    - slider_id: Unique ID for dash callback.
    - min_val, max_val, step, value: numeric slider parameters.
    - marks: marks for the slider.
    - tooltip_text: text for an optional dbc.Tooltip.
    - label_style: style dict for the label (e.g., {'fontSize': '20px'}).
    """
    if label_style is None:
        label_style = {'fontSize': '20px'}

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

###############################################################################
# AptamerBindingModel Class
###############################################################################

# Small epsilon to avoid division-by-zero
EPSILON = 1e-12

class AptamerBindingModel:
    """
    Encapsulates the math/physics/chemistry for an aptamer-based sensor:
    - Binding kinetics
    - Hill–Langmuir model
    - SWV frequency response
    - Butler–Volmer equation
    - Marcus theory (in eV)
    """
    
    # Physical / fundamental constants for non-Marcus calculations
    F = 96485.3329       # Faraday's constant (C/mol)
    R = 8.314            # Universal gas constant (J/(mol·K))
    kB = 1.380649e-23    # Boltzmann constant (J/K) for standard usage in other eqns

    def __init__(self):
        """
        Initialize the model with default parameters. 
        Modify them at runtime via sliders or function calls.
        """
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

        # -------------------------------------------------
        # 4. Marcus Theory Parameters (in eV)
        # -------------------------------------------------
        # All energies are in eV; we use k_B (in eV/K) inside marcus_rate.
        self.A_marcus      = 1e12      # Pre-exponential factor (s⁻¹)
        self.lambda_marcus = 0.2       # Reorganisation energy (eV)
        self.deltaG0       = -0.1      # Standard free energy change (eV)

        # Epsilon for safe divisions
        self.eps = EPSILON

    ###########################################################################
    # 1. Aptamer Binding Kinetics
    ###########################################################################
    def _binding_ode(self, AT, t):
        """
        ODE describing the rate of change of bound aptamer concentration.
        AT: bound aptamer concentration
        t : time
        """
        A_free = self.A_total - AT
        return self.k_on * A_free * self.T_conc - self.k_off * AT

    def compute_binding_curve(self, t_end=4, num_points=200):
        """
        Solve the ODE for a given time range to get fraction bound & sensor current.
        Returns time array, fraction_bound, current.
        """
        t = np.linspace(0, t_end, num_points)
        AT = odeint(self._binding_ode, self.AT0, t).flatten()
        fraction_bound = AT / (self.A_total + self.eps)
        I = self.I_min + (self.I_max - self.I_min) * fraction_bound

        return t, fraction_bound, I

    ###########################################################################
    # 2. Hill–Langmuir Isotherm
    ###########################################################################
    def hill_langmuir(self, T, Kd, n):
        """
        Hill–Langmuir equation for fraction bound:
        fraction_bound = T^n / (Kd^n + T^n)
        """
        return np.power(T, n) / (np.power(Kd, n) + np.power(T, n) + self.eps)

    def compute_hill_langmuir_curve(self, t_min=0.01, t_max=10, num_points=1000):
        """
        Compute fraction bound & current over a range of target concentrations.
        Uses the Hill–Langmuir equation with the current outputs.
        """
        Kd = self.k_off / (self.k_on + self.eps)
        T_conc_range = np.logspace(np.log10(t_min), np.log10(t_max), num_points)

        fraction_bound = self.hill_langmuir(T_conc_range, Kd, self.n)
        I_hl = self.I_min + (self.I_max - self.I_min) * fraction_bound

        return T_conc_range, fraction_bound, I_hl

    ###########################################################################
    # 3. Square Wave Voltammetry (SWV)
    ###########################################################################
    def compute_swv_curve(self, f_min=1, f_max=200, num_points=200):
        """
        Compute SWV frequency response for unbound vs. bound states.
        Returns frequency array, unbound/bound raw curves, normalised curves, and crossover freq.
        """
        frequencies = np.linspace(f_min, f_max, num_points)
        
        # Use a small epsilon to avoid dividing by zero
        denom_unbound = (2 * frequencies + self.eps)
        denom_bound   = (2 * frequencies + self.eps)

        I_net_unbound = self.I0 * np.exp(-self.k_et_unbound / denom_unbound) \
                              * np.exp(-frequencies / (self.f_decay + self.eps))
        I_net_bound   = self.I0 * np.exp(-self.k_et_bound   / denom_bound)   \
                              * np.exp(-frequencies / (self.f_decay + self.eps))

        I_net_unbound_norm = I_net_unbound / (np.max(I_net_unbound) + self.eps)
        I_net_bound_norm   = I_net_bound   / (np.max(I_net_bound)   + self.eps)

        # Find frequency where the two normalised curves are closest
        diff_array = np.abs(I_net_unbound_norm - I_net_bound_norm)
        crossover_index = np.argmin(diff_array)
        crossover_freq = frequencies[crossover_index]

        return (frequencies, 
                I_net_unbound, 
                I_net_bound, 
                I_net_unbound_norm, 
                I_net_bound_norm, 
                crossover_freq)

    ###########################################################################
    # 4. Butler–Volmer Equation
    ###########################################################################
    def butler_volmer(self, eta):
        """
        Butler–Volmer equation:
        j = j0 [exp((1-alpha) * F * eta / (R * T)) - exp(-alpha * F * eta / (R * T))]
        """
        term_forward = np.exp((1 - self.alpha) * self.F * eta / (self.R * self.T))
        term_reverse = np.exp(-self.alpha * self.F * eta / (self.R * self.T))
        return self.j0 * (term_forward - term_reverse)

    def compute_butler_volmer_curve(self, eta_min=-0.5, eta_max=0.5, num_points=200):
        """
        Returns an array of overpotentials (eta) and the corresponding current densities j.
        """
        eta_range = np.linspace(eta_min, eta_max, num_points)
        j_values = self.butler_volmer(eta_range)
        return eta_range, j_values

    ###########################################################################
    # 5. Marcus Theory (All in eV)
    ###########################################################################
    def marcus_rate(self, deltaG_eV: float) -> float:
        """
        Compute the electron-transfer rate using Marcus theory in eV:
        
        k_et = A_marcus * exp[ - ((lambda + deltaG)^2 / (4 * lambda * kB_eV * T)) ]
        
        where:
          - lambda, deltaG are in eV
          - k_B_eV = 8.617333262145e-5 eV/K
          - T is in K
        """
        # Boltzmann constant in eV/K
        k_B_eV = 8.617333262145e-5
        
        # Avoid dividing by zero if lambda is very small
        denom = 4.0 * max(self.lambda_marcus, self.eps) * k_B_eV * self.T + self.eps
        
        numerator = (self.lambda_marcus + deltaG_eV)**2
        exponent = - (numerator / denom)

        return self.A_marcus * np.exp(exponent)

    def compute_marcus_curve(self, deltaG0_eV: float = None, window_eV=0.5, num_points=200):
        """
        Compute Marcus-theory electron-transfer rates for a range of deltaG (in eV) 
        around deltaG0_eV (± window_eV).
        """
        if deltaG0_eV is None:
            # Use the model's stored deltaG0 in eV
            deltaG0_eV = self.deltaG0

        # Build an array around deltaG0: [deltaG0 - window_eV, deltaG0 + window_eV]
        deltaG_min = deltaG0_eV - window_eV
        deltaG_max = deltaG0_eV + window_eV

        deltaG_range = np.linspace(deltaG_min, deltaG_max, num_points)

        # Compute rates for each value
        k_et = np.array([self.marcus_rate(dG) for dG in deltaG_range])

        return deltaG_range, k_et


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
        dbc.Col([
            dcc.Graph(id='kinetics-plot')], width=9, style={'flex': '0 0 75%', 
                                                            'maxWidth': '75%',
                                                            'marginTop': '5px',
                                                            'marginBottom': '5px',
                                                            'padding': '5px'
                                                            }),

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
        ], width=3, style={'flex': '0 0 25%', 'maxWidth': '25%'})
    ], style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}),

    html.Hr(),

    # -------------------------------
    # 2. Hill-Langmuir Layout
    # -------------------------------
    dbc.Row([
        dbc.Col([
            dcc.Graph(id='hill-langmuir-plot')], width=9, style={'flex': '0 0 75%',
                                                                  'maxWidth': '75%'
                                                                  }),

        dbc.Col([
            html.H2("Hill-Langmuir Parameters", style={'color': 'orange'}),
            slider_with_label("Hill Coefficient (n)",
                              'hill-n-slider', 0.5, 4.0, 0.01, model.n,
                              marks={0.5: "0.5", 4.0: "4.0"},
                              tooltip_text="Controls cooperativity in binding.")
        ], width=3, style={'flex': '0 0 25%', 'maxWidth': '25%'})
    ], style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}),

    html.Hr(),

    # -------------------------------
    # 3. Butler-Volmer Layout
    # -------------------------------
    dbc.Row([
        dbc.Col([
            dcc.Graph(id='butler-volmer-plot')], width=9,style={'flex': '0 0 75%',
                                                                 'maxWidth': '75%'}),

        dbc.Col([
            html.H2("Butler-Volmer Parameters", style={'color': 'orange'}),
            slider_with_label("Exchange Current Density (j0)", 'j0-slider',
                              1e-4, 1e-2, 1e-4, model.j0,
                              marks={1e-4: "1e-4", 1e-3: "1e-3", 1e-2: "1e-2"}),
            slider_with_label("Charge Transfer Coefficient (alpha)", 'alpha-slider',
                              0.1, 1.0, 0.01, model.alpha,
                              marks = {round(x,1): str(round(x,1)) for x in [0.1 * i for i in range(1, 11)]}),
        ], width=3, style={'flex': '0 0 25%', 'maxWidth': '25%'})
    ], style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}),

    html.Hr(),

    # -------------------------------
    # 4. Marcus Theory Layout
    # -------------------------------
    dbc.Row([
        dbc.Col([
            dcc.Graph(id='marcus-plot')
        ], width=9, style={'flex': '0 0 75%', 'maxWidth': '75%'}),

        dbc.Col([
            html.H2("Marcus Theory Parameters", style={'color': 'orange'}),
            slider_with_label("Reorganisation Energy (eV)", 'lambda-slider',
                                min_val=0.1,
                                max_val=1.0,
                                step=0.01,
                                value=model.lambda_marcus,
                                marks = {round(x,1): str(round(x,1)) for x in [0.1 * i for i in range(1, 11)]}),
            slider_with_label("Standard Free Energy (eV)", 'deltaG-slider',
                                min_val=-0.5,
                                max_val=0.5,
                                step=0.01,
                                value=model.deltaG0,
                                marks = {round(x,1): str(round(x,1)) for x in [0.1 * i for i in range(-5, 6)]}),
        ], width=3)
    ], style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}),

    html.Hr(),

    # -------------------------------
    # 5. SWV Layout
    # -------------------------------
    dbc.Row([
        dbc.Col([
            dcc.Graph(id='swv-plot')
        ], width=9, style={'flex': '0 0 75%', 'maxWidth': '75%'}),

        dbc.Col([
            html.H2("SWV Simulation Parameters", style={'color': 'orange'}),
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

            html.Div(id='crossover-info', style={'fontSize': 18, 
                                                 'fontWeight': 'bold', 
                                                 'padding': '10px'})
        ], width=3, style={'flex': '0 0 25%', 'maxWidth': '25%'})
    ], style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'})
], fluid=True)

###############################################################################
# Dash Callback
###############################################################################
@app.callback(
    [Output('kinetics-plot', 'figure'),
     Output('hill-langmuir-plot', 'figure'),
     Output('butler-volmer-plot', 'figure'),
     Output('marcus-plot', 'figure'),
     Output('swv-plot', 'figure'),
     Output('crossover-info', 'children')],
    [Input('kon-slider', 'value'),
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
     Input(ThemeSwitchAIO.ids.switch("theme"), "value")]
)
def update_plots(k_on, k_off, T_conc, I_min, I_max, hill_n,
                 j0, alpha, lambda_eV, deltaG_eV,
                 k_et_unbound, k_et_bound, f_decay, I0,
                 theme_switch):

    # Light or dark template?
    plotly_template = "plotly_white" if theme_switch else "plotly_dark"

    # Update model parameters
    model.k_on        = k_on
    model.k_off       = k_off
    model.T_conc      = T_conc
    model.I_min       = I_min
    model.I_max       = I_max
    model.n           = hill_n
    model.j0          = j0
    model.alpha       = alpha
    model.lambda_marcus = lambda_eV
    model.deltaG0       = deltaG_eV
    model.k_et_unbound  = k_et_unbound
    model.k_et_bound    = k_et_bound
    model.f_decay       = f_decay
    model.I0            = I0

    # 1. Kinetics Plot
    t, frac_bound, current = model.compute_binding_curve()
    fig_kinetics = go.Figure()
    fig_kinetics.add_trace(go.Scatter(x=t, y=frac_bound, mode='lines', name='Fraction Bound'))
    fig_kinetics.add_trace(go.Scatter(x=t, y=current,      mode='lines', name='Sensor Current'))
    fig_kinetics.update_layout(
        title='Aptamer Binding Kinetics',
        xaxis_title='Time (s)',
        yaxis_title='Fraction Bound / Current',
        template=plotly_template,
        font_size=16
    )

    # 2. Hill-Langmuir Plot
    T_range, frac_bound_hl, I_hl = model.compute_hill_langmuir_curve()
    fig_hl = go.Figure()
    fig_hl.add_trace(go.Scatter(x=T_range, y=frac_bound_hl, mode='lines', name='Fraction Bound (HL)'))
    fig_hl.add_trace(go.Scatter(x=T_range, y=I_hl,          mode='lines', name='Sensor Current (HL)'))
    fig_hl.update_layout(
        title='Hill-Langmuir Isotherm',
        xaxis_title='Target Concentration',
        yaxis_title='Fraction Bound / Current',
        xaxis_type='log',
        template=plotly_template,
        font_size=16
    )

    # 3. Butler-Volmer Plot
    eta_range, j_values = model.compute_butler_volmer_curve()
    fig_bv = go.Figure()
    fig_bv.add_trace(go.Scatter(x=eta_range, y=j_values, mode='lines', name='Current Density'))
    fig_bv.update_layout(
        title='Butler-Volmer Equation',
        xaxis_title='Overpotential (V)',
        yaxis_title='Current Density (A/cm²)',
        template=plotly_template,
        font_size=16
    )

    # 4. Marcus Theory Plot
    deltaG_range_eV, k_et_vals = model.compute_marcus_curve()
    fig_marcus = go.Figure()
    fig_marcus.add_trace(go.Scatter(x=deltaG_range_eV, y=k_et_vals, mode='lines', name='Electron Transfer Rate'))
    fig_marcus.update_layout(
        title='Marcus Theory Electron Transfer',
        xaxis_title='ΔG (eV)',
        yaxis_title='Rate Constant (s⁻¹)',
        yaxis_type='log',
        yaxis=dict(
        exponentformat="power",  # e.g. 1e+12
        showexponent="all"
                    ),
        template=plotly_template,
        font_size=16
    )
    # fig_marcus.update_yaxes(type='log')

    # 5. SWV Frequency Response
    (freqs, I_u, I_b, I_u_norm, I_b_norm, crossover_freq) = model.compute_swv_curve()
    fig_swv = go.Figure()
    fig_swv.add_trace(go.Scatter(x=freqs, y=I_u_norm, mode='lines', name='Unbound (Normalised)', line=dict(color='blue')))
    fig_swv.add_trace(go.Scatter(x=freqs, y=I_b_norm, mode='lines', name='Bound (Normalised)',   line=dict(color='red')))
    fig_swv.add_trace(go.Scatter(x=[crossover_freq, crossover_freq], y=[0, 1],
                                 mode='lines',
                                 name=f'Crossover = {crossover_freq:.1f} Hz',
                                 line=dict(color='black', dash='dash')))
    fig_swv.update_layout(
        title='SWV Frequency Response',
        xaxis_title='Frequency (Hz)',
        yaxis_title='Normalised SWV Current',
        template=plotly_template,
        font_size=16
    )

    crossover_text = f"Crossover Frequency: {crossover_freq:.1f} Hz"

    return fig_kinetics, fig_hl, fig_bv, fig_marcus, fig_swv, crossover_text

###############################################################################
# Run the App
###############################################################################
if __name__ == '__main__':
    app.run_server(debug=True, host='127.0.0.1', port=8050)
