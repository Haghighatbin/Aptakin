import numpy as np
import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
import plotly.graph_objects as go
from scipy.integrate import odeint
from dash_bootstrap_templates import ThemeSwitchAIO

# =============================================================================
# Model Class Definition
# =============================================================================
THEMES = [dbc.themes.COSMO, dbc.themes.DARKLY]

class AptamerBindingModel:
    def __init__(self):
        # Binding kinetics parameters
        self.A_total = 1.0                      # Total aptamer concentration (normalized)
        self.T_conc = 10.0                      # Target concentration
        self.k_on = 1.0                         # Binding rate constant (1/s)
        self.k_off = 0.5                        # Unbinding rate constant (1/s)
        self.I_min = 0.1                        # Minimum current
        self.I_max = 1.0                        # Maximum current
        self.AT0 = 0.0                          # Initial bound aptamer concentration

        # Hill–Langmuir parameter
        self.n = 1.0                            # Hill coefficient

        # SWV Simulation parameters
        self.k_et_unbound = 10                  # Electron-transfer rate constant for unbound state (s⁻¹)
        self.k_et_bound = 50                    # Electron-transfer rate constant for bound state (s⁻¹)
        self.f_decay = 100                      # Frequency decay constant (Hz)
        self.I0 = 1.0                           # Initial current amplitude (a.u.)

        # Butler–Volmer parameters
        self.j0 = 1e-3                          # Exchange current density (A/cm²)
        self.alpha = 0.5                        # Charge transfer coefficient
        self.F = 96485.3329                     # Faraday's constant (C/mol)
        self.R = 8.314                          # Universal gas constant (J/mol·K)
        self.T = 298.15                         # Temperature (K)

        # Marcus Theory parameters
        self.A_marcus = 1e12                    # Pre-exponential factor (s⁻¹)
        
        self.lambda_marcus = 0.2 * 1.60218e-19  # Reorganisation energy (lambda) in Joules; for example, 0.2 eV ≈ 0.2*1.60218e-19 J
        self.kB = 1.380649e-23                  # Boltzmann constant (J/K)
        
        self.deltaG0 = -0.1 * 1.60218e-19       # Standard free energy change (Delta G⁰) in Joules; for example, -0.1 eV
    # -------------------------------
    # Aptamer Binding Kinetics (Time Domain)
    # -------------------------------
    def binding_kinetics(self, AT, t):
        A_free = self.A_total - AT
        return self.k_on * A_free * self.T_conc - self.k_off * AT

    def compute_binding_curve(self):
        t = np.linspace(0, 4, 200)
        AT = odeint(self.binding_kinetics, self.AT0, t).flatten()
        fraction_bound = AT / self.A_total
        I = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return t, fraction_bound, I

    # -------------------------------
    # Hill–Langmuir Equilibrium Curve
    # -------------------------------
    def hill_langmuir(self, T, Kd, n):
        return T**n / (Kd**n + T**n)

    def compute_hill_langmuir_curve(self):
        Kd = self.k_off / self.k_on
        T_conc_range = np.linspace(0.01, 10, 1000)
        fraction_bound = self.hill_langmuir(T_conc_range, Kd, self.n)
        I_hl = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return T_conc_range, fraction_bound, I_hl

    # -------------------------------
    # SWV Frequency Response Simulation
    # -------------------------------
    def compute_swv_curve(self):
        frequencies = np.linspace(1, 200, 200)
        # SWV model: I_net(f) = I0 * exp[-(k_et/(2f))] * exp[-(f/f_decay)]
        I_net_unbound = self.I0 * np.exp(-self.k_et_unbound / (2 * frequencies)) * np.exp(-frequencies / self.f_decay)
        I_net_bound   = self.I0 * np.exp(-self.k_et_bound / (2 * frequencies)) * np.exp(-frequencies / self.f_decay)
        # Normalize each response by its maximum value
        I_net_unbound_norm = I_net_unbound / np.max(I_net_unbound)
        I_net_bound_norm   = I_net_bound / np.max(I_net_bound)
        # Find the crossover frequency where the two normalized curves are closest
        crossover_index = np.argmin(np.abs(I_net_unbound_norm - I_net_bound_norm))
        crossover_freq = frequencies[crossover_index]
        return frequencies, I_net_unbound, I_net_bound, I_net_unbound_norm, I_net_bound_norm, crossover_freq

    # -------------------------------
    # Butler–Volmer Theory Methods
    # -------------------------------
    def butler_volmer(self, eta):
        term_forward = np.exp((1 - self.alpha) * self.F * eta / (self.R * self.T))
        term_reverse = np.exp(-self.alpha * self.F * eta / (self.R * self.T))
        return self.j0 * (term_forward - term_reverse)

    def compute_butler_volmer_curve(self):
        eta_range = np.linspace(-0.5, 0.5, 200)  # Example range of overpotentials in Volts
        j = self.butler_volmer(eta_range)
        return eta_range, j
    
    # -------------------------------
    # Marcus Theory Methods
    # -------------------------------
    def marcus_rate(self, deltaG):
        exponent = -((self.lambda_marcus + deltaG)**2) / (4 * self.lambda_marcus * self.kB * self.T)
        return self.A_marcus * np.exp(exponent)

    def compute_marcus_curve(self, deltaG0):
        # Define a range for deltaG in Joules (e.g., from -0.5 eV to 0.5 eV)
        deltaG_range = np.linspace(deltaG0 - 0.5, deltaG0 + 0.5, 200)  # convert eV to Joules
        k_et = self.marcus_rate(deltaG_range)
        # Convert deltaG back to eV for plotting
        deltaG_range_eV = deltaG_range / 1.60218e-19
        return deltaG_range_eV, k_et
    
# =============================================================================
# Dash App Layout and Callbacks
# =============================================================================
app = dash.Dash(__name__, external_stylesheets=THEMES)
model = AptamerBindingModel()

app.layout = dbc.Container([
    html.H1("Electrochemical Aptamer-Based (E-AB) Sensor Simulation", style={'marginTop': '5px',
                                                                             'fontSize': '40px',
                                                                             'marginBottom': '5px',
                                                                             'padding': '5px'
                                                                             }),

    dbc.Row([
        dbc.Col(ThemeSwitchAIO(aio_id="theme", themes=THEMES), width=12)
        ]),
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
            html.Label("k_on (Binding Rate Constant, 1/s)", style={'fontSize': '20px'}),
            dcc.Slider(id='kon-slider', min=0.1, max=5.0, step=0.01, value=model.k_on,
                    marks={int(i): str(int(i)) for i in np.arange(1, 6)}),
            html.Label("k_off (Unbinding Rate Constant, 1/s)", style={'fontSize': '20px'}),
            dcc.Slider(id='koff-slider', min=0.1, max=5.0, step=0.01, value=model.k_off,
                    marks={int(i): str(int(i)) for i in np.arange(1, 6)}),
            html.Label("Target Concentration", style={'fontSize': '20px'}),
            dcc.Slider(id='Tconc-slider', min=1, max=50, step=0.01, value=model.T_conc,
                    marks={int(i): str(int(i)) for i in range(1, 51, 10)}),
            html.Label("I_min (Minimum Current)", style={'fontSize': '20px'}),
            dcc.Slider(id='Imin-slider', min=0.0, max=0.5, step=0.01, value=model.I_min,
                    marks={0.0: "0.0", 0.25: "0.25", 0.5: "0.5"}),
            html.Label("I_max (Maximum Current)", style={'fontSize': '20px'}),
            dcc.Slider(id='Imax-slider', min=0.5, max=2.0, step=0.01, value=model.I_max,
                    marks={0.5: "0.5", 1.0: "1.0", 1.5: "1.5", 2.0: "2.0"}),
                ], width=3, style={'flex': '0 0 25%', 'maxWidth': '25%'})
            ],
            style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}
            ),
    html.Hr(),
    dbc.Row([
        dbc.Col([
            dcc.Graph(id='hill-langmuir-plot')
            ], width=9, style={'flex': '0 0 75%', 'maxWidth': '75%'}),
        dbc.Col([
            html.H2("Hill-Langmuir Parameters", style={'color': 'orange'}),
            html.Label("Hill Coefficient (n)", style={'fontSize': '20px'}),
            dcc.Slider(id='hill-n-slider', min=0.5, max=4.0, step=0.01, value=model.n,
                    marks={0.5: "0.5", 1.0: "1.0", 1.5: "1.5", 2.0: "2.0", 2.5: "2.5", 3.0: "3.0", 3.5: "3.5", 4.0: "4.0"})
                ], width=3, style={'flex': '0 0 25%', 'maxWidth': '25%'})
            ],
            style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}
            ),
    html.Hr(),
    # Butler-Volmer Plot
    dbc.Row([
        dbc.Col([dcc.Graph(id='butler-volmer-plot')], width=9),
        dbc.Col([
            html.H2("Butler-Volmer Parameters", style={'color': 'orange'}),
            html.Label("Exchange Current Density (j0)"),
            dcc.Slider(id='j0-slider', min=1e-4, max=1e-2, step=1e-4, value=model.j0,
                       marks={1e-4: "1e-4", 5e-3: "5e-3", 1e-2: "1e-2"}),
            html.Label("Charge Transfer Coefficient (alpha)"),
            dcc.Slider(id='alpha-slider', min=0.1, max=1.0, step=0.01, value=model.alpha,
                       marks={0.1: "0.1", 0.5: "0.5", 1.0: "1.0"})
                ], width=3, style={'flex': '0 0 25%', 'maxWidth': '25%'})
            ],            
            style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}
            ),
    html.Hr(),
    # Marcus Theory Plot
    dbc.Row([
        dbc.Col([dcc.Graph(id='marcus-plot')], width=9),
        dbc.Col([
            html.H2("Marcus Theory Parameters", style={'color': 'orange'}),
            html.Label("Reorganization Energy (eV)"),
            dcc.Slider(id='lambda-slider', min=0.1, max=1.0, step=0.01, value=model.lambda_marcus / 1.60218e-19,
                       marks={0.1: "0.1", 0.5: "0.5", 1.0: "1.0"}),
            html.Label("Standard Free Energy Change (eV)"),
            dcc.Slider(id='deltaG-slider', min=-0.5, max=0.5, step=0.01, value=model.deltaG0 / 1.60218e-19,
                       marks={-0.5: "-0.5", 0.0: "0.0", 0.5: "0.5"})
                ], width=3, style={'flex': '0 0 25%', 'maxWidth': '25%'})
            ],
            style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}
            ),
    dbc.Row([
        dbc.Col([
            dcc.Graph(id='swv-plot')], width=9, style={'flex': '0 0 75%', 'maxWidth': '75%'}),
        dbc.Col([
            html.H2("SWV Simulation Parameters", style={'color': 'orange'}),
            html.Label("k_et Unbound (s⁻¹)", style={'fontSize': '20px'}),
            dcc.Slider(id='ket-unbound-slider', min=1, max=100, step=0.01, value=model.k_et_unbound,
                    marks={1: "1", 25: "25", 50: "50", 75: "75", 100: "100"}),
            html.Label("k_et Bound (s⁻¹)", style={'fontSize': '20px'}),
            dcc.Slider(id='ket-bound-slider', min=1, max=100, step=0.01, value=model.k_et_bound,
                    marks={1: "1", 25: "25", 50: "50", 75: "75", 100: "100"}),
            html.Label("f_decay (Hz)", style={'fontSize': '20px'}),
            dcc.Slider(id='fdecay-slider', min=50, max=300, step=1, value=model.f_decay,
                    marks={50: "50", 100: "100", 150: "150", 200: "200", 250: "250", 300: "300"}),
            html.Label("I0 (Initial Current)", style={'fontSize': '20px'}),
            dcc.Slider(id='I0-slider', min=0.1, max=2.0, step=0.01, value=model.I0,
                    marks={0.1: "0.1", 1.0: "1.0", 2.0: "2.0"}),
            html.Div(id='crossover-info', style={'fontSize': 22, 'fontWeight': 'bold', 'padding': '10px'})
                ], width=3, style={'flex': '0 0 25%', 'maxWidth': '25%'})
            ],
            style={'display': 'flex', 'flexDirection': 'row', 'width': '100%'}
            )
        ], fluid=True)

@app.callback(
    [Output('kinetics-plot', 'figure'),
     Output('hill-langmuir-plot', 'figure'),
     Output('swv-plot', 'figure'),
     Output('crossover-info', 'children'),
     Output('butler-volmer-plot', 'figure'),
     Output('marcus-plot', 'figure')],
    [Input('kon-slider', 'value'),
     Input('koff-slider', 'value'),
     Input('Tconc-slider', 'value'),
     Input('Imin-slider', 'value'),
     Input('Imax-slider', 'value'),
     Input('hill-n-slider', 'value'),
     Input('ket-unbound-slider', 'value'),
     Input('ket-bound-slider', 'value'),
     Input('fdecay-slider', 'value'),
     Input('I0-slider', 'value'),
     Input('j0-slider', 'value'),
     Input('alpha-slider', 'value'),
     Input('lambda-slider', 'value'),
     Input('deltaG-slider', 'value'),
     Input(ThemeSwitchAIO.ids.switch("theme"), "value")]
)
def update_plots(k_on, k_off, T_conc, I_min, I_max, hill_n, 
                 k_et_unbound, k_et_bound, f_decay, I0, j0, alpha, lambda_eV, deltaG_eV, theme_value):
    
    plotly_template = "plotly_white" if theme_value else "plotly_dark"

    # Update model parameters
    model.k_on = k_on
    model.k_off = k_off
    model.T_conc = T_conc
    model.I_min = I_min
    model.I_max = I_max
    model.n = hill_n
    model.k_et_unbound = k_et_unbound
    model.k_et_bound = k_et_bound
    model.f_decay = f_decay
    model.I0 = I0
    model.j0 = j0
    model.alpha = alpha
    model.lambda_marcus = lambda_eV * 1.60218e-19
    model.deltaG0 = deltaG_eV * 1.60218e-19

    # ---- Aptamer Binding Kinetics ----
    t, fraction_bound, I = model.compute_binding_curve()
    kinetics_fig = go.Figure()
    kinetics_fig.add_trace(go.Scatter(
                                x=t, y=fraction_bound, 
                                mode='lines', 
                                name='Fraction Bound')
                                )
    kinetics_fig.add_trace(go.Scatter(
                                x=t,
                                y=I,
                                mode='lines',
                                name='Sensor Current')
                                )
    kinetics_fig.update_layout(
                                title='Aptamer Binding Kinetics', 
                                xaxis_title='Time (s)', 
                                yaxis_title='Fraction Bound / Current',
                                template=plotly_template, 
                                font_size=16, 
                                margin=dict(l=30, r=20, t=50, b=20)
                                )

    # ---- Hill-Langmuir Isotherm ----
    T_conc_range, fraction_bound_hl, I_hl = model.compute_hill_langmuir_curve()
    hill_fig = go.Figure()
    hill_fig.add_trace(go.Scatter(
                                x=T_conc_range, 
                                y=fraction_bound_hl,
                                mode='lines', 
                                name='Fraction Bound (HL)')
                                )
    hill_fig.add_trace(go.Scatter(
                                x=T_conc_range, 
                                y=I_hl, 
                                mode='lines', 
                                name='Sensor Current (HL)')
                                )
    hill_fig.update_layout(
                                title='Hill-Langmuir Isotherm', 
                                xaxis_title='Target Concentration', 
                                yaxis_title='Fraction Bound / Current',
                                template=plotly_template, 
                                font_size=16,
                                margin=dict(l=30, r=20, t=50, b=20),
                                xaxis_type='log'
                           )

    # ---- Butler-Volmer Plot ----
    eta_range, j = model.compute_butler_volmer_curve()
    bv_fig = go.Figure()
    bv_fig.add_trace(go.Scatter(
                                x=eta_range,
                                y=j, 
                                mode='lines', 
                                name='Current Density')
                                )
    bv_fig.update_layout(
                                title='Butler-Volmer Equation',
                                xaxis_title='Overpotential (V)',
                                yaxis_title='Current Density (A/cm²)',
                                template=plotly_template,
                                font_size=16,
                                margin=dict(l=30, r=20, t=50, b=20),
                                )
    
    # ---- Marcus Theory Plot ----
    deltaG_range_eV, k_et = model.compute_marcus_curve(deltaG0=model.deltaG0)
    marcus_fig = go.Figure()
    marcus_fig.add_trace(go.Scatter(
                                x=deltaG_range_eV,
                                y=k_et, 
                                mode='lines', 
                                name='Electron Transfer Rate')
                                )
    marcus_fig.update_layout(
                                title='Marcus Theory Electron Transfer',
                                xaxis_title='ΔG (eV)', 
                                yaxis_title='Rate Constant (s⁻¹)',
                                template=plotly_template,
                                font_size=16,
                                margin=dict(l=30, r=20, t=50, b=20),
                                )

    # ---- SWV Frequency Response ----
    (frequencies, I_net_unbound, I_net_bound, 
     I_net_unbound_norm, I_net_bound_norm, 
     crossover_freq) = model.compute_swv_curve()
    
    swv_fig = go.Figure()
    swv_fig.add_trace(go.Scatter(
                                x=frequencies, 
                                y=I_net_unbound_norm,
                                mode='lines',
                                name='Unbound (Normalised)',
                                line=dict(color='blue'))
                                )
    swv_fig.add_trace(go.Scatter(
                                x=frequencies, 
                                y=I_net_bound_norm, 
                                mode='lines',
                                name='Bound (Normalised)', line=dict(color='red'))
                                )
    swv_fig.add_trace(go.Scatter(
                                x=[crossover_freq,
                                crossover_freq],
                                y=[0, 1],
                                mode='lines',
                                name=f'Crossover ({crossover_freq:.1f} Hz)',
                                line=dict(dash='dash', 
                                color='black'))
                                )
    swv_fig.update_layout(
                            title='SWV Frequency Response',
                            xaxis_title='Frequency (Hz)',
                            yaxis_title='Normalised SWV Current',
                            template=plotly_template, 
                            font_size=16,
                            margin=dict(l=30, r=20, t=50, b=20)
                          )

    # Create a text message with the crossover frequency information
    crossover_text = f"Crossover Frequency: {crossover_freq:.1f} Hz (at which unbound and bound responses are equal)"

    return kinetics_fig, hill_fig, swv_fig, crossover_text, bv_fig, marcus_fig

if __name__ == '__main__':
    app.run_server(debug=True, host='127.0.0.1', port=8050)
