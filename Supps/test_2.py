import numpy as np
import dash
from dash import dcc, html, callback, Output, Input
import dash_bootstrap_components as dbc
from dash_bootstrap_templates import ThemeSwitchAIO
import plotly.graph_objects as go
from scipy.integrate import odeint

# Define themes for dark/light mode switching
THEMES = [dbc.themes.CYBORG, dbc.themes.SLATE]

# Aptamer Binding Model
class AptamerBindingModel:
    def __init__(self):
        self.A_total = 1.0
        self.T_conc = 10.0
        self.k_on = 1.0
        self.k_off = 0.5
        self.I_min = 0.1
        self.I_max = 1.0
        self.AT0 = 0.0
        
        self.n = 1.0
        
        # SWV Simulation parameters
        self.k_et_unbound = 10
        self.k_et_bound = 50
        self.f_decay = 100
        self.I0 = 1.0

    def binding_kinetics(self, AT, t):
        A_free = self.A_total - AT
        return self.k_on * A_free * self.T_conc - self.k_off * AT

    def compute_binding_curve(self):
        t = np.linspace(0, 1, 200)
        AT = odeint(self.binding_kinetics, self.AT0, t).flatten()
        fraction_bound = AT / self.A_total
        I = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return t, fraction_bound, I

    def hill_langmuir(self, T, Kd, n):
        return T**n / (Kd**n + T**n)

    def compute_hill_langmuir_curve(self):
        Kd = self.k_off / self.k_on
        T_conc_range = np.linspace(0.0001, 2, 2000)
        fraction_bound = self.hill_langmuir(T_conc_range, Kd, self.n)
        I_hl = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return T_conc_range, fraction_bound, I_hl
    
    def compute_swv_curve(self):
        frequencies = np.linspace(1, 1000, 200)
        I_net_unbound = self.I0 * np.exp(-self.k_et_unbound / (2 * frequencies)) * np.exp(-frequencies / self.f_decay)
        I_net_bound   = self.I0 * np.exp(-self.k_et_bound / (2 * frequencies)) * np.exp(-frequencies / self.f_decay)
        I_net_unbound_norm = I_net_unbound / np.max(I_net_unbound)
        I_net_bound_norm   = I_net_bound / np.max(I_net_bound)
        crossover_index = np.argmin(np.abs(I_net_unbound_norm - I_net_bound_norm))
        crossover_freq = frequencies[crossover_index]
        return frequencies, I_net_unbound, I_net_bound, I_net_unbound_norm, I_net_bound_norm, crossover_freq

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.COSMO])
model = AptamerBindingModel()

app.layout = dbc.Container([
    ThemeSwitchAIO(aio_id="theme", themes=THEMES),
    
    html.H1("Aptamer Sensor Simulation", style={'fontSize': '40px'}),
    
    html.Div([
        html.H2("Binding Kinetics Parameters", style={'fontSize': '30px'}),
        html.Label("k_on (Binding Rate Constant, 1/s)"),
        dcc.Slider(id='kon-slider', min=0.1, max=5.0, step=0.1, value=model.k_on,
                   marks={i: str(i) for i in range(1, 6)}),

        html.Label("k_off (Unbinding Rate Constant, 1/s)"),
        dcc.Slider(id='koff-slider', min=0.1, max=5.0, step=0.1, value=model.k_off,
                   marks={i: str(i) for i in range(1, 6)}),
    ]),
    dcc.Graph(id='kinetics-plot'),

    html.Div([
        html.H2("Hill-Langmuir Parameters", style={'fontSize': '30px'}),
        html.Label("Hill Coefficient (n)"),
        dcc.Slider(id='hill-n-slider', min=0.5, max=4.0, step=0.1, value=model.n,
                   marks={i: str(i) for i in range(1, 5)}),
    ]),
    dcc.Graph(id='hill-langmuir-plot'),

    html.Div([
        html.H2("SWV Simulation", style={'fontSize': '30px'}),
        html.Label("k_et Unbound (s⁻¹)"),
        dcc.Slider(id='ket-unbound-slider', min=1, max=100, step=1, value=model.k_et_unbound,
                   marks={1: "1", 50: "50", 100: "100"}),
        html.Label("k_et Bound (s⁻¹)"),
        dcc.Slider(id='ket-bound-slider', min=1, max=100, step=1, value=model.k_et_bound,
                   marks={1: "1", 50: "50", 100: "100"}),
    ]),
    html.Div(id='crossover-info', style={'fontSize': '18px', 'fontWeight': 'bold'}),
    dcc.Graph(id='swv-plot')

    ], fluid=True)

@app.callback(
    [Output('kinetics-plot', 'figure'),
     Output('hill-langmuir-plot', 'figure'),
     Output('swv-plot', 'figure'),
     Output('crossover-info', 'children')],
    [Input('kon-slider', 'value'),
     Input('koff-slider', 'value'),
     Input('hill-n-slider', 'value'),
     Input('ket-unbound-slider', 'value'),
     Input('ket-bound-slider', 'value')
     ]
)
def update_plots(k_on, k_off, hill_n, k_et_unbound, k_et_bound):
    model.k_on = k_on
    model.k_off = k_off
    model.n = hill_n
    model.k_et_unbound = k_et_unbound
    model.k_et_bound = k_et_bound


    t, fraction_bound, I = model.compute_binding_curve()
    kinetics_fig = go.Figure()
    kinetics_fig.add_trace(go.Scatter(x=t, y=fraction_bound, mode='lines', name='Fraction Bound'))
    kinetics_fig.update_layout(title='Aptamer Binding Kinetics')

    T_conc_range, fraction_bound_hl, I_hl = model.compute_hill_langmuir_curve()
    hill_fig = go.Figure()
    hill_fig.add_trace(go.Scatter(x=T_conc_range, y=fraction_bound_hl, mode='lines', name='Fraction Bound (HL)'))
    hill_fig.update_layout(title='Hill-Langmuir Isotherm')

    frequencies, _, _, I_net_unbound_norm, I_net_bound_norm, crossover_freq = model.compute_swv_curve()
    swv_fig = go.Figure()
    swv_fig.add_trace(go.Scatter(x=frequencies, y=I_net_unbound_norm, mode='lines', name='Unbound (Norm)'))
    swv_fig.add_trace(go.Scatter(x=frequencies, y=I_net_bound_norm, mode='lines', name='Bound (Norm)'))
    swv_fig.update_layout(title='SWV Frequency Response')

    crossover_text = f"Crossover Frequency: {crossover_freq:.1f} Hz"

    return kinetics_fig, hill_fig, swv_fig, crossover_text

if __name__ == '__main__':
    app.run_server(debug=True, host='127.0.0.1', port=8050)
