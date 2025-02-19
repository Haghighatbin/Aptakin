import numpy as np
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objects as go
from scipy.integrate import odeint

class AptamerBindingModel:
    def __init__(self):
        self.A_total = 1.0  # Total aptamer concentration (normalized)
        self.T_conc = 10.0  # Target concentration
        self.k_on = 1.0  # Binding rate constant
        self.k_off = 0.5  # Unbinding rate constant
        self.I_min = 0.1  # Minimum current
        self.I_max = 1.0  # Maximum current
        self.AT0 = 0.0  # Initial bound aptamer concentration

    def binding_kinetics(self, AT, t):
        A_free = self.A_total - AT
        return self.k_on * A_free * self.T_conc - self.k_off * AT

    def compute_binding_curve(self):
        t = np.linspace(0, 10, 200)
        AT = odeint(self.binding_kinetics, self.AT0, t).flatten()
        fraction_bound = AT / self.A_total
        I = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return t, fraction_bound, I

    def hill_langmuir(self, T, Kd, n):
        return T**n / (Kd**n + T**n)

    def compute_hill_langmuir_curve(self):
        Kd = self.k_off / self.k_on
        n = 1.0  # Hill coefficient
        T_conc_range = np.linspace(0, 20, 200)
        fraction_bound = self.hill_langmuir(T_conc_range, Kd, n)
        I_hl = self.I_min + (self.I_max - self.I_min) * fraction_bound
        return T_conc_range, fraction_bound, I_hl

app = dash.Dash(__name__)
model = AptamerBindingModel()

app.layout = html.Div([
    html.H1("Aptamer Binding Kinetics"),
    
    dcc.Slider(id='kon-slider', min=0.1, max=5.0, step=0.1, value=1.0, marks={i: str(i) for i in range(1, 6)}),
    dcc.Slider(id='koff-slider', min=0.1, max=5.0, step=0.1, value=0.5, marks={i: str(i) for i in range(1, 6)}),
    dcc.Graph(id='kinetics-plot'),
    dcc.Graph(id='hill-langmuir-plot')
])

@app.callback(
    [Output('kinetics-plot', 'figure'), Output('hill-langmuir-plot', 'figure')],
    [Input('kon-slider', 'value'), Input('koff-slider', 'value')]
)
def update_plots(k_on, k_off):
    model.k_on = k_on
    model.k_off = k_off
    
    t, fraction_bound, I = model.compute_binding_curve()
    T_conc_range, fraction_bound_hl, I_hl = model.compute_hill_langmuir_curve()
    
    kinetics_fig = go.Figure()
    kinetics_fig.add_trace(go.Scatter(x=t, y=fraction_bound, mode='lines', name='Fraction Bound'))
    kinetics_fig.add_trace(go.Scatter(x=t, y=I, mode='lines', name='Sensor Current'))
    kinetics_fig.update_layout(title='Aptamer Binding Kinetics', xaxis_title='Time (s)', yaxis_title='Fraction Bound / Current')
    
    hill_fig = go.Figure()
    hill_fig.add_trace(go.Scatter(x=T_conc_range, y=fraction_bound_hl, mode='lines', name='Fraction Bound (HL)'))
    hill_fig.add_trace(go.Scatter(x=T_conc_range, y=I_hl, mode='lines', name='Sensor Current (HL)'))
    hill_fig.update_layout(title='Hill-Langmuir Isotherm', xaxis_title='Target Concentration', yaxis_title='Fraction Bound / Current')
    
    return kinetics_fig, hill_fig

if __name__ == '__main__':

    app.run_server(debug=True, host='127.0.0.1', port=8050)

