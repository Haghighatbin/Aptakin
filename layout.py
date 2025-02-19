import dash
from dash import dcc, html
import dash_bootstrap_components as dbc
from model import AptamerBindingModel
from dash_bootstrap_templates import ThemeSwitchAIO

from utils import slider_with_label, slider_input_with_label

THEMES = [dbc.themes.COSMO, dbc.themes.DARKLY]
app = dash.Dash(__name__, external_stylesheets=THEMES)
model = AptamerBindingModel()

main_layout = dbc.Container([
    html.H1("Electrochemical Aptamer-Based (E-AB) Sensor Simulation",
            style={'marginTop': '5px', 'fontSize': '50px', 
                   'marginBottom': '5px', 'padding': '5px'}),

    dbc.Row([
        dbc.Col(ThemeSwitchAIO(aio_id="theme", themes=THEMES), width=12)
    ]),

    # -------------------------------
    # 1. Binding Kinetics Layout
    # -------------------------------
    html.H2("Binding Kinetics", style={'color': 'green'}),
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
    html.H2("Integrated Model", style={'color': 'green'}),

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
    html.H2("Butler-Volmer vs. Marcus Theory", style={'color': 'green'}),

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

    html.H2("SWV Simulation", style={'color': 'green'}),
    dbc.Row([
        dbc.Col([dcc.Graph(id='swv-plot', style={'height': '100vh'})], width=9, 
                style={'flex': '0 0 80%', 'maxWidth': '80%'}),
        dbc.Col([
            html.H2("Experimental Parameters", style={'color': 'orange'}),
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
