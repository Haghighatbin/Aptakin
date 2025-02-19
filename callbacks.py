import dash
from dash.dependencies import Input, Output, State
from dash_bootstrap_templates import ThemeSwitchAIO
from dash import callback_context
import plotly.graph_objects as go
from config import Config

from app import app

from model import AptamerBindingModel
model = AptamerBindingModel()


# Synchronising Slider & Input Callbacks

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
    fig_kin.add_trace(go.Scatter(x=t, y=frac_bound, mode='lines', name='Fraction Bound'))
    fig_kin.add_trace(go.Scatter(x=t, y=I_kin, mode='lines', name='Current (I)'))
    fig_kin.update_layout(
        title=dict(text='Binding Kinetics', font=dict(size=Config.PLOT_TITLE_FONTSIZE)), 
        xaxis=dict(
            title=dict(text='Time (s)', font=dict(size=Config.PLOT_XAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE)
        ),
        yaxis=dict(
            title=dict(text='Fraction Bound / Current', font=dict(size=Config.PLOT_YAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE)
        ),
        template=plotly_template
    )

    # -------------------------------------------------
    # 2) Hill–Langmuir Plot
    # -------------------------------------------------
    T_range, frac_HL, I_HL = model.compute_hill_langmuir_curve()
    fig_hl = go.Figure()
    fig_hl.add_trace(go.Scatter(x=T_range, y=frac_HL, mode='lines', name='Fraction Bound'))
    fig_hl.add_trace(go.Scatter(x=T_range, y=I_HL, mode='lines', name='Current (I)', yaxis='y2'))
    fig_hl.update_layout(
        title=dict(text='Hill-Langmuir Isotherm', font=dict(size=int(Config.PLOT_TITLE_FONTSIZE))),
        xaxis=dict(
            title=dict(text='[Target] (µM, log scale)', font=dict(size=int(Config.PLOT_XAXIS_FONTSIZE))),
            tickfont=dict(size=int(Config.PLOT_TICK_FONTSIZE)),
            type='log'
        ),
        yaxis=dict(
            title=dict(text='Fraction Bound', font=dict(size=int(Config.PLOT_YAXIS_FONTSIZE))),
            tickfont=dict(size=int(Config.PLOT_TICK_FONTSIZE)),
            side='left'
        ),
        yaxis2=dict(
            title=dict(text='Current (I)', font=dict(size=int(Config.PLOT_YAXIS_FONTSIZE))),
            tickfont=dict(size=int(Config.PLOT_TICK_FONTSIZE)),
            overlaying='y',
            side='right',
            showgrid=False
        ),
        template=plotly_template
    )

    # -------------------------------------------------
    # 3) Integrated Model Plot (Butler-Volmer or Marcus)
    # -------------------------------------------------
    T_int, fbound_int, I_int = model.compute_integrated_current_curve(method=integrated_method)
    fig_int = go.Figure()
    fig_int.add_trace(go.Scatter(x=T_int, y=fbound_int, mode='lines', name='Fraction Bound'))
    fig_int.add_trace(go.Scatter(x=T_int, y=I_int, mode='lines', name='Integrated Current', yaxis='y2'))
    fig_int.update_layout(
        title=dict(text=f'Integrated Model ({integrated_method})', font=dict(size=Config.PLOT_TITLE_FONTSIZE)),
        xaxis=dict(
            title=dict(text='[Target] (µM, log scale)', font=dict(size=Config.PLOT_XAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE)
        ),
        yaxis=dict(
            title=dict(text='Fraction Bound', font=dict(size=Config.PLOT_YAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE),
            side='left'
        ),
        yaxis2=dict(
            title=dict(text='Current (I)', font=dict(size=Config.PLOT_YAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE),
            overlaying='y',
            side='right',
            showgrid=False
        ),
        xaxis_type='log',
        template=plotly_template
    )

    # -------------------------------------------------
    # 4) Butler–Volmer Plot
    # -------------------------------------------------
    eta_range, j_values = model.compute_butler_volmer_curve()
    fig_bv = go.Figure()
    fig_bv.add_trace(go.Scatter(x=eta_range, y=j_values, mode='lines', name='j(η)'))
    fig_bv.update_layout(
        title=dict(text='Butler-Volmer', font=dict(size=Config.PLOT_TITLE_FONTSIZE)),
        xaxis=dict(
            title=dict(text='Overpotential (V)', font=dict(size=Config.PLOT_XAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE)
        ),
        yaxis=dict(
            title=dict(text='Current Density (A/cm²)', font=dict(size=Config.PLOT_YAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE)
        ),
        template=plotly_template
    )

    # -------------------------------------------------
    # 5) Marcus Plot
    # -------------------------------------------------
    dG_range, k_et_vals = model.compute_marcus_curve()
    fig_marcus = go.Figure()
    fig_marcus.add_trace(go.Scatter(x=dG_range, y=k_et_vals, mode='lines', name='k_et(ΔG)'))
    fig_marcus.update_layout(
        title=dict(text='Marcus Theory', font=dict(size=Config.PLOT_TITLE_FONTSIZE)),
        xaxis=dict(
            title=dict(text='ΔG (eV)', font=dict(size=Config.PLOT_XAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE)
        ),
        yaxis=dict(
            title=dict(text='ET Rate (s⁻¹)', font=dict(size=Config.PLOT_YAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE),
            type='log',
            exponentformat='power',
            showexponent='all'
        ),
        template=plotly_template
    )

    # -------------------------------------------------
    # 6) SWV Plot + crossover info
    # -------------------------------------------------
    inv_freqs, I_u_F, I_b_F, I_u_norm, I_b_norm, crossover_freq, crossover_inv_freq = model.compute_swv_curve()
    fig_swv = go.Figure()
    fig_swv.add_trace(go.Scatter(x=inv_freqs, y=I_u_norm, mode='lines', name='Unbound (Norm)'))
    fig_swv.add_trace(go.Scatter(x=inv_freqs, y=I_b_norm, mode='lines', name='Bound (Norm)'))
    fig_swv.add_trace(go.Scatter(
        x=[crossover_inv_freq, crossover_inv_freq], 
        y=[0, 1], 
        mode='lines',
        line=dict(color='black', dash='dash'),
        name=f'Crossover ~ {crossover_freq:.2f} Hz'
    ))
    fig_swv.update_layout(
        title=dict(text='SWV Frequency Response', font=dict(size=Config.PLOT_TITLE_FONTSIZE)),
        xaxis=dict(
            title=dict(text='1 / Frequency (s)', font=dict(size=Config.PLOT_XAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE)
        ),
        yaxis=dict(
            title=dict(text='(Peak Current)/freq (normalized)', font=dict(size=Config.PLOT_YAXIS_FONTSIZE)),
            tickfont=dict(size=Config.PLOT_TICK_FONTSIZE)
        ),
        template=plotly_template
    )
    crossover_text = f"Crossover Frequency: {crossover_freq:.2f} Hz"

    # Return all 7 outputs
    return fig_kin, fig_hl, fig_int, fig_bv, fig_marcus, fig_swv, crossover_text
