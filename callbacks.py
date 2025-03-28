# callbacks.py

import numpy as np
from scipy.optimize import curve_fit
from typing import Tuple, List, Dict, Optional, Any

from dash import Input, Output, State
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from dash_bootstrap_templates import ThemeSwitchAIO

from app import app
from config import Config
from models import SWVModelEvaluator
from plot_utils import get_equation_image

# Instantiate the SWV model evaluator class
model_evaluator = SWVModelEvaluator()

@app.callback(
    Output('swv-table', 'data'),
    Input('add-row-button', 'n_clicks'),
    State('swv-table', 'data'),
    prevent_initial_call=True
)
def add_row(
        n_clicks: int,
        rows: list[dict]
        ) -> list[dict]:
    """
    Adds a new empty row to the SWV data table upon button click.

    Args:
        n_clicks (int): Number of times the 'Add Row' button has been clicked.
        rows (list[dict]): Current list of rows in the table.

    Returns:
        list[dict]: Updated list of rows with an additional blank entry.
    """
    if rows is None:
        rows = []
    rows.append({"Frequency": 0, "Unbound": 0, "Bound": 0})
    return rows

@app.callback(
    [Output('results-table', 'style_cell'),
     Output('results-table', 'style_header'),
     Output('swv-table', 'style_cell'),
     Output('swv-table', 'style_header')],
    Input(ThemeSwitchAIO.ids.switch("theme"), "value")
)
def update_table_styles(
        current_theme: bool
        ) -> tuple[dict, dict, dict, dict]:
    """
    Updates the style of data tables based on the selected UI theme.

    Args:
        current_theme (bool): True if dark mode is enabled, False otherwise.

    Returns:
        tuple: Four dictionaries defining the style of cells and headers
               for both results and SWV data tables.
    """
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

    return style_cell, style_header, style_cell, style_header

#  ------------------  #
# | Helper Functions | #
#  ------------------  #

def _normalise_data(
    data: np.ndarray, 
    zero_protection_val: float = Config.PLT_NODEVZERO_VAL
) -> np.ndarray:
    """
    Normalize data by dividing by the maximum value with zero protection.

    Args:
        data (np.ndarray): Input data array to normalize.
        zero_protection_val (float, optional): Small value to prevent division by zero.

    Returns:
        np.ndarray: Normalized data array.
    """
    return data / (np.max(data) + zero_protection_val)

def _add_data_traces(
    fig: go.Figure, 
    freqs: np.ndarray, 
    unbound_vals: np.ndarray, 
    bound_vals: np.ndarray
) -> None:
    """
    Add data point traces to the Plotly figure.

    Args:
        fig (go.Figure): Plotly figure to add traces to.
        freqs (np.ndarray): Frequency values.
        unbound_vals (np.ndarray): Unbound measurement values.
        bound_vals (np.ndarray): Bound measurement values.
    """
    # Normalize data
    unbound_norm = _normalise_data(unbound_vals)
    bound_norm = _normalise_data(bound_vals)

    # Add data point traces
    fig.add_trace(go.Scatter(
        x=freqs, 
        y=unbound_norm, 
        mode='markers', 
        name='Unbound (Data)', 
        marker=dict(color='blue')
    ))
    fig.add_trace(go.Scatter(
        x=freqs, 
        y=bound_norm, 
        mode='markers', 
        name='Bound (Data)', 
        marker=dict(color='red')
    ))

def _add_fit_traces(
    fig: go.Figure, 
    model_type: str, 
    f_plot: np.ndarray, 
    unbound_params: Optional[np.ndarray], 
    bound_params: Optional[np.ndarray]
) -> None:
    """
    Add fit traces to the Plotly figure.

    Args:
        fig (go.Figure): Plotly figure to add traces to.
        model_type (str): Type of SWV model.
        f_plot (np.ndarray): Frequency array for plotting.
        unbound_params (Optional[np.ndarray]): Fitted parameters for unbound data.
        bound_params (Optional[np.ndarray]): Fitted parameters for bound data.
    """
    # Prepare extrapolation frequency array
    f_min, f_max = f_plot.min(), f_plot.max()
    step_size = (f_max - f_min) / (len(f_plot) - 1)
    num_new_points = int(round((Config.PLT_EXTP_MAXFREQ - Config.PLT_EXTP_MINFREQ) / step_size)) + 1
    f_extp = np.linspace(Config.PLT_EXTP_MINFREQ, Config.PLT_EXTP_MAXFREQ, num_new_points)

    # Add unbound fit traces if parameters exist
    if unbound_params is not None:
        # Model fit within original frequency range
        fit_u_plot = model_evaluator.evaluate(f_plot, model_type, unbound_params)
        fit_u_norm = _normalise_data(fit_u_plot)
        fig.add_trace(go.Scatter(
            x=f_plot, 
            y=fit_u_norm, 
            mode='lines', 
            name='Unbound (Fit)', 
            line=dict(color='blue')
        ))

        # Extrapolated fit trace
        extp_fit_u_plot = model_evaluator.evaluate(f_extp, model_type, unbound_params)
        extp_fit_u_norm = _normalise_data(extp_fit_u_plot)
        fig.add_trace(go.Scatter(
            x=f_extp, 
            y=extp_fit_u_norm, 
            mode='lines', 
            name='EXTP Unbound (Fit)', 
            line=dict(color='cyan', dash='dash')
        ))

    # Add bound fit traces if parameters exist
    if bound_params is not None:
        # Model fit within original frequency range
        fit_b_plot = model_evaluator.evaluate(f_plot, model_type, bound_params)
        fit_b_norm = _normalise_data(fit_b_plot)
        fig.add_trace(go.Scatter(
            x=f_plot, 
            y=fit_b_norm, 
            mode='lines', 
            name='Bound (Fit)', 
            line=dict(color='red')
        ))

        # Extrapolated fit trace
        extp_fit_b_plot = model_evaluator.evaluate(f_extp, model_type, bound_params)
        extp_fit_b_norm = _normalise_data(extp_fit_b_plot)
        fig.add_trace(go.Scatter(
            x=f_extp, 
            y=extp_fit_b_norm, 
            mode='lines', 
            name='EXTP Bound (Fit)', 
            line=dict(color='orange', dash='dash')
        ))

def _generate_results_table(
    model_type: str, 
    unbound_fit: Dict[str, Any], 
    bound_fit: Dict[str, Any]
) -> List[Dict[str, str]]:
    """
    Generate results table based on the model type and fit results.

    Args:
        model_type (str): Type of SWV model.
        unbound_fit (Dict[str, Any]): Fit results for unbound data.
        bound_fit (Dict[str, Any]): Fit results for bound data.

    Returns:
        List[Dict[str, str]]: Formatted results table.
    """
    # Retrieve parameters or use zeros if fit failed
    def _safe_get_params(fit_result, num_params):
        return fit_result['parameters'] if fit_result['parameters'] is not None else [0] * num_params

    # Model-specific parameter mappings
    model_param_mappings = {
        'exponential': {
            'param_names': ['I₀', 'kₑₜ', 'f_conf', 'f_others'],
            'num_params': 4
        },
        'butler-volmer': {
            'param_names': ['I₀', 'J₀', 'α', 'η₀', 'k_freq', 'f_conf', 'f_others'],
            'num_params': 7
        },
        'marcus': {
            'param_names': ['I₀', 'A', 'λ', 'ΔG⁰', 'k_freq', 'f_conf', 'f_others'],
            'num_params': 7
        }
    }

    # Get mapping for the current model
    mapping = model_param_mappings.get(model_type, {})
    param_names = mapping.get('param_names', [])
    num_params = mapping.get('num_params', 0)

    # Retrieve parameters
    unbound_params = _safe_get_params(unbound_fit, num_params)
    bound_params = _safe_get_params(bound_fit, num_params)

    # Generate results table
    results_table = []
    
    # Add parameter values
    for name, u_val, b_val in zip(param_names, unbound_params, bound_params):
        results_table.append({
            "Parameters": name,
            "Unbound": f"{u_val:.4g}",
            "Bound": f"{b_val:.4g}"
        })

    # Add statistical metrics
    metrics = [
        ("SSR", 'ss_res'),
        ("R²", 'r2'),
        ("σ (error STDV)", None),
        ("Reduced-χ²", 'reduced_chi2'),
        ("χ²", 'chi2')
    ]

    for metric_name, metric_key in metrics:
        if metric_key is None:
            # Special handling for fixed sigma values
            results_table.append({
                "Parameters": metric_name,
                "Unbound": f"{Config.CF_UNBOUND_PCT_SIGMA}",
                "Bound": f"{Config.CF_BOUND_PCT_SIGMA}"
            })
        else:
            results_table.append({
                "Parameters": metric_name,
                "Unbound": f"{unbound_fit.get(metric_key, 0):.6g}",
                "Bound": f"{bound_fit.get(metric_key, 0):.6g}"
            })

    return results_table

def parse_table_data(
        table_data: List[Dict]
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parse and validate input table data for curve fitting.

    Args:
        table_data (List[Dict]): Raw input data from the SWV table.

    Returns:
        Tuple of numpy arrays containing frequencies, unbound, and bound values.
    """
    freqs, unbound_vals, bound_vals = [], [], []
    
    for row in table_data:
        try:
            freq = float(row["Frequency"])
            unbound = float(row["Unbound"])
            bound = float(row["Bound"])
            
            # Validate data points
            if freq > 0 and unbound >= 0 and bound >= 0:
                freqs.append(freq)
                unbound_vals.append(unbound)
                bound_vals.append(bound)
        except (ValueError, KeyError) as err:
            print(f"Data parsing error: {err}")
            continue
    
    # Convert to numpy arrays
    return (
        np.array(freqs), 
        np.array(unbound_vals), 
        np.array(bound_vals)
    )

def get_initial_parameters(
    model_type: str, 
    unbound_vals: np.ndarray, 
    bound_vals: np.ndarray
) -> Tuple[List[float], List[float]]:
    """
    Generate initial parameter guesses based on the selected model type.

    Args:
        model_type (str): The type of SWV model ('exponential', 'butler-volmer', 'marcus').
        unbound_vals (np.ndarray): Unbound values for initial parameter estimation.
        bound_vals (np.ndarray): Bound values for initial parameter estimation.

    Returns:
        Tuple of initial parameter guesses for unbound and bound fits.
    """
    # Define parameter configurations for different models
    model_configs = {
        'exponential': {
            'unbound_multiplier': Config.CF_UNBOUND_I0_MUL,
            'unbound_params': [
                Config.UNBOUND_K_ET,
                Config.UNBOUND_F_DECAY1,
                Config.UNBOUND_F_DECAY2,
            ],
            'bound_multiplier': Config.CF_BOUND_I0_MUL,
            'bound_params': [
                Config.BOUND_K_ET,
                Config.BOUND_F_DECAY1,
                Config.BOUND_F_DECAY2,
            ]
        },
        'butler-volmer': {
            'unbound_params': [
                Config.BV_CURR_DENS_J0,
                Config.BV_ALPHA,
                Config.BV_ETA_0,
                Config.BV_UNBOUND_K_EFF,
                Config.BV_UNBOUND_F_DECAY1,
                Config.BV_UNBOUND_F_DECAY2,
            ],
            'bound_params': [
                Config.BV_CURR_DENS_J0,
                Config.BV_ALPHA,
                Config.BV_ETA_0,
                Config.BV_BOUND_K_EFF,
                Config.BV_BOUND_F_DECAY1,
                Config.BV_BOUND_F_DECAY2,
            ]
        },
        'marcus': {
            'unbound_params': [
                Config.MARCUS_PREEXP_A,
                Config.MARCUS_REORG_E_LAMBDA,
                Config.MARCUS_DELTAG,
                Config.MARCUS_UNBOUND_K_FREQ,
                Config.MARCUS_UNBOUND_F_DECAY1,
                Config.MARCUS_UNBOUND_F_DECAY2,
            ],
            'bound_params': [
                Config.MARCUS_PREEXP_A,
                Config.MARCUS_REORG_E_LAMBDA,
                Config.MARCUS_DELTAG,
                Config.MARCUS_BOUND_K_FREQ,
                Config.MARCUS_BOUND_F_DECAY1,
                Config.MARCUS_BOUND_F_DECAY2,
            ]
        }
    }

    # Retrieve configuration for the specific model
    config = model_configs.get(model_type, {})
    
    # Prepare initial parameter guesses
    p0_unbound = [
        np.max(unbound_vals) * config.get('unbound_multiplier', 1)
    ] + config.get('unbound_params', [])
    
    p0_bound = [
        np.max(bound_vals) * config.get('bound_multiplier', 1)
    ] + config.get('bound_params', [])
    
    return p0_unbound, p0_bound

def perform_curve_fit(
    model_type: str, 
    freqs: np.ndarray, 
    values: np.ndarray, 
    initial_params: List[float]
) -> Dict[str, float]:
    """
    Perform curve fitting for a given dataset.

    Args:
        model_type (str): The SWV model type.
        freqs (np.ndarray): Frequency data.
        values (np.ndarray): Measurement values.
        initial_params (List[float]): Initial parameter guesses.

    Returns:
        Dict containing fit results and statistical metrics.
    """
    def model_wrapper(f, *params):
        return model_evaluator.evaluate(f, model_type, params)

    try:
        # Perform curve fitting
        popt, _ = curve_fit(model_wrapper, freqs, values, p0=initial_params, maxfev=10000)
        
        # Calculate fit metrics
        fit_values = model_evaluator.evaluate(freqs, model_type, popt)
        residuals = values - fit_values
        
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((values - np.mean(values))**2)
        r2 = 1 - (ss_res / ss_tot)
        
        # Calculate chi-square
        sigma_array = Config.CF_UNBOUND_PCT_SIGMA * values
        dof = len(freqs) - len(popt)
        chi2 = np.sum((residuals / sigma_array)**2)
        reduced_chi2 = chi2 / dof
        
        return {
            'parameters': popt,
            'ss_res': ss_res,
            'r2': r2,
            'chi2': chi2,
            'reduced_chi2': reduced_chi2
        }
    
    except Exception as e:
        print(f"Curve fit failed: {e}")
        return {
            'parameters': None,
            'ss_res': 0,
            'r2': 0,
            'chi2': 0,
            'reduced_chi2': 0
        }

@app.callback(
    [Output('swv-graph', 'figure'),
     Output('results-table', 'data')],
    [Input('fit-button', 'n_clicks'),
     Input('swv-model-dropdown', 'value'),
     Input(ThemeSwitchAIO.ids.switch("theme"), "value")],
    State('swv-table', 'data')
)
def perform_fit(
    n_clicks: int, 
    model_type: str, 
    current_theme: bool, 
    table_data: List[Dict]
) -> Tuple[go.Figure, List[Dict]]:
    """
    Performs curve fitting using the selected model and updates the graph and results table.

    Args:
        n_clicks (int): Number of times the 'Fit' button has been clicked.
        model_type (str): The selected SWV model type.
        current_theme (bool): Indicates whether the dark theme is enabled.
        table_data (List[Dict]): User-input data from the SWV table.

    Returns:
        Tuple containing Plotly figure and results table data.
    """
    # Set template based on theme
    template = "plotly_dark" if current_theme else "plotly_white"
    
    # Return empty figure if no button click
    if not n_clicks:
        fig = go.Figure()
        fig.update_layout(template=template)
        return fig, []
    
    # Parse and validate input data
    freqs, unbound_vals, bound_vals = parse_table_data(table_data)
    
    # Validate data points
    if len(freqs) < 2:
        print("Not enough data points to perform fit.")
        return go.Figure(template=template), []
    
    # Get initial parameter guesses
    p0_unbound, p0_bound = get_initial_parameters(model_type, unbound_vals, bound_vals)
    
    # Perform curve fitting for unbound and bound data
    unbound_fit = perform_curve_fit(model_type, freqs, unbound_vals, p0_unbound)
    bound_fit = perform_curve_fit(model_type, freqs, bound_vals, p0_bound)
    
    # Prepare plotting data
    f_min, f_max = np.min(freqs), np.max(freqs)
    f_plot = np.linspace(f_min, f_max, Config.PLT_LNSPC_NUM)
    
    # Create figure and add data traces
    fig = go.Figure()
    _add_data_traces(fig, freqs, unbound_vals, bound_vals)
    
    # Add fit traces
    _add_fit_traces(
        fig, 
        model_type, 
        f_plot, 
        unbound_fit['parameters'], 
        bound_fit['parameters']
    )
    
    # Update figure layout
    fig.update_layout(
        title=f"{model_type.capitalize()} Model Fit",
        xaxis_title="Frequency / (Hz)",
        yaxis_title="Normalised Charge / (C)",
        template=template
    )
    
    # Add equation image
    fig.update_layout(images=get_equation_image(model_type))
    
    # Generate results table
    results_table = _generate_results_table(
        model_type, 
        unbound_fit, 
        bound_fit
    )
    
    return fig, results_table

