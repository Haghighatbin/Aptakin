# plot_utils.py

from typing import List, Dict
from app import app

def get_equation_image(model_name: str) -> List[Dict]:
    """
    Returns a list containing a Plotly-compatible dictionary for overlaying
    a mathematical equation image corresponding to the selected model.

    Args:
        model_name (str): The name of the model ('exponential', 'butler-volmer', 'marcus').

    Returns:
        List[Dict]: A list containing one image dictionary for Plotly layout,
                    or an empty list if the model name is unrecognised.
    """
    image_paths = {
        "exponential": "exponential.png",
        "butler-volmer": "butler-volmer.png",
        "marcus": "marcus.png"
    }

    path = image_paths.get(model_name)
    if not path:
        return []

    return [
        dict(
            source=app.get_asset_url(path),
            xref="paper",
            yref="paper",
            x=1.00,
            y=1.00,
            sizex=1.0,
            sizey=1.0,
            xanchor="right",
            yanchor="top",
            opacity=0.99,
            layer="above"
        )
    ]

