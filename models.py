# models.py

import numpy as np
from config import Config


class SWVModelEvaluator:
    """
    A class to evaluate Square Wave Voltammetry (SWV) models based on 
    different theoretical frameworks: empirical double exponential, 
    Butler-Volmer, and Marcus theory.

    Attributes:
        config (Config): Configuration object containing physical constants.
    """

    def __init__(self, config: Config = Config):
        self.config = config

    def double_exponential(
                        self, 
                        f: np.ndarray,
                        I0: float, k_et: float,
                        f_decay1: float,
                        f_decay2: float
                        ) -> np.ndarray:
        """
        Evaluates a double exponential SWV model.

        Args:
            f (np.ndarray): Frequencies.
            I0 (float): Baseline amplitude.
            k_et (float): Electron transfer rate.
            f_decay1 (float): Characteristic decay frequency 1.
            f_decay2 (float): Characteristic decay frequency 2.

        Returns:
            np.ndarray: Modelled SWV current.
        """
        try:
            return (
                I0
                * np.exp(-k_et / (2.0 * f))
                * np.exp(-f / f_decay1)
                * np.exp(-f / f_decay2)
            )
        except Exception as e:
            print(f"Error in double_exponential: {e}")
            return np.zeros_like(f)

    def butler_volmer(
                    self,
                    f: np.ndarray,
                    I0: float,
                    j0: float,
                    alpha: float,
                    eta0: float,
                    k_freq: float,
                    f_decay1: float,
                    f_decay2: float
                    ) -> np.ndarray:
        """
        Evaluates a Butler-Volmer-based SWV model with frequency-dependent effective overpotential.

        Args:
            f (np.ndarray): Frequencies.
            I0 (float): Baseline amplitude.
            j0 (float): Exchange current density.
            alpha (float): Charge transfer coefficient.
            eta0 (float): Max overpotential.
            k_freq (float): Frequency scaling constant for η_eff.
            f_decay1 (float): Characteristic decay frequency 1.
            f_decay2 (float): Characteristic decay frequency 2.

        Returns:
            np.ndarray: Modelled SWV current.
        """
        try:
            F = self.config.FARADAY_CONST
            R = self.config.UNI_GAS_CONST
            T = self.config.STRD_TEMP

            eta_eff = eta0 * (1.0 - np.exp(-k_freq / f))
            term_fwd = np.exp((1 - alpha) * F * eta_eff / (R * T))
            term_rev = np.exp(-alpha * F * eta_eff / (R * T))
            j_bv = j0 * (term_fwd - term_rev)

            return (
                I0
                * j_bv
                * np.exp(-f / f_decay1)
                * np.exp(-f / f_decay2)
            )
        except Exception as e:
            print(f"Error in butler_volmer: {e}")
            return np.zeros_like(f)

    def marcus(
            self,
            f: np.ndarray,
            I0: float,
            A_marcus: float,
            lambda_: float,
            deltaG0: float,
            k_freq: float,
            f_decay1: float,
            f_decay2: float
            ) -> np.ndarray:
        """
        Evaluates a Marcus theory-based SWV model with frequency modulation.

        Args:
            f (np.ndarray): Frequencies.
            I0 (float): Baseline amplitude.
            A_marcus (float): Pre-exponential factor.
            lambda_ (float): Reorganisation energy (λ).
            deltaG0 (float): Gibbs free energy change (ΔG⁰).
            k_freq (float): Frequency scaling constant.
            f_decay1 (float): Characteristic decay frequency 1.
            f_decay2 (float): Characteristic decay frequency 2.

        Returns:
            np.ndarray: Modelled SWV current.
        """
        try:
            k_B = self.config.BOLTZMANN_EVK_CONST
            T = self.config.STRD_TEMP

            exponent = -((lambda_ + deltaG0) ** 2) / (4.0 * lambda_ * k_B * T)
            k_et_base = A_marcus * np.exp(exponent)
            k_et = k_et_base * (1.0 - np.exp(-f / k_freq))

            return (
                I0
                * k_et
                * np.exp(-f / f_decay1)
                * np.exp(-f / f_decay2)
            )
        except Exception as e:
            print(f"Error in marcus: {e}")
            return np.zeros_like(f)

    def evaluate(
                self,
                f: np.ndarray,
                model_type: str,
                params: list[float]
                ) -> np.ndarray:
        """
        General evaluation dispatcher for any supported SWV model.

        Args:
            f (np.ndarray): Frequencies.
            model_type (str): Model type identifier ('exponential', 'butler-volmer', 'marcus').
            params (list): Parameters required by the selected model.

        Returns:
            np.ndarray: Modelled SWV current or charge response.
        """
        try:
            model_type = model_type.lower()
            if model_type == 'exponential':
                return self.double_exponential(f, *params)
            elif model_type == 'butler-volmer':
                return self.butler_volmer(f, *params)
            elif model_type == 'marcus':
                return self.marcus(f, *params)
            else:
                raise ValueError(f"Unsupported model type: '{model_type}'")
        except Exception as e:
            print(f"Error in evaluate(): {e}")
            return np.zeros_like(f)

