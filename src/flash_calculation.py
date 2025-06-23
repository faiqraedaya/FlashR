import numpy as np
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve
import warnings
warnings.filterwarnings('ignore')

class FlashCalculation:
    """Core thermodynamic calculations for flash separation"""
    
    def __init__(self):
        self.components = []
        self.mole_fractions = []
        self.temperature = 298.15  # K
        self.pressure = 101325     # Pa
        
    def antoine_equation(self, component, T):
        """Calculate vapor pressure using Antoine equation"""
        antoine_params = {
            'Water': (8.07131, 1730.63, 233.426),
            'Methane': (6.61184, 389.93, 266.69),
            'Ethane': (6.80266, 663.70, 256.68),
            'Propane': (6.82973, 803.81, 246.99),
            'n-Butane': (6.83029, 945.90, 240.00),
            'Benzene': (6.90565, 1211.03, 220.79),
            'Toluene': (6.95464, 1344.80, 219.48),
            'CO2': (6.81228, 1301.679, 3.494),
            'H2S': (6.99052, 768.12, 273.16),
            'Nitrogen': (6.49457, 255.68, 266.55)
        }
        
        if component in antoine_params:
            A, B, C = antoine_params[component]
            return 10**(A - B/(T - 273.15 + C)) * 133.322  # Convert mmHg to Pa
        else:
            try:
                return CP.PropsSI('P', 'T', T, 'Q', 1, component)
            except:
                return 101325  # Default to atmospheric pressure

    def wilson_k_values(self, T, P):
        """Calculate K-values using Wilson correlation"""
        k_values = []
        for component in self.components:
            try:
                Psat = self.antoine_equation(component, T)
                Tc = CP.PropsSI('Tcrit', component)
                Pc = CP.PropsSI('Pcrit', component)
                w = CP.PropsSI('acentric', component)
                
                # Wilson correlation
                Tr = T / Tc
                Pr = P / Pc
                
                if Tr < 1.0:
                    ln_Psat_Pc = 5.37 * (1 + w) * (1 - 1/Tr)
                    K = (Pc/P) * np.exp(ln_Psat_Pc)
                else:
                    K = Psat / P
                    
                k_values.append(max(K, 1e-6))
            except:
                # Fallback calculation
                k_values.append(1.0)
                
        return np.array(k_values)

    def rachford_rice_equation(self, beta, z, K):
        """Rachford-Rice equation for flash calculation"""
        return np.sum(z * (K - 1) / (1 + beta * (K - 1)))

    def flash_calculation(self, T, P, z):
        """Perform flash calculation"""
        self.temperature = T
        self.pressure = P
        
        K = self.wilson_k_values(T, P)
        
        # Check if all liquid or all vapor
        if np.all(K <= 1.0):
            return 0.0, z, np.zeros_like(z), K  # All liquid
        elif np.all(K >= 1.0):
            return 1.0, np.zeros_like(z), z, K  # All vapor
        
        # Solve Rachford-Rice equation
        try:
            beta_min = 1 / (1 - np.max(K))
            beta_max = 1 / (1 - np.min(K))
            beta_guess = 0.5
            
            if beta_min < beta_max:
                beta = fsolve(self.rachford_rice_equation, beta_guess, 
                            args=(z, K))[0]
                beta = np.clip(beta, 0, 1)
            else:
                beta = 0.5
                
        except:
            beta = 0.5
        
        # Calculate liquid and vapor compositions
        x = z / (1 + beta * (K - 1))
        y = K * x
        
        # Normalize compositions
        x = x / np.sum(x) if np.sum(x) > 0 else x
        y = y / np.sum(y) if np.sum(y) > 0 else y
        
        return beta, x, y, K 