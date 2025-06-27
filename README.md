# FlashR

## Overview
FlashR is a simple desktop application for simulating flash drums for separating vapour/liquid mixtures.

## Features
- Interactive GUI built with PyQt5
- Supports user-defined mixtures using CoolProp
- User input for mixture definition, temperature, pressure, and feed rate
- Calculates vapor and liquid compositions, flow rates, recovery factors, K-values (equilibrium ratios)
- Visualizes composition data, property and phase plots

## Requirements
- Python 3.7+
- [CoolProp](http://www.coolprop.org/)
- [PyQt5](https://pypi.org/project/PyQt5/)
- [matplotlib](https://matplotlib.org/)
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)

## Installation
1. **Clone the repository:**
   ```bash
   git clone https://github.com/faiqraedaya/FlashR
   cd "FlashR"
   ```
2. **Install dependencies:**
   ```bash
   pip install CoolProp PyQt5 matplotlib numpy scipy 
   ```

## Usage
1. To launch the application, run:
  ```bash
  python main.py
  ```
2. Set operating conditions (temperature, pressure, feed rate).
3. Define a mixture by selecting components from CoolProp's database.
4. Specify feed composition (mole/mass fractions)
5. Run the simulation
6. View calculated phase splits, compositions, K-values, and property plots

## License
This project is provided under the MIT License.
