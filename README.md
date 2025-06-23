# FlashR

## Overview
A desktop application for simulating a flash drum separator, a device used to separate a mixture of liquid and vapor phases. Built with PyQt5 for the GUI and CoolProp for thermodynamic property calculations.

## Features
- Input parameters: Temperature, Pressure, Feed rate, Feed composition, Separator design
- Calculates:
  - Vapor and liquid compositions
  - Vapor and liquid flow rates
  - Recovery factors
  - K-values (equilibrium ratios)
  - Property and phase plots
- Interactive GUI for parameter entry and results visualization
- Supports a wide range of chemical components via CoolProp
- Modular codebase for easy extension

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
2. A window will open where you can:
  - Set operating conditions (temperature, pressure, feed rate, units)
  - Select and add components from CoolProp's database
  - Specify feed composition (mole/mass fractions)
  - Run the simulation to view calculated phase splits, compositions, K-values, and property plots

## Project Structure
```
2513P Flash Drum Simulator/
├── main.py                  # Entry point, launches the GUI
├── core/
│   ├── separator_simulator.py  # Main PyQt5 GUI logic and user interface
│   ├── flash_calculation.py    # Core thermodynamic flash calculation logic
│   └── plot_canvas.py          # Matplotlib canvas integration for PyQt5
└── ...
```

## License
This project is provided under the MIT License.
