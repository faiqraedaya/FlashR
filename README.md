# FlashR
## Flash Drum Simulator

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

## Installation

### Prerequisites
- Python 3.7+
- [CoolProp](http://www.coolprop.org/)
- [PyQt5](https://pypi.org/project/PyQt5/)
- [matplotlib](https://matplotlib.org/)
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)

### Install dependencies
It is recommended to use a virtual environment:

```bash
python -m venv venv
# On Windows:
venv\Scripts\activate
# On Unix/Mac:
# source venv/bin/activate

pip install PyQt5 CoolProp matplotlib numpy scipy
```

## Usage

Run the main application:

```bash
python main.py
```

A window will open where you can:
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

## Author
Faiq Raedaya

## License
[Specify your license here, e.g., MIT, GPL, etc.]

## Acknowledgments
- [CoolProp](http://www.coolprop.org/) for property calculations
- [PyQt5](https://riverbankcomputing.com/software/pyqt/intro/) for the GUI framework
- [matplotlib](https://matplotlib.org/) for plotting

---

*For academic or educational use. Contributions and suggestions are welcome!* 
