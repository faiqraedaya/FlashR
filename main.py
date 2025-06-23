"""
Flash Drum Simulator

This program simulates a flash drum separator, which is a device used to separate a mixture of liquid and vapor phases.

The program allows the user to input the following parameters:
- Temperature
- Pressure
- Feed rate
- Feed composition
- Separator design 

The program then calculates the following:
- Vapor and liquid compositions
- Vapor and liquid flow rates
- Recovery factors
- K-values
- Property plots   

Author: Faiq Raedaya
Date: 23/06/2025
Version: 1.0.0

Changelog:
- 1.0.0: 
    - Initial release
- 1.1.0: 
    - Modularized main.py
    - Added more units
"""

import sys
from PyQt5.QtWidgets import QApplication
from core.separator_simulator import SeparatorSimulator

def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    window = SeparatorSimulator()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()