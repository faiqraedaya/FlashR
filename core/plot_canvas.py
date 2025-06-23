from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

class PlotCanvas(FigureCanvas):
    """Custom matplotlib canvas for PyQt5"""
    def __init__(self, parent=None, width=8, height=6, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)
        self.setParent(parent)
        
    def clear_plot(self):
        self.fig.clear()
        self.draw()

    # ... existing code from main.py (PlotCanvas class only) ... 