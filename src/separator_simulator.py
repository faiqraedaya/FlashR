import numpy as np
from PyQt5.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QPushButton, QComboBox, QTableWidget, QTableWidgetItem, QTabWidget, QGroupBox, QDoubleSpinBox, QMessageBox, QSplitter, QTextEdit, QHeaderView)
from PyQt5.QtCore import Qt

from core.flash_calculation import FlashCalculation
from core.plot_canvas import PlotCanvas

class SeparatorSimulator(QMainWindow):
    """Main application window"""
    
    def __init__(self):
        super().__init__()
        self.flash_calc = FlashCalculation()
        self.results = {}
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle("Gas-Liquid Separator Simulator")
        self.setGeometry(100, 100, 1600, 1200)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # Main layout
        main_layout = QHBoxLayout(central_widget)
        
        # Left panel for inputs
        left_panel = self.create_input_panel()
        
        # Right panel for results
        right_panel = self.create_results_panel()
        
        # Splitter
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setSizes([400, 1000])
        
        main_layout.addWidget(splitter)
        
    def create_input_panel(self):
        """Create input panel with all controls"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        
        # Operating conditions
        conditions_group = QGroupBox("Operating Conditions")
        conditions_layout = QGridLayout(conditions_group)
        
        # Temperature
        conditions_layout.addWidget(QLabel("Temperature:"), 0, 0)
        self.temp_input = QDoubleSpinBox()
        self.temp_input.setRange(-200, 2000)
        self.temp_input.setValue(25)
        self.temp_input.setDecimals(2)
        conditions_layout.addWidget(self.temp_input, 0, 1)
        self.temp_unit_combo = QComboBox()
        self.temp_unit_combo.addItems(["°C", "K", "°F"])
        self.temp_unit_combo.setCurrentIndex(0)
        conditions_layout.addWidget(self.temp_unit_combo, 0, 2)
        
        # Pressure
        conditions_layout.addWidget(QLabel("Pressure:"), 1, 0)
        self.pressure_input = QDoubleSpinBox()
        self.pressure_input.setRange(-1, 10000)
        self.pressure_input.setValue(1.013)
        self.pressure_input.setDecimals(3)
        conditions_layout.addWidget(self.pressure_input, 1, 1)
        self.pressure_unit_combo = QComboBox()
        self.pressure_unit_combo.addItems(["bar", "barg", "Pa", "kPa", "MPa", "atm"])
        self.pressure_unit_combo.setCurrentIndex(0)
        conditions_layout.addWidget(self.pressure_unit_combo, 1, 2)
        
        # Feed Rate
        conditions_layout.addWidget(QLabel("Feed Rate:"), 2, 0)
        self.feed_rate_input = QDoubleSpinBox()
        self.feed_rate_input.setRange(0.0001, 100000000)
        self.feed_rate_input.setValue(100)
        self.feed_rate_input.setDecimals(4)
        conditions_layout.addWidget(self.feed_rate_input, 2, 1)
        self.feed_rate_unit_combo = QComboBox()
        self.feed_rate_unit_combo.addItems(["kg/h", "g/s", "kg/s", "t/h", "kmol/h"])
        self.feed_rate_unit_combo.setCurrentIndex(0)
        conditions_layout.addWidget(self.feed_rate_unit_combo, 2, 2)
        
        layout.addWidget(conditions_group)
        
        # Component selection
        comp_group = QGroupBox("Component Selection")
        comp_layout = QVBoxLayout(comp_group)
        
        comp_select_layout = QHBoxLayout()
        # Dynamically get all CoolProp fluids
        import CoolProp.CoolProp as CP
        all_fluids = CP.FluidsList()
        self.component_combo = QComboBox()
        self.component_combo.addItems(sorted(all_fluids))
        comp_select_layout.addWidget(self.component_combo)
        
        add_comp_btn = QPushButton("Add Component")
        add_comp_btn.clicked.connect(self.add_component)
        comp_select_layout.addWidget(add_comp_btn)
        
        comp_layout.addLayout(comp_select_layout)
        
        # Component table with mass and mole fraction columns
        self.component_table = QTableWidget(0, 4)
        self.component_table.setHorizontalHeaderLabels([
            "Component", "Mole Fraction", "Mass Fraction", "Remove"
        ])
        header = self.component_table.horizontalHeader()
        if header is not None:
            try:
                header.setSectionResizeMode(QHeaderView.Stretch)
            except AttributeError:
                header.setStretchLastSection(True)
        comp_layout.addWidget(self.component_table)
        
        # Normalize buttons
        norm_layout = QHBoxLayout()
        self.norm_mole_btn = QPushButton("Normalize Mole Fractions")
        self.norm_mole_btn.clicked.connect(self.normalize_mole_fractions)
        norm_layout.addWidget(self.norm_mole_btn)
        self.norm_mass_btn = QPushButton("Normalize Mass Fractions")
        self.norm_mass_btn.clicked.connect(self.normalize_mass_fractions)
        norm_layout.addWidget(self.norm_mass_btn)
        comp_layout.addLayout(norm_layout)
        
        layout.addWidget(comp_group)
        
        # Calculation button
        calc_btn = QPushButton("Run Simulation")
        calc_btn.clicked.connect(self.run_simulation)
        calc_btn.setMinimumHeight(40)
        layout.addWidget(calc_btn)
        
        # Add default components
        self.add_default_components()
        
        layout.addStretch()
        return panel
        
    def create_results_panel(self):
        """Create results panel with tabs"""
        panel = QWidget()
        layout = QVBoxLayout(panel)
        
        self.results_tabs = QTabWidget()
        
        # Summary tab
        self.summary_tab = QTextEdit()
        self.summary_tab.setReadOnly(True)
        self.results_tabs.addTab(self.summary_tab, "Summary")
        
        # Composition plots
        self.comp_plot = PlotCanvas(width=10, height=6)
        self.results_tabs.addTab(self.comp_plot, "Compositions")
        
        # Phase envelope
        self.phase_plot = PlotCanvas(width=10, height=6)
        self.results_tabs.addTab(self.phase_plot, "Phase Envelope")
        
        # K-values plot
        self.k_plot = PlotCanvas(width=10, height=6)
        self.results_tabs.addTab(self.k_plot, "K-Values")
        
        # Property plots
        self.prop_plot = PlotCanvas(width=10, height=6)
        self.results_tabs.addTab(self.prop_plot, "Properties")
        
        layout.addWidget(self.results_tabs)
        return panel
        
    def add_default_components(self):
        """Add default component mixture"""
        default_components = [
            ('Methane', 0.85),
            ('Ethane', 0.08),
            ('Propane', 0.04),
            ('Water', 0.03)
        ]
        
        for comp, frac in default_components:
            self.add_component_to_table(comp, frac)
            
    def add_component(self):
        """Add selected component to mixture"""
        component = self.component_combo.currentText()
        self.add_component_to_table(component, 0.1)
        
    def add_component_to_table(self, component, mole_fraction):
        """Add component to the table, calculate mass fraction automatically"""
        import CoolProp.CoolProp as CP
        row = self.component_table.rowCount()
        self.component_table.insertRow(row)
        self.component_table.setItem(row, 0, QTableWidgetItem(component))
        frac_item = QTableWidgetItem(f"{mole_fraction:.4f}")
        self.component_table.setItem(row, 1, frac_item)
        # Calculate mass fraction
        total_moles = 0.0
        mass_fractions = []
        for r in range(self.component_table.rowCount()):
            comp_item = self.component_table.item(r, 0)
            mole_item = self.component_table.item(r, 1)
            if comp_item is not None and mole_item is not None:
                comp = comp_item.text()
                try:
                    mole = float(mole_item.text())
                except (ValueError, TypeError):
                    mole = 0.0
                M_i = CP.PropsSI('M', comp) / 1000
                total_moles += mole
                mass_fractions.append(mole * M_i)
            else:
                mass_fractions.append(0.0)
        total_mass = sum(mass_fractions) + 1e-12
        for r in range(self.component_table.rowCount()):
            comp_item = self.component_table.item(r, 0)
            mole_item = self.component_table.item(r, 1)
            if comp_item is not None and mole_item is not None:
                comp = comp_item.text()
                try:
                    mole = float(mole_item.text())
                except (ValueError, TypeError):
                    mole = 0.0
                M_i = CP.PropsSI('M', comp) / 1000
                mass_frac = (mole * M_i) / total_mass if total_mass > 0 else 0.0
                self.component_table.setItem(r, 2, QTableWidgetItem(f"{mass_frac:.4f}"))
            else:
                self.component_table.setItem(r, 2, QTableWidgetItem("0.0000"))
        remove_btn = QPushButton("Remove")
        remove_btn.clicked.connect(lambda: self.remove_component(row))
        self.component_table.setCellWidget(row, 3, remove_btn)
        
    def remove_component(self, row):
        """Remove component from table"""
        self.component_table.removeRow(row)
        
    def normalize_mole_fractions(self):
        """Normalize mole fractions and update mass fractions"""
        import CoolProp.CoolProp as CP
        n = self.component_table.rowCount()
        moles = []
        for r in range(n):
            mole_item = self.component_table.item(r, 1)
            try:
                moles.append(float(mole_item.text()) if mole_item is not None else 0.0)
            except (ValueError, TypeError):
                moles.append(0.0)
        total = sum(moles) + 1e-12
        moles = [x/total for x in moles]
        for r in range(n):
            self.component_table.setItem(r, 1, QTableWidgetItem(f"{moles[r]:.4f}"))
        # Update mass fractions
        mass_fractions = []
        for r in range(n):
            comp_item = self.component_table.item(r, 0)
            if comp_item is not None:
                comp = comp_item.text()
                M_i = CP.PropsSI('M', comp) / 1000
                mass_fractions.append(moles[r] * M_i)
            else:
                mass_fractions.append(0.0)
        total_mass = sum(mass_fractions) + 1e-12
        for r in range(n):
            comp_item = self.component_table.item(r, 0)
            if comp_item is not None:
                comp = comp_item.text()
                M_i = CP.PropsSI('M', comp) / 1000
                mass_frac = (moles[r] * M_i) / total_mass if total_mass > 0 else 0.0
                self.component_table.setItem(r, 2, QTableWidgetItem(f"{mass_frac:.4f}"))
            else:
                self.component_table.setItem(r, 2, QTableWidgetItem("0.0000"))

    def normalize_mass_fractions(self):
        """Normalize mass fractions and update mole fractions"""
        import CoolProp.CoolProp as CP
        n = self.component_table.rowCount()
        mass_fracs = []
        for r in range(n):
            mass_item = self.component_table.item(r, 2)
            try:
                mass_fracs.append(float(mass_item.text()) if mass_item is not None else 0.0)
            except (ValueError, TypeError):
                mass_fracs.append(0.0)
        total = sum(mass_fracs) + 1e-12
        mass_fracs = [x/total for x in mass_fracs]
        for r in range(n):
            self.component_table.setItem(r, 2, QTableWidgetItem(f"{mass_fracs[r]:.4f}"))
        # Update mole fractions
        moles = []
        for r in range(n):
            comp_item = self.component_table.item(r, 0)
            if comp_item is not None:
                comp = comp_item.text()
                M_i = CP.PropsSI('M', comp) / 1000
                moles.append(mass_fracs[r] / M_i if M_i > 0 else 0.0)
            else:
                moles.append(0.0)
        total_moles = sum(moles) + 1e-12
        moles = [x/total_moles for x in moles]
        for r in range(n):
            self.component_table.setItem(r, 1, QTableWidgetItem(f"{moles[r]:.4f}"))

    def get_mixture_data(self):
        """Extract mixture data from table (returns components and mole fractions)"""
        components = []
        mole_fractions = []
        for row in range(self.component_table.rowCount()):
            comp_item = self.component_table.item(row, 0)
            mole_item = self.component_table.item(row, 1)
            if comp_item is not None and mole_item is not None:
                comp = comp_item.text()
                try:
                    frac = float(mole_item.text())
                except (ValueError, TypeError):
                    frac = 0.0
                components.append(comp)
                mole_fractions.append(frac)
        total = sum(mole_fractions)
        if total > 0:
            mole_fractions = [f/total for f in mole_fractions]
        return components, np.array(mole_fractions)
        
    def run_simulation(self):
        """Run the flash separation simulation"""
        try:
            # Get input data
            components, z = self.get_mixture_data()
            if len(components) == 0:
                QMessageBox.warning(self, "Warning", "Please add components!")
                return
            # Temperature
            T = self.temp_input.value()
            temp_unit = self.temp_unit_combo.currentText()
            if temp_unit == "°C":
                T = T + 273.15
            elif temp_unit == "°F":
                T = (T - 32) * 5.0/9.0 + 273.15
            # else assume K
            # Pressure
            P = self.pressure_input.value()
            pressure_unit = self.pressure_unit_combo.currentText()
            if pressure_unit == "bar":
                P = P * 1e5
            elif pressure_unit == "barg":
                P = (P + 1) * 1e5  # gauge to absolute
            elif pressure_unit == "kPa":
                P = P * 1e3
            elif pressure_unit == "MPa":
                P = P * 1e6
            elif pressure_unit == "atm":
                P = P * 101325
            # else assume Pa
            # Feed rate
            F = self.feed_rate_input.value()
            feed_unit = self.feed_rate_unit_combo.currentText()
            import CoolProp.CoolProp as CP
            M = sum([CP.PropsSI('M', comp) / 1000 * z[i] for i, comp in enumerate(components)])
            if feed_unit == "kg/h":
                F = F / M
            elif feed_unit == "g/s":
                F = (F * 3600) / (M * 1000)
            elif feed_unit == "kg/s":
                F = (F * 3600) / M
            elif feed_unit == "t/h":
                F = (F * 1000) / M
            # else assume kmol/h
            # Update flash calculator
            self.flash_calc.components = components
            self.flash_calc.mole_fractions = list(z)
            # Perform flash calculation
            beta, x, y, K = self.flash_calc.flash_calculation(T, P, z)
            # Calculate stream properties
            L = F * (1 - beta)  # Liquid flow rate
            V = F * beta        # Vapor flow rate
            # Store results
            self.results = {
                'components': components,
                'feed_composition': z,
                'liquid_composition': x,
                'vapor_composition': y,
                'k_values': K,
                'vapor_fraction': beta,
                'liquid_fraction': 1 - beta,
                'feed_rate': F,
                'liquid_rate': L,
                'vapor_rate': V,
                'temperature': T,
                'pressure': P
            }
            # Update displays
            self.update_summary()
            self.plot_compositions()
            self.plot_phase_envelope()
            self.plot_k_values()
            self.plot_properties()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Simulation failed: {str(e)}")
            
    def update_summary(self):
        """Update summary text"""
        if not self.results:
            return
        summary = f"""
FLASH SEPARATION RESULTS
========================

Operating Conditions:
- Temperature: {self.results['temperature']:.2f} K ({self.results['temperature']-273.15:.2f} °C)
- Pressure: {self.results['pressure']:.0f} Pa ({self.results['pressure']/1000:.1f} kPa)

Feed Stream:
- Flow Rate: {self.results['feed_rate']:.1f} kmol/h
- Composition:
"""
        for i, comp in enumerate(self.results['components']):
            summary += f"  {comp}: {self.results['feed_composition'][i]:.4f}\n"
        summary += f"""
Separation Results:
- Vapor Fraction: {self.results['vapor_fraction']:.4f}
- Liquid Fraction: {self.results['liquid_fraction']:.4f}

Vapor Stream:
- Flow Rate: {self.results['vapor_rate']:.1f} kmol/h
- Composition:
"""
        for i, comp in enumerate(self.results['components']):
            summary += f"  {comp}: {self.results['vapor_composition'][i]:.4f}\n"
        summary += f"""
Liquid Stream:
- Flow Rate: {self.results['liquid_rate']:.1f} kmol/h
- Composition:
"""
        for i, comp in enumerate(self.results['components']):
            summary += f"  {comp}: {self.results['liquid_composition'][i]:.4f}\n"
        summary += f"""
K-Values:
"""
        for i, comp in enumerate(self.results['components']):
            summary += f"  {comp}: {self.results['k_values'][i]:.4f}\n"
        self.summary_tab.setText(summary)
        
    def plot_compositions(self):
        """Plot feed, liquid, and vapor compositions"""
        if not self.results:
            return
            
        self.comp_plot.clear_plot()
        
        fig = self.comp_plot.fig
        ax = fig.add_subplot(111)
        
        components = self.results['components']
        x_pos = np.arange(len(components))
        width = 0.25
        
        ax.bar(x_pos - width, self.results['feed_composition'], width, 
               label='Feed', alpha=0.8, color='blue')
        ax.bar(x_pos, self.results['liquid_composition'], width, 
               label='Liquid', alpha=0.8, color='red')
        ax.bar(x_pos + width, self.results['vapor_composition'], width, 
               label='Vapor', alpha=0.8, color='green')
        
        ax.set_xlabel('Components')
        ax.set_ylabel('Mole Fraction')
        ax.set_title('Stream Compositions')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(components, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        fig.tight_layout()
        self.comp_plot.draw()
        
    def plot_phase_envelope(self):
        """Plot phase envelope and operating point"""
        if not self.results:
            return
            
        self.phase_plot.clear_plot()
        
        fig = self.phase_plot.fig
        ax = fig.add_subplot(111)
        
        # Generate temperature range for phase envelope
        T_range = np.linspace(200, 500, 50)
        P_bubble = []
        P_dew = []
        
        z = self.results['feed_composition']
        
        for T in T_range:
            try:
                # Bubble point calculation
                K_vals = self.flash_calc.wilson_k_values(T, 101325)
                P_bub = 101325 / np.sum(z * K_vals)
                P_bubble.append(P_bub)
                
                # Dew point calculation  
                P_dew_calc = 101325 * np.sum(z / K_vals)
                P_dew.append(P_dew_calc)
                
            except:
                P_bubble.append(np.nan)
                P_dew.append(np.nan)
        
        # Plot phase envelope
        ax.plot(T_range, np.array(P_bubble)/1000, 'b-', label='Bubble Point', linewidth=2)
        ax.plot(T_range, np.array(P_dew)/1000, 'r-', label='Dew Point', linewidth=2)
        
        # Plot operating point
        ax.plot(self.results['temperature'], self.results['pressure']/1000, 
                'ko', markersize=10, label='Operating Point')
        
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('Pressure (kPa)')
        ax.set_title('Phase Envelope')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        fig.tight_layout()
        self.phase_plot.draw()
        
    def plot_k_values(self):
        """Plot K-values for all components"""
        if not self.results:
            return
            
        self.k_plot.clear_plot()
        
        fig = self.k_plot.fig
        ax = fig.add_subplot(111)
        
        components = self.results['components']
        k_values = self.results['k_values']
        
        bars = ax.bar(components, k_values, alpha=0.8, color='purple')
        
        # Add horizontal line at K=1
        ax.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='K = 1')
        
        ax.set_xlabel('Components')
        ax.set_ylabel('K-Value')
        ax.set_title('Equilibrium K-Values')
        ax.set_xticklabels(components, rotation=45, ha='right')
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, k_val in zip(bars, k_values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                   f'{k_val:.3f}', ha='center', va='bottom')
        
        fig.tight_layout()
        self.k_plot.draw()
        
    def plot_properties(self):
        """Plot various thermodynamic properties"""
        if not self.results:
            return
            
        self.prop_plot.clear_plot()
        
        fig = self.prop_plot.fig
        
        # Create subplots
        ax1 = fig.add_subplot(221)
        ax2 = fig.add_subplot(222)
        ax3 = fig.add_subplot(223)
        ax4 = fig.add_subplot(224)
        
        components = self.results['components']
        
        # Plot 1: Vapor pressures
        T = self.results['temperature']
        vapor_pressures = []
        for comp in components:
            try:
                Psat = self.flash_calc.antoine_equation(comp, T)
                vapor_pressures.append(Psat/1000)  # Convert to kPa
            except:
                vapor_pressures.append(0)
                
        ax1.bar(components, vapor_pressures, alpha=0.8, color='orange')
        ax1.set_ylabel('Vapor Pressure (kPa)')
        ax1.set_title('Component Vapor Pressures')
        ax1.tick_params(axis='x', rotation=45)
        
        # Plot 2: Relative volatility
        if len(self.results['k_values']) > 1:
            k_ref = self.results['k_values'][-1]  # Use last component as reference
            rel_vol = self.results['k_values'] / k_ref
            ax2.bar(components, rel_vol, alpha=0.8, color='cyan')
            ax2.set_ylabel('Relative Volatility')
            ax2.set_title('Relative Volatility (vs. last component)')
            ax2.tick_params(axis='x', rotation=45)
        
        # Plot 3: Flow rate distribution
        liquid_flows = self.results['liquid_rate'] * self.results['liquid_composition']
        vapor_flows = self.results['vapor_rate'] * self.results['vapor_composition']
        
        x_pos = np.arange(len(components))
        width = 0.35
        
        ax3.bar(x_pos - width/2, liquid_flows, width, label='Liquid', alpha=0.8, color='blue')
        ax3.bar(x_pos + width/2, vapor_flows, width, label='Vapor', alpha=0.8, color='red')
        ax3.set_ylabel('Flow Rate (kmol/h)')
        ax3.set_title('Component Flow Rates')
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(components, rotation=45)
        ax3.legend()
        
        # Plot 4: Recovery factors
        recovery_liquid = (liquid_flows / (self.results['feed_rate'] * self.results['feed_composition'])) * 100
        recovery_vapor = (vapor_flows / (self.results['feed_rate'] * self.results['feed_composition'])) * 100
        
        ax4.bar(x_pos - width/2, recovery_liquid, width, label='Liquid Recovery', alpha=0.8, color='green')
        ax4.bar(x_pos + width/2, recovery_vapor, width, label='Vapor Recovery', alpha=0.8, color='yellow')
        ax4.set_ylabel('Recovery (%)')
        ax4.set_title('Component Recovery')
        ax4.set_xticks(x_pos)
        ax4.set_xticklabels(components, rotation=45)
        ax4.legend()
        
        fig.tight_layout()
        self.prop_plot.draw() 