from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from labjack import ljm
from backend.AdaptiveData import Buffer

class PressureTransducer:
    def __init__(self, name, input_channel_1, max_voltage, max_psi, scalar_offset, linear_offset, x, y, parent):
        self.name = name
        self.input_channel_1 = input_channel_1
        self.max_voltage = max_voltage
        self.max_psi = max_psi
        self.label = QtWidgets.QLabel("0.0 psi", parent)
        self.label.setGeometry(x, y, 100, 25)
        self.label.setStyleSheet("background-color: #00baff; color: white; font-size: 12pt; font-weight: bold; font-family: 'Courier New', monospace;")
        self.label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.pressure = 0.0
        self.data = Buffer(5)
        self.redline = None
        self.linear_offset = linear_offset
        self.scalar_offset = scalar_offset

    def update_pressure(self, handle):
        try:
            voltage_1 = ljm.eReadName(handle, self.input_channel_1)
            # print(voltage_1)
            self.pressure = self.scalar_offset * voltage_1 / self.max_voltage * self.max_psi - self.linear_offset
            self.label.setText(f"{self.pressure:.1f} psi")
            self.data.append(self.pressure)
        except Exception as e:
            self.pressure = float('nan')
            self.data.append(self.pressure)
            raise Exception(f"pressure_transducer.py: Error reading pressure from {self.input_channel_1}: {e}")
