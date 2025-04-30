from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from labjack import ljm
from AdaptiveData import Buffer

class PressureTransducer:
    def __init__(self, name, input_channel_1, input_channel_2, x, y, parent):
        self.name = name
        self.input_channel_1 = input_channel_1
        self.input_channel_2 = input_channel_2
        self.label = QtWidgets.QLabel("0.00 psi", parent)
        self.label.setGeometry(x, y, 100, 25)
        self.label.setStyleSheet("background-color: #00baff; color: white; font-size: 12pt; font-weight: bold;")
        self.label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.pressure = 0.0
        self.data = Buffer(100)
        self.redline = None

    def update_pressure(self, handle):
        try:
            voltage_1 = ljm.eReadName(handle, self.input_channel_1)
            # print(voltage_1)
            self.pressure = voltage_1
            # voltage_2 = ljm.eReadName(handle, self.input_channel_2)
            # excitation_voltage = 12.0
            # full_scale_output = 3.00 - 0.015 * 1e-3 * excitation_voltage
            # pressure_full_scale = 10000  # 10,000 PSI max
            # self.pressure = ((voltage_1 - voltage_2) / full_scale_output) * pressure_full_scale
            self.label.setText(f"{self.pressure:.2f} psi")
            # self.pressure = 500 * ((abs((voltage_1/1000)-.5))/4)
            self.data.append(self.pressure)
        except Exception as e:
            self.pressure = float('nan')
            self.data.append(self.pressure)
            raise Exception(f"Error reading pressure from {self.input_channel_1}, {self.input_channel_2}: {e}")
