from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from labjack import ljm

class PressureTransducer:
    def __init__(self, name, input_channel, x, y, parent):
        self.name = name
        self.input_channel = input_channel
        self.label = QtWidgets.QLabel("0.00 psi", parent)
        self.label.setGeometry(x, y, 100, 25)
        self.label.setStyleSheet("background-color: #00baff; color: white; font-size: 12pt; font-weight: bold;")
        self.label.setAlignment(Qt.AlignCenter)
        self.pressure = 0.0
        self.data = []

    def update_pressure(self, handle):
        try:
            voltage = ljm.eReadName(handle, self.input_channel)
            self.pressure = max(0.0, min(500.0, ((voltage - 0.5) / 4.0) * 500))
            self.label.setText(f"{self.pressure:.2f} psi")
            self.data.append(self.pressure)
        except Exception as e:
            print(f"Error reading pressure from {self.input_channel}: {e}")
