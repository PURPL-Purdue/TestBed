#Startup
# 1.	Open SN-N2-01
# 2.	Close SN-N2-01
# 3.	Immediately Open SN-O2-01
# 4.	Start the spark.
# 5.	Open SN-H2-01 fuel 
# 6.	Stop the spark

# Shutdown
# 1.	Open SN-N2-01
# 2.	Close SN-H2-01
# 3.	Close SN-O2-01
# 4.	Close SN-N2-01

from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt, QTimer
import sys
from labjack import ljm
from pyqtgraph import PlotWidget, plot

class ValveControl(QtWidgets.QPushButton):
    def __init__(self, name, labjack_output, x, y, parent=None):
        super().__init__(parent)
        self.name = name
        self.labjack_output = labjack_output
        self.valve_open = False  # Start with valve closed
        self.device_connected = False  # Track connection state
        self.handle = None  # LabJack handle
        self.setText("Valve Closed")  # Initial text
        self.update_button_style()
        self.clicked.connect(self.confirm_toggle_valve)
        self.move(x, y)
        self.adjustSize()  # Adjust button size to fit text exactly

    def set_labjack_connection(self, handle, connected):
        """Sets the LabJack handle and connection status."""
        self.handle = handle
        self.device_connected = connected
        self.setEnabled(connected)  # Enable or disable button based on connection status

    def confirm_toggle_valve(self):
        if not self.device_connected:
            # Do nothing if device is not connected
            return

        # Create a confirmation message box with dark mode styling
        msg = QtWidgets.QMessageBox(self)
        msg.setWindowTitle("Confirmation")
        msg.setText(f"Are you certain you want to {'close' if self.valve_open else 'open'} {self.name}?")
        
        # Apply dark mode styling specifically to the dialog background and text color only
        msg.setStyleSheet(
            "QMessageBox { background-color: #2e2e2e; }"  # Dark gray background
            "QLabel { color: white; }"  # White text for labels, no background
            "QPushButton { background-color: #4a4a4a; color: white; font-size: 10pt; }"  # Dark buttons with white text
        )

        # Add Yes and No buttons
        yes_button = msg.addButton("Yes", QtWidgets.QMessageBox.AcceptRole)
        no_button = msg.addButton("No", QtWidgets.QMessageBox.RejectRole)
        msg.setDefaultButton(no_button)
        msg.exec_()

        if msg.clickedButton() == yes_button:
            self.toggle_valve()

    def toggle_valve(self):
        if not self.device_connected:
            # Prevent toggling if device is not connected
            return

        self.valve_open = not self.valve_open
        self.update_button_style()
        self.update_labjack_output()

    def update_button_style(self):
        text, color = ("Valve Open", "green") if self.valve_open else ("Valve Closed", "red")
        self.setText(text)
        self.setStyleSheet(f"background-color: {color}; color: white; font-size: 12pt; font-weight: bold;")
        self.adjustSize()  # Adjust button size after text change

    def update_labjack_output(self):
        """Send the valve state to the LabJack output."""
        if self.device_connected and self.handle:
            # Set output to low (0) if valve is open, high (1) if closed
            output_value = 0 if self.valve_open else 1
            try:
                ljm.eWriteName(self.handle, self.labjack_output, output_value)
                print(f"{self.labjack_output} set to {output_value} for {'open' if self.valve_open else 'closed'} valve state.")
            except Exception as e:
                print(f"Error writing to {self.labjack_output}: {e}")

class PressureTransducer:
    def __init__(self, name, input_channel, x, y, parent):
        self.name = name
        self.input_channel = input_channel
        self.label = QtWidgets.QLabel("0.00 psi", parent)
        self.label.setGeometry(x, y, 100, 25)
        self.label.setStyleSheet("background-color: #00baff; color: white; font-size: 12pt; font-weight: bold;")
        self.label.setAlignment(Qt.AlignCenter)
        self.pressure = 0.0
        self.data = []  # Store pressure data for plotting

    def update_pressure(self, handle):
        try:
            voltage = ljm.eReadName(handle, self.input_channel)
            self.pressure = max(0.0, min(500.0, ((voltage - 0.5) / 4.0) * 500))
            self.label.setText(f"{self.pressure:.2f} psi")
            self.data.append(self.pressure)
        except Exception as e:
            print(f"Error reading pressure from {self.input_channel}: {e}")

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TeenyK P&ID")
        self.setGeometry(100, 100, 1320, 700)

        # Initialize connection state
        self.device_connected = False
        self.handle = None

        # Main layout
        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)

        # Background image setup
        bg_label = QtWidgets.QLabel(self)
        bg_pixmap = QtGui.QPixmap("TeenyK_BG.png").scaled(1000, 700, Qt.KeepAspectRatioByExpanding)
        bg_label.setPixmap(bg_pixmap)
        bg_label.setGeometry(0, 0, 1000, 700)

        # Connection Status Label
        self.connection_status = QtWidgets.QLabel("LabJack T7: Connection Missing", self)
        self.connection_status.setGeometry(740, 20, 240, 25)
        self.connection_status.setAlignment(Qt.AlignCenter)
        self.connection_status.setStyleSheet("background-color: red; color: white; font-weight: bold;")

        # Valve Control button using ValveControl class
        self.valve_button = ValveControl("SN-H2-01", "CIO3", 570, 438, parent=self)
        self.valve2 = ValveControl("SN-OX-01", "CIO2", 356, 438, parent=self)
        self.valve3 = ValveControl("SN-N2-01", "CIO3", 676, 200, parent=self)

        # Exit button
        self.exit_button = QtWidgets.QPushButton("Exit", self)
        self.exit_button.setStyleSheet("background-color: grey; color: black; font-size: 12pt;")
        self.exit_button.clicked.connect(self.close)
        self.exit_button.setGeometry(880, 55, 100, 40)

        # Pressure Transducers
        self.transducer1 = PressureTransducer("PT-O2-01", "AIN0", 303, 280, self)
        self.transducer2 = PressureTransducer("PT-H2-01", "AIN1", 584, 280, self)
        self.transducer3 = PressureTransducer("PT-TO-01", "AIN2", 457, 392, self)

        # Graphs for pressure readings
        self.graph1 = PlotWidget(self)
        self.graph1.setGeometry(1020, 0, 280, 230)
        self.graph1.setBackground('w')
        self.graph1.setTitle("PT1 Pressure")

        self.graph2 = PlotWidget(self)
        self.graph2.setGeometry(1020, 230, 280, 230)
        self.graph2.setBackground('w')
        self.graph2.setTitle("PT2 Pressure")

        self.graph3 = PlotWidget(self)
        self.graph3.setGeometry(1020, 460, 280, 230)
        self.graph3.setBackground('w')
        self.graph3.setTitle("PT3 Pressure")

        # Timers for updating pressure value and checking connection
        self.pressure_timer = QTimer(self)
        self.pressure_timer.timeout.connect(self.update_pressure)
        self.pressure_timer.start(1000)

        self.reconnect_timer = QTimer(self)
        self.reconnect_timer.timeout.connect(self.connect_to_labjack)
        self.reconnect_timer.start(5000)

        # Attempt initial connection to LabJack
        self.connect_to_labjack()

    def connect_to_labjack(self):
        if not self.device_connected:
            try:
                self.handle = ljm.openS("ANY", "ANY", "ANY")
                self.device_connected = True
                self.update_connection_status(True)
                print("Connection attempt: Success - Connection Established")
            except Exception as e:
                print("Connection attempt: Failed - Could not connect to LabJack:", e)
                self.device_connected = False
                self.update_connection_status(False)

    def update_connection_status(self, connected):
        if connected:
            self.connection_status.setText("LabJack T7: Connection Established")
            self.connection_status.setStyleSheet("background-color: green; color: white; font-weight: bold;")
            self.valve_button.set_labjack_connection(self.handle, True)
            self.valve2.set_labjack_connection(self.handle, True)
        else:
            self.connection_status.setText("LabJack T7: Connection Missing")
            self.connection_status.setStyleSheet("background-color: red; color: white; font-weight: bold;")
            self.valve_button.set_labjack_connection(None, False)
            self.valve2.set_labjack_connection(None, False)

    def update_pressure(self):
        if self.device_connected:
            self.transducer1.update_pressure(self.handle)
            self.transducer2.update_pressure(self.handle)
            self.transducer3.update_pressure(self.handle)

            # Update graphs
            self.graph1.plot(self.transducer1.data, clear=True)
            self.graph2.plot(self.transducer2.data, clear=True)
            self.graph3.plot(self.transducer3.data, clear=True)

    def closeEvent(self, event):
        if self.device_connected and self.handle:
            ljm.close(self.handle)
        event.accept()

# Set up the application
app = QtWidgets.QApplication(sys.argv)
main_window = MainWindow()
main_window.show()
sys.exit(app.exec_())
