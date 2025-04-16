import sys
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt
from valve_control import ValveControl
from pressure_transducer import PressureTransducer
from labjack_connection import LabJackConnection
from transducer_data_logger import TransducerDataLogger
from sequencer import Sequencer
from pyqtgraph import PlotWidget
from PyQt5.QtCore import QTimer

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("TeenyK P&ID")
        self.setGeometry(100, 100, 1200, 900)
        
        self._transducers = []
        self._graphs = []
        self._devices = []

        # Background
        bg_label = QtWidgets.QLabel(self)
        bg_pixmap = QtGui.QPixmap("Torch_Hot_Fire/torch_bk_bg.png")
        scaled_pixmap = bg_pixmap.scaled(1200, 900, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        bg_label.setPixmap(scaled_pixmap)
        bg_label.setGeometry(0, 0, 1200, 900)

        # Connection Status Label
        self.connection_status = QtWidgets.QLabel("LabJack T7: Connection Missing", self)
        self.connection_status.setGeometry(740, 20, 240, 25)
        self.connection_status.setAlignment(Qt.AlignCenter)
        self.connection_status.setStyleSheet("background-color: red; color: white; font-weight: bold;")

        # LabJack Connection
        self.labjack = LabJackConnection(self.connection_status)
        self.labjack.connect_to_labjack()

        # Devices
        self._devices.append(ValveControl("SN-H2-01", "CIO0", 562, 425, parent=self))
        self._devices.append(ValveControl("SN-OX-01", "CIO1", 317, 416, parent=self))
        self._devices.append(ValveControl("SN-N2-01", "EIO7", 676, 170, parent=self))
        self._devices.append(ValveControl("Spark-Plug", "CIO3", 590, 530, parent=self))
        '''
        self._devices.append(ValveControl("SN-H2-01", "CIO0", 562, 425, parent=self))
        self._devices.append(ValveControl("SN-OX-01", "CIO1", 317, 416, parent=self))
        self._devices.append(ValveControl("SN-N2-01", "EIO7", 676, 170, parent=self))
        self._devices.append(ValveControl("Spark-Plug", "CIO3", 590, 530, parent=self))
        '''

        # Device Mapping
        self.device_map = {
            "SN-H2-01": self._devices[0],
            "SN-OX-01": self._devices[1],
            "SN-N2-01": self._devices[2],
            "Spark-Plug": self._devices[3],
        }

        # Create the sequencer with the events and devices
        self.sequencer = Sequencer(self.device_map, 445, 573, parent=self)

        # Pressure Transducers
        self._transducers.append(PressureTransducer("PT-O2-01", "AIN13", "", 332, 318, self))
        # self._transducers.append(PressureTransducer("PT-H2-01", "AIN98", "AIN99", 562, 318, self))
        # self._transducers.append(PressureTransducer("PT-TO-01", "AIN100", "AIN101", 457, 425, self))

        # TODO: Clean up whether or not we want graphs
        # # Graphs for pressure readings
        # self._graphs.append(PlotWidget(self))
        # self._graphs[0].setGeometry(1020, 0, 280, 230)
        # self._graphs[0].setBackground('w')
        # self._graphs[0].setTitle("AIN96 Pressure")

        # self._graphs.append(PlotWidget(self))
        # self._graphs[1].setGeometry(1020, 230, 280, 230)
        # self._graphs[1].setBackground('w')
        # self._graphs[1].setTitle("AIN98 Pressure")

        # self._graphs.append(PlotWidget(self))
        # self._graphs[2].setGeometry(1020, 460, 280, 230)
        # self._graphs[2].setBackground('w')
        # self._graphs[2].setTitle("AIN100 Pressure")

        # Data Logger
        self.data_logger = TransducerDataLogger(self._transducers, self._devices)

        # Timers for updating pressure value and checking connection
        self.pressure_timer = QTimer(self)
        self.pressure_timer.timeout.connect(self.update_pressure)
        self.pressure_timer.start(10)  # Changes Data Timing as Well

        self.reconnect_timer = QTimer(self)
        self.reconnect_timer.timeout.connect(self.labjack.connect_to_labjack)
        self.reconnect_timer.start(5000)

        QtWidgets.QApplication.instance().aboutToQuit.connect(self.perform_shutdown)

    def update_pressure(self):
        if self.labjack.connection_status:
            for i in range(len(self._transducers)):
                try:
                    # Update pressure
                    self._transducers[i].update_pressure(self.labjack.handle)

                    # # Update Graphs
                    # self._graphs[i].plot(self._transducers[i].data, clear=True)
                except Exception as e:
                    print(f"Error reading pressure from {self._trandsucers[i].input_channel_1}: {e}")

            self.data_logger.log_data()
    
    def perform_shutdown(self):
        print("Shutting down")
        # Turn off all devices
        for device in self._devices:
            try:
                # Ensure connection to LabJack
                if not device.device_connected:
                    device.connect_to_labjack()
                
                # Force the valve closed regardless of UI state
                if device.device_connected and device.handle:
                    # Update the UI state to match
                    device.valve_open = False
                    device.update_button_style()
                    device.toggle_valve_off()
                    
                    print(f"Closed valve: {device.name}")
                else:
                    print(f"WARNING: Could not close {device.name} - No connection")
            except Exception as e:
                print(f"ERROR closing {device.name}: {e}")
        self.data_logger.stop()
        self.labjack.close_connection()

    def closeEvent(self, event):
        self.perform_shutdown()
        event.accept()

# Main application
app = QtWidgets.QApplication(sys.argv)
main_window = MainWindow()
main_window.show()
sys.exit(app.exec_())