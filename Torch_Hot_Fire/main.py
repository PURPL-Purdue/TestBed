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
        self.setGeometry(100, 100, 1320, 700)

        # Background
        bg_label = QtWidgets.QLabel(self)
        bg_pixmap = QtGui.QPixmap("TeenyK_BG.png").scaled(1000, 700)
        bg_label.setPixmap(bg_pixmap)
        bg_label.setGeometry(0, 0, 1000, 700)

        # Connection Status Label
        self.connection_status = QtWidgets.QLabel("LabJack T7: Connection Missing", self)
        self.connection_status.setGeometry(740, 20, 240, 25)
        self.connection_status.setAlignment(Qt.AlignCenter)
        self.connection_status.setStyleSheet("background-color: red; color: white; font-weight: bold;")

        # LabJack Connection
        self.labjack = LabJackConnection(self.connection_status)
        self.labjack.connect_to_labjack()

        # Devices
        self.device_one = ValveControl("SN-H2-01", "CIO0", 570, 438, parent=self)
        self.device_two = ValveControl("SN-OX-01", "CIO1", 356, 438, parent=self)
        self.device_three = ValveControl("SN-N2-01", "CIO2", 676, 200, parent=self)
        self.device_four = ValveControl("Spark-Plug", "CIO3", 100, 100, parent=self)

        # Device Mapping
        self.device_map = {
            "SN-H2-01": self.device_one,
            "SN-OX-01": self.device_two,
            "SN-N2-01": self.device_three,
            "Spark-Plug": self.device_four,
        }

        # Create the sequencer with the events and devices
        self.sequencer = Sequencer(self.device_map, 500, 500, parent=self)

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

        # Data Logger
        self.data_logger = TransducerDataLogger()

        # Timers for updating pressure value and checking connection
        self.pressure_timer = QTimer(self)
        self.pressure_timer.timeout.connect(self.update_pressure)
        self.pressure_timer.start(100)  # Changes Data Timing as Well

        self.reconnect_timer = QTimer(self)
        self.reconnect_timer.timeout.connect(self.labjack.connect_to_labjack)
        self.reconnect_timer.start(5000)

    def update_pressure(self):
        if self.labjack.connection_status:
            self.transducer1.update_pressure(self.labjack.handle)
            self.transducer2.update_pressure(self.labjack.handle)
            self.transducer3.update_pressure(self.labjack.handle)

            # Update Graphs
            self.graph1.plot(self.transducer1.data, clear=True)
            self.graph2.plot(self.transducer2.data, clear=True)
            self.graph3.plot(self.transducer3.data, clear=True)

            # Log Data
            self.data_logger.log_data(self.transducer1.pressure, self.transducer2.pressure, self.transducer3.pressure)


    def closeEvent(self, event):
        self.data_logger.export_to_csv()
        self.labjack.close_connection()
        print("Shutting down")
        event.accept()

# Main application
app = QtWidgets.QApplication(sys.argv)
main_window = MainWindow()
main_window.show()
sys.exit(app.exec_())
