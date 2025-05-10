from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt
from valve_control import ValveControl
from pressure_transducer import PressureTransducer
from labjack_connection import LabJackConnection
from data_logger import DataLogger
from sequencer import Sequencer
from PyQt5.QtCore import QTimer
import pyqtgraph as pg
import statistics

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        desktop = QtWidgets.QApplication.desktop()
        screen_rect = desktop.availableGeometry()
        screen_width = screen_rect.width()
        screen_height = screen_rect.height()
        # print(screen_height, screen_width)
        # margin = 30
        self.windim_x, self.windim_y = 1728, 973
        self.setWindowTitle("TeenyK P&ID")
        1003, 1728
        window_x = 0
        window_y = 0
        
        self.setGeometry(window_x, window_y, self.windim_x, self.windim_y)
        
        self._transducers = []
        self._graphs = []
        self._solenoids = []
        self.is_closing = False
        
        # # Store a reference to the solenoid window 
        # # (will be set by main.py after both windows are created)
        # self.valve_window = SolenoidWindow(self)
        # self.torch_window = TorchWindow(self)

        # Background
        bg_label = QtWidgets.QLabel(self)
        bg_pixmap = QtGui.QPixmap("GG_Test/P&ID.png")
        scaled_pixmap = bg_pixmap.scaled(bg_pixmap.width(), self.windim_y, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        bg_label.setPixmap(scaled_pixmap)
        bg_label.setGeometry(self.windim_x - scaled_pixmap.width(), 0, scaled_pixmap.width(), self.windim_y)

        # Get dimensions of the scaled pixmap
        self.scaled_width = scaled_pixmap.width()

        # Set border
        self.border_frame = QtWidgets.QFrame(self)
        self.border_frame.setFrameStyle(QtWidgets.QFrame.Box)
        self.border_frame.setLineWidth(5)
        self.border_frame.setGeometry(0, 0, self.windim_x, self.windim_y)

        # Connection Status Label
        self.connection_status = QtWidgets.QLabel("LabJack T7: Connection Missing", self)
        self.connection_status.setGeometry(self.windim_x - 260, 20, 240, 25)
        self.connection_status.setAlignment(Qt.AlignCenter)
        self.connection_status.setStyleSheet("background-color: red; color: white; font-weight: bold;")

        # LabJack Connection
        self.labjack = LabJackConnection(self.connection_status)
        self.labjack.connect_to_labjack()  # Initial connection attempt

        self.shutdown_button = QtWidgets.QPushButton("Emergency Shutdown", self)
        self.shutdown_button.setGeometry(10, self.windim_y-110, self.windim_x - self.scaled_width - 15, 100)  # Position below connection status
        self.shutdown_button.setStyleSheet("background-color: red; color: white; font-weight: bold;")
        self.shutdown_button.clicked.connect(self.perform_shutdown)

        # Hydrogen Valves
        self._solenoids.append(ValveControl("SN-H2-01", "CIO0", 749, 402, parent=self))
        # Oxygen valves
        self._solenoids.append(ValveControl("SN-O2-01", "CIO1", 638, 398, parent=self))
        self._solenoids.append(ValveControl("SN-O2-02", "", 591, 403, parent=self))
        # Nitrogen Valves
        self._solenoids.append(ValveControl("SN-N2-01", "CIO3", 694, 365, parent=self))
        self._solenoids.append(ValveControl("SN-N2-02", "", 525, 365, parent=self))
        self._solenoids.append(ValveControl("SN-N2-03", "", 398, 419, parent=self))
        self._solenoids.append(ValveControl("SN-N2-04", "", 541, 68, horizontal=True, parent=self))
        self._solenoids.append(ValveControl("SN-N2-05", "", 802, 532, horizontal=True, parent=self))
        self._solenoids.append(ValveControl("SN-N2-06", "", 802, 627, horizontal=True, parent=self))
        # Fuel Valves
        self._solenoids.append(ValveControl("SN-FU-01", "", 472, 403, parent=self))
        # Spark plug
        self._solenoids.append(ValveControl("Spark Plug", "EIO4", 540, 594, parent=self))
        # # Label Spark Plug
        # self.label = QtWidgets.QLabel("Spark Plug", self)
        # self.label.setGeometry(445, 210, 100, 25)
        # self.label.setStyleSheet("background-color: #FFFFFF; color: black; font-size: 12pt")
        # self.label.setAlignment(Qt.AlignCenter)

        # Pressure Transducers
        self._transducers.append(PressureTransducer("PT-TI-01", "AIN90", 10, 1500, 0, 528, 553, self))
        self._transducers.append(PressureTransducer("PT-GG-01", "", 10, 1500, 0, 528, 627, self))
        self._transducers.append(PressureTransducer("PT-H2-01", "", 10, 1500, 0, 1383, 410, self))
        self._transducers.append(PressureTransducer("PT-H2-02", "", 10, 1500, 0, 1123, 480, self))
        self._transducers.append(PressureTransducer("PT-H2-03", "AIN91", 10, 1500, 0, 634, 482, self))
        self._transducers.append(PressureTransducer("PT-O2-01", "", 5, 7500, 0, 1300, 281, self))
        self._transducers.append(PressureTransducer("PT-O2-02", "", 5, 10000, 0, 1317, 373, self))
        self._transducers.append(PressureTransducer("PT-O2-03", "", 5, 7500, 0, 1123, 374, self))
        self._transducers.append(PressureTransducer("PT-O2-04", "", 5, 7500, 0,  654, 515, self))
        self._transducers.append(PressureTransducer("PT-O2-05", "AIN89", 10, 1500, 0, 528, 465, self))
        self._transducers.append(PressureTransducer("PT-N2-01", "", 5, 10000, 0, 1254, 215, self))
        self._transducers.append(PressureTransducer("PT-N2-02", "", 5, 10000, 0, 1392, 228, self))
        self._transducers.append(PressureTransducer("PT-N2-03", "", 5, 10000, 0, 1391, 382, self))
        self._transducers.append(PressureTransducer("PT-N2-04", "", 5, 10000, 0, 1121, 218, self))
        self._transducers.append(PressureTransducer("PT-N2-05", "", 5, 10000, 0, 1121, 238, self))
        self._transducers.append(PressureTransducer("PT-N2-06", "", 5, 10000, 0, 1121, 259, self))
        self._transducers.append(PressureTransducer("PT-N2-07", "", 5, 10000, 0, 414, 76, self))
        self._transducers.append(PressureTransducer("PT-N2-08", "", 5, 10000, 0, 389, 454, self))
        self._transducers.append(PressureTransducer("PT-N2-09", "", 5, 10000, 0, 640, 541, self))
        # self._transducers.append(PressureTransducer("PT-FU-01", "", 5, 10000, 0, 640, 541, self))

        # Device Mapping
        self.device_map = {}

        for i in range(len(self._solenoids)):
            self.device_map[self._solenoids[i].name] = self._solenoids[i]
        for i in range(len(self._transducers)):
            self.device_map[self._transducers[i].name] = self._transducers[i]

        # Data Logger
        self.data_logger = DataLogger(self._transducers, self._solenoids, width = self.windim_x - self.scaled_width - 15, height = 100, parent=self)
        self.data_logger.move(10, 10) 

        # Connect the state_changed signal to update the main window border
        self.data_logger.state_changed.connect(self.update_border_color)
        # self.data_logger.state_changed.connect(self.valve_window.update_border_color)
        # self.data_logger.state_changed.connect(self.torch_window.update_border_color)
        
        # Set initial border style
        self.update_border_color(self.data_logger.high_speed_mode)

        # Graphs for pressure readings
        self._graph_number = 3
        self._graph_height = int((self.windim_y - 325 - (self._graph_number + 1) * 5)/3)
        # Torch graph
        self._graphs.append(pg.PlotWidget(self))
        self._graphs[0].setGeometry(10, 220, self.windim_x - self.scaled_width - 15, self._graph_height)
        self._graphs[0].setYRange(0, 200)
        self._graphs[0].setBackground('w')
        self._graphs[0].setTitle("PT-TI-01 Pressure")
        self._graphs[0].showGrid(y=True) 
        # GG graph
        self._graphs.append(pg.PlotWidget(self))
        self._graphs[1].setGeometry(10, 225 + self._graph_height, self.windim_x - self.scaled_width - 15, self._graph_height)
        self._graphs[1].setYRange(0, 200)
        self._graphs[1].setBackground('w')
        self._graphs[1].setTitle("PT-GG-01 Pressure")
        self._graphs[1].showGrid(y=True) 
        # OF graph
        self._graphs.append(pg.PlotWidget(self))
        self._graphs[2].setGeometry(10, 230 + 2 * self._graph_height, self.windim_x - self.scaled_width - 15, self._graph_height)
        self._graphs[2].setYRange(0, 200)
        self._graphs[2].setBackground('w')
        self._graphs[2].setTitle("OF Graph")
        self._graphs[2].showGrid(y=True) 


        # Timer for updating pressure value
        self.pressure_timer = QTimer(self)
        self.pressure_timer.timeout.connect(self.update_pressure)
        self.pressure_timer.start(500)  # Initial timing at 500ms

        # Give the data logger a reference to the timer
        self.data_logger.set_timer(self.pressure_timer)

        # Create the sequencer with the events and devices
        self.sequencer = Sequencer(self.device_map, self.data_logger, width = self.windim_x - self.scaled_width - 15, height = 100, parent=self)
        self.sequencer.move(10, 115)

        # self.valve_window.show()
        # self.torch_window.show()

    def update_pressure(self):
        """Update pressure readings from all transducers if LabJack is connected"""
        if not self.labjack.connection_status:
            print("Warning: Cannot update pressure - LabJack not connected")
            return
            
        for i in range(len(self._transducers)):
            try:
                # Update pressure
                self._transducers[i].update_pressure(self.labjack.handle)
                if self._transducers[i].redline is not None:
                    if statistics.median(self._transducers[i].data) > self._transducers[i].redline:
                        print(f"CRITICAL: {self._transducers[i].name} exceeded redline value! Initiating shutdown.")
                        print(self._transducers[i].data)
                        self.perform_shutdown()

                if self._transducers[i].name == "PT-TI-01":
                    if self.sequencer.running:
                        # Update Graphs
                        self.sequencer.pressure_data.append(self._transducers[i].pressure)
                        self._graphs[0].plot(self.sequencer.pressure_data, pen=pg.mkPen(color='b', width=3), clear=True)
                elif self._transducers[i].name == "PT-GG-01":
                    if self.sequencer.running:
                        # Update Graphs
                        self.sequencer.pressure_data.append(self._transducers[i].pressure)
                        self._graphs[1].plot(self.sequencer.pressure_data, pen=pg.mkPen(color='r', width=3), clear=True)
            except Exception as e:
                print(f"CRITICAL ERROR: Failed reading pressure from {self._transducers[i].name} ({self._transducers[i].input_channel_1}): {e}")
                # Log the error more prominently for rocket test infrastructure
        self.data_logger.log_data()

    def update_border_color(self, high_speed_mode):
        """Update the border color based on the button state"""
        if high_speed_mode:
            # Green border for high speed mode
            border_style = "border: 5px solid #4CAF50;"
        else:
            # Red border for normal speed mode
            border_style = "border: 5px solid #f44336;"
        
        # Update main window border
        self.border_frame.setStyleSheet(border_style)

        # # Update solenoid window border if it exists
        # if hasattr(self, 'valve_window') and self.valve_window:
        #     self.valve_window.border_frame.setStyleSheet(border_style)
        #     print(f"Updated solenoid window border: {'green' if high_speed_mode else 'red'}")
    
    def perform_shutdown(self):
        print("Shutting down")
        # Stop sequencer
        if self.sequencer.running:
            self.sequencer.stop_sequencer()
        # Turn off all devices
        for device in self._solenoids:
            try:
                # Force the valve closed regardless of UI state
                if self.labjack.connection_status and self.labjack.handle:
                    # Update the UI state to match
                    device.update_button_style()
                    device.toggle_valve_off()
                    
                    print(f"Closed valve: {device.name}")
                else:
                    print(f"WARNING: Could not close {device.name} - No connection")
            except Exception as e:
                print(f"ERROR closing {device.name}: {e}")

    def closeEvent(self, event):
        # Prevent multiple shutdown executions
        if not self.is_closing:
            self.is_closing = True
            
            # If we have a valve window and it's not already closing itself,
            # close it without triggering its closeEvent to close the main window again
            if hasattr(self, 'valve_window') and self.valve_window and not self.valve_window.is_closing:
                self.valve_window.is_closing = True  # Prevent it from closing main window again
                self.valve_window.close()
            
            # If we have a torch window and it's not already closing itself,
            # close it without triggering its closeEvent to close the main window again
            if hasattr(self, 'torch_window') and self.torch_window and not self.torch_window.is_closing:
                self.torch_window.is_closing = True  # Prevent it from closing main window again
                self.torch_window.close()

            # Perform shutdown tasks - only do this from the main window
            self.perform_shutdown()
            self.data_logger.stop()
            self.labjack.close_connection()
        event.accept()
