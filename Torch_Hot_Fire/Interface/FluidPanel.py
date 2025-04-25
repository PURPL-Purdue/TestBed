from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt
from valve_control import ValveControl
from pressure_transducer import PressureTransducer
from labjack_connection import LabJackConnection
from Torch_Hot_Fire.data_logger import DataLogger
from sequencer import Sequencer
from PyQt5.QtCore import QTimer
from Interface.SolenoidPanel import SolenoidWindow

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.windim_x, self.windim_y = 900, 680
        self.setWindowTitle("TeenyK P&ID")
        desktop = QtWidgets.QApplication.desktop()
        screen_rect = desktop.screenGeometry()
        screen_width = screen_rect.width()
        
        window_x = screen_width - self.windim_x - 100
        window_y = 100 
        
        self.setGeometry(window_x, window_y, self.windim_x, self.windim_y)
        
        self._transducers = []
        self._graphs = []
        self._devices = []
        self.is_closing = False
        
        # Store a reference to the solenoid window 
        # (will be set by main.py after both windows are created)
        self.valve_window = SolenoidWindow(self)

        # Background
        bg_label = QtWidgets.QLabel(self)
        bg_pixmap = QtGui.QPixmap("Torch_Hot_Fire/torch_bk_bg.png")
        scaled_pixmap = bg_pixmap.scaled(bg_pixmap.width(), self.windim_y, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        bg_label.setPixmap(scaled_pixmap)
        bg_label.setGeometry(self.windim_x - scaled_pixmap.width(), 0, scaled_pixmap.width(), self.windim_y)

        # Set border
        self.border_frame = QtWidgets.QFrame(self)
        self.border_frame.setFrameStyle(QtWidgets.QFrame.Box)
        self.border_frame.setLineWidth(5)
        self.border_frame.setGeometry(0, 0, self.windim_x, self.windim_y)

        # Connection Status Label
        self.connection_status = QtWidgets.QLabel("LabJack T7: Connection Missing", self)
        self.connection_status.setGeometry(640, 20, 240, 25)
        self.connection_status.setAlignment(Qt.AlignCenter)
        self.connection_status.setStyleSheet("background-color: red; color: white; font-weight: bold;")

        # LabJack Connection
        self.labjack = LabJackConnection(self.connection_status)
        self.labjack.connect_to_labjack()

        self.shutdown_button = QtWidgets.QPushButton("Emergency Shutdown", self)
        self.shutdown_button.setGeometry(10, self.windim_y-110, 194, 100)  # Position below connection status
        self.shutdown_button.setStyleSheet("background-color: red; color: white; font-weight: bold;")
        self.shutdown_button.clicked.connect(self.perform_shutdown)

        # Create the valve controls in the main window but don't display them
        # These will serve as the "backend" for the solenoid panel
        self._devices.append(ValveControl("SN-H2-01", "CIO0", 465, 117, parent=self.valve_window))
        self._devices.append(ValveControl("SN-O2-01", "CIO1", 271, 117, parent=self.valve_window))
        self._devices.append(ValveControl("SN-N2-01", "EIO7", 369, 52, parent=self.valve_window))
        self.valve_window.show()

        # Device Mapping
        self.device_map = {}

        for i in range(len(self._devices)):
            self.device_map[self._devices[i].name] = self._devices[i]

        # Pressure Transducers
        self._transducers.append(PressureTransducer("PT-O2-01", "AIN0", "", 540, 332, self))
        self._transducers.append(PressureTransducer("PT-O2-03", "AIN96", "", 210, 463, self))
        self._transducers.append(PressureTransducer("PT-N2-01", "AIN0", "", 439, 188, self))
        self._transducers.append(PressureTransducer("PT-N2-04", "AIN72", "", 210, 231, self))
        self._transducers.append(PressureTransducer("PT-H2-01", "AIN0", "", 607, 530, self))
        self._transducers.append(PressureTransducer("PT-H2-02", "AIN72", "", 210, 623, self))

        # Data Logger
        self.data_logger = DataLogger(self._transducers, self._devices, parent=self)
        self.data_logger.move(10, 10) 

        # Connect the state_changed signal to update the main window border
        self.data_logger.state_changed.connect(self.update_border_color)
        
        # Set initial border style
        self.update_border_color(self.data_logger.high_speed_mode)

        # Timers for updating pressure value and checking connection
        self.pressure_timer = QTimer(self)
        self.pressure_timer.timeout.connect(self.update_pressure)
        self.pressure_timer.start(500)  # Initial timing at 500ms

        # Give the data logger a reference to the timer
        self.data_logger.set_timer(self.pressure_timer)

        # Create the sequencer with the events and devices
        self.sequencer = Sequencer(self.device_map, self.data_logger, parent=self)
        self.sequencer.move(10, 115)

        self.reconnect_timer = QTimer(self)
        self.reconnect_timer.timeout.connect(self.labjack.connect_to_labjack)
        self.reconnect_timer.start(5000)

        # QtWidgets.QApplication.instance().aboutToQuit.connect(self.perform_shutdown)

    def update_pressure(self):
        if self.labjack.connection_status:
            for i in range(len(self._transducers)):
                try:
                    # Update pressure
                    self._transducers[i].update_pressure(self.labjack.handle)
                except Exception as e:
                    print(f"Error reading pressure from {self._transducers[i].input_channel_1}: {e}")

            self.data_logger.log_data()

    def update_border_color(self, high_speed_mode):
        """Update the border color based on the button state"""
        if high_speed_mode:
            # Green border for high speed mode
            self.border_frame.setStyleSheet("border: 5px solid #4CAF50;")
        else:
            # Red border for normal speed mode
            self.border_frame.setStyleSheet("border: 5px solid #f44336;")
    
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

    def closeEvent(self, event):
        # Prevent multiple shutdown executions
        if not self.is_closing:
            self.is_closing = True
            
            # If we have a valve window and it's not already closing itself,
            # close it without triggering its closeEvent to close the main window again
            if hasattr(self, 'valve_window') and self.valve_window and not self.valve_window.is_closing:
                self.valve_window.is_closing = True  # Prevent it from closing main window again
                self.valve_window.close()
            
            # Perform shutdown tasks - only do this from the main window
            self.perform_shutdown()
            self.data_logger.stop()
            self.labjack.close_connection()
        event.accept()