from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt
from valve_control import ValveControl
from PyQt5.QtCore import QTimer

class SolenoidWindow(QtWidgets.QWidget):
    def __init__(self, main_window):
        super().__init__()  # No parent - independent window
        self.windim_x, self.windim_y = 400, 300
        self.setWindowTitle("Valve Control Panel")
        self.setGeometry(1000, 100, self.windim_x, self.windim_y)  # Position to the right of main window
        
        # Store reference to main window to access its data
        self.main_window = main_window
        self.is_closing = False
        
        # Background
        bg_label = QtWidgets.QLabel(self)
        bg_pixmap = QtGui.QPixmap("Torch_Hot_Fire/solenoid_bg.png")
        scaled_pixmap = bg_pixmap.scaled(self.windim_x, self.windim_y, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        bg_label.setPixmap(scaled_pixmap)
        bg_label.setGeometry(0, 0, self.windim_x, self.windim_y)
        
        # Set border - matching main window style
        self.border_frame = QtWidgets.QFrame(self)
        self.border_frame.setFrameStyle(QtWidgets.QFrame.Box)
        self.border_frame.setLineWidth(5)
        self.border_frame.setGeometry(0, 0, self.windim_x, self.windim_y)
        self.border_frame.setStyleSheet("border: 5px solid #f44336;")  # Red border initially
        
        # Title label - positioned at top like in MainWindow
        title = QtWidgets.QLabel("Valve Control Panel", self)
        title.setGeometry(100, 10, 200, 30)
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("font-size: 14pt; font-weight: bold; color: white; background-color: rgba(0, 0, 0, 120);")
        
        # Connection Status Label - similar to MainWindow
        self.connection_status = QtWidgets.QLabel("LabJack Status: Disconnected", self)
        self.connection_status.setGeometry(80, self.windim_y - 40, 240, 25)
        self.connection_status.setAlignment(Qt.AlignCenter)
        self.connection_status.setStyleSheet("background-color: red; color: white; font-weight: bold;")
        
        # Create the valve controls directly on the window (like in MainWindow)
        self.create_valve_controls()
        
        # Timer to update connection status
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_connection_status)
        self.timer.start(1000)  # Update every second
    
    def create_valve_controls(self):
        # Clear existing devices in main window
        self.main_window._devices = []
        self.main_window.device_map = {}
        
        # Create new valve controls with absolute positioning (like in MainWindow)
        valve_configs = [
            ("SN-H2-01", "CIO0", 80, 80),
            ("SN-OX-01", "CIO1", 240, 80),
            ("SN-N2-01", "EIO7", 80, 170),
            ("Spark-Plug", "CIO3", 240, 170)
        ]
        
        for name, labjack_output, x_pos, y_pos in valve_configs:
            # Create the valve control - using absolute positioning like in MainWindow
            valve = ValveControl(name, labjack_output, x_pos, y_pos, parent=self)
            
            # Add a label above the valve
            label = QtWidgets.QLabel(name, self)
            label.setGeometry(x_pos, y_pos - 25, 80, 20)
            label.setAlignment(Qt.AlignCenter)
            label.setStyleSheet("font-weight: bold; color: white; background-color: rgba(0, 0, 0, 120);")
            
            # Add to main window's device collections
            self.main_window._devices.append(valve)
            self.main_window.device_map[name] = valve
    
    def update_connection_status(self):
        if self.main_window.labjack.connection_status:
            self.connection_status.setText("LabJack Status: Connected")
            self.connection_status.setStyleSheet("background-color: green; color: white; padding: 5px; font-weight: bold;")
            # Update border to match data logger state
            if hasattr(self.main_window, 'data_logger') and self.main_window.data_logger.high_speed_mode:
                self.border_frame.setStyleSheet("border: 5px solid #4CAF50;")  # Green border for high speed
            else:
                self.border_frame.setStyleSheet("border: 5px solid #f44336;")  # Red border for normal speed
        else:
            self.connection_status.setText("LabJack Status: Disconnected")
            self.connection_status.setStyleSheet("background-color: red; color: white; padding: 5px; font-weight: bold;")
    
    # Override closeEvent to close the main window when this window is closed
    def closeEvent(self, event):
        if not self.is_closing:
            self.is_closing = True
            self.main_window.close()
        event.accept()