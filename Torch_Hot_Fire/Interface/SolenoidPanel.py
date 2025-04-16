from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt
from valve_control import ValveControl
from PyQt5.QtCore import QTimer

class SolenoidWindow(QtWidgets.QWidget):
    def __init__(self, main_window):
        super().__init__()  # No parent - independent window
        self.setWindowTitle("Valve Control Panel")
        self.setGeometry(1000, 100, 400, 300)  # Position to the right of main window
        
        # Store reference to main window to access its data
        self.main_window = main_window
        self.is_closing = False
        
        # Create layout
        layout = QtWidgets.QVBoxLayout(self)
        
        # Title label
        title = QtWidgets.QLabel("Valve Control Panel")
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("font-size: 14pt; font-weight: bold;")
        layout.addWidget(title)
        
        # Create a grid layout for valve controls
        grid = QtWidgets.QGridLayout()
        layout.addLayout(grid)
        
        # Create the valve controls in this window
        self.create_valve_controls(grid)
        
        # Add connection status label
        self.connection_status = QtWidgets.QLabel()
        layout.addWidget(self.connection_status)
        self.update_connection_status()
        
        # Timer to update connection status
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_connection_status)
        self.timer.start(1000)  # Update every second
    
    def create_valve_controls(self, grid):
        # Create valve controls and add them to the grid
        row, col = 0, 0
        
        # Clear existing devices in main window
        self.main_window._devices = []
        self.main_window.device_map = {}
        
        # Create new valve controls
        valve_configs = [
            ("SN-H2-01", "CIO0", row, col),
            ("SN-OX-01", "CIO1", row, col+1),
            ("SN-N2-01", "EIO7", row+1, col),
            ("Spark-Plug", "CIO3", row+1, col+1)
        ]
        
        for name, labjack_output, row_pos, col_pos in valve_configs:
            # Create a frame for each valve
            frame = QtWidgets.QFrame()
            frame.setFrameStyle(QtWidgets.QFrame.Box | QtWidgets.QFrame.Raised)
            frame_layout = QtWidgets.QVBoxLayout(frame)
            
            # Create label for valve name
            label = QtWidgets.QLabel(name)
            label.setAlignment(Qt.AlignCenter)
            label.setStyleSheet("font-weight: bold;")
            frame_layout.addWidget(label)
            
            # Create the valve control button
            valve = ValveControl(name, labjack_output, 0, 0, parent=frame)
            valve.setFixedHeight(50)  # Make buttons bigger in this window
            valve.setFixedWidth(100)
            
            # Update styling and remove positioning (since we're using layout)
            valve.clicked.disconnect(valve.confirm_toggle_valve)  # Remove old connection
            valve.clicked.connect(valve.confirm_toggle_valve)  # Add back to maintain functionality
            
            # Add to layout in a centered way
            button_container = QtWidgets.QHBoxLayout()
            button_container.addStretch()
            button_container.addWidget(valve)
            button_container.addStretch()
            frame_layout.addLayout(button_container)
            
            # Add frame to grid
            grid.addWidget(frame, row_pos, col_pos)
            
            # Add to main window's device collections
            self.main_window._devices.append(valve)
            self.main_window.device_map[name] = valve
    
    def update_connection_status(self):
        if self.main_window.labjack.connection_status:
            self.connection_status.setText("LabJack Status: Connected")
            self.connection_status.setStyleSheet("background-color: green; color: white; padding: 5px;")
        else:
            self.connection_status.setText("LabJack Status: Disconnected")
            self.connection_status.setStyleSheet("background-color: red; color: white; padding: 5px;")
    
    # Override closeEvent to close the main window when this window is closed
    def closeEvent(self, event):
        if not self.is_closing:
            self.is_closing = True
            self.main_window.close()
        event.accept()