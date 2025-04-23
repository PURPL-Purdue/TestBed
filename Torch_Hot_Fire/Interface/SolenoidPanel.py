
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt
from valve_control import ValveControl
from PyQt5.QtCore import QTimer

class SolenoidWindow(QtWidgets.QWidget):
    def __init__(self, main_window):
        super().__init__()  # No parent - this creates an independent window
        self.windim_x, self.windim_y = 600, 230
        self.setWindowTitle("Valve Control Panel")
        
        # Calculate position to be on the left side of main window
        main_x = main_window.geometry().x()
        main_y = main_window.geometry().y()
        window_x = main_x - self.windim_x  # 10px gap
        if window_x < 0:
            window_x = 0  # Prevent positioning off-screen
            
        window_y = main_y
        self.setGeometry(window_x, window_y, self.windim_x, self.windim_y)
        
        # Store reference to main window to access its data
        self.main_window = main_window
        self.is_closing = False
        
        # Background
        bg_label = QtWidgets.QLabel(self)
        bg_pixmap = QtGui.QPixmap("Torch_Hot_Fire/solenoid_bg.png")
        
        # Handle image alignment to right side of window
        if bg_pixmap.width() > self.windim_x:
            # If image is wider than window, crop it from the right side
            cropped_pixmap = bg_pixmap.copy(bg_pixmap.width() - self.windim_x, 0, self.windim_x, bg_pixmap.height())
            bg_label.setPixmap(cropped_pixmap)
        else:
            # If image is not wider, just scale it and align to right
            scaled_pixmap = bg_pixmap.scaled(bg_pixmap.width(), self.windim_y, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            bg_label.setPixmap(scaled_pixmap)
            bg_label.setGeometry(self.windim_x - scaled_pixmap.width(), 0, scaled_pixmap.width(), self.windim_y)
        
        # Set border - matching main window style
        self.border_frame = QtWidgets.QFrame(self)
        self.border_frame.setFrameStyle(QtWidgets.QFrame.Box)
        self.border_frame.setLineWidth(5)
        self.border_frame.setGeometry(0, 0, self.windim_x, self.windim_y)
        self.border_frame.setStyleSheet("border: 5px solid #f44336;")  # Red border initially
        
        # # Title label - positioned at top like in MainWindow
        # title = QtWidgets.QLabel("Valve Control Panel", self)
        # title.setGeometry(215, 10, 200, 30)  # Centered in the wider window
        # title.setAlignment(Qt.AlignCenter)
        # title.setStyleSheet("font-size: 14pt; font-weight: bold; color: white; background-color: rgba(0, 0, 0, 120);")
        
        # # Connection Status Label - similar to MainWindow
        # self.connection_status = QtWidgets.QLabel("LabJack Status: Disconnected", self)
        # self.connection_status.setGeometry(195, self.windim_y - 40, 240, 25)  # Centered at bottom
        # self.connection_status.setAlignment(Qt.AlignCenter)
        # self.connection_status.setStyleSheet("background-color: red; color: white; font-weight: bold;")
        
        # # Create the valve controls directly on the window (like in MainWindow)
        # self.create_valve_controls()
        
        # # Timer to update connection status
        # self.timer = QTimer(self)
        # self.timer.timeout.connect(self.update_connection_status)
        # self.timer.start(1000)  # Update every second
    
    # def create_valve_controls(self):
        # # Clear existing devices in main window
        # self.main_window._devices = []
        # self.main_window.device_map = {}
        
        # # Create new valve controls with absolute positioning (like in MainWindow)
        # Spread across the wider 630px window
        # valve_configs = [
        #     ("SN-H2-01", "CIO0", 120, 90),
        #     ("SN-OX-01", "CIO1", 320, 90),
        #     ("SN-N2-01", "EIO7", 230, 160),
        #     ("Spark-Plug", "CIO3", 430, 90)
        # ]
        
        # for name, labjack_output, x_pos, y_pos in valve_configs:
        #     # Create the valve control - using absolute positioning like in MainWindow
        #     valve = ValveControl(name, labjack_output, x_pos, y_pos, parent=self)
            
        #     # Add to main window's device collections
        #     self.main_window._devices.append(valve)
        #     self.main_window.device_map[name] = valve
    
    # def update_connection_status(self):
    #     if self.main_window.labjack.connection_status:
    #         self.connection_status.setText("LabJack Status: Connected")
    #         self.connection_status.setStyleSheet("background-color: green; color: white; padding: 5px; font-weight: bold;")
    #         # Update border to match data logger state
    #         if hasattr(self.main_window, 'data_logger') and self.main_window.data_logger.high_speed_mode:
    #             self.border_frame.setStyleSheet("border: 5px solid #4CAF50;")  # Green border for high speed
    #         else:
    #             self.border_frame.setStyleSheet("border: 5px solid #f44336;")  # Red border for normal speed
    #     else:
    #         self.connection_status.setText("LabJack Status: Disconnected")
    #         self.connection_status.setStyleSheet("background-color: red; color: white; padding: 5px; font-weight: bold;")
    
    # # Override closeEvent to avoid closing loops
    # def closeEvent(self, event):
    #     if not self.is_closing:
    #         self.is_closing = True
    #         # Don't close the main window from here
    #         # Just accept the close event for this window
    #     event.accept()