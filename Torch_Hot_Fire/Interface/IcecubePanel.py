from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt
# from valve_control import ValveControl
# from PyQt5.QtCore import QTimer

class TorchWindow(QtWidgets.QWidget):
    def __init__(self, main_window):
        super().__init__()  # No parent - this creates an independent window
        self.windim_x, self.windim_y = 600, 350
        self.setWindowTitle("Torch Panel")
        
        # Calculate position to be on the left side of main window
        main_x = main_window.geometry().x()
        main_y = main_window.geometry().y()
        window_x = main_x - self.windim_x 
        if window_x < 0:
            window_x = 0  # Prevent positioning off-screen
            
        window_y = main_y + 260
        self.setGeometry(window_x, window_y, self.windim_x, self.windim_y)
        
        
        # Store reference to main window to access its data
        self.main_window = main_window
        self.is_closing = False
        
        # Background
        bg_label = QtWidgets.QLabel(self)
        bg_pixmap = QtGui.QPixmap("Torch_Hot_Fire/icecubebg.png")
        scaled_pixmap = bg_pixmap.scaled(self.windim_x, self.windim_y, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        bg_label.setPixmap(scaled_pixmap)
        bg_label.setGeometry(self.windim_x - scaled_pixmap.width(), 0, scaled_pixmap.width(), self.windim_y)
        
        # Set border - matching main window style
        self.border_frame = QtWidgets.QFrame(self)
        self.border_frame.setFrameStyle(QtWidgets.QFrame.Box)
        self.border_frame.setLineWidth(5)
        self.border_frame.setGeometry(0, 0, self.windim_x, self.windim_y)
        self.border_frame.setStyleSheet("border: 5px solid #f44336;")  # Red border initially

    def update_border_color(self, high_speed_mode):
        """Update the border color to match the main window"""
        if high_speed_mode:
            # Green border for high speed mode
            self.border_frame.setStyleSheet("border: 5px solid #4CAF50;")
        else:
            # Red border for normal speed mode
            self.border_frame.setStyleSheet("border: 5px solid #f44336;")