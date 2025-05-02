from PyQt5.QtCore import QObject
from PyQt5.QtWidgets import QLabel
from labjack import ljm

class LabJackConnection(QObject):
    def __init__(self, status_label: QLabel):
        super().__init__()
        self.handle = None
        self.connection_status = False
        self.status_label = status_label

    def connect_to_labjack(self):
        if not self.connection_status:
            try:                
                self.handle = ljm.openS("ANY", "ANY", "ANY")
                self.connection_status = True
                print("Connection attempt: Success - Connection Established")
            except Exception as e:
                print("Connection attempt: Failed - Could not connect to LabJack:", e)
                self.connection_status = False
        # Update Status
        self.update_connection_status(self.connection_status)

    def update_connection_status(self, connected):
        """Update the QLabel with the connection status."""
        if self.status_label:
            self.connection_status = connected
            if connected:
                self.status_label.setText("LabJack T7: Connection Established")
                self.status_label.setStyleSheet("background-color: green; color: white; font-weight: bold;")
            else:
                self.status_label.setText("LabJack T7: Connection Missing")
                self.status_label.setStyleSheet("background-color: red; color: white; font-weight: bold;")

    def close_connection(self):
        """Close the LabJack connection."""
        if self.handle:
            ljm.close(self.handle)
            self.handle = None
            self.connection_status = False
            print("Connection closed.")
        self.update_connection_status(self.connection_status)
