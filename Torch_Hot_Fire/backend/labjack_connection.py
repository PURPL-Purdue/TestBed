from PyQt5.QtCore import QObject, QTimer
from PyQt5.QtWidgets import QLabel
from labjack import ljm

class LabJackConnection(QObject):
    def __init__(self, status_label: QLabel):
        super().__init__()
        self.handle = None
        self.connection_status = False
        self.status_label = status_label
        self.consecutive_failures = 0
        self.max_failures_before_disconnect = 1  # Require 3 consecutive failures before disconnecting

        # Heartbeat Timer
        self.heartbeat_timer = QTimer(self)
        self.heartbeat_timer.timeout.connect(self.heartbeat_check)
        self.heartbeat_interval_ms = 200  

    def connect_to_labjack(self):
        # Always try to close any existing handle first
        if self.handle is not None:
            try:
                ljm.close(self.handle)
            except Exception:
                pass  # Ignore errors on close
            self.handle = None
        
        try:
            # Completely fresh connection attempt
            self.handle = ljm.openS("ANY", "ANY", "ANY")
            # ljm.writeLibraryConfigS(ljm.constants.LJME_OPEN_TIMEOUT_MS, 100)  # 0.5 seconds
            ljm.writeLibraryConfigS("LJM_OPEN_TCP_DEVICE_TIMEOUT_MS", 100)  # Reduce TCP timeout to 500 ms
            ljm.writeLibraryConfigS("LJM_SEND_RECEIVE_TIMEOUT_MS", 100)
            
            # Test the connection immediately
            serial = ljm.eReadName(self.handle, "SERIAL_NUMBER")
            print(f"Connection successful - Device serial: {serial}")
            
            self.consecutive_failures = 0
            self.connection_status = True
            
            # Start heartbeat timer if not already running
            if not self.heartbeat_timer.isActive():
                self.heartbeat_timer.start(self.heartbeat_interval_ms)
                
        except Exception as e:
            print(f"Connection attempt failed: {e}")
            self.connection_status = False
            self.handle = None
            
            # Start heartbeat for reconnection attempts if not running
            if not self.heartbeat_timer.isActive():
                self.heartbeat_timer.start(self.heartbeat_interval_ms)

        self.update_connection_status(self.connection_status)

    def update_connection_status(self, connected):
        """Update the QLabel with the connection status."""
        self.connection_status = connected
        if self.status_label:
            if connected:
                self.status_label.setText("LabJack T7: Connection Established")
                self.status_label.setStyleSheet("background-color: green; color: white; font-weight: bold;")
            else:
                self.status_label.setText("LabJack T7: Connection Missing")
                self.status_label.setStyleSheet("background-color: red; color: white; font-weight: bold;")

    def heartbeat_check(self):
        """Periodically called to verify connection is still alive and to reconnect if needed."""
        # If we think we're connected and have a handle
        if self.handle is not None:
            try:
                # Try reading the serial number to test connection
                # print("Achieved1")
                ljm.eReadName(self.handle, "SERIAL_NUMBER")
                # print("Achieved2")
                # Connection is good
                self.consecutive_failures = 0
                
                if not self.connection_status:
                    print("Heartbeat: Connection restored")
                    self.connection_status = True
                    self.update_connection_status(True)
                    
            except Exception as e:
                # Increment failure counter
                self.consecutive_failures += 1
                
                if self.consecutive_failures >= self.max_failures_before_disconnect:
                    print(f"Heartbeat: Connection lost after {self.consecutive_failures} failures - {e}")
                    
                    # Clean up the handle
                    try:
                        ljm.close(self.handle)
                    except Exception:
                        pass  # Ignore errors on close
                    
                    self.handle = None
                    if self.connection_status != False:
                        self.connection_status = False
                        self.update_connection_status(False)
                else:
                    print(f"Heartbeat: Potential connection issue (failure {self.consecutive_failures}/{self.max_failures_before_disconnect})")
                    
        # If we don't have a handle, try to establish a new connection
        elif not self.handle:
            # Reset consecutive failures counter since we're trying fresh
            self.consecutive_failures = 0
            
            try:
                # Fresh connection attempt
                # print("Ping 1")
                self.handle = ljm.openS("ANY", "ANY", "ANY")
                # print("Ping2")
                # Test the connection immediately
                serial = ljm.eReadName(self.handle, "SERIAL_NUMBER")
                print(f"Heartbeat: Connected to device with serial: {serial}")
                
                self.connection_status = True
                self.update_connection_status(True)
                
            except Exception as e:
                # Only log failures occasionally to avoid spamming the console
                # print("Ping3")
                if self.connection_status: # or self.consecutive_failures % 10 == 0
                    print(f"Heartbeat: Connection attempt failed - {e}")
                
                # Clean up any partially created handle
                if self.handle is not None:
                    try:
                        ljm.close(self.handle)
                    except Exception:
                        pass
                    self.handle = None
                
                self.connection_status = False
                self.update_connection_status(False)

    def close_connection(self):
        """Close the LabJack connection."""
        if self.handle is not None:
            try:
                ljm.close(self.handle)
            except Exception as e:
                print(f"Error closing connection: {e}")
            
            self.handle = None
            
        self.connection_status = False
        self.heartbeat_timer.stop()  # Stop the heartbeat when connection is closed
        self.update_connection_status(False)
