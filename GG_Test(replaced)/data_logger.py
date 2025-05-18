import csv
import threading
import queue
import os
from datetime import datetime
import time
from PyQt5 import QtWidgets
from PyQt5.QtCore import pyqtSignal, Qt

class DataLogger(QtWidgets.QPushButton):
    state_changed = pyqtSignal(bool)
    
    def __init__(self, transducers_list, devices_list, width = 194, height = 100, parent=None, path="GG_Test/data"):
        super().__init__(parent)  # Initialize the QPushButton parent class
        self.path = path
        self.log_queue = queue.Queue()
        self.running = True
        self._transducers = transducers_list
        self._devices = devices_list
        self.high_speed_mode = False  # Track current mode
        self.base_name = "start"

        # Create layout for the button and the filename textbox
        self.button_layout = QtWidgets.QVBoxLayout(self)
        
        # Create label for the textbox
        self.filename_label = QtWidgets.QLabel("File Name:", self)
        self.filename_label.setStyleSheet("color: white; font-weight: bold;")
        
        # Create the filename textbox
        self.filename_textbox = QtWidgets.QLineEdit(self)
        self.filename_textbox.setText(self.base_name)
        self.filename_textbox.setStyleSheet("background-color: white; color: black;")

        # Configure Button Display
        self.setFixedHeight(height) 
        self.setFixedWidth(width) 

        # Create a widget for the button text
        self.button_text = QtWidgets.QLabel("LOGGING SPEED: \n Low Speed", self)
        self.button_text.setAlignment(Qt.AlignCenter)
        self.button_text.setStyleSheet("color: white; font-weight: bold;")
        
        # Add widgets to layout
        self.button_layout.addWidget(self.filename_label)
        self.button_layout.addWidget(self.filename_textbox)
        self.button_layout.addWidget(self.button_text)
        self.button_layout.setContentsMargins(10, 10, 10, 10)

        # Set initial color
        self.update_button_style()

        # Connect the button click to toggle function
        self.clicked.connect(self.toggle_sample_rate)

        # Create a new log file with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if self.high_speed_mode:
            mode_str = "high"
        else:
            mode_str = "low"
        self.filename = f"{self.path}/{self.base_name}_{mode_str}_{timestamp}.csv"
        print(f"Creating new log file: {self.filename}")
        
        # Reference to the timer (to be set from main window)
        self.pressure_timer = None

        # Create CSV file and write header if it doesn't exist
        if not os.path.exists(self.filename):
            os.makedirs(os.path.dirname(self.filename), exist_ok=True)  # Create directory if needed
            with open(self.filename, mode='w', newline='') as file:
                entry = ["Timestamp"]
                for transducer in self._transducers:
                    entry.append(f"{transducer.name} Pressure")
                for device in self._devices:
                    entry.append(f"{device.name} State")
                writer = csv.writer(file)
                writer.writerow(entry)

        # Start the background thread to handle writing to the CSV
        self.thread = threading.Thread(target=self._process_queue)
        self.thread.daemon = True  # Make thread daemon so it exits when main program exits
        self.thread.start()
    
    def set_timer(self, timer):
        """Set reference to the pressure timer"""
        self.pressure_timer = timer
    
    def toggle_sample_rate(self):
        """Toggle between high speed (10ms) and low speed (500ms) and create a new log file"""
        if not self.pressure_timer:
            return
        
        # Get the base name from the textbox
        self.base_name = self.filename_textbox.text()
        if not self.base_name:  # If empty, use default
            self.base_name = "log"
            self.filename_textbox.setText(self.base_name)
        
        # Create a new log file with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if self.high_speed_mode:
            # Switch to low speed
            mode_str = "low"
            self.pressure_timer.stop()
            print("Setting timer interval to 500ms")
            self.pressure_timer.setInterval(500)
            self.pressure_timer.start()
            self.button_text.setText("LOGGING SPEED: \n Low Speed")
            self.high_speed_mode = False
        else:
            # Switch to high speed
            mode_str = "high"
            self.pressure_timer.stop()
            print("Setting timer interval to 10ms")
            self.pressure_timer.setInterval(10)
            self.pressure_timer.start()
            self.button_text.setText("LOGGING SPEED: \n High Speed")
            self.high_speed_mode = True

        # Create new filename with timestamp and mode
        new_filename = f"{self.path}/{self.base_name}_{mode_str}_{timestamp}.csv"
        print(f"Creating new log file: {new_filename}")
        
        # Create the new file with headers
        if not os.path.exists(self.path):
            os.makedirs(self.path, exist_ok=True)
        
        with open(new_filename, mode='w', newline='') as file:
            entry = ["Timestamp"]
            for transducer in self._transducers:
                entry.append(f"{transducer.name} Pressure")
            for device in self._devices:
                entry.append(f"{device.name} State")
            writer = csv.writer(file)
            writer.writerow(entry)
        
        # Update the filename to use for future logging
        self.filename = new_filename
        
        # Update button color
        self.update_button_style()
        
        # Emit signal with current state
        self.state_changed.emit(self.high_speed_mode)

    def log_data(self):
        """Put the log data in the queue for the background thread to process."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        entry = [timestamp]
        for transducer in self._transducers:
            entry.append(transducer.pressure)
        for device in self._devices:
            entry.append(device.valve_open)
        self.log_queue.put(entry)

    def _process_queue(self):
        """Background thread function to process the queue and write to the CSV file."""
        while self.running or not self.log_queue.empty():
            try:
                # Write the data to the CSV
                with open(self.filename, mode='a', newline='') as file:
                    writer = csv.writer(file)
                    while not self.log_queue.empty():
                        # Get data with a timeout to prevent hanging
                        data = self.log_queue.get(timeout=1)
                        writer.writerow(data)
                        self.log_queue.task_done()
            except queue.Empty:
                # If timeout occurs file will close and reopen, continue the loop
                continue
            except Exception as e:
                print(f"Error writing to log file: {e}")
                time.sleep(1)  # Prevent CPU spinning on persistent errors

    def stop(self):
        """Stop the logging and gracefully shut down the background thread."""
        print("Terminating Data Logger")
        self.running = False
        if self.thread.is_alive():
            self.thread.join(timeout=2)  # Wait for the thread with timeout

    def update_button_style(self):
        """Update button color based on current mode"""
        if self.high_speed_mode:
            # Green for high speed
            self.setStyleSheet("""
                QPushButton {
                    background-color: #4CAF50;  /* Green */
                    border: none;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #45a049;
                }
                QPushButton:pressed {
                    background-color: #3e8e41;
                }
            """)
        else:
            # Red for low speed
            self.setStyleSheet("""
                QPushButton {
                    background-color: #f44336;  /* Red */
                    border: none;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #d32f2f;
                }
                QPushButton:pressed {
                    background-color: #b71c1c;
                }
            """)