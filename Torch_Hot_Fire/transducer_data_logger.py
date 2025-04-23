import csv
import threading
import queue
import os
from datetime import datetime
import time
from PyQt5 import QtWidgets
from PyQt5.QtCore import pyqtSignal

class TransducerDataLogger(QtWidgets.QPushButton):
    state_changed = pyqtSignal(bool)
    
    def __init__(self, transducers_list, devices_list, parent=None, filename="Torch_Hot_Fire/transducer_log.csv"):
        super().__init__(parent)  # Initialize the QPushButton parent class
        self.filename = filename
        self.log_queue = queue.Queue()
        self.running = True
        self._transducers = transducers_list
        self._devices = devices_list
        self.high_speed_mode = False  # Track current mode
        
        # Configure Button Display
        self.setFixedHeight(100) 
        self.setFixedWidth(194) 
        self.setText("LOGGING SPEED: \n Low Speed")  # Initial text
        
        # Set initial color
        self.update_button_style()
        
        # Connect the button click to toggle function
        self.clicked.connect(self.toggle_sample_rate)

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
        """Toggle between high speed (10ms) and low speed (500ms)"""
        if not self.pressure_timer:
            return
            
        if self.high_speed_mode:
            # Switch to low speed
            self.pressure_timer.stop()
            print("Setting timer interval to 500ms")
            self.pressure_timer.setInterval(500)
            self.pressure_timer.start()
            self.setText("LOGGING SPEED: \n Low Speed")
            self.high_speed_mode = False
        else:
            # Switch to high speed
            self.pressure_timer.stop()
            print("Setting timer interval to 10ms")
            self.pressure_timer.setInterval(10)
            self.pressure_timer.start()
            self.setText("LOGGING SPEED: \n High Speed")
            self.high_speed_mode = True
        
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
                    color: white;
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
                    color: white;
                    font-weight: bold;
                }
                QPushButton:hover {
                    background-color: #d32f2f;
                }
                QPushButton:pressed {
                    background-color: #b71c1c;
                }
            """)