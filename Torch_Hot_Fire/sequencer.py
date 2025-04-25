from sequence_reader import load_sequence_from_csv
from PyQt5.QtWidgets import QPushButton, QMessageBox
from PyQt5.QtCore import QTimer, Qt

class Sequencer(QPushButton):
    def __init__(self, device_map, x=0, y=0, parent=None):
        super(Sequencer, self).__init__(parent)
        self.current_event_index = 0
        self.running = False
        self.timer = None  # Store the timer reference so it can be stopped
        
        # Configure Button Display
        self.setFixedHeight(100) 
        self.setFixedWidth(194)
        self.setText("Start Sequencer")
        
        # Position the button if coordinates provided
        if x and y:
            self.move(x, y)
        
        # Set initial color
        self.update_button_style()

        # Store device map reference
        self.device_map = device_map
        
        # Load sequence data
        self.devices, self.events = load_sequence_from_csv(device_map)

        # Connect the button click to toggle sequencer state
        self.clicked.connect(self.toggle_sequencer)
    
    def toggle_sequencer(self):
        """Toggle between starting and stopping the sequencer."""
        if not self.running:
            self.confirm_start_sequencer()
        else:
            self.stop_sequencer()
    
    def confirm_start_sequencer(self):
        """Display confirmation dialog before starting the sequencer."""
        msg = QMessageBox(self)
        msg.setWindowTitle("Confirmation")
        msg.setText("Are you sure you want to start the sequencer?")
        msg.setStyleSheet("QMessageBox { background-color: #2e2e2e; } "
                         "QLabel { color: white; } "
                         "QPushButton { background-color: #4a4a4a; color: white; font-size: 10pt; }")
        yes_button = msg.addButton("Yes", QMessageBox.AcceptRole)
        no_button = msg.addButton("No", QMessageBox.RejectRole)
        msg.setDefaultButton(no_button)
        msg.exec_()

        if msg.clickedButton() == yes_button:
            self.start_sequencer()

    def start_sequencer(self):
        """Start the sequencing process."""
        if not self.devices or not self.events:
            QMessageBox.warning(self, "Sequencer Error", 
                               "No sequence data loaded. Please check the sequence file.")
            return
            
        self.running = True
        self.current_event_index = 0
        self.setText("Stop Sequencer")
        self.update_button_style()
        self._trigger_event()
    
    def stop_sequencer(self):
        """Stop the sequencing process."""
        # Display confirmation dialog
        msg = QMessageBox(self)
        msg.setWindowTitle("Confirmation")
        msg.setText("Are you sure you want to stop the sequencer?")
        msg.setStyleSheet("QMessageBox { background-color: #2e2e2e; } "
                         "QLabel { color: white; } "
                         "QPushButton { background-color: #4a4a4a; color: white; font-size: 10pt; }")
        yes_button = msg.addButton("Yes", QMessageBox.AcceptRole)
        no_button = msg.addButton("No", QMessageBox.RejectRole)
        msg.setDefaultButton(no_button)
        msg.exec_()

        if msg.clickedButton() == yes_button:
            self.running = False
            if self.timer:
                self.timer.stop()
            self.setText("Start Sequencer")
            self.update_button_style()
            print("Sequencer stopped by user.")
            
    def _trigger_event(self):
        """Trigger events at the specified intervals."""
        if not self.running or self.current_event_index >= len(self.events):
            if self.running:  # We've finished all events
                self.running = False
                self.setText("Start Sequencer")
                self.update_button_style()
                print("Sequencer completed all events.")
            return

        # Get the current time delay
        current_delay = self.events[self.current_event_index]

        # Collect all devices with the same delay
        devices_to_trigger = []
        while (self.current_event_index < len(self.events) and 
               self.events[self.current_event_index] == current_delay):
            devices_to_trigger.append(self.devices[self.current_event_index])
            self.current_event_index += 1

        # Trigger all devices with the same delay
        for device in devices_to_trigger:
            try:
                if not device.device_connected:
                    device.connect_to_labjack()  # Try to connect if not connected

                if device.device_connected:
                    device.toggle_valve()
                    print(f"Triggered device: {device.name} at time {current_delay}ms")
                else:
                    QMessageBox.warning(self, "Connection Error",
                                       f"Failed to connect to {device.name}. Please check the connection.")
            except Exception as e:
                print(f"Error triggering device {device.name}: {e}")

        # Schedule the next event if there are more events remaining
        if self.current_event_index < len(self.events):
            next_delay = self.events[self.current_event_index]
            delay_time = next_delay - current_delay
            self.timer = QTimer.singleShot(delay_time, self._trigger_event)
            print(f"Next event scheduled in {delay_time}ms")

    def update_button_style(self):
        """Update button color based on current mode"""
        if self.running:
            # Green for running
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
            # Red for stopped
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