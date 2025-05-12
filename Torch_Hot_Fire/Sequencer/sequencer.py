# from Sequencer.sequence_reader import load_data_from_csv
from PyQt5.QtWidgets import QPushButton, QMessageBox
from PyQt5.QtCore import QTimer, Qt
import csv

class Sequencer(QPushButton):
    def __init__(self, device_map, data_logger, x=0, y=0, parent=None):
        super(Sequencer, self).__init__(parent)
        self.current_event_index = 0
        self.running = False
        self.timer = None  # Store the timer reference so it can be stopped
        self.data_logger = data_logger
        self.pressure_data = []
        self.current_event_index = 0
        # Configure Button Display
        self.setFixedHeight(100) 
        self.setFixedWidth(194)
        self.setText("Start Sequencer")
        
        self.faulty_sequencer = True
        
        # Position the button if coordinates provided
        if x and y:
            self.move(x, y)
        
        # Set initial color
        self.update_button_style()

        # Store device map reference
        self.device_map = device_map
        
        # Load sequence data
        self.devices, self.events = self.load_data_from_csv(device_map)
        print(f"Loaded sequence with {len(self.devices) - 1} devices and {len(self.events)} events")

        # Connect the button click to toggle sequencer state
        self.clicked.connect(self.toggle_sequencer)
    
    def toggle_sequencer(self):
        """Toggle between starting and stopping the sequencer."""
        if not self.running:
            self.confirm_start_sequencer()
        else:
            self.confirm_stop_sequencer()
    
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
        if not self.devices or not self.events or self.faulty_sequencer:
            print("Error: SEQUENCER ERROR")
            QMessageBox.warning(self, "Sequencer Error", 
                               "INVALID SEQUENCE. Please check the sequence file.")
            return
        self.running = True
        print(f"Starting sequencer - setting running=True")
        if not self.data_logger.high_speed_mode:
            self.data_logger.toggle_sample_rate()
        self.current_event_index = 0
        self.setText("Stop Sequencer")
        self.update_button_style()
        self._trigger_event()
    
    def confirm_stop_sequencer(self):
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
            self.stop_sequencer()
    
    def stop_sequencer(self):
        """Stop the sequencing process."""
        print(f"Stopping sequencer - setting running=False")
        self.pressure_data = []
        self.running = False
        self.setText("Start Sequencer")
        self.update_button_style()
        
        # Check all valves
        print("Checking valve states:")
        valves_closed = 0
        for device in self.devices:
            if device == "Timestamp (ms)":
                continue
            if self.device_map[device].valve_open:
                print(f"Closing valve: {device}")
                self.device_map[device].toggle_valve()
                valves_closed += 1
            else:
                print(f"Valve already closed: {device}")
        print(f"Closed {valves_closed} open valves")
        

    def _trigger_event(self):
        """Recursive function to trigger events at the specified intervals."""
        # Check if we should still be running
        if not self.running:
            print("Trigger cancelled: sequencer not running")
            return
        
        if self.current_event_index >= len(self.events):
            if self.running:  # We've finished all events
                print(f"Completed all {len(self.events)} events")
                self.running = False
                self.setText("Start Sequencer")
                self.update_button_style()
            return
        
        if self.current_event_index == 0:
            # Verify intial state matches expected - extra layer of safety between Software Handler and Sequencer
            initial_state = self.events[0]
            if initial_state[0] != 0:
                raise ValueError(f"Incorrect initial timestamp in CSV")
            for i in range(1, len(initial_state)):
                if initial_state[i] == 0:
                    if self.device_map[self.devices[i]].valve_open:
                        self.stop_sequencer()
                        print("Error: Initial State Error")
                        QMessageBox.warning(self, "Initial State Error", 
                               "INVALID STARTING STATE. Please check the sequence file.")
                        raise ValueError(f"Sequencer is not prepared to start: Mismatch for {self.devices[i]}")
                if initial_state[i] == 1:
                    if not self.device_map[self.devices[i]].valve_open:
                        self.stop_sequencer()
                        print("Error: Initial State Error")
                        QMessageBox.warning(self, "Initial State Error", 
                               "INVALID STARTING STATE. Please check the sequence file.")
                        
                        raise ValueError(f"Sequencer is not prepared to start: Mismatch for {self.devices[i]}")
        # Initial state approved, continue

        # # Get the current time delay
        # current_delay = self.events[self.current_event_index]

        # Current event
        event = self.events[self.current_event_index]
        self.current_event_index += 1

        print("TRIGGERING EVENT")
        for i in range(1, len(event)):
            if (event[i] == 0): 
                self.device_map[self.devices[i]].toggle_valve_off()
            elif (event[i] == 1):
                self.device_map[self.devices[i]].toggle_valve_on()
            else:
                raise ValueError(f"Faulty input for {self.devices[i]} state for event: {event}")

        # Schedule the next event if there are more events remaining
        if self.current_event_index < len(self.events) and self.running:
            delay_time = self.events[self.current_event_index][0]
            print(f"Next event scheduled in {delay_time}ms")
            QTimer.singleShot(delay_time, self._trigger_event)
        else:
            # Cooldown timer
            QTimer.singleShot(2000, lambda: (
                self.data_logger.toggle_sample_rate(),
                self.stop_sequencer()
            ))
            if not self.running:
                print("Not scheduling next event - sequencer stopped")

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

    def load_data_from_csv(self, device_map):
        """
        Load both limits and sequence data from a CSV file.
        Processes limits section before "Sequence" line and device timing after.

        Args:
            device_map (dict): Dictionary mapping device names to device objects.

        Returns:
            tuple: A tuple containing (limits_dict, devices, events)
            - limits_dict: Dictionary mapping device names to pressure limits
            - devices: List of device objects in sequence
            - events: List of timing events in milliseconds
        """
        eventDevices = []
        events = []
        
        filename = 'Torch_Hot_Fire/Sequencer/Sequencer_Info/new_sequence_test2.csv'
        
        try:
            with open(filename, "r") as file:
                reader = csv.reader(file)
                
                # First section should be limits
                header = next(reader)
                if header[0] != "Limits":
                    print(header)
                    raise ValueError("Sequencer file is missing Limits header")
                # read in the redline devices
                redlineDevices = next(reader)
                redlines = next(reader)
                if len(redlines) != len(redlineDevices):
                    raise ValueError("Mismatched redline header and value rows")
                for i in range (0, len(redlines)):
                    device_name = redlineDevices[i]
                    limit = redlines[i]
                    if int(limit) == -1:
                        print(f"No redline added for {device_name}")
                        continue
                    try:
                        device_map[device_name].redline = float(limit)
                        print(f"Redline of {limit} psi added for {device_name}")
                    except ValueError:
                        print(f"Warning: Invalid inputs for {device_name}: {limit}")
                    except Exception as e:
                        print(f"Random Error for device: {device_name}, limit: {limit}, error: {e}")
                
                # Read in events header
                if next(reader)[0] != "Sequence":
                    raise ValueError("Sequence header not found")
                eventDevices = next(reader)
                prev_time = 0                    
                for event in reader:
                    if len(event) != len(eventDevices):
                        raise ValueError(f"Incorrect event row format: {event}")
                    delay = int(event[0]) - prev_time  # Convert timing to integer
                    if delay < 0:
                        raise ValueError(f"Negative delay detected - faulty timestamp for: {event}")
                    prev_time = int(event[0])
                    event[0] = delay
                    for i in range(1, len(event)):
                        event[i] = int(event[i])
                    events.append(event)
            if not eventDevices:
                print("Warning: No sequence data found.")
            
            for i in range(1, len(events[-1])):
                if events[-1][i] != 0:
                    raise(ValueError("Turn everything off in your sequencer shithead"))
            self.faulty_sequencer = False
        except Exception as e:
            print(f"Error loading data from {filename}: {e}")

        return eventDevices, events