from sequence_reader import load_sequence_from_csv
from PyQt5.QtWidgets import QPushButton, QMessageBox
from PyQt5.QtCore import QTimer

from sequence_reader import load_sequence_from_csv


class Sequencer(QPushButton):
    def __init__(self, device_map, x, y, parent=None):
        super(Sequencer, self).__init__(parent)
        self.current_event_index = 0
        self.running = False
        self.setText("Start Sequencer")
        self.move(x, y)
        self.adjustSize()
        self.devices, self.events = load_sequence_from_csv(device_map)

        # Connect the button click to start the sequencer
        self.clicked.connect(self.confirm_start_sequencer)

    def confirm_start_sequencer(self):
        # Display confirmation dialog before starting the sequencer.
        msg = QMessageBox(self)
        msg.setWindowTitle("Confirmation")
        msg.setText("Are you sure you want to start the sequencer?")
        msg.setStyleSheet("QMessageBox { background-color: #2e2e2e; } QLabel { color: white; } QPushButton { background-color: #4a4a4a; color: white; font-size: 10pt; }")
        yes_button = msg.addButton("Yes", QMessageBox.AcceptRole)
        no_button = msg.addButton("No", QMessageBox.RejectRole)
        msg.setDefaultButton(no_button)
        msg.exec_()

        if msg.clickedButton() == yes_button:
            self.start_sequencer()

    def start_sequencer(self):
        # Start the sequencing process.
        self.running = True
        self.current_event_index = 0
        self._trigger_event()
        print("Sequencer ended.")

    def _trigger_event(self):
        """Trigger events at the specified intervals."""
        if not self.running or self.current_event_index >= len(self.events):
            return

        # Get the current time delay
        current_delay = self.events[self.current_event_index]

        # Collect all devices with the same delay
        devices_to_trigger = []
        while self.current_event_index < len(self.events) and self.events[self.current_event_index] == current_delay:
            devices_to_trigger.append(self.devices[self.current_event_index])
            self.current_event_index += 1

        # Trigger all devices with the same delay
        for device in devices_to_trigger:
            if not device.device_connected:
                device.connect_to_labjack()  # Try to connect if not connected

            if device.device_connected:
                device.toggle_valve()
            else:
                QMessageBox.warning(self, "Connection Error",
                                    f"Failed to connect to {device.name}. Please check the connection.")

        # Schedule the next event if there are more events remaining
        if self.current_event_index < len(self.events):
            next_delay = self.events[self.current_event_index]
            QTimer.singleShot(next_delay - current_delay, self._trigger_event)

