from PyQt5 import QtWidgets
from labjack import ljm
import datetime as dt

class ValveControl(QtWidgets.QPushButton):
    def __init__(self, name, labjack_output, x, y, parent=None):
        super(ValveControl, self).__init__(parent)
        self.name = name
        self.labjack_output = labjack_output
        self.valve_open = False
        self.device_connected = False
        self.handle = None
        self.setText("Valve Closed")
        self.update_button_style()
        self.clicked.connect(self.confirm_toggle_valve)
        self.move(x, y)
        self.adjustSize()

    def connect_to_labjack(self):
        """Establish connection to LabJack if not already connected."""
        if not self.device_connected:
            try:
                # Attempt to open LabJack device
                self.handle = ljm.openS("T7", "ANY", "ANY")
                self.device_connected = True
            except Exception as e:
                print(f"Failed to connect to LabJack: {e}")
                self.device_connected = False

    def confirm_toggle_valve(self):
        """Display confirmation dialog before toggling the valve."""
        # Ensure connection is established before showing confirmation
        if not self.device_connected:
            self.connect_to_labjack()

        if not self.device_connected:
            QtWidgets.QMessageBox.warning(self, "Connection Error", "Failed to connect to LabJack. Please try again.")
            return

        # Create a confirmation message box
        msg = QtWidgets.QMessageBox(self)
        msg.setWindowTitle("Confirmation")
        msg.setText(f"Are you certain you want to {'close' if self.valve_open else 'open'} {self.name}?")

        # Apply dark mode styling
        msg.setStyleSheet(
            "QMessageBox { background-color: #2e2e2e; }"
            "QLabel { color: white; }"
            "QPushButton { background-color: #4a4a4a; color: white; font-size: 10pt; }"
        )

        # Add Yes and No buttons
        yes_button = msg.addButton("Yes", QtWidgets.QMessageBox.AcceptRole)
        no_button = msg.addButton("No", QtWidgets.QMessageBox.RejectRole)
        msg.setDefaultButton(no_button)
        msg.exec_()

        if msg.clickedButton() == yes_button:
            self.toggle_valve()

    def toggle_valve(self):
        """Toggle the valve state and update LabJack output."""
        if not self.device_connected:
            self.connect_to_labjack()
            print("Failed to toggle valve because the device is not connected.")
            return

        self.valve_open = not self.valve_open
        self.update_button_style()
        self.update_labjack_output()

    def update_labjack_output(self):
        """Send the valve state to the LabJack output."""
        if self.device_connected and self.handle:
            output_value = 0 if self.valve_open else 1
            try:
                ljm.eWriteName(self.handle, self.labjack_output, output_value)
                print(f"{self.labjack_output} set to {output_value} for {'open' if self.valve_open else 'closed'} valve state at {dt.datetime.now()}")
            except Exception as e:
                print(f"Error writing to {self.labjack_output}: {e}")

    def update_button_style(self):
        """Update the button's text and style based on the valve state."""
        text, color = ("Valve Open", "green") if self.valve_open else ("Valve Closed", "red")
        self.setText(text)
        self.setStyleSheet(f"background-color: {color}; color: white; font-size: 12pt; font-weight: bold;")
