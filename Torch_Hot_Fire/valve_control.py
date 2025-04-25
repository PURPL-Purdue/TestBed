from PyQt5 import QtWidgets
from labjack import ljm
import datetime as dt

# TODO: Fix connect_to_labjack redundancies. Many try catches
class ValveControl(QtWidgets.QPushButton):
    def __init__(self, name, labjack_output, x, y, parent=None):
        super(ValveControl, self).__init__(parent)
        print(f"Creating valve {name} with parent: {parent}")
        self.name = name
        self.labjack_output = labjack_output
        self.valve_open = False
        self.device_connected = False
        self.handle = None
        # self.setText("C")
        # self.setText("Valve Closed")
        self.update_button_style()
        self.clicked.connect(self.toggle_valve)
        self.move(x, y)
        # self.adjustSize()
        self.setFixedHeight(30)
        self.setFixedWidth(12)

        # Initialize a connection and read initial values
        self.connect_to_labjack()
        if self.device_connected:
            try:
                # Read the current state from the LabJack
                current_state = ljm.eReadName(self.handle, self.labjack_output)
                print(f"This is the current state: {current_state}")
                # Update our internal state based on actual hardware state (0=open, 1=closed)
                self.valve_open = (current_state == 0)

                if (self.valve_open):
                    # Update button appearance to match actual state
                    self.update_button_style()
                    # eReadName resets the switch to off, this will switch it back to on
                    self.toggle_valve_on()

                print(f"Initialized {self.name} - Current state: {'OPEN' if self.valve_open else 'CLOSED'}")
            except Exception as e:
                print(f"Error reading initial state of {self.name}: {e}")

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
        else: 
            print("Device already connected")

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
            if (self.valve_open == True):
                self.toggle_valve_off()
            elif (self.valve_open == False):
                self.toggle_valve_on()

    def toggle_valve(self):
        # Ensure connection is established before showing confirmation
        if not self.device_connected:
            self.connect_to_labjack()

        if not self.device_connected:
            QtWidgets.QMessageBox.warning(self, "Connection Error", "Failed to connect to LabJack. Please try again.")
            return
        
        if self.valve_open:
            self.toggle_valve_off
        else:
            self.toggle_valve_on

    def toggle_valve_on(self):
        """Toggle on the valve and update LabJack output"""
        self.valve_open = True
        self.update_button_style()
        self.update_labjack_output()

    def toggle_valve_off(self):
        """Toggle off the valve and update LabJack output"""
        self.valve_open = False
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
        # text, color = ("Valve Open", "green") if self.valve_open else ("Valve Closed", "red")
        text, color = ("O", "green") if self.valve_open else ("C", "red")
        self.setText(text)
        self.setStyleSheet(f"background-color: {color}; color: white; font-size: 12pt; font-weight: bold;")
