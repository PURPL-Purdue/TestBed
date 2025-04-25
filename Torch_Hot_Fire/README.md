# TeenyK P&ID Control System

A PyQt5-based application for controlling and monitoring a fluid system with LabJack interface for valve control, pressure readings, data logging, and automated sequencing.

## Overview

This application provides a graphical interface for controlling a fluid control system with the following features:

- Real-time pressure monitoring from multiple transducers
- Solenoid valve control with visual state indicators
- Data logging with adjustable sampling rates
- Automated sequencing capability for valve operations
- Emergency shutdown functionality

## System Requirements

- Python 3.x
- PyQt5
- LabJack T7 device and LJM library
- Additional Python libraries: csv, threading, queue, datetime

## Installation

1. Ensure you have Python 3.x installed
2. Install required packages:
   ```
   pip install pyqt5 labjack-ljm
   ```
3. Clone or download this repository to your local machine

## Running the Application

The application is started by running the main.py file:

```
python main.py
```

This will open the fluid panel interface, which is the main entry point for the program. The fluid panel will open all other necessary panels automatically.

## Application Structure

- **main.py**: Application entry point
- **FluidPanel.py**: Main control window showing system P&ID
- **SolenoidPanel.py**: Dedicated panel for valve control
- **data_logger.py**: Handles data logging functionality
- **pressure_transducer.py**: Interface for pressure sensors
- **valve_control.py**: Controls for solenoid valves
- **sequencer.py**: Handles automated valve sequencing

## Usage Instructions

### Main Interface

The application opens with two windows:

1. **TeenyK P&ID** - Main interface showing system diagram with pressure readings
2. **Valve Control Panel** - Dedicated interface for controlling solenoid valves

### Data Logging

- The **LOGGING SPEED** button toggles between high-speed (10ms) and low-speed (500ms) data collection
- Enter a base filename in the text field before starting logging
- Log files are saved in the "Torch_Hot_Fire/data" directory with timestamps
- The border color of both panels indicates the current logging speed (red = low, green = high)

### Valve Control

- Valves are shown in the Valve Control Panel
- Click a valve to change its state (open/closed)
- "O" (green) indicates an open valve, "C" (red) indicates a closed valve

### Sequencer

- The **Start Sequencer** button initiates an automated valve sequence
- Sequences are loaded from a CSV file via the sequence_reader module
- The sequencer automatically switches to high-speed logging when started

### Emergency Shutdown

- The **Emergency Shutdown** button will close all valves and stop any active sequence
- Use this instead of the sequencer's stop button due to known issues

## File Formats

### Data Log Files

Log files are CSV format with the following columns:

- Timestamp
- Pressure readings for each transducer
- State of each valve (True/False)

## Closing the Application

Close the Fluid Panel to close the entire program. This will perform a proper shutdown of all components.

## Known Issues

- **Sequencer Stop**: The "Stop Sequencer" button doesn't reliably stop the sequencer. Use the Emergency Shutdown button instead.

## Future Development

- Consolidate all panels into one unified P&ID Panel
- Fix sequencer stop functionality
- Add additional error handling and connection recovery

## Directory Structure

The application expects the following directory structure:

```
/
├── main.py
├── Interface/
│   ├── FluidPanel.py
│   └── SolenoidPanel.py
├── valve_control.py
├── pressure_transducer.py
├── data_logger.py
├── sequencer.py
├── sequence_reader.py
├── labjack_connection.py
└── Torch_Hot_Fire/
    ├── data/
    ├── torch_bk_bg.png
    └── solenoid_bg.png
```

## Permissions and Hardware Access

This application requires appropriate permissions to access the LabJack hardware. Ensure the user running the application has the required system permissions.
