import csv

def load_data_from_csv(device_map):
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
    devices = []
    events = []
    filename = 'Torch_Hot_Fire/Sequencer_Info/sequence.csv'
    
    try:
        with open(filename, "r") as file:
            reader = csv.reader(file)
            
            # First section should be limits
            current_section = None
            for row in reader:
                # Skip empty rows
                if not row:
                    continue
                    
                # Check if this is a section header
                if len(row) == 1:
                    current_section = row[0]
                    # Skip the column headers row after a section marker
                    next(reader)
                    continue
                    
                # Process based on current section
                if current_section == "Limits":
                    if len(row) >= 2:
                        device_name, limit = row
                        try:
                            device_map[device_name].redline = float(limit)
                            print(f"Redline of {limit} psi added for {device_name}")
                        except ValueError:
                            print(f"Warning: Invalid inputs for {device_name}: {limit}")
                            
                elif current_section == "Sequence":
                    if len(row) >= 2:
                        device_name, timing = row
                        try:
                            timing = int(timing)  # Convert timing to integer
                            
                            if device_name in device_map:
                                devices.append(device_map[device_name])
                                events.append(timing)
                            else:
                                print(f"Warning: Device '{device_name}' not found in device map.")
                        except ValueError:
                            print(f"Warning: Invalid timing value: {timing}")

        if not devices:
            print("Warning: No sequence data found.")
            
    except Exception as e:
        print(f"Error loading data from {filename}: {e}")

    return devices, events