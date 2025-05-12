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
    
    filename = 'Torch_Hot_Fire/Sequencer/Sequencer_Info/torch_sequence8.csv'
    
    try:
        with open(filename, "r") as file:
            reader = csv.reader(file)
            
            # First section should be limits
            if next(reader) != "Limits":
                raise ValueError("Sequencer file is missing Limits header")
            # read in the redline devices
            redlineDevices = next(reader)
            redlines = next(reader)
            if len(redlines) != len(redlineDevices):
                raise ValueError("Mismatched redline header and value rows")
            for i in range (0, len(redlines)):
                device_name = redlineDevices[i]
                limit = redlines[i]
                if limit == -1:
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
            if next(reader) != "Sequence":
                raise ValueError("Sequence header not found")
            eventDevices = next(reader)
            prev_time = 0
            for event in reader:
                if len(event) != len(eventDevices):
                    raise ValueError(f"Incorrect event row format: {event}")
                try:
                    delay = int(event[0]) - prev_time  # Convert timing to integer
                    if delay < 0:
                        raise ValueError(f"Negative delay detected - faulty timestamp for: {event}")

                    for i in range(1, len(event)):
                        device_name = event[i]
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