import csv

def load_sequence_from_csv(device_map):
    """
    Load device sequence and timing data from a CSV file.

    Args:
        device_map (dict): Dictionary mapping device names to device objects.

    Returns:
        tuple: A tuple of two lists (devices, events).
    """
    devices = []
    events = []
    filename = 'Sequencer_Info/sequence.csv'
    try:
        with open(filename, "r") as file:
            reader = csv.reader(file)
            next(reader)  # Skip the header row

            for row in reader:
                if len(row) < 2:
                    print(f"Warning: Skipping invalid row {row}")
                    continue

                device_name, timing = row
                timing = int(timing)  # Convert the timing to an integer

                if device_name in device_map:
                    devices.append(device_map[device_name])
                    events.append(timing)
                else:
                    print(f"Warning: Device '{device_name}' not found in device map.")
    except Exception as e:
        print(f"Error loading sequence from {filename}: {e}")

    return devices, events
