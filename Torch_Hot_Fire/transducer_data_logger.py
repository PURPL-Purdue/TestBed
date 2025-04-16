import csv
import threading
import queue
import os
from datetime import datetime

class TransducerDataLogger:
    def __init__(self, transducers_list, devices_list, filename="Torch_Hot_Fire/transducer_log.csv"):
        self.filename = filename
        self.log_queue = queue.Queue()
        self.running = True
        self._transducers = transducers_list
        self._devices = devices_list

        # Create CSV file and write header if it doesn't exist
        if not os.path.exists(self.filename):
            with open(self.filename, mode='w', newline='') as file:
                entry = ["Timestamp"]
                for transducer in self._transducers:
                    entry.append(f"{transducer.input_channel_1} Pressure")
                for device in self._devices:
                    entry.append(f"{device.name} State")
                writer = csv.writer(file)
                writer.writerow(entry)

        # Start the background thread to handle writing to the CSV
        self.thread = threading.Thread(target=self._process_queue)
        self.thread.start()

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

    def stop(self):
        """Stop the logging and gracefully shut down the background thread."""
        self.running = False
        self.thread.join()  # Wait for the thread to finish processing the remaining logs
