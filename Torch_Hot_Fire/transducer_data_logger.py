import csv
import threading
import queue
import os
from datetime import datetime

class TransducerDataLogger:
    def __init__(self, filename="transducer_log.csv"):
        self.filename = filename
        self.log_queue = queue.Queue()
        self.running = True

        # Create CSV file and write header if it doesn't exist
        if not os.path.exists(self.filename):
            with open(self.filename, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(["Timestamp", "Transducer1 Pressure", "Transducer2 Pressure", "Transducer3 Pressure"])

        # Start the background thread to handle writing to the CSV
        self.thread = threading.Thread(target=self._process_queue)
        self.thread.start()

    def log_data(self, pressure1, pressure2, pressure3):
        """Put the log data in the queue for the background thread to process."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.log_queue.put((timestamp, pressure1, pressure2, pressure3))

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
                # If timeout occurs, just continue the loop
                continue

    def stop(self):
        """Stop the logging and gracefully shut down the background thread."""
        self.running = False
        self.thread.join()  # Wait for the thread to finish processing the remaining logs
