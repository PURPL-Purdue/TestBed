import pandas as pd
import datetime as dt

class TransducerDataLogger:
    def __init__(self):
        self.df = pd.DataFrame(columns=["Time", "Transducer 1", "Transducer 2", "Transducer 3"])

    def log_data(self, transducer1, transducer2, transducer3):
        new_row = pd.DataFrame(
            [[dt.datetime.now(), transducer1, transducer2, transducer3]],
            columns=self.df.columns
        )
        self.df = pd.concat([self.df, new_row], ignore_index=True)

    def export_to_csv(self, filename="Transducer Data.csv"):
        self.df.to_csv(filename, index=False)
