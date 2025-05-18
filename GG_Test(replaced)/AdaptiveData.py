# Adaptive Memory Buffer
# Ryan Wans

# Dynamically limit the size of arrays stored in memory

class Buffer:
    def __init__(self, max_size):
        self.max_size = max_size
        self.data = []

    def append(self, value):
        if len(self.data) >= self.max_size:
            self.data.pop(0)
        self.data.append(value)

    def get(self):
        return self.data

    def clear(self):
        self.data = []

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        self.data[index] = value

    def __delitem__(self, index):
        del self.data[index]

    def __len__(self):
        return len(self.data)
    
    def __iter__(self):
        return iter(self.data)

    def __repr__(self):
        return f"Buffer(max_size={self.max_size}, data={self.data})"
