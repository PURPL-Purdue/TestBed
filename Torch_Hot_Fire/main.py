import sys
from PyQt5 import QtWidgets
from Interface.FluidPanel import MainWindow

def main():
    # Main application
    app = QtWidgets.QApplication(sys.argv)
    
    # Create the main window (fluid panel)
    fluid_window = MainWindow()
    
    # # Link the windows properly for shutdown
    # fluid_window.valve_window = solenoid_window
    
    # This ensures that the application won't quit until all windows are closed
    app.setQuitOnLastWindowClosed(True)
    
    # Show both windows
    fluid_window.show()
    # solenoid_window.show()
    
    return app.exec_()

if __name__ == "__main__":
    sys.exit(main())