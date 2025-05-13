import sys
from PyQt5 import QtWidgets
from Interface.MainPanel import MainWindow

def main():
    # Main application
    app = QtWidgets.QApplication(sys.argv)
    
    # Create the main window (fluid panel)
    main_window = MainWindow()
    
    # This ensures that the application won't quit until all windows are closed
    app.setQuitOnLastWindowClosed(True)
    
    main_window.show()
    
    return app.exec_()

if __name__ == "__main__":
    sys.exit(main())