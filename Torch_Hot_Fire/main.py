import sys
from PyQt5 import QtWidgets
from Interface.FluidPanel import MainWindow
from Interface.SolenoidPanel import SolenoidWindow

# Main application
app = QtWidgets.QApplication(sys.argv)
main_window = MainWindow()
main_window.show()
solenoid_window = SolenoidWindow(main_window)
solenoid_window.show()
sys.exit(app.exec_())