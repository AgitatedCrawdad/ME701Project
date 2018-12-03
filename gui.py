import sys
import platform
import numpy as np
from temp import massflow, walltemp
#import matplotlib
#matplotlib.use('Qt5Agg')

from PyQt5.QtWidgets import (QMainWindow, QApplication, QDialog, QLineEdit, 
                             QVBoxLayout, QAction, QMessageBox,QFileDialog,
                             QSizePolicy, QWidget, QComboBox, QInputDialog,
                             QPushButton,QHBoxLayout,QLabel,)
from PyQt5.QtCore import QT_VERSION_STR, PYQT_VERSION_STR
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure
import math

class MainWindow(QMainWindow) :
    
    def __init__(self, parent=None) :
        super(MainWindow, self).__init__(parent)

        ########################################################################

        # ADD MENU ITEMS
        
        ########################################################################
        
        # Create the File menu
        self.menuFile = self.menuBar().addMenu("&File")
        self.actionSaveAs = QAction("&Save As", self)
        self.actionSaveAs.triggered.connect(self.saveas)
        self.actionQuit = QAction("&Quit", self)
        self.actionQuit.triggered.connect(self.close)
        self.menuFile.addActions([self.actionSaveAs, self.actionQuit])
    
        # Create the Help menu
        self.menuHelp = self.menuBar().addMenu("&Help")
        self.actionAbout = QAction("&About",self)
        self.actionAbout.triggered.connect(self.about)
        self.menuHelp.addActions([self.actionAbout])
        
        ########################################################################
        # CREATE CENTRAL WIDGET
        ########################################################################

        self.widget = QDialog()
        self.plot = MatplotlibCanvas()
        
        self.edit1 = QLineEdit("99")
        self.edit2 = QLineEdit("1.4")
        
        self.edit3 = QLineEdit("5")
        self.edit4 = QLineEdit("0.0005")

        self.edit5 = QLineEdit("0")
        self.edit6 = QLineEdit("5")
        
        self.edit7 = QLineEdit("0.7")
        self.edit8 = QLineEdit("3")
        
        self.runbutton = QPushButton("Run")
        self.cancelbutton = QPushButton("Clear")
        hbox0 = QHBoxLayout()
        label = QLabel("Mass Flow Rate Settings",self)
        label.setAlignment(Qt.AlignCenter)
        hbox0.addWidget(label)
        
        hbox1 = QHBoxLayout()
        label = QLabel("Inlet Temperature [C]:",self)
        hbox1.addWidget(label)
        hbox1.addWidget(self.edit1)
        label = QLabel("Riser Height [m]:",self)
        hbox1.addWidget(label)
        hbox1.addWidget(self.edit2)

        hbox2 = QHBoxLayout()
        label = QLabel("Loss Coefficient []:",self)
        hbox2.addWidget(label)
        hbox2.addWidget(self.edit3)
        label = QLabel("Flow Area [m^2]:",self)
        hbox2.addWidget(label)
        hbox2.addWidget(self.edit4)
        
        hbox3 = QHBoxLayout()
        label = QLabel("Starting Power [kW]:",self)
        hbox3.addWidget(label)
        hbox3.addWidget(self.edit5)
        label = QLabel("Ending Power [kW]:",self)
        hbox3.addWidget(label)
        hbox3.addWidget(self.edit6)
        
        hbox4 = QHBoxLayout()
        label = QLabel("Wall Temperature Settings",self)
        label.setAlignment(Qt.AlignCenter)
        hbox4.addWidget(label)       
        
        hbox5 = QHBoxLayout()
        label = QLabel("Heated Length [m]:",self)
        hbox5.addWidget(label)
        hbox5.addWidget(self.edit7)
        label = QLabel("Heat in [kW]:",self)
        hbox5.addWidget(label)
        hbox5.addWidget(self.edit8)
        
        
        hbox = QHBoxLayout()
        hbox.addWidget(self.runbutton)
#        hbox.addWidget(self.cancelbutton)
        # signals + slots ()
        
        self.runbutton.clicked.connect(self.update)
#        self.cancelbutton.clicked.connect(self.clear)
        self.edit1.returnPressed.connect(self.update)
        
#        self.comboBox = QComboBox(self)
#        self.comboBox.addItems(['x**(1/2)','x**2','np.sin(x)',''])
#        self.comboBox.setEditable(True)
        

        layout = QVBoxLayout()
        
        layout.addWidget(self.plot)
#        layout.addWidget(self.comboBox)
#        layout.addWidget(self.edit1)
#        layout.addWidget(self.edit2) 
        layout.addLayout(hbox0)
        layout.addLayout(hbox1)
        layout.addLayout(hbox2)
        layout.addLayout(hbox3)
        layout.addLayout(hbox4)
        layout.addLayout(hbox5)
        
        layout.addLayout(hbox)
        
        
        self.widget.setLayout(layout)        
        self.setCentralWidget(self.widget) 
    def saveas(self):
        """Save something
        
        Hint: look up QFileDialog.getSaveFileName.
        """
        name = QFileDialog.getSaveFileName(self, 'Save File')
        tempstr = name[0].split('/')
#        print(str(tempstr)c)
        np.savetxt(str(tempstr[-1]),np.c_[self.x,self.y])
#        file = open(str(tempstr[-1]),'w')
##        text = self.textEdit.toPlainText()
#        file.write(str(self.x.T))
#        file.close()
                
    def about(self):
        QMessageBox.about(self, 
            "About Function Evaluator",
            """<b>Function Evaluator</b>
               <p>Copyright &copy; 2016 Jeremy Roberts, All Rights Reserved.
               <p>Python %s -- Qt %s -- PyQt %s on %s""" %
            (platform.python_version(),
             QT_VERSION_STR, PYQT_VERSION_STR, platform.system()))

    def update(self):
        """Update the figure title.
        
        Of course, this is trivial(ish), but it serves as a guid for how
        to update other parts of the figure (see the reference for more).
        """
#        title = str(self.edit1.text())
#        self.plot.axes.set_title(title)
        Temp_in = float(self.edit1.text())
        Length = float(self.edit2.text())
        K_loss = float(self.edit3.text())
        Area = float(self.edit4.text())
        q_in = float(self.edit5.text())
        q_out = float(self.edit6.text())
        heated_length = float(self.edit7.text())
        power_in = float(self.edit8.text())
#        var = 'self.x='
#        if span[0].isdigit():
#            totalspan = var+'np.array(['+span+'])'
#        else:
#            totalspan = var+span
#        
#        exec(totalspan)
        self.x, self.y = massflow(Tin=Temp_in,L=Length,K=K_loss,A=Area,q_start=q_in,q_end=q_out)
        
        fraction = int(((power_in-q_in)/(q_out-q_in))*1000)
        
        self.Twall, self.Tbulk, self.z = walltemp(q_in = power_in, z=heated_length, m_dot = self.y[fraction])
        
#        fun = str(self.comboBox.currentText())
      
#        fun = fun.replace('x','self.x')
        
#        self.y = eval(fun)
#        x = np.linspace(0,10)
#        y = x**2
#        self.plot.draw()
        
        self.plot.redraw(self.x,self.y, self.Twall, self.Tbulk, self.z)
#        self.edit2.setText(str(self.y))

    def clear(self):
        self.plot.redraw(0,0)
        self.edit2.setText(str([])) 
    

        
class Form(QDialog) :

    def __init__(self, parent=None) :
        super(Form, self).__init__(parent)
        self.function_edit = QLineEdit("x**2")
        self.function_edit.selectAll()
        self.parameter_edit = QLineEdit("np.linspace(0,1,4)")
        self.parameter_edit.selectAll()
        self.output_edit = QLineEdit(" ")
        self.output_edit.selectAll()
        self.plot = MatplotlibCanvas()
        layout = QVBoxLayout()
        layout.addWidget(self.plot)
        layout.addWidget(self.function_edit)
        layout.addWidget(self.parameter_edit)     
        layout.addWidget(self.output_edit)  
        self.setLayout(layout)
        self.function_edit.setFocus()
        self.output_edit.returnPressed.connect(self.updateUi)
        self.setWindowTitle("Function Evaluator")
        self.x = None
        self.f = None

    def updateUi(self) :
        #try : 
            self.x = str(self.parameter_edit.text())
            x = eval(self.x) 
            if len(x) > 1 :
                x = np.array(x)
            # Is there a cleaner way?
            f = eval(str(self.function_edit.text()))
            self.f = str(f)
            self.f = self.f.replace("[","").replace("]","")
            self.f = ",".join(self.f.split())
            self.output_edit.setText(self.f)
            
            self.plot.redraw(x, f)
            #self.plot.axes.plot(x, f)
            #self.plot.draw()
        #except :
        #    self.output_edit.setText("error! check function or parameter.")
                

class MatplotlibCanvas(FigureCanvas) :
    """ This is borrowed heavily from the matplotlib documentation;
        specifically, see:
        
    """
    def __init__(self):
        
        # Initialize the figure and axes
        self.fig = Figure()
        self.axes = self.fig.add_subplot(221)
        self.axes2 = self.fig.add_subplot(222)
        # Give it some default plot (if desired).  
#        x = np.arange(0.0, 3.0, 0.01)
#        y = np.sin(2*np.pi*x)
#        self.axes.plot(x, y)
#        self.axes.set_xlabel('Heat Input [kW]')
#        self.axes.set_ylabel('Mass Flow Rate [kg/s]')   
        
        # Now do the initialization of the super class
        FigureCanvas.__init__(self, self.fig)
        #self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
         
        
    def redraw(self, x, y,Twall, Tbulk,z) :
        """ Redraw the figure with new x and y values.
        """
        # clear the old image (axes.hold is deprecated)
        self.axes.clear()
        self.axes.plot(x, y)
        self.axes.set_xlabel('Heat Input [kW]')
        self.axes.set_ylabel('Mass Flow Rate [kg/s]')
        self.axes.grid(which='both',axis='both')

        self.axes2.clear()
        self.axes2.plot(z,Tbulk)
        self.axes2.plot(z,Twall)
        self.axes2.set_xlabel("Wall Distance [m]")
        self.axes2.set_ylabel("Temperature [C]")
        self.axes2.grid(which='both',axis='both')
        
        self.draw()    
    def styleChoice(self,text):
        self.stylechoice.setText(text)
        QApplication.setStyle(QStyleFactory.create(text))
        
        
app = QApplication(sys.argv)
form = MainWindow()
form.show()
app.exec_()
