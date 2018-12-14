import sys
import platform
import numpy as np
from temp import massflow, massflow2, massflow3, massflow4, walltemp
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
        self.edit4 = QLineEdit("0.0254")

        self.edit5 = QLineEdit("0")
        self.edit6 = QLineEdit("5")
        
        self.edit7 = QLineEdit("0.7")
        self.edit8 = QLineEdit("3")
        
        self.edit9 = QLineEdit("0.101325")
        
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
        label = QLabel("Flow Diameter [m]:",self)
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
  
        hbox6 = QHBoxLayout()
        label = QLabel("Operating Pressure [MPa]:",self)
        hbox6.addWidget(label)
        hbox6.addWidget(self.edit9)
    
        
        hbox = QHBoxLayout()
        hbox.addWidget(self.runbutton)
        
        self.runbutton.clicked.connect(self.update)
        self.edit1.returnPressed.connect(self.update)
       

        layout = QVBoxLayout()
        
        layout.addWidget(self.plot)

        layout.addLayout(hbox0)
        layout.addLayout(hbox1)
        layout.addLayout(hbox2)
        layout.addLayout(hbox3)
        layout.addLayout(hbox4)
        layout.addLayout(hbox5)
        layout.addLayout(hbox6)
        
        layout.addLayout(hbox)
        
        
        self.widget.setLayout(layout)     
        self.widget.setFixedSize(1100,850)
        self.setCentralWidget(self.widget) 
        
    def saveas(self):
        """Save something
        
        Hint: look up QFileDialog.getSaveFileName.
        """
        name = QFileDialog.getSaveFileName(self, 'Save File')
        tempstr = name[0].split('/')
        np.savetxt(str(tempstr[-1]),np.c_[self.x,self.y, self.x2, self.y2, self.x3, self.y3, self.x4, self.y4, self.Twall, self.Tbulk, self.z, self.Twall_F, self.Tbulk_F
                         ,self.x_thermo, self.x_levy,self.x_thermo_F, self.x_levy_F, self.alpha_levy, self.alpha_levy_F, self.alpha_levy_M])
                
    def about(self):
        QMessageBox.about(self, 
            "Two-Phase Flow Models",
            """<b>Two-Phase Flow Models</b>
               <p>Python %s -- Qt %s -- PyQt %s on %s""" %
            (platform.python_version(),
             QT_VERSION_STR, PYQT_VERSION_STR, platform.system()))

    def update(self):
        """Update the figure title.
        
        Of course, this is trivial(ish), but it serves as a guid for how
        to update other parts of the figure (see the reference for more).
        """
        Temp_in = float(self.edit1.text())
        Length = float(self.edit2.text())
        K_loss = float(self.edit3.text())
        Diam = float(self.edit4.text())
        q_in = float(self.edit5.text())
        q_out = float(self.edit6.text())
        heated_length = float(self.edit7.text())
        power_in = float(self.edit8.text())
        press = float(self.edit9.text())
        self.x, self.y = massflow(Tin=Temp_in,L=Length,K=K_loss,d=Diam,q_start=q_in,q_end=q_out,Pin=press)
        self.x2, self.y2= massflow2(Tin=Temp_in,L=Length,K=K_loss,d=Diam,q_start=q_in,q_end=q_out,Pin=press)
        self.x3, self.y3= massflow3(Tin=Temp_in,L=Length,K=K_loss,d=Diam,q_start=q_in,q_end=q_out,Pin=press)
        self.x4, self.y4= massflow4(Tin=Temp_in,L=Length,K=K_loss,d=Diam,q_start=q_in,q_end=q_out,Pin=press)
        fraction = int(((power_in-q_in)/(q_out-q_in))*100)-1
        
        self.Twall,     self.Tbulk, self.x_thermo, self.x_levy, self.alpha_levy ,self.z = walltemp(Tin=Temp_in, q_in = power_in, z=heated_length, m_dot = self.y[fraction])
        self.Twall_F, self.Tbulk_F, self.x_thermo_F, self.x_levy_F, self.alpha_levy_F,self.z = walltemp(Tin=Temp_in, q_in = power_in, z=heated_length, m_dot = self.y3[fraction])
        self.Twall_M, self.Tbulk_M, self.x_thermo_M, self.x_levy_M, self.alpha_levy_M,self.z = walltemp(Tin=Temp_in, q_in = power_in, z=heated_length, m_dot = self.y4[fraction])
        
        
        self.plot.redraw(self.x,self.y, self.x2, self.y2, self.x3, self.y3, self.x4, self.y4, self.Twall, self.Tbulk, self.z, self.Twall_F, self.Tbulk_F
                         ,self.x_thermo, self.x_levy,self.x_thermo_F, self.x_levy_F, self.alpha_levy, self.alpha_levy_F, self.alpha_levy_M)

    def clear(self):
        self.plot.redraw(0,0)
        self.edit2.setText(str([])) 
    

        

class MatplotlibCanvas(FigureCanvas) :
    """ This is borrowed heavily from the matplotlib documentation;
        specifically, see:
        
    """
    def __init__(self):
        
        # Initialize the figure and axes
        self.fig = Figure()
        self.axes = self.fig.add_subplot(221)
        self.axes2 = self.fig.add_subplot(222)
        self.axes3 = self.fig.add_subplot(224)
        self.axes4 = self.fig.add_subplot(223)
        
        # Now do the initialization of the super class
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
         
        
    def redraw(self, x, y, x2, y2, x3, y3, x4, y4, Twall, Tbulk, z, Twall_F, Tbulk_F, 
               x_thermo,x_levy, x_thermo_F, x_levy_F, alpha_levy, alpha_levy_F, alpha_levy_M):
        """ Redraw the figure with new x and y values.
        """
        # clear the old image (axes.hold is deprecated)
        self.axes.clear()
        self.axes.plot(x, y,'c')
        self.axes.plot(x2,y2,'k')
        self.axes.plot(x3,y3,'g')
        self.axes.plot(x4,y4,'r')
        self.axes.set_xlabel('Heat Input [kW]')
        self.axes.set_ylabel('Mass Flow Rate [kg/s]')
        self.axes.grid(which='both',axis='both')
        self.axes.legend(['Hom. 1','Hom. 2','Friedel','Lockhart & Martinelli'])

        self.axes2.clear()
        self.axes2.plot(z,Tbulk,'b')
        self.axes2.plot(z,Twall,'r--')
        self.axes2.plot(z,Tbulk_F,'r')
        self.axes2.plot(z,Twall_F,'b--') 
        self.axes2.set_xlabel("Wall Distance [m]")
        self.axes2.set_ylabel("Temperature [C]")
        self.axes2.yaxis.set_label_position("right")
        self.axes2.yaxis.tick_right()
        self.axes2.grid(which='both',axis='both')
        self.axes2.legend(['Bulk Temp.-Hom. 1','Wall Temp.-Hom. 1', 'Bulk Temp.-Friedel','Wall Temp.-Friedel'])
        
        self.axes3.clear()
        self.axes3.plot(z,x_thermo,'r')
        self.axes3.plot(z,x_levy,'r--')
        self.axes3.plot(z,x_thermo_F,'b')
        self.axes3.plot(z,x_levy_F,'b--')
        self.axes3.set_xlabel("Wall Distance [m]")
        self.axes3.set_ylabel('Quality')
        self.axes3.yaxis.set_label_position("right")
        self.axes3.yaxis.tick_right()
        self.axes3.grid(which='both',axis='both')
        self.axes3.legend(['thermo. x-Hom. 1','Levy x-Hom. 1','thermo. x-Friedel','Levy x-Friedel'])
  
        self.axes4.clear()
        self.axes4.plot(z,alpha_levy,'c')
        self.axes4.plot(z,alpha_levy_F,'g')
        self.axes4.plot(z,alpha_levy_M,'r')        
        self.axes4.set_xlabel("Wall Distance [m]")
        self.axes4.set_ylabel('Void Fraction')
        self.axes4.grid(which='both',axis='both')
        self.axes4.legend(['Levy-Hom. 1','Levy-Friedel','Levy-Lockhart & Martinelli'])
    
        self.draw()    
    def styleChoice(self,text):
        self.stylechoice.setText(text)
        QApplication.setStyle(QStyleFactory.create(text))
        
        
app = QApplication(sys.argv)
form = MainWindow()
form.show()
app.exec_()
