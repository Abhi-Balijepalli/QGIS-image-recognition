# -*- coding: utf-8 -*-
"""
/***************************************************************************
 ReedCanaryProject
                                 A QGIS plugin
 2021-2022 CS46x Capstone project
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2022-02-01
        git sha              : $Format:%H$
        copyright            : (C) 2022 by Logan, Abhi, Nuocheng
        email                : Logan@kleditz.com
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
# ----===USER SETTINGS===----

#PATH FOR THE TEMPORARY LAYERS (string (must be a valid directory with a '/' at the end; example: "c:/rcg/"))
 #rC:\Users\logan\Desktop\Stuff_for_Capstone" + "/" # add an input that is the file directory 
#print(directory)
#This directory will determine where the temporary layers are placed.
#THIS MUST BE SET TO A VALID DIRECTORY for this script to run.
#However, the directory does not need to already exist, as the script will create a new directory if it's not already there

#ROI CLASS IDENTIFIER (string)
roi_class_identifier = "class-id"
#This string is used to determine the class of each ROI, the name of the field (attribute), whose value is set to a number for the class of the ROI
#This will be set to the attribute name that identifies the roi feature within the mle-roi shapefile

#TREATMENT AREA IDENTIFIER (string)
treatment_area_identifier = "PlotName"
#This script was made to iterate over treatment areas and provide results for each of these treatments,
#This will be set to the attribute name that identifies the treatment area feature within the treatment-areas shapefile

#ROI VALIDIATION PROPORTION (min = 0, max = 1)
roi_validation_proportion = 0.33
#Determines the proption of ROIs in each class that will be selected for validation
#Example: If roi_validation_proportion == 0.33, then 1/3 of the ROIs will be used for
#validation, and the remaining 2/3 will be used for classifiction

#ROI VALIDATION SELECTION PSEUDO-RANDOM (True/False)
roi_validation_pseudo_random = True
#If True, the random selection will only be re shuffled when ROIs are added or removed
#If False, every time the program is run (even if the ROIs are not changed) the selection
#of which are used for validation/classification will be shuffled again, which will yield different results

#MINIMUM PROBABILITY PER BAND (min = 0, max = 1)
mle_minimum_probability = 0.001
#Upon iteration of the Normal Distribution Probability Density Function (NDPDF) at all the bands within 1 pixel, the determined
#probability (equal to (the result of NDPDF) ^ (the number of bands)) must be greater than (mle_minimum_probability) ^ (the number of bands)
#otherwise, the pixel will be classified as "unclassified"

temp = []

#--=ADDITIONAL DEBUG OPTIONS=--

#PIXEL READ LOOP MAX (min = 0, max = infinity)
pixel_read_loop_max = 10000
#IF NOT DEBUGGING, LEAVE THIS AT 0
#If this set to 0, all pixels will be read as normal
#If this is set to above 0, only that many pixels will be read before the algorithm stops early.

#PIXEL READ PROPORTION (min = 0, max = 1)
pixel_read_proportion = 1
#IF NOT DEBUGGING, LEAVE THIS AT 1
#If this is set to a value below 1, only that proportion of pixels will be read,
#When this is below 1, pixels can also be skipped for validation, so this really messes with the validation results


#--==end of user settings==--

from qgis.gui import QgsMessageBar
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction, QFileDialog, QDockWidget
from qgis.core import * #QgsProject, Qgis, QgsFields, QgsField
from qgis import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from qgis.utils import iface
from qgis.analysis import QgsZonalStatistics
from PIL import Image as im
import numpy as np
from qgis.core import *
from matplotlib.pyplot import imshow
import os
from qgis.PyQt.QtCore import QCoreApplication, QLocale, QThread, qDebug
from qgis.PyQt.QtWidgets import QPushButton, QApplication
from qgis.core import *
from qgis.gui import QgsMessageBar
from qgis.utils import iface
from qgis.analysis import QgsZonalStatistics
from osgeo import gdal
import collections
import sys
# Initialize Qt resources from file resources.py
from .resources import *
# Import the code for the dialog
from .reed_canary_project_dialog import ReedCanaryProjectDialog
from .reed_canary_project_dialog import Test1Dialog
from .reed_canary_project_dialog import AlgoDialog
from .reed_canary_project_dialog import Image_SelectDialog
import os.path

#temp_list = ["T:\\Teach\\Classes\\CS461\\ReedCanary\\Images\\pix4D_2019-04-03_transparent_reflectance_w550nm.tif"]
items_list = [] #Probably can get rid of this since this feature was basically scrapped
items_name = [] #Same here


#filename = ""
class ReedCanaryProject:
    """QGIS Plugin Implementation."""

    Algorithm_use_list = []

    Data_file_list = []

    Training_file_list = []

    directory_save_location = ""

    location_Data = ""


    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'ReedCanaryProject_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)
            QCoreApplication.installTranslator(self.translator)

        self.dlg = ReedCanaryProjectDialog()
        self.dlg2 = Test1Dialog()
        self.dlg3 = AlgoDialog()
        self.dlg4 = Image_SelectDialog()

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Reed Canary Project')

        # Check if plugin was started the first time in current QGIS session
        # Must be set in initGui() to survive plugin reloads
        self.first_start = None

        self.toolbar = self.iface.addToolBar(u'&Reed Canary Project')
        self.toolbar.setObjectName(u'&Reed Canary Project')


    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('ReedCanaryProject', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            # Adds plugin icon to Plugins toolbar
            self.iface.addToolBarIcon(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/reed_canary_project/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'ReedCanaryProject2021'),
            callback=self.run,
            parent=self.iface.mainWindow())

        # will be set False in run()
        self.first_start = True


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&Reed Canary Project'),
                action)
            self.iface.removeToolBarIcon(action)

    #def select_file4(self): #selects a file to open
        #filename, _filer = QFileDialog.getOpenFileName(self.dlg, "Open File", "nm", '*.tif')
        #self.dlg.lineEdit_4.setText(filename)

    def refresh_listView(self): #takes all the layers  and places them into the listView and are checked (checks will be used to make vraster) also some other stuff maybe IDFK
        qListView = QDockWidget()
        model = QStandardItemModel()
        temp_array = []

        for layer in QgsProject.instance().layerTreeRoot().children():  #gets the layers
            item = QStandardItem(layer.name())
            check = Qt.Checked
            item.setCheckState(check)
            item.setCheckable(True)
            model.appendRow(item)
            temp_array.append(layer.name())
            items_name.append(layer.name())
        for i in temp_array: #gets their path, this fucking sucks btw
            temp = QgsProject.instance().mapLayersByName(i)
            items_list.append(temp[0].dataProvider().dataSourceUri())


        self.dlg2.listView.setModel(model)
        self.dlg2.listView.show()


        # show the dialog 

    def check_array(self, algo_array):
        
        if collections.Counter(algo_array) == collections.Counter(self.Algorithm_use_list):
            return True
        else:
            return False

    def open_ui_2(self):

        self.dlg2.show() #shows the second UI

    def image_select_open(self): #allows the user to select the file that is the actual data, checks the location data
        #global Data_file_list 
        #global location_Data
        if len(self.Data_file_list) != 0:
            self.Data_file_list.clear()
        self.dlg4.show()
        qListView = QDockWidget()
        temp_2 = []
        model = QStandardItemModel()
        if len(self.Data_file_list) != 0:
            self.Data_file_list.clear
        for layer in QgsProject.instance().layerTreeRoot().children():  #gets the layers
            item = QStandardItem(layer.name())
            check = Qt.Checked
            #item.setCheckState(check)
            item.setCheckable(True)
            model.appendRow(item)
            temp_2.append(layer.name())

        self.dlg4.listView.setModel(model)
        self.dlg4.listView.show()
        result = self.dlg4.exec_()
        if len(temp_2) != 0:
            check = []
            for i in range(0, len(temp_2)):
                item = model.item(i)
                if item.checkState() == QtCore.Qt.Unchecked:
                    check.append(0)
                else:
                    check.append(1)
            for i in range(0, len(check)):
                if check[i] == 1:
                    self.Data_file_list.append(temp_2[i])
                else:
                    continue
            if len(self.Data_file_list) == 0:          
                print("ERROR: no data selected")
                self.Data_file_list.clear()
                return
        
        if len(self.Data_file_list) != 0:
            cnt = 0
            for i in self.Data_file_list: #gets the selected layers path
                temp = QgsProject.instance().mapLayersByName(i)
                if self.location_Data == "":
                    self.location_Data = temp[0].crs()
                else:
                    if self.location_Data != temp[0].crs():
                        self.iface.messageBar().pushMessage("Error", "Selected Data ESPG does not match", Qgis.Critical)
                        self.Data_file_list.clear()
                        return
                self.Data_file_list[cnt] = temp[0].dataProvider().dataSourceUri()
                cnt += 1
        else: #just clear the list and return as if nothing happened (I know this is shoddy code)
            self.Data_file_list.clear()
            return

       # if len(Data_file_list) != 0: 

        print("Data_FIle_LIST", self.Data_file_list)

    def trainingdata_select_open(self): #opens the dialog box that allows the user to select the training data file, and when they click ok it adds it to a list, and checks the location data 
       # global Training_file_list
        #global location_Data
        if len(self.Training_file_list) != 0:
            self.Training_file_list.clear()
        self.dlg4.show()
        qListView = QDockWidget()
        temp_2 = []
        model = QStandardItemModel()
        if len(self.Training_file_list) != 0:
           self.Training_file_list.clear
        for layer in QgsProject.instance().layerTreeRoot().children():  #gets the layers
            item = QStandardItem(layer.name())
            check = Qt.Checked
            #item.setCheckState(check)
            item.setCheckable(True)
            model.appendRow(item)
            temp_2.append(layer.name())

        self.dlg4.listView.setModel(model)
        self.dlg4.listView.show()
        result = self.dlg4.exec_()
        if len(temp_2) != 0:
            check = []
            for i in range(0, len(temp_2)):
                item = model.item(i)
                if item.checkState() == QtCore.Qt.Unchecked:
                    check.append(0)
                else:
                    check.append(1)
            for i in range(0, len(check)):
                if check[i] == 1:
                    self.Training_file_list.append(temp_2[i])
                else:
                    continue
            if len(self.Training_file_list) == 0:          
                print("ERROR: no training data selected")
                self.Training_file_list.clear()
                return
        
        if len(self.Training_file_list) != 0:
            cnt = 0
            for i in self.Training_file_list: #gets their path, this fucking sucks btw
                temp = QgsProject.instance().mapLayersByName(i)
                if self.location_Data == "":
                    self.location_Data = temp[0].crs()
                if self.location_Data != temp[0].crs():
                    self.iface.messageBar().pushMessage("Error", "Selected Data and Training Data ESPG does not match", Qgis.Critical)
                    self.Training_file_list.clear()
                    return
                self.Training_file_list[cnt] = temp[0].dataProvider().dataSourceUri()
                cnt += 1
        else:
            self.Training_file_list.clear()
            return
        print("Training_file_list", self.Training_file_list)

    def openAlgoPage(self):
        self.dlg3.show()
        qListView = QDockWidget()
        model = QStandardItemModel()
        algorithm_array = ["Maximum likelihood", "Random Forest"]
        if len(self.Algorithm_use_list) == 0:
            for i in range(0, len(algorithm_array)):
                item = QStandardItem(algorithm_array[i])
                check = Qt.Checked
                #item.setCheckState(check)
                item.setCheckable(True)
                model.appendRow(item)
                self.Algorithm_use_list.append(algorithm_array[i])
        else:
            self.Algorithm_use_list.clear()
            #print(Algorithm_use_list)
            for i in range(0, len(algorithm_array)):
                item = QStandardItem(algorithm_array[i])
                check = Qt.Checked
                #item.setCheckState(check)
                item.setCheckable(True)
                model.appendRow(item)
                self.Algorithm_use_list.append(algorithm_array[i])
   
        self.dlg3.listView.setModel(model)
        self.dlg3.listView.show()
        result = self.dlg3.exec_()
        if len(self.Algorithm_use_list) != 0:
            for i in range(0, len(self.Algorithm_use_list)):
                item = model.item(i)
                if item.checkState() == QtCore.Qt.Unchecked:
                    self.Algorithm_use_list.remove(self.Algorithm_use_list[i])
                    i = 0
        print("299",self.Algorithm_use_list)

    # def vraster(self): #makes the new vrt and loads it into the project
    #     model = self.dlg2.listView.model()
    #     print(items_list, "before removing things")
    #     for i in range(0,len(items_list)):
    #         print(i, items_list[i])
    #     temp = 0
    #     count = 0
    #     for index in range(model.rowCount()): #checks to see if any of the items have been unchecked, if they are unchecked remove them from the list that will be used
    #         item = model.item(index)
    #         print(item, index)
    #         if item.checkState() == QtCore.Qt.Unchecked:
    #                     print("item at index ", temp, items_list[temp])
    #                     items_list.remove(items_list[temp])
    #                     #items_name.remove(temp)
    #                     index = 0

            
        
    #     print("items that are checked", len(items_list), items_list)
    #     outfile = self.dlg2.lineEdit.text() + ".tif"
    #     processing.run("gdal:buildvirtualraster", 
    #     {'INPUT':items_list,'RESOLUTION':0,'SEPARATE':True,'PROJ_DIFFERENCE':False,'ADD_ALPHA':False,'ASSIGN_CRS':None,'RESAMPLING':0,'SRC_NODATA':'',
    #     'OUTPUT':outfile})
    #     time.sleep(1)
    #     #self.iface.addRasterLayer(outfile, "Virtual") #load the virtual raster into the project under the name Virtual
    #     rlayer = QgsRasterLayer(outfile, "data")
    #     if not rlayer.isValid():
    #         print("Layer failed to load!")

    def dir_path(self):
        #global directory_save_location
        if self.directory_save_location != "":
            self.directory_save_location = ""
        #if button == 2: #if the button is being pushed to retrieve a file (I do not think this does anything....)
        print("pressed button 2")
        filename, _filter = QFileDialog.getSaveFileName(  
        self.dlg2, "Select output location", "", '*.txt')  
        self.directory_save_location = filename
        temp = self.directory_save_location
        self.directory_save_location = ""
        temp_2 = ""
        temp = temp[::-1]
        #print("after making the new file and flip", temp)
        for elements in temp:
            
            if elements == "/":
                #print("found the first /", temp_2, temp)
                #print("temp_2", temp_2, "temp", temp)
                temp = temp[::-1]
                temp_2 = temp_2[::-1]
                #print("temp and temp_2 before the replace after changing to right side around", temp, temp_2)
                #print(len(temp), len(temp_2))
                temp.replace(temp_2, '', 1)
                #print("temp after replace", temp)
                try:
                    self.directory_save_location = temp
                    self.directory_save_location = self.directory_save_location[:len(self.directory_save_location) - len(temp_2)]
                except:
                    print("did not work")
                
                #print("end result", self.directory_save_location)
                break

            temp_2 = temp_2 + elements
        print(self.directory_save_location)




    def newshp(self):
        #global location_Data
        if self.location_Data == "":
            self.iface.messageBar().pushMessage("Error", "No data selected", Qgis.Critical)
            return 
        #This allows the user to search for a place to save a file, it also opens the directory traversal window
        if len(self.Data_file_list) == 0:
            self.iface.messageBar().pushMessage("Error", "No data selected", Qgis.Critical)
            return

        filename2, _filter = QFileDialog.getSaveFileName(  
        self.dlg4, "Select output filename and destination","mle-roi", '*.shp') 

            

        # create fields
        #Creates the field 'class-id' inside the new .shp file
        layerFields = QgsFields()
        layerFields.append(QgsField("class-id", QVariant.Int, "Integer", 20))
        print("layerfields", layerFields[0])

        #this actuall makes the new .shp layer with the Oregon projection, UTF-8  and with polygons as the roi

        writer = QgsVectorFileWriter(filename2, "UTF-8", layerFields, QgsWkbTypes.Polygon, QgsCoordinateReferenceSystem(self.location_Data), 'ESRI Shapefile')
        if writer.hasError() != QgsVectorFileWriter.NoError:
            print("Error when creating shapefile: ", writer.errorMessage())

        #adds new shp file to the layers
        vlayer = QgsVectorLayer(filename2, "mle-roi", "ogr")
        vlayer.dataProvider().addAttributes(layerFields)
        vlayer.updateFields()
        QgsProject.instance().addMapLayer(vlayer)
        QgsProject.instance().removeMapLayer(vlayer)
        vlayer = QgsVectorLayer(filename2, "mle-roi", "ogr")
        QgsProject.instance().addMapLayer(vlayer)
        del(writer)

    def run_algo(self):
        try:
            print(self.directory_save_location, self.Data_file_list, self.Training_file_list)
            print("run algo", self.Algorithm_use_list.index('Maximum likelihood', 0))
            #open(os.path.join(sys.path[0], "Some file.txt"), "r")
            #cwd = os.getcwd()
            #print()
            f = open(os.path.join(os.path.dirname(os.path.realpath(__file__)), "test.txt"), "a")
            print(f.name, os.getcwd())
            f.write(self.directory_save_location + "\n")
            for i in range(0, len(self.Data_file_list)):
                f.write(self.Data_file_list[i] + "\n")
            for i in range(0, len(self.Training_file_list)):
                f.write(self.Training_file_list[i]+ "\n")
            #f.close()
            #if Algorithm_use_list.index('Maximum likelihood', 0) >= 0:
            exec(open(self.plugin_dir + "/readmleroisours.py").read())
            #os.system('python readmleroisours.py')
            #os.system("python3 ./readmleroisours.py {self.directory_save_location}")
            #print("afyer the os.system")
            #else:
                #print("Did not have Maximum likelihood algorthm selected, go to the algorithm tab and select it now")
        except Exception as e:
            print("Error occured when trying to run the Maximum likelihood algorithm", e)
        

    def run(self,checked):
        global directory_save_location
        """Run method that performs all the real work"""

        # Create the dialog with elements (after translation) and keep reference
        # Only create GUI ONCE in callback, so that it will only load when the plugin is started
        if self.first_start == True:
            self.first_start = False
            self.dlg = ReedCanaryProjectDialog()
            self.dlg.pushButton.clicked.connect(self.open_ui_2)
            # self.dlg2.pushButton_3.clicked.connect(self.vraster)
            self.dlg.pushButton_3.clicked.connect(self.run_algo)
            self.dlg.pushButton_2.clicked.connect(self.openAlgoPage)
            self.dlg2.pushButton_4.clicked.connect(self.newshp)
            self.dlg2.pushButton.clicked.connect(self.image_select_open) #makes the new shp
            self.dlg2.pushButton_2.clicked.connect(self.trainingdata_select_open)
            self.dlg2.pushButton_5.clicked.connect(self.dir_path)

        # items = self.refresh_listView()

        
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            if directory_save_location == "":
                print("Did not select a directory to save the output from this plugin")
                #sys.exit()
            filename = directory_save_location
            pass