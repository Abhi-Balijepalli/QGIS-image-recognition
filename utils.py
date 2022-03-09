# This file is the pre-requirements for running the real-ml-rois.py file

# ----===USER SETTINGS===----

#PATH FOR THE TEMPORARY LAYERS (string (must be a valid directory with a '/' at the end; example: "c:/rcg/"))
directory = "T:\Teach\Classes\CS461\ReedCanary" + "/" # add an input that is the file directory 
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
pixel_read_loop_max = 100
#IF NOT DEBUGGING, LEAVE THIS AT 0
#If this set to 0, all pixels will be read as normal
#If this is set to above 0, only that many pixels will be read before the algorithm stops early.
#PIXEL READ PROPORTION (min = 0, max = 1)
pixel_read_proportion = 1
#IF NOT DEBUGGING, LEAVE THIS AT 1
#If this is set to a value below 1, only that proportion of pixels will be read,
#When this is below 1, pixels can also be skipped for validation, so this really messes with the validation results

#--==end of user settings==--

#Imports, etc.
from PIL import Image as im
import urllib.request
import random
import numpy as np
from qgis.core import *
from matplotlib.pyplot import imshow
import os
import math
import time
from qgis.PyQt.QtCore import QCoreApplication, QLocale, QThread, qDebug
from qgis.PyQt.QtWidgets import QPushButton, QApplication
from qgis.core import *
from qgis.gui import QgsMessageBar
from qgis.utils import iface
from qgis.analysis import QgsZonalStatistics
from osgeo import gdal
gdal.AllRegister()

#Create the directory if it doesn't already exist
if not os.path.exists(directory):
    os.mkdir(directory)
if not os.path.exists(directory+"temp/"):
    os.mkdir(directory+"temp/")
if not os.path.exists(directory+"temp/shapes"):
    os.mkdir(directory+"temp/shapes")
if not os.path.exists(directory+"temp/visual_output"):
    os.mkdir(directory+"temp/visual_output")

#These 4 variables will all later be set the their apporopriate shapefiles
#An empty string is only used for a placeholder until the script can find the layers
roi_shapefile = ""
merged_roi_shapefile_classification = ""
merged_roi_shapefile_validation = ""
treatment_areas = ""

#this will be used to store the statistics calculated on the "merged-rois-classification" shapefile, for faster access in the MLE calculation
class_statistics_dict = {}
#These indeces are used as keys in the class_statistics_dict dictionary
index_mean = 0
index_stdev = 1


#This dictonary will be used to sort the ROIs from the "mle-roi" shapefile,
#so that some can be used for training, and some for validation
class_roi_dict = {}
#These indeces are used as keys in the class_roi_dict dictionary
index_training_roi = 1
index_validation_roi = 2

#The numbers will be used in the Normal_Distribution_Probability_Density function
euler = math.e
pi = math.pi

#Make a string to track all of the missing layers as an error message. The program will only run if this reamins empty
missing_layers_error_message = ""

#Read and set variables that were unset above
#the layers variable will be used to access all the map layers in the current QGIS project
layers = QgsProject.instance().mapLayers().values()

#This function will return a list of unique colors from a number n.
#The colors returned will be as different as possible, given n
#This function could definetley be shortned to a fraction of the lines
unique_roi_ids = []
class_colors = None
mle_adjusted_minimum_probability = None
num_rois = 0
raster_pixel_scale_x = None
raster_pixel_scale_y = None

classified_pixel_dictionary = {}
output_image_x_min = 9999999999
output_image_x_max = 0
output_image_y_min = 9999999999
output_image_y_max = 0
cm_characters_per_entry = 5
confusion_matrix = {}

def SpaceText(val,num_chars,first_char,last_char):
    text = str(val)
    for i in range(len(text),num_chars):
        if i%2 == 0:
            text = text+" "
        else:
            text = " "+text
    return first_char+text+last_char

def Normal_Distribution_Probability_Density(x,mean,sd):
    return (1/(sd*(math.pow(2*pi,0.5))))*math.pow(euler,(-1/2)*math.pow(((x-mean)/(sd)),2))