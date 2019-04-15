# -*- coding: utf-8 -*-
"""
/***************************************************************************
 LosAnalyzer
                                 A QGIS plugin
 Analyzes the Visibility of UAV from a given GCS
                              -------------------
        begin                : 2019-01-02
        git sha              : $Format:%H$
        copyright            : (C) 2019 by Jonas Hener - Avy BV
        email                : jonas@avy.eu
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
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from PyQt4.QtGui import QAction, QIcon, QFileDialog
from qgis.core import *
from qgis.core import QgsMapLayerRegistry
from osgeo import gdal, ogr
from qgis.gui import QgsMessageBar
import struct
import sys
import processing
import numpy as np
from math import *
from osgeo.gdalconst import *
# Initialize Qt resources from file resources.py
import resources
# Import the code for the dialog
from los_analyzer_dialog import LosAnalyzerDialog
import os.path


class LosAnalyzer:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):


        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgisInterface
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
            'LosAnalyzer_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        self.dlg = LosAnalyzerDialog()
        
        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&LOS Analyzer ')
        # TODO: We are going to let the user set this up in a future iteration
        self.toolbar = self.iface.addToolBar(u'LosAnalyzer')
        self.toolbar.setObjectName(u'LosAnalyzer')
        #Clear Output File text and connect PushButton
        self.dlg.OutputText.clear()
        self.dlg.OutputButton.clicked.connect(self.select_output_file)

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
        return QCoreApplication.translate('LosAnalyzer', message)


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
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/LosAnalyzer/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Line of Sight Analysis'),
            callback=self.run,
            parent=self.iface.mainWindow())


    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&LOS Analyzer '),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    def select_output_file(self):
        
        filename = QFileDialog.getSaveFileName(self.dlg, "Select output file ","", '*.tif')
        self.dlg.OutputText.setText(filename)

    def getPointGeo(self, myLayers):
        #assign Vector Layer from dlg
        vectorLayerIndex = self.dlg.VectorBox.currentIndex()
        vectorLayer = myLayers[vectorLayerIndex]
        features = vectorLayer.getFeatures()
        for feature in features:
            geom = feature.geometry()
            # show some information about the feature
            if geom.type() == QGis.Point:
                x = geom.asPoint()
                geoX = int(x[0])
                geoY = int(x[1])
                return geoX, geoY
            else:
                print "No Point Geometry Feature"

    def geoToInd(self, geoX, geoY, rasterGT):
        #only works for geotransforms with no rotation.
        indX = int((geoX - rasterGT[0]) / rasterGT[1]) #x pixel
        indY = int((geoY - rasterGT[3]) / rasterGT[5]) #y pixel

        return indY, indX #reverse order intentional, don't ask why...
            
            
    def cropInds(self, rasterGT, indX, indY):
        # crop to box of dimesions radius
        r = int(self.dlg.RadiusBox.text())*1000  #input from dlg in km
        rPix = r/rasterGT[1]
        loX = indX-rPix
        hiX = indX+rPix+1
        loY = indY-rPix
        hiY = indY+rPix+1
            
        return loX, hiX, loY, hiY
    
    def drawFrame(self, fs):
        edgeT = []
        edgeR = []
        edgeB = []
        edgeL = []
        
 

        for i in range(fs):
            edgeT.append([0, i])
            edgeR.append([i, fs-1])
            edgeB.append([fs-1, i])
            edgeL.append([i, 0])
                
        frame = [edgeT, edgeR, edgeB, edgeL]
        
    
        return frame


    def rasterised_line(self, x, y, x2, y2, interpolation = True, crop=0):

        dx = abs(x2 - x); dy = abs(y2 - y)
        steep = (dy > dx)
        #direction of movement : plus or minus in the coord. system
        sx = 1 if (x2 - x) > 0 else -1
        sy = 1 if (y2 - y) > 0 else -1
    
        if steep: # if the line is steep: change x and y
            #x,y = y,x they are the same !!
        
            dx,dy = dy,dx
            sx,sy = sy,sx

        D = 0
    
        #for interpolation
        # slope = dy / dx *sx *sy #!!
        #the operators for negative quadrants (we do need dx, dy as absolute to verify steepness, to search distances...)
    
        dx_short = dx-crop
    
        #store indices 1) los, 2) neighbours, 3) error
        mx_line = np.zeros((dx_short, 2), dtype=int)
    
        if interpolation:
            mx_neighbours = np.zeros((dx_short,2), dtype=int)
            mx_err = np.zeros((dx_short))
            msk = np.ones((dx_short),dtype=bool)


        for i in range (0, dx_short):
        
        # ---- Bresenham's algorithm (understandable variant)
        # http://www.cs.helsinki.fi/group/goa/mallinnus/lines/bresenh.html
            x += sx
            if 2* (D + dy) < dx:
                D += dy # y remains
            else:
                y += sy
                D += dy - dx

        #unfortunately, np allows for negative values...
        # when x or y are too large, the break is on try-except below


        # --------- coordinates not unpacked ! ------------

            mx_line[i, :] = [y, x] if not steep else [x, y]

            if interpolation:
    
                if D:
                    sD = -1 if D < 0 else 1
                    interp = y + sy *sD
                    
                    if steep:
                        mx_neighbours[i, :] = x, interp
                    else:
                        mx_neighbours[i, :] = interp, x
                    
                    mx_err [i]=D

                else:   msk[i]=False

        if interpolation:
        #give zero-error points themselves as neighbours
        # NB. this is not needed because error is 0; the result *= 0,
        # but it's more clear this way, and will eliminate the possibility of stepping out of matrix
            mx_neighbours[~msk, :]= mx_line[~msk, :]
    
        #error is actually D / dx !
            mx_err[msk]/= dx # zero values will give nans on division!
            return mx_line
        # also outpurs in original function: mx_neighbours, abs(mx_err)
    
        else: return mx_line


    def getDistance(self, x, y, pos, res): #get distance along LoS
        dX = abs(x-pos[0])*res
        dY = abs(y-pos[1])*res
        d = sqrt(dX**2+dY**2)

        return d
    
    # actual line of sight analysis
    def analysis(self, obsX, obsY, DEM, obsH, obsA, res, clear, curvature, refraction, rCoeff, RF, n, p, freq): #auxfile for debugging
        
        if RF:
            #convert frequency to wavelength for Fresnel calculation
            c = 0.299792458 #Gm/s sorry
            lam = c/freq
        
        #compute map of terrain elevation difference w.r.t. obs
        deltaH = DEM - (obsH * np.ones(DEM.shape))
        
        #empty array to write results to
        blank = np.zeros(DEM.shape)
        
        #draw targets for LoS
        frame = self.drawFrame(DEM.shape[0])

        for i in range(len(frame)):
            for j in range(len(frame[1])):
                #call precalculated los indices
                los = self.rasterised_line(obsX, obsY, frame[i][j][0], frame[i][j][1])
                
                #create/reset los specific values
                peak = [[0,0]]
                aMax = int(0)
                curv = 0
                precurv = 0
                
                if curvature:
                    #calculate earth curvature contribution:
                    R = 6371000 #m not best, in advanced version this could be function of location
                    #with refraction:
                    if refraction:
                        k = rCoeff
                        den = float((1-k))
                        R = R/den
            
                #now moving along los
                #purposely skipps first 3 pixels (i.e. assumes you position yourself at highest point in nearest surrounding)
                for k in range(3, len(los)):
                    #compute distance from obs to current step
                    d_i = self.getDistance(obsX, obsY, los[k], res)
                    d_i = float(d_i)
                    #store earth curvature contribution of previous step on los
                    precurv = curv
                    #compute earth curvature contribution for current step
                    if curvature:
                        curv = d_i**2/2./R
                    #get terrain elevation on current step
                    DEMi = DEM[los[k,0], los[k,1]]
                    #get height difference of obs to current step terrain elevation
                    dH = deltaH[los[k,0], los[k,1]]
                    #visual dH = vdH, i.e. with curvature term subtracted
                    #subtract optional earth curvature term from peaks
                    vdH = dH - curv
                    
                    #bool for protection against Nans (doesn't do anything rn I think, will be implemented later)
                    number = True
                    
                    #position is lower and there has not been a peak yet
                    if vdH <= peak[-1][0] and aMax == int(0):
                        
                        new = clear+curv
                        #auxfile.write('newA: ' + str(new) + '\n')
                    
                    
                    #position is lower and there has been a peak
                    elif vdH <= peak[-1][0] and aMax != int(0):
                        #compute new LoS based on distance increase
                        new = d_i*aMax + curv
                        #auxfile.write('newB: ' + str(new) + '\n')
                        
                        if RF and len(peak) > 1:
                            #check fresnel zones on critical spots on LoS (namely peaks), if peaks exist already
                            #avoid pair of zeros to start peak list
                            for l in range(len(peak)-1):
                                #go over list of peaks, from most recent to first
                                d_p = peak[-1-l][1]
                                h_p = peak[-1-l][0]
                                LoSHP = d_p*new/d_i - h_p    #height over peak on direct LoS
                                #actual Fresnel radius formula at critical spot
                                F1 = sqrt((n*lam*d_p*(d_i-d_p))/d_i)
                                F_eff = F1*sqrt(p)      #protected fraction area-->radius
                    
                                if F_eff >= LoSHP:
                                    #adjust new flight height if violation has been detected
                                    new = d_i * (F_eff+h_p)/d_p
                                else:
                                    pass
                        else:
                            pass
                
                
                    #position is higher
                    elif vdH/d_i > aMax:
                    
                        peak.append([vdH, d_i])
                        aMax = vdH/d_i
                        
                        new = d_i*aMax + curv
                        #auxfile.write('newC: ' + str(new) + '\n')
            
                    else:
                        #auxfile.write('d')
                        number = False
                        new = d_i*aMax + curv
                    
                        if RF and len(peak) > 1:
                            #check fresnel zones on critical spots on LoS (namely peaks), if peaks exist already
                            #avoid pair of zeros to start peak list
                            for l in range(len(peak)-1):
                                #go over list of peaks, from most recent to first
                                d_p = peak[-1-l][1]
                                h_p = peak[-1-l][0]
                                LoSHP = d_p*new/d_i - h_p    #height over peak on direct LoS
                                #actual Fresnel radius formula at critical spot
                                F1 = sqrt((n*lam*d_p*(d_i-d_p))/d_i)
                                F_eff = F1*sqrt(p)      #protected fraction area-->radius
                                    
                                if F_eff >= LoSHP:
                                    #adjust new flight height if violation has been detected
                                    new = d_i * (F_eff+h_p)/d_p
                                else:
                                    pass
                                        
                        #auxfile.write('newD: ' + str(new) + '\n')
                    
                    
                    #make sure there is a clearance
                    if new + obsH <= DEMi + clear and number:
                        new = DEMi + clear - obsH
                        #auxfile.write('values: ' + str(new+obsH) + ' < ' + str(DEMi + clear) + '\n' + 'newE: ' + str(new) + '\n')
                                
                    #make sure you don't go lower than before
                    if blank[los[k-1,0], los[k-1,1]] > new:
                        new = blank[los[k-1,0], los[k-1,1]]-precurv+curv
                        #auxfile.write('corr1: ' + str(new))
                    
                    ##and you don't overwrite higher value
                    ##this line had to be erased since it weirdly propagated high elevations to the side at angles(xy plane of close to 45 deg... at the cost of possibly overwriting a higher value for one LoS)
                    #if blank[los[k,0], los[k,1]] > new:
                        #new = blank[los[k,0], los[k,1]]
                        ##auxfile.write('corr2: ' + str(new))

                    blank[los[k,0], los[k,1]] = new

                        
        blank = blank + obsA*np.ones(blank.shape)
        return blank

        
    def writeOut(self, openRaster, filename, npOut, frameBand, loX, hiX, loY, hiY, noData = -99):
        # create the output image1 from frame of open raster DEM and output array
        driver = openRaster.GetDriver()
        rasterGT = openRaster.GetGeoTransform()
        rows = openRaster.RasterYSize
        cols = openRaster.RasterXSize
        
        #print driver
        outRaster = driver.Create(filename, cols, rows, 1, GDT_Int32)
        if outRaster is None:
            print 'Could not create Tif File'
            sys.exit(1)
        
        outBand = outRaster.GetRasterBand(1)
        outData = noData*np.ones(frameBand.shape) #for when cropping DEM is allowed
        outData[loX:hiX, loY:hiY] = npOut   #for when cropping DEM is allowed
        
        
        # write the data
        outBand.WriteArray(outData, 0, 0)
        
        # flush data to disk, set the NoData value and calculate stats
        outBand.FlushCache()
        outBand.SetNoDataValue(noData)
        
        # georeference the image and set the projection
        outRaster.SetGeoTransform(openRaster.GetGeoTransform())
        outRaster.SetProjection(openRaster.GetProjection())
        
        del outData

    
            
    def run(self):

        #gets the layers loaded in QGIS and adds it to the comboBox object from the plugin dialog
        iface = self.iface
        #clear combos        
        self.dlg.RasterBox.clear();self.dlg.VectorBox.clear()
        
        #take layers from canvas and sort
        myLayers = iface.mapCanvas().layers()
        
        for i in range(len(myLayers)):
            myLayer = myLayers[i]
        
            #add all to RasterBox
            self.dlg.RasterBox.addItem(myLayer.name(),myLayer.id())

            #add all to VectorBox
            self.dlg.VectorBox.addItem(myLayer.name(),myLayer.id())
        
            # TODO: Sort for valid inputs and sync ComboBox index to mapCanvas().layers()
            
        #-----------------------------------------------------------------------------#
    
        # show the dialog
        self.dlg.show()
        # Run the dialog event loop
        result = self.dlg.exec_()
        # See if OK was pressed
        if result:
            
            #for debugging and logging
            filename = self.dlg.OutputText.text()
            parts = filename.split(".")
            #auxfilename = parts[0]+"-aux.txt"
            #auxfile = open(auxfilename, 'w')
            
            #for logging metadata, mainly analysis parameters
            metafilename = parts[0] + "-0.txt"
            metafile = open(metafilename, 'w')
            
            #assign Raster Layer from dlg
            rasterLayerIndex = self.dlg.RasterBox.currentIndex()
            rasterLayer = myLayers[rasterLayerIndex]
            
            #-----------------------------------------------------------------------------#
            
            #open raster, get dims and band1
            gdal.AllRegister()
            provider = rasterLayer.dataProvider()
            openRaster = gdal.Open(str(provider.dataSourceUri()), gdal.GA_Update)
            
            if openRaster is None:
                print 'Could not open image file'
                sys.exit(1)
            
            rasterGT = openRaster.GetGeoTransform()
            rows = openRaster.RasterYSize
            cols = openRaster.RasterXSize
            band1 = openRaster.GetRasterBand(1)
            npBand1 = band1.ReadAsArray(0,0,cols,rows)
            
            #-----------------------------------------------------------------------------#
            
            #get Geo Coordinates of Observer
            geoX, geoY = self.getPointGeo(myLayers)
            #Convert from map to pixel coordinates
            indX, indY = self.geoToInd(geoX, geoY, rasterGT)
            
            #determine Observer Height
            obsZ = npBand1[indX, indY]
            obsA = int(self.dlg.HeightBox.text())
            obsH = obsZ + obsA
            
            #get min ground clearance:
            clear = int(self.dlg.ClearanceBox.text())
            
            loX, hiX, loY, hiY = self.cropInds(rasterGT, indX, indY)
            croppedBand = npBand1[loX:hiX, loY:hiY]
            newX = int(indX - loX)
            newY = int(indY - loY)
            
            #-----------------------------------------------------------------------------#
            
            #check if curvature is selected:
            curvature = self.dlg.CurvatureCheck.isChecked()
            #check if refraction is selected:
            refraction = self.dlg.RefractionCheck.isChecked()
            #set refraction coefficient
            if refraction:
                rCoeff = float(self.dlg.kBox.text())
            else:
                rCoeff = None
            
            #-----------------------------------------------------------------------------#
            
            #check if Fresnel is selected:
            RF = self.dlg.FresnelCheck.isChecked()
            #set Fresnel Details
            if RF:
                n = int(self.dlg.nBox.text())
                pp = int(self.dlg.pBox.text())
                p = pp/100.
                freq = float(self.dlg.fBox.text())
            else:
                n = None
                p = None
                pp = None
                freq = None
            
            #-----------------------------------------------------------------------------#
            
            #actual LoS Analysis:
            out = self.analysis(newX, newY, croppedBand, obsH, obsA, rasterGT[1], clear, curvature, refraction, rCoeff, RF, n, p, freq)
            
            #for debugging
            #np.savetxt(str(parts[0]) + '.csv', out, delimiter=",")
            
            #-----------------------------------------------------------------------------#
            
            #save LoS map
            filename1 = parts[0] + "-1.tif"
            self.writeOut(openRaster, filename1, out, npBand1, loX, hiX, loY, hiY)
            
            #-----------------------------------------------------------------------------#
            
            #make and save a 'obstruction' map
            obstructions =  croppedBand - obsH*np.ones(croppedBand.shape)
            obstructions[obstructions < 0] = -99
            filename2 = parts[0] + "-2.tif"
            self.writeOut(openRaster, filename2, obstructions, npBand1, loX, hiX, loY, hiY)
            
            #-----------------------------------------------------------------------------#
            
            

            metafile.write("input raster:  " + str(self.dlg.RasterBox.currentText()) + "\n" + "input vector:  " + str(self.dlg.VectorBox.currentText()) + "\n" + "observer height:  " + str(obsA) + " m" + "\n" + "radius of analysis:  " + str(self.dlg.RadiusBox.text()) + " km" + "\n" + "minimum clearance:  " + str(clear) + " m" + "\n" + "Earth Curvature:  " + str(curvature) + "\n" + "Refraction:  " + str(refraction) + ", k = " + str(rCoeff) + "\n" + "Fresnel Clearance: " + str(RF) + ", Zone " + str(n) + ", Clearance: " + str(pp) + "%" + ", Frequency: " + str(freq) + "GHz")
            
            #-----------------------------------------------------------------------------#
            '''
            #talk to user about clipping
            commands = ">>> import processing" + '\n' ">>>processing.runandload('gdalogr:cliprasterbymasklayer', '" + filename1 + "', subsitute path/to/shapefile.shp as string here, '0', False, False, False, 0, 0, 1, 1, 1, False, 0, False, '', None)"
            
            
            instructions = "If you want to clip the results with a shapefile and have them appear in the canvas temporarily, I suggest to adjust and run the commands in " + auxfilename +  " in the QGIS Python console."
            
            iface.messageBar().pushMessage("HINT!", instructions, level=QgsMessageBar.INFO, duration = 30)
            
            auxfile.write(commands)
            
            #-----------------------------------------------------------------------------#
            '''



