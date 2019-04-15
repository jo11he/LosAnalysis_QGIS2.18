# -*- coding: utf-8 -*-
"""
/***************************************************************************
 LosAnalyzer
                                 A QGIS plugin
 Analyzes the Visibility of UAV from a given GCS
                             -------------------
        begin                : 2019-01-02
        copyright            : (C) 2019 by Jonas Hener - Avy BV
        email                : jonas@avy.eu
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load LosAnalyzer class from file LosAnalyzer.

    :param iface: A QGIS interface instance.
    :type iface: QgisInterface
    """
    #
    from .los_analyzer import LosAnalyzer
    return LosAnalyzer(iface)
