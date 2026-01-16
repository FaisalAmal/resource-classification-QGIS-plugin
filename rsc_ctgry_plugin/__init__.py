# -*- coding: utf-8 -*-

def classFactory(iface):
    """Load ResourceClassification class from file main.
    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    from .main import ResourceClassificationPlugin
    return ResourceClassificationPlugin(iface)
