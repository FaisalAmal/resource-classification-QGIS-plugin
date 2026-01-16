# -*- coding: utf-8 -*-

import os
from qgis.PyQt.QtWidgets import (QAction, QDialog, QVBoxLayout, QFormLayout, QLabel, QLineEdit, 
                                 QPushButton, QFileDialog, QDoubleSpinBox, QSpinBox, QGroupBox,
                                 QDialogButtonBox, QComboBox, QRadioButton, QButtonGroup, QCheckBox)
from qgis.PyQt.QtGui import QIcon
from qgis.core import QgsProject, QgsVectorLayer, Qgis, QgsApplication, QgsMapLayerProxyModel
from qgis.gui import QgsMapLayerComboBox
from . import core
import geopandas as gpd

class ResourceClassificationDialog(QDialog):
    def __init__(self, iface, parent=None):
        super(ResourceClassificationDialog, self).__init__(parent)
        self.iface = iface
        self.setWindowTitle("Resource Classification")

        layout = QVBoxLayout()
        form_layout = QFormLayout()

        # Input Selection Group
        self.input_group = QGroupBox("Input Selection")
        input_group_layout = QVBoxLayout()
        
        self.radio_active_layer = QRadioButton("Active Legend Layer")
        self.radio_file = QRadioButton("File System Path")
        self.radio_active_layer.setChecked(True)
        
        input_group_layout.addWidget(self.radio_active_layer)
        
        # Active Layer Combo
        self.layer_combo = QgsMapLayerComboBox()
        self.layer_combo.setFilters(QgsMapLayerProxyModel.PointLayer)
        input_group_layout.addWidget(self.layer_combo)
        
        input_group_layout.addWidget(self.radio_file)
        
        # File Input
        self.input_file_widget = QGroupBox()
        file_layout = QVBoxLayout()
        self.input_layer_edit = QLineEdit()
        browse_input_button = QPushButton("Browse")
        browse_input_button.clicked.connect(self.browse_input_layer)
        file_layout.addWidget(self.input_layer_edit)
        file_layout.addWidget(browse_input_button)
        self.input_file_widget.setLayout(file_layout)
        input_group_layout.addWidget(self.input_file_widget)
        
        self.input_group.setLayout(input_group_layout)
        layout.addWidget(self.input_group)
        
        # Connect radios to visibility toggles
        self.radio_active_layer.toggled.connect(self.toggle_input_mode)
        self.toggle_input_mode() # Initial state

        # Output Layer
        self.output_layer_edit = QLineEdit()
        browse_output_button = QPushButton("Browse")
        browse_output_button.clicked.connect(self.browse_output_layer)
        output_layout = QVBoxLayout()
        output_layout.addWidget(self.output_layer_edit)
        output_layout.addWidget(browse_output_button)
        form_layout.addRow("Output Polygon Layer:", output_layout)
        
        # Parameters
        self.measured_distance = QDoubleSpinBox()
        self.measured_distance.setRange(0, 10000)
        self.measured_distance.setValue(71)
        form_layout.addRow("Measured Distance (m):", self.measured_distance)

        self.indicated_distance = QDoubleSpinBox()
        self.indicated_distance.setRange(0, 10000)
        self.indicated_distance.setValue(121)
        form_layout.addRow("Indicated Distance (m):", self.indicated_distance)

        self.measured_buffer = QDoubleSpinBox()
        self.measured_buffer.setRange(0, 1000)
        self.measured_buffer.setValue(25)
        form_layout.addRow("Measured Buffer (m):", self.measured_buffer)

        self.indicated_buffer = QDoubleSpinBox()
        self.indicated_buffer.setRange(0, 1000)
        self.indicated_buffer.setValue(50)
        form_layout.addRow("Indicated Buffer (m):", self.indicated_buffer)
        
        # Alpha Measures
        # Measured Alpha
        self.measured_alpha = QDoubleSpinBox()
        self.measured_alpha.setRange(0, 10000)
        self.measured_alpha.setValue(0) # Default 0 means Auto
        self.measured_alpha_auto = QCheckBox("Auto")
        self.measured_alpha_auto.setChecked(True)
        self.measured_alpha.setEnabled(False)
        self.measured_alpha_auto.toggled.connect(lambda chk: self.measured_alpha.setEnabled(not chk))
        
        alpha_m_layout = QVBoxLayout()
        alpha_m_layout.addWidget(self.measured_alpha_auto)
        alpha_m_layout.addWidget(self.measured_alpha)
        form_layout.addRow("Measured Alpha Radius:", alpha_m_layout)

        # Indicated Alpha
        self.indicated_alpha = QDoubleSpinBox()
        self.indicated_alpha.setRange(0, 10000)
        self.indicated_alpha.setValue(0) # Default 0 means Auto
        self.indicated_alpha_auto = QCheckBox("Auto")
        self.indicated_alpha_auto.setChecked(True)
        self.indicated_alpha.setEnabled(False)
        self.indicated_alpha_auto.toggled.connect(lambda chk: self.indicated_alpha.setEnabled(not chk))
        
        alpha_i_layout = QVBoxLayout()
        alpha_i_layout.addWidget(self.indicated_alpha_auto)
        alpha_i_layout.addWidget(self.indicated_alpha)
        form_layout.addRow("Indicated Alpha Radius:", alpha_i_layout)

        self.min_drillholes = QSpinBox()
        self.min_drillholes.setRange(1, 100)
        self.min_drillholes.setValue(4)
        form_layout.addRow("Min Drillholes per Cluster:", self.min_drillholes)

        self.area_factor = QDoubleSpinBox()
        self.area_factor.setRange(0.1, 1.0)
        self.area_factor.setSingleStep(0.05)
        self.area_factor.setValue(0.6)
        form_layout.addRow("Area Factor:", self.area_factor)

        self.edge_tolerance = QDoubleSpinBox()
        self.edge_tolerance.setRange(1.0, 3.0)
        self.edge_tolerance.setSingleStep(0.1)
        self.edge_tolerance.setValue(1.5)
        form_layout.addRow("Edge Tolerance:", self.edge_tolerance)

        layout.addLayout(form_layout)

        # Buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.run_processing)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

        self.setLayout(layout)

    def toggle_input_mode(self):
        is_active = self.radio_active_layer.isChecked()
        self.layer_combo.setEnabled(is_active)
        self.input_file_widget.setEnabled(not is_active)

    def browse_input_layer(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select Drillhole Layer", "", "GeoPackage (*.gpkg);;All Files (*)")
        if path:
            self.input_layer_edit.setText(path)

    def browse_output_layer(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save Output Layer", "", "GeoPackage (*.gpkg)")
        if path:
            self.output_layer_edit.setText(path)

    def run_processing(self):
        output_path = self.output_layer_edit.text()

        if not output_path:
            self.iface.messageBar().pushMessage("Error", "Output path must be specified.", level=Qgis.Critical)
            return
            
        # Prepare Input Data
        input_data = None
        
        try:
            if self.radio_active_layer.isChecked():
                layer = self.layer_combo.currentLayer()
                if not layer:
                    self.iface.messageBar().pushMessage("Error", "No active layer selected.", level=Qgis.Critical)
                    return
                # Convert QGIS layer to GeoDataFrame
                # This could be slow for huge layers, but fine for typical drillholes
                import pandas as pd
                
                # Get all features
                feats = [f for f in layer.getFeatures()]
                if not feats:
                     self.iface.messageBar().pushMessage("Error", "Selected layer has no features.", level=Qgis.Critical)
                     return

                # Quick conversion (attributes + geometry)
                # Note: This is a simple conversion. For complex types, one might save to tmp file.
                data = []
                for f in feats:
                    attrs = f.attributes()
                    # Mapping fields is complex without knowing schema, but we mainly need geometry
                    # For simplicty, let's try reading from the source if it's a file
                    # If not possible (memory layer), we construct manually
                    row = dict(zip([f.name() for f in layer.fields()], attrs))
                    geom = f.geometry()
                    if geom:
                        wkt = geom.asWkt()
                        from shapely import wkt as shapely_wkt
                        row['geometry'] = shapely_wkt.loads(wkt)
                    data.append(row)
                
                gdf = gpd.GeoDataFrame(data)
                # Set CRS
                gdf.crs = layer.crs().toWkt()
                input_data = gdf
                
            else:
                input_path = self.input_layer_edit.text()
                if not input_path or not os.path.exists(input_path):
                     self.iface.messageBar().pushMessage("Error", "Valid input file path required.", level=Qgis.Critical)
                     return
                input_data = input_path

            # Prepare Parameters
            measured_alpha = None if self.measured_alpha_auto.isChecked() else self.measured_alpha.value()
            indicated_alpha = None if self.indicated_alpha_auto.isChecked() else self.indicated_alpha.value()

            self.iface.messageBar().pushMessage("Info", "Processing started...", level=Qgis.Info)
            
            result_gdf, _ = core.generate_resource_polygons(
                input_data=input_data,
                output_path=output_path,
                measured_distance=self.measured_distance.value(),
                indicated_distance=self.indicated_distance.value(),
                measured_buffer=self.measured_buffer.value(),
                indicated_buffer=self.indicated_buffer.value(),
                min_drillholes=self.min_drillholes.value(),
                area_factor=self.area_factor.value(),
                edge_tolerance=self.edge_tolerance.value(),
                measured_alpha=measured_alpha,
                indicated_alpha=indicated_alpha
            )

            if not result_gdf.empty:
                self.iface.messageBar().pushMessage("Success", "Processing complete.", level=Qgis.Success)
                # Load the result layer
                layer_name = os.path.splitext(os.path.basename(output_path))[0]
                result_layer = QgsVectorLayer(output_path, layer_name, "ogr")
                if result_layer.isValid():
                    QgsProject.instance().addMapLayer(result_layer)
                else:
                    self.iface.messageBar().pushMessage("Error", "Could not load result layer.", level=Qgis.Critical)
            else:
                self.iface.messageBar().pushMessage("Warning", "Processing resulted in no polygons.", level=Qgis.Warning)
            
            self.accept()

        except Exception as e:
            self.iface.messageBar().pushMessage("Error", f"An error occurred: {e}", level=Qgis.Critical)
            # self.reject() # Don't close on error so user can adjust

class ResourceClassificationPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.plugin_dir = os.path.dirname(__file__)
        self.actions = []
        self.menu = u'&Resource Classification'
        self.toolbar = self.iface.addToolBar(u'ResourceClassification')
        self.toolbar.setObjectName(u'ResourceClassification')

    def initGui(self):
        icon_path = os.path.join(self.plugin_dir, 'icon.png')
        action = QAction(QIcon(icon_path), u'Run Resource Classification', self.iface.mainWindow())
        action.triggered.connect(self.run)
        self.iface.addPluginToMenu(self.menu, action)
        self.toolbar.addAction(action)
        self.actions.append(action)

    def unload(self):
        for action in self.actions:
            self.iface.removePluginMenu(u'&Resource Classification', action)
        
        if self.toolbar:
            self.iface.mainWindow().removeToolBar(self.toolbar)
        del self.toolbar

    def run(self):
        dialog = ResourceClassificationDialog(self.iface)
        dialog.exec_()

def classFactory(iface):
    return ResourceClassificationPlugin(iface)
