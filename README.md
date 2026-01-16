# Resource Classification QGIS Plugin

This plugin creates resource category polygons based on drillhole spacing.

## Installation

1.  **Compile Resources:**
    *   Open a shell or command prompt that has access to your QGIS Python environment (or an environment where `PyQt5-tools` is installed).
    *   Navigate to this plugin directory (`qgis_plugin`).
    *   Run the following command to compile the resources file:
        ```bash
        pyrcc5 -o resources.py resources.qrc
        ```
    *   Alternatively, you can run the provided `pb_tool.py` script with the QGIS python interpreter.

2.  **Install Dependencies:**
    *   This plugin requires the `geopandas` and `scipy` libraries.
    *   You need to install them into the Python environment used by QGIS. You can do this from the QGIS Python console, or by using `pip` from an OSGeo4W shell or similar.
        ```
        pip install geopandas scipy
        ```

3.  **Install the Plugin in QGIS:**
    *   Copy the entire `qgis_plugin` directory to your QGIS plugins directory. You can find this directory in QGIS under `Settings > User Profiles > Open Active Profile Folder`, then navigate to `python/plugins`.
    *   Restart QGIS.

4.  **Enable the Plugin:**
    *   In QGIS, go to `Plugins > Manage and Install Plugins...`.
    *   Find "Resource Classification" in the list and check the box to enable it.
    *   You should see a new icon in the toolbar and a new menu under `Plugins > Resource Classification`.

## Usage

1.  Click the "Resource Classification" icon in the toolbar or select it from the plugins menu.
2.  In the dialog, select your input drillhole layer (GeoPackage).
3.  Specify the path for the output polygon layer (GeoPackage).
4.  Adjust the classification parameters as needed.
5.  Click "OK" to run the classification.
6.  The resulting polygon layer will be added to your QGIS project.
