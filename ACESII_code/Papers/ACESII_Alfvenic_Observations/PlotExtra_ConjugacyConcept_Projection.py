# --- PlotExtra.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: read in NOM_predict trajectories in order to roughly plot the ACESII rockets


import pygmt
from ACESII_code.myImports import *








# ----------------------
# Load the attitude Data
# ----------------------
attitudeFolderPath = f'{ACES_data_folder}\\attitude'
inputFilesTraj = [glob(f'{attitudeFolderPath}\\{fliers[0]}\\*.cdf*')[0],
                  glob(f'{attitudeFolderPath}\\{fliers[1]}\\*.cdf*')[0]]
data_dicts_attitude = [loadDictFromFile(inputFilesTraj[0]),loadDictFromFile(inputFilesTraj[1])]



fig = pygmt.Figure()
figure_height = 20
figure_width = 20
fig.set_figwidth(figure_width)
fig.set_figheight(figure_height)

# Format: 'GLon/Lat/width
projectType = 'G' # or 'g'
lon = 4
lat = 52
width = 12
azimuth = 30 # the "a" must be there
tilt = 45
vwidth= 60
vheight = 60
twist = 0
altitude = 250


# --- OLD CODE
fig = pygmt.Figure()
region_rocket = [8, 20, 65, 76, 0, 400] #
perspective = [-120, 35]  # azimuth, elevation (in deg)
resolution = "15s"
# projection = "G16.020833/69.294167/12c+a0+t45+v60/80+w0+z400"
projection = "M15c"
styleR1 = "3p,red,-"
styleR2 = "3p,blue,-"
zscale = 0.05
frame =['+t"Andoya, Norway"','xa2g','ya2g']
cmap = "geo"
registration = "gridline"
registration = None
# ----




projectionString = f"{projectType}{lon}/{lat}/{width}c+a{azimuth}+t{tilt}+v{vwidth}/{vheight}+w{twist}+z{altitude}"
fig.coast(
    projection=projectionString,
    region="g",
    frame=["x10g10", "y5g5"],
    land="gray",
)
fig.show()

#
# # Toggles
# BigAllSky_costLineSize = 2
# lonW = -10
# lonE = 40
# latS = 60
# latN = 80
# BigAllSky_GridSize = 1
# res = '50m'
#
# # ----------------
# # --- PLOTTING ---
# # ----------------
# # projProjection = ccrs.Orthographic(central_longitude=15, central_latitude=70)
# projProjection = ccrs.NearsidePerspective(central_longitude=0.0,
#                                           central_latitude=0.0,
#                                           satellite_height=35785831,
#                                           false_easting=10000000,
#                                           false_northing=0, globe=None)
# projTransform = ccrs.PlateCarree()
# fig, axBigAllSky = plt.subplots(1,subplot_kw=dict(projection=projProjection))
# figure_height = 20
# figure_width = 20
# fig.set_figwidth(figure_width)
# fig.set_figheight(figure_height)
#
#
#
# # axBigAllSky.set_extent([lonW, lonE, latS, latN])  # controls lat/long axes display
# gl = axBigAllSky.gridlines(draw_labels=True, linewidth=BigAllSky_GridSize,
#                                    alpha=0.4,
#                                    linestyle='--',
#                                    color='black')
#
#
# # coastlines
# axBigAllSky.coastlines(resolution=res, color='black',  alpha=1,linewidth=BigAllSky_costLineSize)  # adds coastlines with resolution
# plt.show()