import os

import pydeck as pdk

token = 'sk.eyJ1IjoicnVuZHVzIiwiYSI6ImNsdzlyYmpwczA1NnQybXFrNHhiY3Iyd2EifQ.1hJLbkFo-aODO22MOiV-XA'

UK_ACCIDENTS_DATA = 'https://raw.githubusercontent.com/visgl/deck.gl-data/master/examples/3d-heatmap/heatmap-data.csv'


# rocket data
from ACESII_code.myImports import *
attitudeFolderPath = f'{ACES_data_folder}\\attitude'
inputFiles = [glob(f'{attitudeFolderPath}\\{fliers[0]}\\*.cdf*')[0],
                  glob(f'{attitudeFolderPath}\\{fliers[1]}\\*.cdf*')[0]]

data_dict_attitude_high = loadDictFromFile(inputFiles[0])
data_dict_attitude_low = loadDictFromFile(inputFiles[1])

os.environ['MAPBOX_API_KEY'] = token
# print(os.getenv('MAPBOX_API_KEY'))

layer = pdk.Layer(
    'HexagonLayer',  # `type` positional argument is here
    UK_ACCIDENTS_DATA,
    get_position=['lng', 'lat'],
    auto_highlight=True,
    elevation_scale=50,
    pickable=True,
    elevation_range=[0, 100],
    extruded=True,
    coverage=1)

# Set the viewport location
view_state = pdk.ViewState(
    longitude=-1.415,
    latitude=52.2323,
    zoom=6,
    min_zoom=5,
    max_zoom=500,
    pitch=90,
    bearing=-50)

# Combined all of it and render a viewport
r = pdk.Deck(layers=[layer], initial_view_state=view_state)
r.show()
r.to_html(r'C:\Users\cfelt\OneDrive\Desktop\hexagon-example.html')