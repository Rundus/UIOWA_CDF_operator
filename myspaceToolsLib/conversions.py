# --- conversions.py ---
# --- Author: C. Feltman ---
# DESCRIPTION: Place to store all the physics conversions I often use

# IMPORTS
from numpy import cos,radians,pi

# Variables
lat_to_meter = 111.319488  # 1 deg latitude to kilometers on Earth
Re = 6357 # radius of earth in kilometer

def long_to_meter(long, lat):
    return long*(lat_to_meter * cos(radians(lat)))

def meter_to_long(long_km, lat_km):
    latDeg = lat_km/lat_to_meter
    return (long_km/lat_to_meter) * 1/(cos(radians(latDeg)))

def calculateLong_to_meter(Lat): # determines the meters/long conversion for each input lattitude using earth's radius as a perfect sphere
    return (pi/180) * Re * cos(radians(Lat))