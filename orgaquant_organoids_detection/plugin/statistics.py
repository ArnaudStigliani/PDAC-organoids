import numpy as np

def area_ellipse(diameter1, diameter2):
    return np.around(1/4 * diameter1 * diameter2 * np.pi, decimals=1)
