import numpy as np

def period2freq(period):
    return 1.0 / period * 2.0 * np.pi

def freq2period(freq):
    return 1.0 / freq * 2.0 * np.pi

