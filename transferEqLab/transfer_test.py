import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

a =  np.linspace(0, 1, 11)
b =  np.linspace(0, 1, 41)
i, = np.where(np.isclose(b, 0.1))
print(i)
