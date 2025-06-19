import numpy as np
import matplotlib.pyplot as plt

u = np.zeros((10,10))
dx = 0.1
dy = 0.1
for i in range(10):
    for j in range(10):
        u[i,j] = np.exp(-np.power((dx*i * dy*j) - 5.0, 2) / 0.5)

plt.imshow(u.T,origin='lower')
plt.show()        
