import numpy as np
import matplotlib as plt

#%% Plot the fields
X,Y = np.meshgrid( np.arange(-4,4,.2), np.arange(-4,4,.2) )
Ex = (X + 1)/((X+1)**2 + Y**2) - (X - 1)/((X-1)**2 + Y**2)
Ey = Y/((X+1)**2 + Y**2) - Y/((X-1)**2 + Y**2)

plt.figure()
plt.streamplot(X,Y,Ex,Ey)
plt.title('streamplot')

plt.figure()
plt.quiver(X,Y,Ex,Ey,scale=50)
plt.title('quiver')