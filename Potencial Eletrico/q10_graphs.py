import numpy as np
import matplotlib.pyplot as plt

def sph2cart(az, el, r):
    rcos_theta = r * np.cos(el)
    x = rcos_theta * np.cos(az)
    y = rcos_theta * np.sin(az)
    z = r * np.sin(el)
    return x, y, z

num = 50
xv = np.linspace(-20,20,num)
yv = np.linspace(-20,20,num)
zv = np.linspace(-20,20,num)

raio = 8.5

rrl= np.linspace(0,raio,num)
thetal= np.linspace(0,2*np.pi,num)
phil= np.linspace(0,np.pi,num)

X,Y = np.meshgrid(xv,yv)


# frist X,Y

pv = 10e-12
a = 8.5
intervalo  = 10
valor = 100
ke = 1/(4*np.pi*8.85418e-12)

V = np.zeros((num,num, num))
print(V)

# x -> i
# y -> j
# x = a + intervalo/2 
i = j = 0
for xi in xv:
    print(xi)
    for yj in yv:
        for zk in zv:

            for r in rrl:
                for t in thetal:
                    for p in phil:
                        x = a + intervalo/2 
                        # calcula o valor da carga do intevalor de linha
                        r[x,y,z] = [xi,yj,zk] 
                        rl[xl, yl, zl] = sph2cart(t, p, r)

                        Q = ( (kl*x) / ((x**2) + (a**2)) ) * intervalo  
                        d = np.sqrt(r - rl)

                        dV = (r^2)*sin(t)*(np.pi/num)*(pi/180)*(raio/num)*(2*np.pi/num)*(np.pi/180);
                        if d<0.01:
                            d == 0.01
                        V[j][i][k] += ke*pv*(Q/d)
                          

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_wireframe(X, Y, V, color='black')  
plt.show()

print(V[0][0][0])
print(V[0][1][0])
