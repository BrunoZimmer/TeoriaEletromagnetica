# Dado o campo de potencial, no espaço livre, V = 100xz/(x2 + 4)
# V, considerando que a superfície z = 0 seja a superfície de um
# condutor, calcule numericamente a distribuição superficial de
# cargas para qualquer lugar no espaço, e determine a carga
# total na porção do condutor definida por 0 < x < 8,8 m ; -6,9 m < y < 0.

import numpy as np


# apos deducoes matematicas e fisicas ps = Ds

# Ds = -e * d(v)dz = ps

# ps = -e100x/(x^2 + 4) 

# Q = integral(ps) in x and y
# q = 9.2 e-9

e = 8.85418e-12
num = 1000

xi = 0
xf = 6.9
dx = abs(xi-xf)/num

yi = -1.2
yf = 0
dy = abs(yi-yf)/num

dA = dx*dy

x = np.linspace(0,8.8,num)
y = np.linspace(-6.9,0,num)

Q = 0 

for X in x:
    for Y in y:
        xv = X + dx/2 
        dq =  (-e*100 * xv) / ( xv**2 + 4)
        dq = dq * dA 
        Q += dq
print(Q)




