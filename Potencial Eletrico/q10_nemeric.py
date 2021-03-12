# Calcule numericamente o potencial produzido por uma linha de cargas 
# ρL = k.x/(x² +a²) que se estende ao longo do eixo x, 
# de x = a até +∞, com a = 8,2 m e k = 1,0 pC, em todo o espaço. 
# Considere o zero de referência no infinito. 
# Trace o gráfico das linhas equipotenciais nos planos yz e xz. 
# Para validar sua resposta, calcule V na origem.


# importante

# ρL = k.x/(x² +a²)
# a = 8,2
# k = 1,0 pC
# v = 0,86e-3
import numpy as np

a = 8.2
k = 1e-12
ke = 1/(4*np.pi*8.85418e-12)
intervalo = 0.001 # 1mm



valor = 100
x = a + intervalo/2 
V = 0
while valor > 0.000000000001:
    
    # calcula o valor da carga do intevalor de linha
    Q = ( (k*x) / ((x**2) + (a**2)) ) * intervalo
    # valor do potencial dado por aquela carga
    valor = ke*(Q/x) #k * Q * dx 
    # somoo valor 
    V += valor
    # avanco X
    x = x + intervalo
print(V)

