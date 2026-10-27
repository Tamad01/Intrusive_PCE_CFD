import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return 2

def solex(x):
    return (x*(2-x))/k0

a = 0
b = 2
ne = 20
ns = ne+1
h = (b-a)/ne
Nodes = np.linspace(a,b,ns)

ua = 0
ub = 0

F_base = np.zeros(ns)
for i in range(ne):
    Fki = (h/2)*np.array([f(Nodes[i]), f(Nodes[i+1])])
    F_base[i] += Fki[0]
    F_base[i+1] += Fki[1]

A_base = np.zeros((ns,ns))
for i in range(ne):
    Aki = (1/h)*np.array([[1,-1],[-1,1]])
    A_base[i, i] += Aki[0,0]
    A_base[i,i+1] += Aki[0,1]
    A_base[i+1, i] += Aki[1,0]
    A_base[i+1, i+1] += Aki[1,1]

k0 = 1.0  # Moyenne du matériau
k1 = 0.2  # Fluctuation 20%

#SUPER-MATRICE
A_super = np.zeros((2*ns, 2*ns))
F_super = np.zeros(2*ns)

# Remplissage par blocs
#d00 = k0, d01 = k1, d10 = k1, d11 = k0
A_super[0:ns, 0:ns]       = k0 * A_base  # Bloc d00
A_super[0:ns, ns:2*ns]    = k1 * A_base  # Bloc d01
A_super[ns:2*ns, 0:ns]    = k1 * A_base  # Bloc d10
A_super[ns:2*ns, ns:2*ns] = k0 * A_base  # Bloc d11

F_super[0:ns] = F_base

#CONDITIONS AUX LIMITES

A_super[0, :] = 0
A_super[0, 0] = 1
F_super[0] = ua

A_super[ns-1, :] = 0
A_super[ns-1, ns-1] = 1
F_super[ns-1] = ub

A_super[ns, :] = 0
A_super[ns, ns] = 1
F_super[ns] = 0

A_super[2*ns-1, :] = 0
A_super[2*ns-1, 2*ns-1] = 1
F_super[2*ns-1] = 0

#RÉSOLUTION
U_super = np.linalg.solve(A_super, F_super)

U0 = U_super[0:ns]
U1 = U_super[ns:2*ns]

#Calcul de l'écart-type (pour P=1 avec Hermite, Variance = U1^2*1)
std_dev = np.abs(U1)

plt.figure(figsize=(10, 5))
plt.plot(Nodes, U0, 'b-', label='Moyenne (U0)', linewidth=2)
plt.fill_between(Nodes, U0 - 3*std_dev, U0 + 3*std_dev, color='blue', alpha=0.2, label=r'Intervalle de confiance ($\pm 3\sigma$)')
plt.plot(Nodes, solex(Nodes), 'k--', label='Solution Exacte (k=1)')
plt.title("Solution avec Polynôme de Chaos Intrusif (P=1)")
plt.xlabel("x")
plt.ylabel("u(x)")
plt.legend()
plt.grid(True)
plt.show()