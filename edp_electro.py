import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import *
from matplotlib import cm
import mpl_toolkits.mplot3d as axes3D
import seaborn as sns


# constante de discretización
h = 0.05

# se definen las constantes iniciales
epsilon_o = 8.85418782*(10**-12)    # C**2 / N m**2
epsilon = 10*epsilon_o
rho = 1.0 * 10**(-12)  # C /m**2
Vo = 4.7
cte1 = (rho/epsilon) * (h**2)
cte2 = 0
LEFT = RIGHT = UP = DOWN = 0 # bordes del diagrama


# se definen las dimensiones
[a, b, t, o, s] = [15, 20, 4, 5.5, 8]

# se escalan las dimensiones respecto al valor de la discretización
a = int(a/h)  # 150
b = int(b/h)  # 200
t = int(t/h)
o = int(o/h)
s = int(s/h)

# dimensiones de columnas y filas respectivamente
hi = b - 1
hj = a - 1

# matriz de resultados
mr = np.zeros((hj, hi))

# se desea resolver la ecuacion matricial Ax = b
N = hi*hj
A = np.zeros((N, N), dtype=int)
b = np.zeros((N,), dtype=int)


def getValue(i, j):
    return j*hi + i


def getPosition(k):
    return (k % hi, k//hi)


# se asocian los valores geometricamente al observar el diagrama
j1, j2, j3, j4 = t-1, o-1, o+t-1, a-t-1
i1, i2 = s-1, s+t-1


for i in range(0, hi):
    for j in range(0, hj):
        # numero asociado a un punto en la matriz
        k = getValue(i, j)

        # valor asociado en sus cuatro lados
        k_up = getValue(i, j+1)
        k_down = getValue(i, j-1)
        k_right = getValue(i+1, j)
        k_left = getValue(i-1, j)


# dielectrico:
        # puntos internos del dielectrico izquierda
        if 1 <= i <= i1-2 and j1+1 <= j <= j4-1:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1
        # puntos internos del dielectrico
        elif i2+2 <= i <= hi-2 and j1+1 <= j <= j4-1:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1
        # minisecciones dielectrico arriba
        elif j1+1 <= j <= j2-2 and i1-1 <= i <= i2+1:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1
        # miniseccion dielectrico DOWN
        elif j3+2 <= j <= j4-1 and i1-1 <= i <= i2+1:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1
        # dirichlet centro LEFT
        elif i == i1-1 and j2 <= j <= j3:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1 - Vo
        # dirichlet centro RIGHT
        elif i == i2+1 and j2 <= j <= j3:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k] = -4
            b[k] = cte1 - Vo
        # dirchlet centro UP
        elif j == j3+1 and i1 <= i <= i2:
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1 - Vo

        # dirchlet centro DOWN
        elif j == j2-1 and i1 <= i <= i2:
            A[k, k_up] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1 - Vo

        # condicion interior
        elif i1 <= i <= i2 and j2 <= j <= j3: 
            A[k, k] = 1
            b[k] = Vo

        # dirichlet borde izquierda
        elif i == 0 and j1 <= j <= j4:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1 - LEFT

        # dirichlet borde derecha
        elif i == hi-1 and j1 <= j <= j4:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k] = -4
            b[k] = cte1 - RIGHT

        # condicion UP de Neumann
        elif j == j4 and 1 <= i <= hi-2:
            A[k, k_down] = 2
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1 - 2*h*rho

        # condicion DOWN de Neumann
        elif j == j1 and 1 <= i <= hi-2:
            A[k, k_up] = 2
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1 + 2*h*rho
# caso vacio
        # parte interior arriba
        elif 1 <= i <= hi-2 and 1 <= j <= j1-1:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte2

        # parte interior abajo
        elif 1 <= i <= hi-2 and j4+1 <= j <= hj-2:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte2

        # dirchlet RIGHT arriba
        elif i == 0 and 1 <= j <= j1-1:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k] = -4
            b[k] = cte2 - RIGHT
        # dirchlet RIGHT abajo
        elif i == 0 and j3+2 <= j <= hj-2:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k] = -4
            b[k] = cte2 - RIGHT
        # dirchlet LEFT
        elif i == hi-1 and 1 <= j <= j1-1:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte2 - LEFT

        # dirchlet LEFT
        elif i == hi-1 and j3+2 <= j <= hj-2:
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte2 - LEFT
        # dirchlet UP
        elif j == hj-1 and 1 <= i <= hi-2:
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte2 - UP

        # dirchlet DOWN
        elif j == 0 and 1 <= i <= hi-2:
            A[k, k_up] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte2 - DOWN
# casos problematicos
        elif (j, i) == (j2-1, i1-1):
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1
        elif (j, i) == (j3+1, i1-1):
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1
        elif (j, i) == (j2-1, i2+1):
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1
        elif (j, i) == (j3+1, i2+1):
            A[k, k_up] = 1
            A[k, k_down] = 1
            A[k, k_left] = 1
            A[k, k_right] = 1
            A[k, k] = -4
            b[k] = cte1
            
# esquinas culias
        elif (j, i) == (hj-1, 0):
            A[k, k] = 1
            b[k] = -DOWN-LEFT+cte2
        elif (j, i) == (0, 0):
            A[k, k] = 1
            b[k] = -UP-LEFT+cte2
        elif (j, i) == (0, hi-1):
            A[k, k] = 1
            b[k] = -RIGHT-UP+cte2
        elif (j, i) == (hj-1, hi-1):
            A[k, k] = 1
            b[k] = -RIGHT-DOWN+cte2

        else:
            print("Point (" + str(j) + ", " + str(i) + ") missed!")
            raise Exception
# test.
print(A)


# se resuelve la ecuacion
x = np.linalg.solve(A, b)

for k in range(0, N):
    i, j = getPosition(k)
    mr[j, i] = x[k]

# se agregan los valores de las condiciones de borde
ub = np.zeros((hj+2, hi+2))
ub[1:hj+1, 1:hi+1] = mr[:]

# valor en el borde
# Dirichlet boundary condition
# right
ub[0: hj + 2, hi + 1] = RIGHT
#  left
ub[0:hj + 2, 0] = LEFT
# top
ub[0, 1:hi + 1] = LEFT
# buttom
ub[hj + 1, 1:hi + 1] = DOWN




plt.contour(ub, 100)


plt.figure(figsize=(9, 6))
sns.heatmap(ub, center=1, vmin=0,vmax = 4.0,xticklabels=20, yticklabels=20)


# Make data.
X = np.arange(0, ub.shape[0], 1, dtype=int)
Y = np.arange(0, ub.shape[1], 1, dtype=int)
X, Y = np.meshgrid(X, Y)
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_title('Solución de la Ecuación de Laplace')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('valor de \phi')

# Plot the surface.
surf = ax.plot_surface(X, Y, ub.T, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=6)

plt.show()
