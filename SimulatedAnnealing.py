"""
Simulated Annealing Algorithm for wildfire spread on 3D surfaces
Written by J. Diego Guerrero Morales
Contact: diegoguerrero@comunidad.unam.mx
Modified by Gabriel Peytral Borja
Inputs: A vectorized Matrix representation of a three-dimensional function
Output: Maxima of fire distribution
"""

import numpy as np
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def nextPoint(lowBound, uppBound):
    # Chooses a next point randomly
    nextP = np.random.randint(low=lowBound, high=uppBound)
    return nextP


def probFun(deltaE, T):
    # Probability function
    prob = np.exp(deltaE / T)
    return prob


def vectorize(matrix):
    # Vectorize a matrix
    n = len(matrix)
    matrixB = np.zeros_like(matrix)
    for i in range(n):
        for k in range(n):
            matrixB[n * i + k] = matrix[i][k]
    return matrixB


def repmat(item, f):
    # Create a matrix filled with an item
    matOut = np.zeros(f)
    for i in range(f):
        matOut[i] = item
    return matOut


def costFunction(wind, Tmin, Tmax, size):
    # Cost function definition, uses speed wind an assumes a parabolic growth of temperature
    costF = np.zeros(size)
    temp = costF
    for T in range(Tmin, Tmax):
        temp[T] = 2*T  # Parabolic temperature

        costF = temp * wind[T] / wind[0]+0.001  # Adding contribution of wind speed
    return costF


def simAnn(funVector, Tmin, Tmax, TempFunc, cInit, maxIter):
    # Definition of the simulated annealing method
    c = cInit  # Initial point
    iter = 0  # Iteration variable
    history = np.zeros(maxIter)  # Array for storing candidates
    while iter < maxIter:
        for T in range(Tmin, Tmax):
            Ec = funVector[c]  # Starting value
            n = nextPoint(0, Tmax)  # New random position
            En = funVector[n]  # New random value
            deltaE = En - Ec  # Difference between the initial and the random value
            if deltaE < 0:  # If maxima: > 0 , If minima: < 0
                c = n  # The random point becomes our initial point
            elif probFun(deltaE, TempFunc[T]) > np.random.random():
                c = n  # Statistical consideration
        history[iter] = c
        iter += 1

    return c, history


def simAnnParallel(funVector, Tmin, Tmax, cInit, maxIter):
    c = cInit
    iter = 0
    while (iter < maxIter):
        for k in range(len(c)):
            for T in range(Tmin, Tmax):
                Ec = funVector[c[k]]
                n = nextPoint(0, len(funVector))  # New random position
                En = funVector[n]
                deltaE = En - Ec  # Difference between the initial and the random position
                if deltaE < 0:  # If maxima: > 0 , If minima: < 0
                    c[k] = n  # The random point becomes our initial point
                elif probFun(deltaE, T) > np.random.random():
                    c[k] = n  # Statistical consideration.
            iter += 1
    return c, 10


def neighAnn(z, Tmax, Tmin, tempfunc, cinit, maxI, X, Y, Z, x, y):
    # Simulated Annealing with multiple initial points
    coors_x = []
    coors_y = []
    coors_z = []
    # Inicializo contador
    # count = 0
    print('Número de puntos: \n')
    npoint = int(input())

    for i in range(npoint):
        print("Da un punto inicial 0 a ", len(z))
        valor_k = int(input())
        min, hist = simAnn(z, Tmin, Tmax, tempfunc, cinit, maxI)
        print('Primer punto: z[', min, ']=', z[min])
        coors_x.append(x[int(min)])
        coors_y.append(y[int(min)])
        coors_z.append(z[int(min)])
        print("arreglos en interación", valor_k, "\n", coors_x, coors_y, coors_z)

    # Matplotlib visuals
    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.gist_stern
                    , linewidth=1, antialiased=True)
    ax.scatter(coors_x, coors_y, coors_z, c='orchid', marker="D", linewidth=40)
    plt.show()


def addHist(hist, X, Y, Z, x, y, z):
    # Show the record of candidates of a simulated annealing result
    coors_x = []
    coors_y = []
    coors_z = []
    for i in range(len(hist)):
        print('Primer punto: z[', hist[i], ']=', z[int(hist[i])])
        coors_x.append(x[int(hist[i])])
        coors_y.append(y[int(hist[i])])
        coors_z.append(z[int(hist[i])])
        print("arreglos en interación", hist[i], "\n", coors_x, coors_y, coors_z)

    # Matplotlib visuals
    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.gist_stern
                    , linewidth=1, antialiased=True)
    # Plotting first point
    ax.scatter(coors_x[0], coors_y[0], coors_z[0], c='cyan', marker="d", linewidth=60)

    # Plotting path taken
    ax.scatter(coors_x, coors_y, coors_z, c='red', marker="d", linewidth=5)
    # Plotting final point
    ax.scatter(coors_x[-1], coors_y[-1], coors_z[-1], c='green', marker="d", linewidth=60)
    print('Minimo python:', np.min(z), 'minimo recocido', coors_z[-1])
    plt.show()  # Para que la ventana de visualización permanezca estática


def main():
    # Basic variables
    l = 100  # length
    x = np.random.randint(0, 100, l)
    y = np.random.randint(0, 100, l)
    z = np.random.randint(0, 100, l)

    # Constraint variables
    Tmax = l
    Tmin = 0
    tol = 1000
    cinit = 3

    # Cost variables
    velViento = np.random.randint(4, 16, size=(Tmax))
    costfn = costFunction(velViento, Tmin, Tmax, l)

    ### Starts simulated Annealing

    # Single evaluation
    cpoint, hist = simAnn(z, Tmin, Tmax, costfn, cinit, tol)

    ## Gabo graph
    # construcción de la grilla 2D
    xi = np.linspace(np.amin(x), np.amax(x))
    print('Gridd :', np.amin(x), ',')
    yi = np.linspace(np.amin(y), np.amax(y))
    print('Gridd:', np.amin(y), ',')
    X, Y = np.meshgrid(xi, yi)

    # interpolación
    # Puntos, valores, xi, método
    Z = griddata((x, y,), z, (xi[None, :], yi[:, None]), method='cubic')

    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.gist_stern
                    , linewidth=1, antialiased=True)
    ptoMin = [x[cinit], y[cinit], z[cinit]]
    plt.show()

    print('pintos minimos en :', ptoMin)
    ax.scatter3D(ptoMin[0], ptoMin[1], ptoMin[2] + 20, c='red', marker="X", linewidth=30)

    ### Additional Features
    ## Multiple search
    neighAnn(z, Tmax, Tmin, costfn, cinit, tol, X, Y, Z, x, y)

    ## Show History
    addHist(hist, X, Y, Z, x, y, z)


if __name__ == "__main__":
    main()
