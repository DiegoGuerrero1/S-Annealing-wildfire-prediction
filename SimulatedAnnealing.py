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
# import numba
from matplotlib import cm, animation
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



def simAnn(funVector, TempFunc, cInit, maxIter):
    velViento = np.random.randint(4, 16, size=(len(TempFunc)))
    # Definition of the simulated annealing method
    c = cInit  # Initial point
    iter = 0  # Iteration variable
    history = np.zeros(maxIter)  # Array for storing candidates
    alpha = 10                   # Amplification parameter for wind (How much does it affect)

    while iter < maxIter:
        for T in range(1,len(funVector)):
            Ec = funVector[c]  # Starting value
            n = nextPoint(0, len(funVector))  # New random position
            En = funVector[n]  # New random value
            deltaE = En - Ec  # Difference between the initial and the random value
            if deltaE < 0:  # If maxima: > 0 , If minima: < 0
                c = n  # The random point becomes our initial point
            elif probFun(deltaE, TempFunc[T]*alpha*(velViento[T]-velViento[T - 1])+1) > np.random.random():
                c = n  # Statistical consideration
        history[iter] = c
        iter += 1

    return c, history



def neighAnn(z, x, y, tempfunc, cinit, maxI, X, Y, Z, ):
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
        min, hist = simAnn(z, tempfunc, cinit, maxI)
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
    ax.scatter(coors_x, coors_y, coors_z, c='green', marker="X", linewidth=40)
    plt.show()


def plot3D(z, x, y, Z, X, Y, cinit):
    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.gist_stern
                    , linewidth=1, antialiased=True)
    ptoMin = [x[cinit], y[cinit], z[cinit]]
    print('punto minimo en :', ptoMin)
    ax.scatter3D(ptoMin[0], ptoMin[1], ptoMin[2] + 20, c='green', marker="X", linewidth=30)
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
    ax.scatter(coors_x[0], coors_y[0], coors_z[0], c='orange', marker="d", linewidth=60)

    # Plotting path taken
    ax.scatter(coors_x, coors_y, coors_z, c='red', marker="d", linewidth=5)
    # Plotting final point
    ax.scatter(coors_x[-1], coors_y[-1], coors_z[-1], c='green', marker="d", linewidth=60)
    print('Minimo python:', np.min(z), 'minimo recocido', coors_z[-1])
    plt.show()  # Para que la ventana de visualización permanezca estática


# Principal method. arguments: initial point and size of the arrays (10)
def heat_equation_vector(k, n):
    # Uses k and n to build matrix representing the initial point and its surroundings
    def create_area_pi():
        area_pto_ini = np.zeros(n * n)
        area_pto_ini[k] = 1

        if k > 0:
            area_pto_ini[k - 1] = 1
        if k < len(area_pto_ini):
            area_pto_ini[k + 1] = 1
        if k >= n:
            area_pto_ini[k - n] = 1
        if k >= n - 1:
            area_pto_ini[k - n - 1] = 1
        if k >= n + 1:
            area_pto_ini[k - n + 1] = 1
        if k < len(area_pto_ini) - n - 1:
            area_pto_ini[k + n - 1] = 1
        if k < len(area_pto_ini) - n + 1:
            area_pto_ini[k + n + 1] = 1
        if k < len(area_pto_ini) - n:
            area_pto_ini[k + n] = 1
        mat_area_pi = np.reshape(area_pto_ini, (n, n))
        return mat_area_pi

    print('Initial temperature')

    # Recibes area of initial point matrix and substitutes with fire and surroundings temperatures
    def sust_temps(area_pi, temp_fuego, temp_ext):
        area_bool = area_pi > 0.9
        temp_pt_ini = temp_fuego
        tmp_area = temp_ext
        T = area_pi + tmp_area
        T[area_bool] = temp_pt_ini
        return T

    # Turns a matrix into a vector. Necessary for method output
    def vectorize_matrix(A):
        AvSize = len(A) * len(A)
        Av = np.zeros(AvSize)
        for i in range(n):
            for j in range(n):
                Av[i * n + j] = A[i][j]
        return Av

    # For loops solving numerically for heat equation
    def heat_equation(T, maxIter):
        delta = 1

        for iteration in range(0, maxIter):
            for i in range(1, n - 1, delta):
                for j in range(1, n - 1, delta):
                    T[i, j] = 0.25 * (T[i + 1][j] + T[i - 1][j] + T[i][j + 1] + T[i][j - 1])
        return vectorize_matrix(T)

    fire_temperature = 1000
    surroundings_temperature = 30
    iterations = 15
    area_initial_point = create_area_pi()
    temperature_matrix = sust_temps(area_initial_point, fire_temperature, surroundings_temperature)
    heat_vector = heat_equation(temperature_matrix, iterations)
    print('Heat vector:', heat_vector)

    return heat_vector


def heatAnn(funVector, cInit, maxIter, lMat):           # Simulated annealing with temperatures defined by heat equation
    print('Please write the number of deltas of time')
    timeIterator = int(input())
    initPoint = cInit
    maximaPoints = []
    for i in range(timeIterator):
        newTemp = heat_equation_vector(initPoint, lMat)
        cNew, history = simAnn(funVector, newTemp, cInit, maxIter)
        initPoint = cNew

        maximaPoints.append(cNew)

    return initPoint, maximaPoints


def main():
    # Basic variables
    l = 100  # length
    x = np.random.randint(0, 100, l) # Array of x values
    y = np.random.randint(0, 100, l) # Array of z values
    z = np.random.randint(0, 100, l) # Array of y values

    # Constraint variables
    maxIterations = 1000
    print('Write the initial point:')
    cinit = int(input()) # Initial point

    ### Single simulated Annealing

    # Single evaluation
    cpoint, hist = simAnn(z, heat_equation_vector(cinit, l), cinit, maxIterations) # Returns the maximum point and the history of candidates

    ## Graphing
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
    ax.scatter3D(ptoMin[0], ptoMin[1], ptoMin[2] + 20, c='green', marker="d", linewidth=60)
    plt.show()

    print('pintos minimos en :', ptoMin)


    ### Additional Features

    ## Multiple search
    neighAnn(z, x, y, heat_equation_vector(cinit,l), cinit, maxIterations, X, Y, Z)

    ## Show History
    addHist(hist, X, Y, Z, x, y, z)

    ## Alongside heat equation ***** Experimental ********

    cmax, heatAhist = heatAnn(z, cinit, 1000, l)
    print('Maxima index using heat equation:', heatAhist)
    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.gist_stern
                    , linewidth=1, antialiased=True)
    ptoMin = [x[cmax], y[cmax], z[cmax]]
    for i in range(len(heatAhist)-1):
        lcolor = np.linspace(0,1,len(heatAhist))
        baseColor  = (lcolor[i], 0.2, 0.5)
        ax.scatter3D(x[heatAhist[i]], y[heatAhist[i]], z[heatAhist[i]] + 20, c=baseColor, marker="v", linewidth=60)
    ax.scatter3D(x[heatAhist[-1]], y[heatAhist[-1]], z[heatAhist[-1]] + 20, c='green', marker="*", linewidth=70)
    plt.show()
    # Done

if __name__ == "__main__":
    main()
