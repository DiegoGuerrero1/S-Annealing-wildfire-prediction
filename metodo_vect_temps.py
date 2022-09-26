import numpy as np


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

    return heat_vector
