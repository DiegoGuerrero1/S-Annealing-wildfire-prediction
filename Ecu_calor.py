import numpy as np


def ecu_calor(areaPI, coefDT, tempiOBJ, tmpEXT, n):
    # Indexación booleana de matriz de área asociada a punto inicial
    area_bool = areaPI < 0.9

    # Coeficiente de difusividad térmica para madera
    dt_madera = coefDT

    # Temperatura inicial del ambiente y el punto inicial (madera quemándose)
    temp_pto_ini = 273.15 + tempiOBJ
    tmp_area = 273.15 + tmpEXT

    # Matriz de ceros en la que se sustituye con la temperatura inicial del área cerca del PI
    temp_ini = np.zeros([n, n]) + tmp_area

    # Usamos la indexación booleana para colocar el punto de inicio de incendio dentro de su área
    temp_ini[area_bool] = temp_pto_ini

    # Tomamos frames de como se propaga la temperatura
    times = 36000  # Número de iteraciones
    time_sp = 3600  # Número de frame
    f = int(times / time_sp)

    # Creamos arreglo donde se va a almacenar la evolución de la temperatura
    heat_frames = np.zeros([time_sp, n, n])
    heat_frames[0] = temp_ini

    # Verificamos parámetros del método numérico
    x = 0.5
    dx = 0.5 / 100
    dt = 1

    # La siguiente relación tiene que estar por debajo del 0.25 para que funcione el método
    dt_madera * dt / dx ** 2

    # Función para el método numérico


    # LLamamos la función que aprovecha indexación booleana para distinguir
    # sobre que partes del arreglo se tiene que iterar
    heat_frames = solve_heat(heat_frames, area_bool)
    # lo pasamos a grados celsius
    heat_frames -= 273.15

    return heat_frames[0], heat_frames[3599]  # Ajustar


def create_area_pi(k, n):  # Recibe el punto inicial k y el tamaño de los arreglos
    area_pto_ini = np.ones(n*n)
    area_pto_ini[k] = 0  # Asigno un 1 en punto inicial para que en ecu_calor ese 1 se convierta en 300

    if k > 0:
        area_pto_ini[k - 1] = 0
    if k < len(area_pto_ini):
        area_pto_ini[k + 1] = 0
    if k >= n:
        area_pto_ini[k - n] = 0
    if k >= n - 1:
        area_pto_ini[k - n - 1] = 0
    if k >= n + 1:
        area_pto_ini[k - n + 1] = 0
    if k < len(area_pto_ini) - n - 1:
        area_pto_ini[k + n - 1] = 0
    if k < len(area_pto_ini) - n + 1:
        area_pto_ini[k + n + 1] = 0
    if k < len(area_pto_ini) - n:
        area_pto_ini[k + n] = 0

    mat_area_pi = np.reshape(area_pto_ini, (n, n))
    return mat_area_pi


def main():
    n = 10  # Largo de la matriz radio del punto inicial
    k = 55  # Punto inicial
    matriz_pi = create_area_pi(k, n)
    coef = 8.2e-8
    temp_fuego = 300
    temp_amb = 25
    mat_temps, mat_temps2 = ecu_calor(matriz_pi, coef, temp_fuego, temp_amb, n)
    print('\n', mat_temps, '\n', mat_temps2)
    # Vectorizar matriz de temperaturas

if __name__ == "__main__":
    main()
