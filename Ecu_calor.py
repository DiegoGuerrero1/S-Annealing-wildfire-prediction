import numpy as np


def ecu_calor(areaPI, coefDT, tempiOBJ, tmpEXT, n):
    # Indexación booleana de matriz de área asociada a punto inicial
    # Tiene que ser de n*n y entrar en blanco y negro
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
    def solve_heat(heatmap, hay_esf):
        cs = heatmap[0].copy()  # Estado actual (current state)
        length = len(cs[0])
        cf = 0  # frame actual
        for t in range(1, times):
            ns = cs.copy()  # estado nuevo
            for i in range(1, length - 1):
                for j in range(1, length - 1):
                    if hay_esf[j][i]:
                        a = dt_madera
                        ns[j][i] = cs[j][i] + a * dt / dx ** 2 * (cs[j + 1][i] + cs[j - 1][i] + \
                                                                  cs[j][i + 1] + cs[j][i - 1] - \
                                                                  4 * cs[j][i])
            cs = ns.copy()
            if t % f == 0:
                cf = cf + 1
                heatmap[cf] = cs

        return heatmap

    # LLamamos la función que aprovecha indexación booleana para distinguir
    # sobre que partes del arreglo se tiene que iterar
    heat_frames = solve_heat(heat_frames, area_bool)
    # lo pasamos a grados celsius
    heat_frames -= 273.15
    
    return heat_frames[1]  # Ajustar


def main():
    n = 10  # Largo de la matriz radio del punto inicial
    radio_pi = np.random.randint(0, 2, size=(n, n))
    coef = 8.2e-8
    temp_fuego = 300
    temp_amb = 25
    mat_temps = ecu_calor(radio_pi, coef, temp_fuego, temp_amb, n)
    print(mat_temps)

if __name__ == "__main__":
    main()
