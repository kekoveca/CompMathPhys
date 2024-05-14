import numpy as np
import pandas as pd

# Вывод матрицы на экран
def print_arr( string, namevec, a ):
    if (type(a) == int) or (type(a) == float):
        print(a)
    else:
        print( string )
        for k in range(len(a)):   
            print("{}[{}] = {:8.4f}".format(namevec, k, a[k]))

# Проверка 3х-диаг. матрицы коэффициентов на корректность
def isCorrectArray(a):
    n = len(a)
    
    for row in range(0, n):
        if( len(a[row]) != n ):
            print('Не соответствует размерность')
            return False
        
    for row in range(1, n - 1):
        if(abs(a[row][row]) < abs(a[row][row - 1]) + abs(a[row][row + 1])):
            print('Не выполнены условия достаточности')
            return False
 
    if (abs(a[0][0]) < abs(a[0][1]))or(abs(a[n - 1][n - 1]) < abs(a[n - 1][n - 2])):
        print('Не выполнены условия достаточности')
        return False
        
    
    for row in range(0, len(a)):
        if( a[row][row] == 0 ):
            print('Нулевые элементы на главной диагонали')
            return False
    return True

# Процедура нахождения решения 3-х диагональной матрицы
def solution(a, b):
    if( not isCorrectArray(a) ):
        print('Ошибка в исходных данных')
        return -1 
 
    n = len(a)
    x = np.array([.0 for k in range(0, n)]) # обнуление вектора решений
    # print('Размерность матрицы: ',n,'x',n)
    
    # Прямой ход
    v = np.array([.0 for k in range(0, n)])
    u = np.array([.0 for k in range(0, n)])
    # для первой 0-й строки
    v[0] = a[0][1] / (-a[0][0]) 
    u[0] = ( - b[0]) / (-a[0][0]) 
    for i in range(1, n - 1): # заполняем за исключением 1-й и (n-1)-й строк матрицы
        v[i] = a[i][i+1] / ( -a[i][i] - a[i][i-1]*v[i-1] )
        u[i] = ( a[i][i-1]*u[i-1] - b[i] ) / ( -a[i][i] - a[i][i-1]*v[i-1] )
    # для последней (n-1)-й строки
    v[n-1] = 0
    u[n-1] = (a[n-1][n-2]*u[n-2] - b[n-1]) / (-a[n-1][n-1] - a[n-1][n-2]*v[n-2])
    
    # print_arr('Прогоночные коэффициенты v: ','v', v)
    # print_arr('Прогоночные коэффициенты u: ','u', u)
    
    # Обратный ход
    x[n-1] = u[n-1]
    for i in range(n-1, 0, -1):
        x[i-1] = v[i-1] * x[i] + u[i-1]
        
    return x

# параметры
eps = 10 ** -2

# аналитическое решение
def u_analytical(t, x, y):
    return ((1+x+y) ** 2) / ((13 - 12 * t))

# из аналитического решения
# начальное условие
def u_init(x, y):
    return u_analytical(0, x, y)

# граничное условие при x = 0
def u_bound_x_left(t, y):
    return u_analytical(t, 0, y)

# граничное условие при x = 1
def u_bound_x_right(t, y):
    return u_analytical(t, 1, y)

# граничное условие при y = 0
def u_bound_y_left(t, x):
    return u_analytical(t, x, 0)

# граничное условие при y = 1
def u_bound_y_right(t, x):
    return u_analytical(t, x, 1)

# прогонка №1
def a_l(u_prev, u, u_next, x_prev, x, x_next, tau, h_x):
    return -1 * (u_next + u) * tau / (2 * h_x ** 2)

def c_l(u_prev, u, u_next, x_prev, x, x_next, tau, h_x):
    return -1 * (u + u_prev) * tau / (2 * h_x ** 2)

# прогонка №2
def a_m(u_prev, u, u_next, x, tau, h_y):
    return -(u_next + u) * tau / (2 * h_y ** 2)

def c_m(u_prev, u, u_next, x, tau, h_y):
    return -(u + u_prev) * tau / (2 * h_y ** 2)


# сетка
# x часть
L_0 = 6
koef_L = 4
L = koef_L * (L_0 - 1) + 1
x_grid, h_x = np.linspace(0, 1, L, retstep=True)
x_grid_0 = np.linspace(0, 1, L_0)  # сетка для печати

# y часть
M_0 = 6
koef_M = 4
M = koef_M * (M_0 - 1) + 1
y_grid, h_y = np.linspace(0, 1, M, retstep=True)
y_grid_0 = np.linspace(0, 1, M_0)  # сетка для печати

# временнАя часть
T = 1  # окончание
N = 2000  # количество точек по времени
t_grid, tau = np.linspace(0, T, N, retstep=True)

print(f"h_x = {h_x}, h_y = {h_y}, tau = {tau}")

# сеточная функция
u = np.full((L, M), np.nan)
# задание начального условия
for l in range(L):
    for m in range(M):
        u[l, m] = u_init(x_grid[l], y_grid[m])
u_prev = []
u_k    = []
u_tmp  = []
# print(f"u = \n{u}")


# аналитическое решение
u_an    = np.full((L, M), np.nan)
delta_u = np.full((L, M), np.nan)
norma = 0
for l in range(L):
    for m in range(M):
        u_an[l, m]    = u_analytical(0, x_grid[l], y_grid[m])
        delta_u[l, m] = abs(u_an[l, m] - u[l, m])
        norma = max(norma, delta_u[l, m])

# print(f"u = \n{u}")
# print(f"u_an = \n{u_an}")
# print(f"delta_u = \n{delta_u}")
# print(f"norma = {norma}") ## раскомментить

# цикл по времени
counter = 1
for t in t_grid[1:]:
    # print(f"{counter}. t = {t}")## раскомментить
    counter += 1

    # переносим текущий слой в предыдущий, обнуляем текущий слой
    u_prev = u.copy()

    # цикл по k
    flag_k = True
    k = 0
    while flag_k:
        k += 1
        # print(f"    k = {k}")
        u_k   = u.copy()
        u     = np.full((L, M), np.nan)
        u_tmp = np.full((L, M), np.nan)
        # первая прогонка
        for m in range(1, M - 1, 1):
            # print(f"    m = {m}")
            A = np.full((L, L), .0)
            A[ 0,  0] = 1
            A[-1, -1] = 1
            d = np.full(L, .0)
            # print(f"d = \n{d}")
            d[ 0] = u_bound_x_left (t, y_grid[m])
            d[-1] = u_bound_x_right(t, y_grid[m])
            # print(f"d = \n{d}")
            for l in range(1, L - 1, 1):
                x_prev = (x_grid[l - 1] + x_grid[l]) / 2
                x_next = (x_grid[l + 1] + x_grid[l]) / 2
                A[l, l - 1] = c_l(u_k[l - 1, m], u_k[l, m], u_k[l + 1, m], x_prev, x_grid[l], x_next, tau, h_x)
                A[l, l + 1] = a_l(u_k[l - 1, m], u_k[l, m], u_k[l + 1, m], x_prev, x_grid[l], x_next, tau, h_x)
                A[l,     l] = 1 - A[l, l - 1] - A[l, l + 1]
                d[l]        = u_prev[l, m]
            sol = solution(A, d)
            # print(f"sol = {sol}")
            u_tmp[:, m] = sol
            # print(f"A = \n{A}")
            # print(f"d = \n{d}")
        # print(f"u_tmp = \n{u_tmp}")
        u_tmp[:,  0] = u_bound_y_left (t, x_grid)
        u_tmp[:, -1] = u_bound_y_right(t, x_grid)
        # print(f"u_tmp = \n{u_tmp}")

        # вторая прогонка
        for l in range(1, L - 1, 1):
            # print(f"    l = {l}")
            A = np.full((M, M), .0)
            A[ 0,  0] = 1
            A[-1, -1] = 1
            d = np.full(L, .0)
            d[ 0] = u_bound_y_left (t, x_grid[l])
            d[-1] = u_bound_y_right(t, x_grid[l])
            for m in range(1, M - 1, 1):
                A[m, m - 1] = c_m(u_k[l, m - 1], u_k[l, m], u_k[l, m + 1], x_grid[l], tau, h_y)
                A[m, m + 1] = a_m(u_k[l - 1, m], u_k[l, m], u_k[l + 1, m], x_grid[l], tau, h_y)
                A[m,     m] = 1 - A[m, m - 1] - A[m, m + 1]
                d[m]        = u_tmp[l, m]
            u[l, :] = solution(A, d)
        #     print(f"A = \n{A}")
        #     print(f"d = \n{d}")
        # print(f"u = \n{u}")
        u[ 0, :] = u_bound_x_left (t, y_grid)
        u[-1, :] = u_bound_x_right(t, y_grid)
        # print(f"u = \n{u}")

        # проверка условия перехода на следующий временной слой
        cond_max = 0
        for l in range(1, L - 1, 1):
            for m in range(1, M - 1, 1):
                cond_max = max(cond_max, abs((u[l, m] - u_k[l, m]) / u[l, m]))
                if (cond_max>1)and(abs((u[l, m] - u_k[l, m]) / u[l, m])>=cond_max) :
                     print(f"    cond_max = {cond_max}",(u[l, m] - u_k[l, m]),u[l, m],l,m)
        if cond_max < eps:
            flag_k = False

    # аналитическое решение
    u_an    = np.full((L, M), np.nan)
    delta_u = np.full((L, M), np.nan)
    norma = 0
    for l in range(L):
        for m in range(M):
            u_an[l, m]    = u_analytical(t, x_grid[l], y_grid[m])
            delta_u[l, m] = abs(u_an[l, m] - u[l, m])
            norma = max(norma, delta_u[l, m])

    # print(f"u = \n{u}")
    # print(f"u_an = \n{u_an}")
    # print(f"delta_u = \n{delta_u}")
    # print(f"norma = {norma}")

# аналитическое решение
u_an    = np.full((L, M), np.nan)
delta_u = np.full((L, M), np.nan)
norma = 0
for l in range(L):
    for m in range(M):
        u_an[l, m]    = u_analytical(T, x_grid[l], y_grid[m])
        delta_u[l, m] = abs(u_an[l, m] - u[l, m])
        norma = max(norma, delta_u[l, m])

print(f"u = \n{u}")
print(f"u_an = \n{u_an}")
print(f"delta_u = \n{delta_u}")
print(f"norma = {norma}")


# численное решение для печати
u_print = u[::koef_L, ::koef_M]

u_an_print   = np.full((L_0, M_0), np.nan)  # аналитическое решение
u_both_print = [0 for k in range(L_0)]      # оба решения для печати
for l in range(L_0):
    u_both_print[l] = [None for m in range(M_0)]
delta_abs = np.full((L_0, M_0), np.nan)  # абсолютная погрешность
# delta_rel = np.full((L_0, M_0), np.nan)  # относительная погрешность
norma_abs = 0  # норма абсолютной погрешности
# norma_rel = 0  # норма относительной погрешности
for l in range(L_0):
    for m in range(M_0):
        u_an_print[l, m]    = u_analytical(T, x_grid_0[l], y_grid_0[m])
        u_both_print[l][m]  = str(u_print[l, m]) + ' ' + str(u_an_print[l, m])
        delta_abs[l, m]     = abs(u_an_print[l, m] - u_print[l, m])
        norma_abs = max(norma_abs, delta_abs[l, m])
        # if u_an_print[l, m] != 0:
        #     delta_rel[l, m] = abs(delta_abs[l, m] / u_an_print[l, m])
        #     norma_rel = max(norma_rel, delta_rel[l, m])
        # else:
        #     delta_rel[l, m] = None
print(f"norma_abs = {norma_abs}")



print("Значения функции. Первое число - численное решение, второе - аналитическое")
# pd.options.display.float_format ='{:,.4f}'.format
# pd.set_option('display.float_format', '{:.2f}'.format)
pd.DataFrame({**{"y/x": y_grid_0}, **{x_grid_0[i]:u_both_print[i] for i in range(L_0)}})


print("Абсолютная погрешность")
print(f"norma_abs = {norma_abs}")
pd.DataFrame({**{"y/x": y_grid_0}, **{x_grid_0[i]:delta_abs[i] for i in range(L_0)}})