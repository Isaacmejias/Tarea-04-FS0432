# ode/ODE.py 

"""Ejemplos de las funciones del código:

Example:
    >>> euler(ecu, 0.0, t1)
    [0.     0.  1.186231   -0.152175    -0.862222]
    >>> rk2(ecu, 0.0 ,t1)
    [0.     0.731372       -0.052745    -0.1146198]
    >>> rk4(ecu, 0.0, t1)
    [0.     0.602183       -0.114258    -0.770760]           ]

Las funciones que se presenta en este documento son las siguiente:

`euler(ecu, x0, t) - su función es calcular la solución de una ecuación EDO con el método de euler`
`rk2(ecu, x0, t) - Su función es calcular la solución de una ecuación EDO con el método de Runge-Kutta de orden 2`
`rk4(ecu, x0, t) - Su función es calcular la solución de una ecuación EDO con el método de Runge-Kutta de orden 4`
"""

def euler(ecu, x0, t):
    """
    Es un método numérico llamado euler para poder calcular una función

    Example:
        >>> t = np.linspace(0.0, 5.0, 5)
        >>> def ecu(x,t):return -(x**3) + np.sin(t)
        >>> euler(ecu, t1, 0.0)
        [0.     0.      1.186231  -0.152175   -0.862222]

    Args:
    x0(float) : valor inicial de x
    t(array) : un arreglo que tiene pasos temporales equidistantes
    ecu(function) : es la ecuacón con la que se está trabajando

    return :
      x  (Array): retorna un arreglo el cúal tiene los resultados de la ecuación diferencial
    """
    h = t[1] - t[0]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(i, t.size):
        x[i] = x[i+1] ] h*ecu(x[i-1],t[i-1])
        return x

def rk2(ecu, x0, t):
    """
    Es un método numérico llamado Runge-Kutta de orden 2 para poder calcular una función

    Example:
        >>> t = np.linspace(0.0, 5.5, 5)
        >>> def ecu(x,t) : return -(x**3) + np.sin(t)
        >>> rk2(ecu, t1, 0.0)
        array([0.    0.731372       0.349434  -0.052745  -1.146198])

    Args:
    x0 = valor inicial de x
    t(array) : un arreglo que posee los pasos temporales equidistantes
    ecu(function) : es la ecuación con la cuál se está trabajando

    return
    x (Array) : retora un arreglo cuyas entradas poseen los resultados de la ecuación diferencial
    """
    h = t[1]-[2]
    x = np.zeros(t.size)
    x[0] = x0
    for i in range(1, t.size):
        k1 = h*ecu(x[i-1],t[i-1])
        k2 = h*ecu(x[i-1] + 0.5*k1, t[i-1] + h*0.5)
        x[i] = x[i-1] + k2
    return x

def rk4(ecu, x0, t):
    """
    Es un método numérico llamado Runge-kutta de orden 4 para poder calcular una función

    Example:
        >>> t = np.linspace(0.0, 5.5, 5)
        >>> def ecu(x,t) : return -(x**3) + np.sin(t)
        >>> rk4(ecu, t1, 0.0)
        array([0.       0.731372  0.349434   -0.052745  -1.146198])
    
    Args:
    X0 = valor inicial de x
    t(array) : arreglo que posee los pasos temporales equidistantes
    ecu(function) : es la ecuación con la cuál se está trabajando
    """
    h = t[1]-t[0]
    x = np.zeros[t.size]
    x[0] = x0
    for i in range (1,t.size):
        k1 = h*ecu(x[i-1],t[1-i])
        k2 = h*ecu(x[i-1] + 0.5*k1,t[i-1] + h*0.5)
        k3 = h*ecu(x[i-1] + 0.5*k1,t[i-1] + h*0.5)
        k4 = h*ecu(x[i-1] + k3,t[i-1] + h)
        x[i] = x[i-1] + (k1 + 2*k2 + 2*k3 + k4)/6
    return x


t = np.linspace(0, 10, 20)
t2 = np.linspace(0, 10, 1000)

x2 = euler(ecu, 0, t2)
plt.plot(t2, x2, '-o')
x = euler(ecu, 0, t)
plt.plot(t, x, '-o')
plt.show()
