# RUN: python %s

# CHECK: Fail!

# Extracted from: https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

import numpy as np
import chemcomp
import chemcomp.lang as cl
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Example: Pendulum
@chemcomp.jit
def pend(y, t, b, c):
    # Let y be the vector [theta, omega]. We implement this system in Python as:
    # theta'(t) = omega(t)
    # omega'(t) = -b*omega(t) - c*sin(theta(t))
    theta, omega = y
    dydt = [omega, -b*omega - c*np.sin(theta)]
    return dydt

b = 0.25
c = 5.0
y0 = [np.pi - 0.1, 0.0]
t = np.linspace(0, 10, 101)

# Solve
sol = odeint(pend, y0, t, args=(b, c))

# Plot the results
plt.plot(t, sol[:, 0], 'b', label='theta(t)')
plt.plot(t, sol[:, 1], 'g', label='omega(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
