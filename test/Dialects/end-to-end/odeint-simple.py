# RUN: python %s

# CHECK: Fail

# Extracted from: 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html

# Example Simple ODE

import numpy as np
import chemcomp
import chemcomp.lang as cl
from scipy.integrate import odeint
import matplotlib.pyplot as plt

@chemcomp.jit
def model(y, t, k):
    dydt = -k * y
    return dydt

# Initial condition
y0 = 5

# Time points
t = np.linspace(0, 20, 100)

# Solve the ODE for different values of k
k_values = [0.1, 0.5, 1.0]
for k in k_values:
    sol = odeint(model, y0, t, args=(k,))
    plt.plot(t, sol, label=f'k = {k}')

# Plot the results
plt.xlabel('Time')
plt.ylabel('y(t)')
plt.title('Exponential Decay')
plt.legend()
plt.show()