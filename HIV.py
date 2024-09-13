import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#use the same initial value in paper ”MODELLING AND CONTROL OF HIV DYNAMICS”
lambda_ = 1 #rate of generation of new TCD4+ T-cells
d = 0.1 #natural death rate of healthy CD4+ T-cells
beta = 0.5 #infection rate of healthy CD4+ T-cells by free virions
a = 0.2 #death rate of infected CD4+ T-cells
p = 0.1 #rate at which infected CD4+ T-cells are killed by CTLe
k = 0.8 #production rate of free virions from infected CD4+ T-cells
u = 0.3 #death rate of free virions
c = 0.1 #proliferation rate of CTLp due to interactions with infected cells
q = 0.2 #conversion rate of CTLp to CTLe
b = 0.05 #death rate of CTLp
h = 0.05 #death rate of CTLe
x0 = 1000 #initial healthy CD4+ T-cell number
y0 = 100 #initial infected CD4+ T-cell number
v0 = 100 #initial free virion number
w0 = 50 #initial CTLp number
z0 = 20 #initial CTLe number
ic = [x0, y0, v0, w0, z0]
t = np.linspace(0, 200, 1000)

def hiv(y, t, lambda_, d, beta, a, p, k, u, c, q, b, h, s):
    x, y, v, w, z = y
    dx = lambda_ - d * x - beta * x * v
    dy = beta * x * v - a * y - p * y * z
    dv = k * (1 - s) * y - u * v
    dw = c * x * y - q * y * w - b * w
    dz = q * y * w - h * z
    return [dx, dy, dv, dw, dz]

s_list = [0.2, 0.5, 0.8] #simulate with 
results = {}

for s in s_list:
    solution = odeint(hiv, ic, t, args=(lambda_, d, beta, a, p, k, u, c, q, b, h, s))
    results[s] = solution.T

fig, axs = plt.subplots(2, 3, figsize=(15, 10))
for i, s in enumerate(s_list):
    x, y, v, w, z = results[s]
    axs[0, i].plot(t, w, '-m', label='CTLp cells', linewidth=2)
    axs[0, i].set_title(f'Socio-economic status: {s}')
    axs[0, i].set_xlabel('Time')
    axs[0, i].set_ylabel('CTLp cells')
    axs[1, i].plot(t, z, '-k', label='CTLe cells', linewidth=2)
    axs[1, i].set_xlabel('Time')
    axs[1, i].set_ylabel('CTLe cells')
fig.tight_layout()
plt.show()
