# BioStatistics: HIV Dynamics Simulation

This project applies a biostatistical and mathematical modeling approach to study HIV infection dynamics. The simulation is based on a system of ordinary differential equations (ODEs) describing interactions among healthy CD4+ T-cells, infected CD4+ T-cells, free HIV virions, and immune response cells.

The project uses Python to numerically solve the HIV dynamic model and visualize how different socio-economic status assumptions may affect immune system response over time.

## Project Overview

HIV infection is a complex biological process involving viral replication, immune cell depletion, and immune system response. In this project, a mathematical model is used to represent the relationship between:

- Healthy CD4+ T-cells
- Infected CD4+ T-cells
- Free HIV virions
- CTLp cells
- CTLe cells

The model is implemented as a system of differential equations and solved using numerical integration. The goal is to observe how immune response variables change under different socio-economic status values.

## Model Description

The Python script simulates the following biological components:

| Variable | Meaning |
|---|---|
| `x` | Healthy CD4+ T-cells |
| `y` | Infected CD4+ T-cells |
| `v` | Free HIV virions |
| `w` | CTLp cells |
| `z` | CTLe cells |

The model includes parameters related to T-cell generation, natural death rate, infection rate, virion production, virion death, CTL proliferation, and CTL death.

The socio-economic status parameter `s` is introduced into the virion production term. The project compares three values:

```python
s_list = [0.2, 0.5, 0.8]
