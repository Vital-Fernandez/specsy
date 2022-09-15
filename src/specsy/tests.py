import matplotlib.pyplot as plt
import numpy as np

from lmfit import Minimizer, Parameters


def func(pars, x, data=None):
    a, b, c = pars['a'], pars['b'], pars['c']
    model = a * np.exp(-b*x) + c
    if data is None:
        return model
    return model - data


def dfunc(pars, x, data=None):
    a, b = pars['a'], pars['b']
    v = np.exp(-b*x)
    return np.array([v, -a*x*v, np.ones(len(x))])


def f(var, x):
    return var[0] * np.exp(-var[1]*x) + var[2]


params = Parameters()
params.add('[a_01', value=10)
params.add('b', value=10)
params.add('c', value=10)

a, b, c = 2.5, 1.3, 0.8
x = np.linspace(0, 4, 50)
y = f([a, b, c], x)
np.random.seed(2021)
data = y + 0.15*np.random.normal(size=x.size)