import numpy as np


# the base, per-bp costs, time (days)
primer = np.array([0, 0.6, 1.5])
ultramer = np.array([0, 0.8, 2.0])


a1 = np.arange(0, 40).reshape(-1, 1)
a2 = np.arange(40, 80).reshape(-1, 1)
a = np.concatenate([a1, a2])

b = np.ones(a.shape)
x = np.block([b, a, b * 20.0])

p1 = np.broadcast_to(primer, (a1.shape[0], x.shape[1]))
p2 = np.broadcast_to(ultramer, (a2.shape[0], x.shape[1]))
p = np.concatenate([p1, p2])

P = p.sum(axis=1)

_a, _b = np.meshgrid(P, P)
d = np.expand_dims(_a + _b, axis=2)

_a, _b = np.meshgrid(a, a)
ext = _a + _b
# use ext and junction array to compute efficiency


e = np.array([0, 0, 0, 0.1, 0.2, 0.4, 0.2, 0.8, 0.9, 1.0])
f = np.multiply(d, e)
print(f.shape)
g = f.min(axis=(0, 1))

import pylab as plt

plt.plot(g)
