import numpy as np

N = 1000
D = 100
x = np.random.random(size=(N,D))
z = np.random.random(size=(D,))


from util.system import Timer
t = Timer()

t.start()
np.linalg.norm(x - z, axis=1)
t.stop()
print("Norm time:", t.total)

t.start()
np.sum((x - z)**2, axis=1)
t.stop()
print("Sq time:", t.total)

t.start()
np.sum(x**2, axis=1) + np.sum(z**2) - 2 * np.dot(x, z)
t.stop()
print("Dot time:", t.total)
