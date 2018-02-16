import numpy as np

size = (10,)

a = np.random.random(size=size)
b = np.random.random(size=size)

print(a)
print(b)
print(sum(a*b))
print(sum(a*a) - 2*sum(a*b) + sum(b*b))
print(sum((a-b)**2))

c = np.random.random(size=(10,2))
print(c)
c = (c.T / np.sum(c, axis=1)).T
print(c)
print(np.sum(c,axis=1))
