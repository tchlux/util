import numpy as np
from util.plot import Plot
from util.system import Timer

t = Timer()
sizes = []
times = []
for i in range(2,28):
    print("i: ",i)
    data = np.random.random(size=(2**i,))
    t.start()
    np.median(data)
    sizes.append(len(data))
    times.append( t.stop() )
    
p = Plot()
p.add("Runtime", sizes, times, mode="markers+lines")
p.show()
