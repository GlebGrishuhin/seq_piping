from matplotlib import pyplot as plt
from matplotlib import animation


def animate(n):
    plt.cla()
    plt.autoscale(enable=False)
    plt.grid(visible=True)
    ax.set_xlim(limX)
    ax.set_ylim(limY)    
    line = plt.plot(x[n * N: (n + 1) * N], y[n * N: (n + 1) * N], color='b')
    #print('x=', *x[n * N: (n + 1) * N])
    #print('y=', *y[n * N: (n + 1) * N])
    return line


fig = plt.figure()
ax = plt.gca()

with open('data.txt') as fobj:
    t, x, y = zip(*([float(i) for i in line.split(',')] for line in fobj))
    
limX = [min(x) - 0.1 * (max(x) - min(x)), max(x) + 0.1 * (max(x) - min(x))]
limY = [min(y) - 0.1 * (max(y) - min(y)), max(y) + 0.1 * (max(y) - min(y))]

N = len(set(x))

anim = animation.FuncAnimation(fig, animate, frames=len(set(t)), interval=100, repeat=False)
plt.show()