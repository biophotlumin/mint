import matplotlib.animation as animation
import numpy as np
from pylab import *
import imageio
import pandas as pd

dpi = 200

def ani_frame():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

    im = ax.imshow(rand(300,300),cmap='gray',interpolation='nearest')
    im.set_clim([0,1])
    fig.set_size_inches([5,5])


    tight_layout()


    def update_img(n):
        tmp = rand(300,300)
        im.set_data(tmp)
        return im

    #legend(loc=0)
    ani = animation.FuncAnimation(fig,update_img,300,interval=30)
    writer = animation.writers['ffmpeg'](fps=30)

    ani.save('demo.mp4',writer=writer,dpi=dpi)
    return ani

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_aspect('equal')
# ax.get_xaxis().set_visible(False)
# ax.get_yaxis().set_visible(False)
frames = imageio.volread(r'/home/baptiste/Documents/vis/190507_kif5a_nKTP.lif - Series011.tif')
# im = ax.imshow(frames[0])

# # fig.set_size_inches([5,5])


# tight_layout()
plt.show()
data = pd.read_csv(r'/home/baptiste/Documents/vis/190507_kif5a_nKTP.lif - Series011.tif_rejoined.csv',sep='\t')
subdata = data[data.rejoined_particle==18]

subframes = frames[subdata.frame.min():subdata.frame.max()]

fig, ax = plt.subplots()
ax.imshow(frames[0])
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
tight_layout()
# plt.show()

cmap = matplotlib.cm.get_cmap('Spectral')

cmap_list = []

for i in range(len(subdata)):
    print(i)
    cmap_list.append(cmap(i/len(subdata)))

print(cmap_list)

def updt(n):
    print(n)
    ax.cla()
    ax.imshow(subframes[n],cmap='gray')
    x = subdata.x[0:n]
    y = subdata.y[0:n]
    ax.scatter(x,y,color=cmap_list[0:n],s=1)
    return ax

# plt.show(updt(1))

def updt2(n):
    print(n)
    ax.imshow(frames[n])
    x = subdata.x[0:n]
    y = subdata.y[0:n]
    ax.scatter(x,y,color=cmap(n/len(subdata)),s=1)
    return ax
# # fig.show()

ani = animation.FuncAnimation(fig,updt,int(len(subdata)),interval=50,blit=False)
writer = animation.writers['ffmpeg'](fps=30)

ani.save('demo.mp4',writer=writer,dpi=dpi)
