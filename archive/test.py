"""import imageio
import numpy as np

file_path = r'/media/baptiste/SHG_tracking_data/Zebrafish data/paper2/210330_nktp_kif5a.lif - Series041.tif'

frames = imageio.volread(file_path)

print(frames.shape)
print(*frames.shape)
frames_array_init = np.zeros(frames.shape)
print(frames_array_init.shape)
print(len(frames))"""
import pandas as pd
import numpy as np
import scipy
import seaborn as sns
import matplotlib.pyplot as plt


"""data = pd.read_csv(r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20210921_182942 line average 4/124 Results - 20211029_133252 switch/Trajectory average parameters.csv',sep="\t")
wt = data[data.condition=='WT']
het = data[data.condition=='HET']
hom = data[data.condition=='HOM']
print(len(wt))
print(len(het))
print(len(hom))
print(len(data['switch']))
print(np.mean(wt.switch))
print(np.mean(het.switch))
print(np.mean(hom.switch))
print(len(wt.loc[wt.switch==0]))
print(len(het.loc[het.switch==0]))
print(len(hom.loc[hom.switch==0]))

sns.histplot(data=hom, x='switch')
plt.show()"""

"""data = pd.read_csv(r'/media/baptiste/SHG_tracking_data/Zebrafish data/Dyna_tri Results - 20210922_145805 line average 4/Dyna_tri Results - 20211029_151512/Trajectory average parameters.csv',sep="\t")
ctrl = data[data.condition=='CONTROL']
dyna = data[data.condition=='DYNAPYRAZOLE']
print(len(ctrl))
print(len(dyna))

print(len(data['switch']))
print(np.mean(ctrl.switch))
print(np.mean(dyna.switch))

print(np.median(ctrl.switch))
print(np.median(dyna.switch))

print(len(ctrl.loc[ctrl.switch==0]))
print(len(dyna.loc[dyna.switch==0]))

sns.histplot(data=dyna, x='switch')
plt.show()

sns.boxplot(x=data['condition'],
            y=data['switch'],  width=0.35, notch=True,
            data=data, showfliers =False,linewidth=1)
sns.despine(trim=True)
plt.show()"""

"""from utils import *
a = np.array((([2,2],[2,2],[2,2]),([2,2],[2,2],[2,2]),([2,2],[2,2],[2,2]),([2,2],[2,2],[2,2])))
print(a.shape)
b = np.expand_dims(a,axis=0)
print(b.shape)
c = b
np.append(b,c,axis=0)
print(b.shape)
print(line_average(a))"""

import imageio
from scipy.optimize import curve_fit

image = imageio.volread(r'/home/baptiste/1250 - ROI8_LA2.tif')
#from trajectory_calculations import _fit_spot_by_gaussian,_gaussian,_spot_moments 
print(image.shape)
"""print(_fit_spot_by_gaussian(img[0]))
print(np.sum(img[0]))
gaussian = _gaussian(1.06,215.44157434,6.64701936,6.87500787,1.75830112,2.10481685)
print(scipy.integrate.dblquad(gaussian,a=0,b=10,gfun=0,hfun=10))"""
"""img = np.array(image)
x = img[:,1]
y = img[:,0]
print(type(x))
print(y.shape)
popt,pcov = curve_fit(_gaussian,x,y)
print(popt)"""

"""print(_spot_moments(image)) 
intensités_fit=np.linspace(0,10,100)
gauss = _gaussian(*_spot_moments(image))
print(type(gauss))
axis = [0,1,2,3,4,5,6,7,8,9,10,11]
p,c = curve_fit(_gaussian,image[:,5],image[5])
print(p)"""
"""plt.plot(intensités_fit,gauss(intensités_fit,0),linestyle="--",color='black')
print(image[:,1])
#plt.scatter(x = axis,y=image[:,1])
plt.show()
"""
"""# 100 linearly spaced numbers
x = np.linspace(-np.pi,np.pi,100)

# the function, which is y = sin(x) here
y = np.sin(x)
print(y)

# plot the function
plt.plot(x,y, 'b')

# show the plot
plt.show()"""
ax = [0,1,2,3,4,5,6,7,8,9,10,11]
params = [ -3.16165109 ,206.7905355   , 5.58480279  , 6.06519139 ,  2.08386441 , 2.06363984]
gauss = _gaussian(*params)
X,Y = np.mgrid[-10:20.1:0.5, -10:20.1:0.5]
#print(intensités_fit.shape)
print(gauss(5.58480279,6.06519139))
npax = np.linspace(-5,15,1000)
plt.plot(npax,gauss(npax,6.06519139),linestyle="--",color='black')
plt.scatter(ax,(image[:,6]))
plt.show()

"""npax = np.linspace(-5,15,1000)
plt.plot(npax,gauss(npax,6.06519139),linestyle="--",color='black')
plt.scatter(ax,image[5])
plt.show()"""