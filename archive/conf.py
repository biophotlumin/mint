import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from data_extraction import confinement, inst_velocity
from trajectory_calculations import minimization


"""vitesse = r'/media/baptiste/SHG_tracking_data/Zebrafish data/124 Results - 20220214_180447/124 Results - 20220215_113910/Per phase parameters.csv'
df = pd.read_csv(vitesse,sep='\t')

print(len(df.loc[df.phase==2]))

n_points = []

for i in set(df.trajectory.unique()):
    n_points.append(len(df.loc[df.trajectory==i]))

print(np.mean(n_points))
plt.hist(n_points, density=True, bins=30)
plt.show()"""

data = pd.read_csv(r'/home/baptiste/Documents/trajs/4630.txt',sep=',')
print(data)
data = data[['x','y']]
print(data)
#data.columns = ['x','y','Signal','Noise']
conf_threshold = 50*1e-3/0.05
data = minimization(data,{'px':0.173,'sigma':129})
r_conf = confinement(data.x,data.y,{'sliding_window':3})
v_inst = inst_velocity(data.x,data.y,{'dt':0.05})
size = len(data.x)
phase = np.zeros(size)

"""if (size>=2*3+1):
    for i in range(len(r_conf)):
        if (r_conf[i]>0.81):
                phase[i] = 2 #GO phase

    for i in range(size): #STOP phase refinment
        if (phase[i]==2) & (v_inst[i]<conf_threshold):
            print('ok')
            phase[i] = 0
else:
    
    for i in range(size): #STOP phase refinment
        phase[i]==2 
        if (phase[i]==2) & (v_inst[i]<conf_threshold):
            phase[i] = 0"""
for i in range(len(r_conf)):
    if v_inst[i]>=conf_threshold:
        phase[i]=2
    else:
        phase[i]=0

print((phase.shape))
data = pd.concat([data,pd.DataFrame(phase)],axis=1)
data.columns = ['x','y','conf']
print(data)
data_go = data.loc[data['conf']==2]
data_stop = data.loc[data['conf']==0]
#plt.plot(data.x,data.y)
plt.scatter(data_go.x,data_go.y,s=0.5)
plt.scatter(data_stop.x,data_stop.y,s=0.5)
plt.show()



"""df = pd.read_csv(r'/media/baptiste/SHG_tracking_data/Zebrafish data/ANCIEN + DONNÉES BRUTES/KTP-dyneine-mw20/données brutes Results - 20220221_184802/données brutes/200114/WT/larve_15/oeil_gauche/200114_nKTP_dynein.lif - Series105.tif/200114_nKTP_dynein.lif - Series105.tif_rejoined.csv',sep='\t')
for i in set(df.rejoined_particle.unique()):
    subdata = df.loc[df.rejoined_particle==i]
    subdata.to_csv(r'/home/baptiste/Documents/trajs/'+str(i)+'.txt')"""