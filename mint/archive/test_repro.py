import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

# D1 = [5,5,5,10,np.nan,10,15,15,15]
# D2 = [5,5,5,10,10,10,15,15,15]
# D3 = [5,5,5,10,10,10,15,15,15]
# CAT = ['CAT1','CAT1','CAT1','CAT2','CAT2','CAT2','CAT3','CAT3','CAT3']

# data = pd.DataFrame({'D1':D1,
#                      'D2':D2,
#                      'D3':D3,
#                      'CAT':CAT})

# print(data)

# sns.barplot(data=data, y=data['D1'].dropna(), x=data['CAT'],errorbar=None,error_kw={'elinewidth':2,'capsize':4,'capthick':2})
# plt.savefig('test.png')


def kruskal(data,variable):
    """Kruskal-Wallis test.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :return: p-value.
    :rtype: float
    """    

    list_of_arrays = []
    for index in data.condition.unique():
        list_of_arrays.append(data.loc[data.condition==index, variable])

    p = stats.kruskal(*list_of_arrays,nan_policy='omit')[1]
    
    return p

variable = 'pausing_frequency'

data = pd.read_csv(r'/media/lumin/DATA/Demo_BioProbe Results - 20231107_115319 mintj/Demo_BioProbe Results - 20231107_120339/20231107_120339_Trajectory average parameters.csv',sep='\t')
print(kruskal(data,variable))
p_value = kruskal(data,variable)
print(str(round(p_value,6)))

pf = pd.read_csv(r'/home/lumin/Documents/GitHub/mint/mint/pf.csv')
p_valuep = kruskal(pf,variable)
print(str(round(p_valuep,6)))