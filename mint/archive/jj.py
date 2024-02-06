import imagej
import jpype

ij = imagej.init()
jimage = ij.io().open(r'/media/lumin/DATA/Demo_BioProbe/Exp1_20190205_06_kif5a_nKTP/HET/larve3/oeil_gauche/low.odt')
frames = ij.py.from_java(jimage)
try:
    frames.shape
except AttributeError:
    print("No reader found")

def barplot(data,variable,input_folder,p,parameters):
    """Bar plot.

    Generate a bar plot, with SEM as error bars, and saves it as a file under `extension_out`.

    :param data: DataFrame containing trajectory parameters.
    :type data: DataFrame
    :param variable: Current variable.
    :type variable: str
    :param input_folder: Folder where the file will be saved.
    :type input_folder: str or Path
    :param p: p-value.
    :type p: float
    :param parameters: Dictionary containing parameters.
    :type parameters: dict
    """

    error = []
    for i in data['condition'].unique():
        error.append(stats.sem(data[variable].loc[data['condition']==i],nan_policy='omit') if len(data)>3 else 0)
    sns.barplot(y=data[variable].dropna(),x=data['condition'],estimator=mean,yerr=error, errorbar=None,\
            error_kw={'elinewidth':2,'capsize':4,'capthick':2})
    sns.despine(trim=True)
    plt.xlabel("Condition")
    plt.ylabel(stats_vars[variable])
    plt.annotate((f'p-value : {p}'),xy=(195,310),xycoords='figure points')
    plt.savefig(Path(input_folder).joinpath(f'Barplot {str(data["condition"].unique())} {variable}.{parameters["extension_out"]}'),dpi=parameters["dpi"])
    plt.close()