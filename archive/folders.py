import os
from pathlib import Path
root = Path(r'/media/baptiste/SHG_tracking_data/COLOC')
for i in range(15):
    folder = root.joinpath("larve_"+str(i+1))
    os.makedirs(folder)
    os.makedirs(folder.joinpath("oeil_gauche"))
    os.makedirs(folder.joinpath("oeil_droit"))
#for path, subfolder, files in os.walk(root): #Scans entire folder structure for files
#    pathl = Path(path)
"""if path.endswith("gauche"):
        
        os.rename(path,pathl.parent.joinpath("oeil_gauche_stb_d"))
    if path.endswith("droit"):
        
        os.rename(path,pathl.parent.joinpath("oeil_droit_stb_g"))"""
"""    if path.endswith("_stb_d"):
        
        os.rename(path,pathl.parent.joinpath("oeil_droit"))
    if path.endswith("_stb_g"):
        
        os.rename(path,pathl.parent.joinpath("oeil_gauche"))"""

print(os.environ)
print(os.getenv("PROJECT_DIR"))