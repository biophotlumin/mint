import imagej
import jpype

ij = imagej.init()
jimage = ij.io().open(r'/media/lumin/DATA/Demo_BioProbe/Exp1_20190205_06_kif5a_nKTP/HET/larve3/oeil_gauche/low.odt')
frames = ij.py.from_java(jimage)
try:
    frames.shape
except AttributeError:
    print("No reader found")