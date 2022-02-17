from heptools.histplot import PyHist_Class_2D
import matplotlib.pyplot as plt

def quick2D(ROOT_histogram,**kwargs):

    pyhist = PyHist_Class_2D.Histogram_Wrapper(ROOT_histogram,"Name")
    fig,ax,pc=pyhist.plot_2d(normed=True,**kwargs)
    # plt.show()
    # input()
    return pyhist,fig,ax,pc