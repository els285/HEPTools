from heptools.histplot import PyHist_Class_2D
import matplotlib.pyplot as plt

def quick2D(ROOT_histogram):

    pyhist = PyHist_Class_2D.Histogram_Wrapper(ROOT_histogram,"Name")
    fig,a,b=pyhist.plot_2d(normed=True)
    plt.show()
    input()