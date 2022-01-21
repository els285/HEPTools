
import matplotlib.pyplot as plt
from heptools.eft.ParseFromTXT import auto_parse
from heptools.eft.EFTLimitPlotter_MultiExclusion import ATLAS_HiggsStyle_EFT_Plot, LHC_EFT_Plot,EFT_Plot

x = auto_parse("/home/ethan/EFTfitterSpinCorr.jl/results_ptz_laurynas")

# print(x)
# input()


y = LHC_EFT_Plot(df=x,experiment="ATLAS",#to_plot="all",
                to_plot=["Linear+Quadratic (Marg.)","Linear (Marg.)" ,"Linear+Quadratic (Indp.)"],
                orientation="vertical",
                plot_zero_line=False,
                colours=["blue","red","green"])

# To define the same order every time, it is necessary to specify the to_plot field in the above

# fig,ax=y.make_horizontal_plot()
# y.plot_global_mode_horizontal()

fig,ax=y.make_plot(orientation="vertical",plot_global_modes=True)
# print(y.__dict__)
# input()

y.include_metadata("Internal",True,139,2017)
y.additional_label(r"SMEFT $\Lambda = 1$ TeV")

# ax.set_xlim(-15,20)


plt.show()
input()

fig.set_size_inches(10, 6)

# plt.savefig("lhc_plot.png",dpi=300)