from Make1D_Posterior import make_1D_posterior_LHC_plot
import matplotlib.pyplot as plt

exit()
Plot_Object  =  make_1D_posterior_LHC_plot(WC="cHt",experiment="ATLAS",
                            samples_filename="../results_ptz_laurynas_ethantest/BATobjects/samples.jld2",
                            limits_filename="/home/ethan/EFTfitterSpinCorr.jl/nov15_eft_limits.json")

Plot_Object.include_metadata(text="",data=True,lumi=50,year=2017)
Plot_Object.ax.set_xlim(-50,50)
plt.show()
input()
