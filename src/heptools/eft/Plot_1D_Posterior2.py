# For plotting 1D posterior distribution (along with global model and another ancillary data)

import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
import boost_histogram as bh


class _1D_Posterior_Plot:

    bin_width = 0.1

    def make_binning(self,low,high):
        return np.concatenate([np.arange(low,high,self.bin_width)  ,  np.asarray([high])])

    def __init__(self,posterior,limits,**kwargs):

        self.posterior = posterior  # This should be a numpy array 
        self.limits    = limits
     
        mini=round(min(self.posterior),1)
        maxi=round(max(self.posterior),1)

        # self.Primary_Hist = np.histogram(posterior, bins=self.make_binning(mini,maxi),density=True)
        # hep.histplot(self.Primary_Hist,fill=True)
        # plt.show()
        # input()

        hist = bh.Histogram(bh.axis.Regular(200,mini, maxi))

        hist.fill(self.posterior)
        self.Primary_Hist = hist/hist.sum()
 
        self.Region_Hists = {"_3sig_left"   : self.Primary_Hist[ bh.loc(mini)                     : bh.loc(self.limits["2sigma"][0])  ],
                             "_2sig_left"   : self.Primary_Hist[ bh.loc(self.limits["2sigma"][0]) : bh.loc(self.limits["1sigma"][0])  ],
                             "_1sig"        : self.Primary_Hist[ bh.loc(self.limits["1sigma"][0]) : bh.loc(self.limits["1sigma"][1])  ], 
                             "_2sig_right"  : self.Primary_Hist[ bh.loc(self.limits["1sigma"][1]) : bh.loc(self.limits["2sigma"][1])  ],
                             "_3sig_right"  : self.Primary_Hist[ bh.loc(self.limits["2sigma"][1]) : bh.loc(maxi)                      ]
        } 
 
    def make_multi_plot(self):


        # Coloured histogram bounds
        hep.histplot(self.Region_Hists["_3sig_left"],   color="red")
        hep.histplot(self.Region_Hists["_3sig_right"],  color="red")
        hep.histplot(self.Region_Hists["_2sig_left"],   color="yellow")
        hep.histplot(self.Region_Hists["_2sig_right"],  color="yellow")
        hep.histplot(self.Region_Hists["_1sig"],        color="lawngreen")

        # Black outline
        # hep.histplot(self.Primary_Hist,fill=False,color="black",linewidth=0.4)

        self.fig = plt.gcf() 
        self.ax  = plt.gca()



class _1D_LHC_Posterior_Plot(_1D_Posterior_Plot):


    allowed_experiment_styles = ["ATLAS","ALICE","LHCb2","CMS","ATLASTex"]

    def initialise_LHC_plot(self,experiment):

        assert any([experiment == x for x in self.allowed_experiment_styles]), "Experiment style not defined"

        if   experiment=="ATLAS": plt.style.use(hep.style.ATLAS) # This is the correct syntax for Pytohn3.6.9 version (mplhep 0.2.8)
        elif experiment=="CMS"  : plt.style.use(hep.style.CMS)
        elif experiment=="ALICE": plt.style.use(hep.style.ALICE)
        elif experiment=="LHCb2": plt.style.use(hep.style.LHCb2)
        elif experiment=="ATLASTex": 
            print("This is not supported because do not have a working TexLive distribution")
            # plt.style.use(hep.style.ATLASTex) # T

    def include_metadata(self,text,data,lumi,year):
        hep.atlas.label(text, data=data, lumi=lumi, year=year)

    def additional_label(self,text,**kwargs):
        plt.text(0.05,0.81,text,transform=self.ax.transAxes)


    def __init__(self,posterior,limits,experiment,**kwargs,):
        
        super().__init__(posterior,limits,**kwargs)
        self.experiment = experiment

        self.initialise_LHC_plot(experiment)


        
        






        # hist = bh.Histogram(bh.axis.Regular(10,-10.0, 10.0))
        # hist.fill(self.posterior)
        # print(hist)
        # input()
        # print(self.post)


    # def ATLAS_wrap(plot_func):
    #     def wrapper(self):
    #         plt.style.use(hep.style.ATLAS)
    #         plot_func(self)
    #     return wrapper



    # @ATLAS_wrap
    # def make_plot(self):
    #     plt.style.use(hep.style.ATLAS)
    #     N = hep.histplot(self.post_hist,color="black",linewidth=0.5,fill=True)#,self.bins)
    #     # patch = StepPatch(values=self.posterior,
    #     #           edges=range(-11, 11),
    #     #           label=('Patch derived underlying object\n'
    #     #                  'with default edge/facecolor behaviour'))


    #     # fig,ax = plt.subplots()
    #     # ax.add_patch(patch)
    #     print(N[0][0].__dict__)
    #     input()
    #     plt.show()
    #     input()

    # def pure_plt_plot(self,num_bins):

    #     N, bins, patches = plt.hist(self.posterior, num_bins)

    #     for i in range(0,33):
    #         patches[i].set_facecolor("red")
    #     for i in range(33,66):
    #         patches[i].set_facecolor("green")
    #     for i in range(66,num_bins):
    #         patches[i].set_facecolor("blue")

    #     plt.show()
    #     input()

    # # def plot_stairs_fill(self):

    #     plt.stairs(self.posterior,edges=)


        


