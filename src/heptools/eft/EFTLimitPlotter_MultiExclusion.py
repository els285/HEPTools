# EFT plotter

# Design features - swap orientation: vertical or horizontal
# Should be run on Python3.6 / Python3.7

# Pass a pandas dataframe with a set of Wilson coefficients ad corresponding values

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import mplhep as hep


from heptools.eft.EFToperators import dic_of_EFT_operators as dEFT
"""
Structure of the DataFrame
WC: Linear, Linear+Quadratic, ...
Each part should then be a list of tuples: if there is a single region found, more if multiple (per operator)
Each tuple should be of the form (global mode, lower bound, upper bound)
"""
"""
The data format input could be a json file nested in the structure (EFT model,Wilson coefficient,number of bounds), 
    with each entry containing (value, lower bound, upper bound).


Plotting is limited to 4 different types of fits at once, those this can be extended though must be hard-coded.
The MultiExclusion works only for a specific type of DataFrame    
"""



class EFT_Plot:

    """
    Main EFT Plot Class
    Takes a pandas DataFrame as input
    """

    def __init__(self,df,**kwargs):

        self.df = df
        self.orientation = kwargs["orientation"] if "orientation" in kwargs else "vertical"
        self.colours = kwargs["colours"] if "colours" in kwargs else plt.rcParams['axes.prop_cycle'].by_key()['color']
        self.all_fits = set(self.df.columns.get_level_values(0))
        self.plot_zero_line = kwargs["plot_zero_line"] if "plot_zero_line" in kwargs else False

        plot_design_dict = {"fit legend marker":"line",
                        "global mode on legend":True,
                        "legend 68 95" : True   ,
                        "inner capsize": 5,
                        "outer capsize": 7,
                        "outer line style":"dashed",
                        "outer alpha": 0.7,
                        "inner linewidth":2,
                        "outer linewidth":2,
                        "global mode size": 20,
                        "inner capstyle": "butt"}

        self.plot_design_dict = {**plot_design_dict,**kwargs["plot_design_dict"]} if "plot_design_dict" in kwargs else plot_design_dict

        # Assign which columns to plot
        if "to_plot" in kwargs:
            if kwargs["to_plot"]=="all": 
                self.to_plot = self.all_fits
            if isinstance(kwargs["to_plot"],list):# and all([x in self.all_fits for x in kwargs["to_plot"]]):
                self.to_plot = kwargs["to_plot"]
        else:
            self.to_plot = self.all_fits

        self.number2plot = len(self.to_plot)
        self.generate_wc_ordinates()


    def get_points(self,tup):

        """Takes a tuple and computes the midpoint for plt.errorbar, returns bounds+midpoint as np.array"""

        get_value = lambda u,l: 0.5*(u-l)
        lower_bound = tup[0]
        upper_bound = tup[1]
        plot_value  = get_value(upper_bound,lower_bound)

        lower_error = plot_value - lower_bound
        upper_error = upper_bound - plot_value

        error_array = np.array([[lower_error ,upper_error]]).T

        return plot_value,error_array


    def generate_wc_ordinates(self):


        """ For plotting several different fits for each operator, need to introduce vertical separation """
        assert self.number2plot < 5 and isinstance(self.number2plot,int), "Too many branches to plot"
        if self.number2plot==1:
            self.wc_array = [0]
        if self.number2plot==2:
            self.wc_array = [-0.05,0.05]
        if self.number2plot==3:
            self.wc_array = [-0.2,0,0.2]
        if self.number2plot==4:
            self.wc_array = [-0.15,-0.05,0.05,0.15]


    def make_plot(self,**kwargs):

        """Makes plot defined by orientation, which can be pre-set or passed to this function"""

        if "orientation" in kwargs: self.orientation=kwargs["orientation"]
        assert self.orientation, "Orientation of plot not defined"
        if self.orientation=="horizontal":
            fig,ax=self.make_horizontal_plot()
            if "plot_global_modes" in kwargs and kwargs["plot_global_modes"]:
                self.plot_global_mode_horizontal()

        elif self.orientation=="vertical":
            fig,ax=self.make_vertical_plot()
            if "plot_global_modes" in kwargs and kwargs["plot_global_modes"]:
                self.plot_global_mode_vertical()
            
        return fig,ax


    def make_vertical_plot(self):
        
        fig, ax = plt.subplots()
        self.fig = fig
        self.ax = ax
        
        y_points = np.arange(0,self.df.shape[0]+1)

        if self.plot_zero_line:
            ax.vlines(0,-100,y_points[-1]+0.4,linewidth=1,color='darkgrey',zorder=0)        


        for col,xp,color in zip(self.to_plot,self.wc_array,self.colours):       # Looping over fit columns e.g Linear,Independent etc.
            for i,(index,row) in enumerate(self.df.iterrows()):                 # Looping over Wilson coefficients
                _68_data , _95_data = row[col].Bounds["68"] , row[col].Bounds["95"]
                for tup in _68_data:
                    if tup !=None:
                        value,asymmetric_error = self.get_points(tup)
                        p1,p2,p3=plt.errorbar(x=value,y=i+1+xp,xerr=asymmetric_error,
                                    capsize=self.plot_design_dict["inner capsize"],
                                    color=color,
                                    elinewidth=self.plot_design_dict["inner linewidth"])
                        plt.setp(p3[0],capstyle=self.plot_design_dict["inner capstyle"])#, color="orange")        
                                    
                for tup in _95_data:
                    if tup !=None:
                        value,asymmetric_error = self.get_points(tup)
                        p1,p2,p3 = plt.errorbar(x=value,y=i+1+xp,xerr=asymmetric_error,
                                        capsize=self.plot_design_dict["outer capsize"],
                                        color=color,
                                        alpha=self.plot_design_dict["outer alpha"],
                                        elinewidth=self.plot_design_dict["outer linewidth"])
                        p3[0].set_linestyle(self.plot_design_dict["outer line style"])


        ax.set_yticks(y_points[1:])
        operator_names_latex = [dEFT[op] for op in self.df.index]
        ax.set_yticklabels(operator_names_latex)
        ax.set_ylim(0.5,len(y_points)+1)

        ax.tick_params(axis='y', which='minor', left=False, right=False)

        self.make_legend()
        plt.xlabel(r"$C/\Lambda^{2} \, \,  [\mathrm{TeV}^{-2}]$")
        # plt.ylabel("Wilson Coefficients")

        return fig,ax


    def make_horizontal_plot(self):

        # Separate by columns based on fit, then loop over rows of Wilson coefficients
        # Don't loop over different bounds, instead hard-coded _68_data and 95_data bounds
        # Possible extension to do this more flexibly

        fig, ax = plt.subplots()
        self.fig = fig
        self.ax = ax

        x_points = np.arange(0,self.df.shape[0]+1)

        if self.plot_zero_line:
            ax.hlines(0,-100,x_points[-1]+0.4,linewidth=1,color='darkgrey',zorder=0)        

        for col,xp,color in zip(self.to_plot,self.wc_array,self.colours):       # Looping over fit columns e.g Linear,Independent etc.
            for i,(index,row) in enumerate(self.df.iterrows()):                 # Looping over Wilson coefficients
                _68_data , _95_data = row[col].Bounds["68"] , row[col].Bounds["95"]
                for tup in _68_data:                                                # Reuiq
                    if tup !=None:
                        value,asymmetric_error = self.get_points(tup)
                        p1,p2,p3=plt.errorbar(x=i+1+xp,y=value,yerr=asymmetric_error,
                                    capsize=self.plot_design_dict["inner capsize"],
                                    color=color,
                                    elinewidth=self.plot_design_dict["inner linewidth"])
                        plt.setp(p3[0],capstyle=self.plot_design_dict["inner capstyle"])#, color="orange")     

                for tup in _95_data:
                    if tup !=None:
                        value,asymmetric_error = self.get_points(tup)
                        p1,p2,p3 = plt.errorbar(x=i+1+xp,y=value,yerr=asymmetric_error,
                                            capsize=self.plot_design_dict["outer capsize"],
                                            color=color,
                                            alpha=self.plot_design_dict["outer alpha"],
                                            elinewidth=self.plot_design_dict["outer linewidth"])
                                            
                        p3[0].set_linestyle(self.plot_design_dict["outer line style"])


        ax.set_xticks(x_points[1:])
        operator_names_latex = [dEFT[op] for op in self.df.index]
        ax.set_xticklabels(operator_names_latex)
        ax.set_xlim(0.5,len(x_points)+1.6)

        ax.tick_params(axis='x', which='minor', bottom=False, top=False)

        self.make_legend()

        # plt.xlabel("Wilson Coefficients")
        plt.ylabel(r"$C/\Lambda^{2} \, \, [\mathrm{TeV}^{-2}]$")
        return fig,ax



    def make_legend(self):

        """
        Legend can include either coloured lines or errorbars, plus additional global mode dots and information on 68 95 lines
        """


        from matplotlib.lines import Line2D
        from matplotlib.patches import Circle

        handles, labels = plt.gca().get_legend_handles_labels()

        # Straight line legend markers
        if self.plot_design_dict["fit legend marker"]=="line":
            for attrib,colour in zip(self.to_plot,self.colours):
                line = Line2D([0], [0], label=attrib, color=colour)
                handles.extend([line])

        if self.plot_design_dict["fit legend marker"]=="errorbar":
            from matplotlib.container import ErrorbarContainer
            from matplotlib.collections import LineCollection
            for attrib,colour in zip(self.to_plot,self.colours):
                line = Line2D([],[], ls="none",color=colour,label=attrib)
                barline = LineCollection(np.empty((2,2,2)),color=colour)
                if self.orientation=="horizontal":
                    err = ErrorbarContainer((line, [line], [barline]), has_xerr=False, has_yerr=True,label=attrib)
                if self.orientation=="vertical":
                    err = ErrorbarContainer((line, [line], [barline]), has_xerr=True, has_yerr=False,label=attrib)
                handles.extend([err])
        handles.reverse()


        # Global model dot
        if self.plot_design_dict["global mode on legend"]:
            global_dot = Line2D([], [], color='dimgrey', marker='o', linestyle='None',
                            markersize=5, label='Global mode')
            handles.extend([global_dot])# input()        

        if self.plot_design_dict["legend 68 95"]:
            line_inner = Line2D([0], [0], label=r"68% confidence", color="dimgrey",linestyle="solid")
            line_outer = Line2D([0], [0], label=r"95% confidence", color="dimgrey",linestyle=self.plot_design_dict["outer line style"],
                                alpha=0.7)
            handles.extend([line_inner,line_outer])


        if self.orientation   == "vertical":     plt.legend(handles=handles,loc="upper right",ncol=2)

        elif self.orientation == "horizontal" :  plt.legend(handles=handles,loc="center right")


    def plot_global_mode_vertical(self):

        assert self.orientation=="vertical", "Cannot plot vertical global modes when main plot is horizontal"
        for col,yp,color in zip(self.to_plot,self.wc_array,self.colours): # Loop over columns
            for i,(index,row) in enumerate(self.df.iterrows()):
                global_mode_value=row[col]['Global Mode']
                plt.scatter(x=global_mode_value,y=i+1+yp,color=color,s=self.plot_design_dict["global mode size"])


    def plot_global_mode_horizontal(self):

        assert self.orientation=="horizontal", "Cannot plot horizontal global modes when main plot is vertical"
        for col,yp,color in zip(self.to_plot,self.wc_array,self.colours): # Loop over columns
            for i,(index,row) in enumerate(self.df.iterrows()):
                global_mode_value=row[col]['Global Mode']
                plt.scatter(y=global_mode_value,x=i+1+yp,color=color,s=self.plot_design_dict["global mode size"])


    def add_known_bounds(self):
        pass

        


class LHC_EFT_Plot(EFT_Plot):

    """
    A daughter class for making the EFT plots in the style of the LHC experiments
    """

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


    def __init__(self,df,experiment,**kwargs,):
        
        super().__init__(df,**kwargs)
        self.experiment = experiment

        self.initialise_LHC_plot(experiment)



class ATLAS_HiggsStyle_EFT_Plot(LHC_EFT_Plot):

    "A daughter daughter class which makes the EFT plot in the style of ATLAS Higgs group"

    default_colours = ["#2754B5","#EB3223","#EDBE5F"] #blue,red,yellow

    default_AH_plt_dict = {"inner capsize":0,
                           "outer capsize": 0, 
                           "outer line style": "dotted",
                           "outer alpha":1,
                           "inner linewidth": 3,
                           "outer linewidth": 1.5,
                           "global mode size": 120,
                           "inner capstyle":"round"}

    # plot_design_

    def __init__(self,df,**kwargs):
       
        super().__init__(df,experiment="ATLAS",colours=self.default_colours,plot_design_dict=self.default_AH_plt_dict,**kwargs)
        # self.experiment = experiment

        # self.initialise_LHC_plot(experiment)














