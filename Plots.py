"""
Author: tzs820, s240670, Emil Henningsen

Plotting module
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.animation import FuncAnimation
import Lattice
import Hamiltonian
from ExperimentData import ExperimentData
from ParameterSet import ParameterSet
from numpy import ndarray, c_, r_
from cycler import cycler
from typing import Literal, Union
from collections.abc import Callable
from textwrap import wrap

class Plots:

    figures = {}

    def __init__(self, figures : dict = None):
        """
        Initialize plot manager. 
        """

        if not figures == None:
            self.figures = figures

        #Upon init, set matplotlib standards. 
        plt.rcParams.update({
            "text.usetex": True,
            "font.family": "Helvetica",
            "axes.labelsize": 16,
            "axes.labelweight": "bold",
            "axes.titlesize": 16,
            "axes.titlelocation": "center",
            "axes.prop_cycle": cycler(color=[
                "#a30000", #dark red
                "#6204b5", #dark purple
                "#000fb5", #dark blue 
                "#00a5ad", #dark cyan
                "#00ad26", #dark green
                "#ada200", #dark yellow
                "#ad6200", #dark orange
                "#fc5dfa", #pink
                "#fc4c4c", #light red
                "#985dfc", #light purple
                "#5d5dfc", #light blue
                "#5dfcfa", #light cyan
                "#83fc5d", #light green
                "#fcf45d", #light yellow
                "#fcb85d", #light orange
            ]),
            "axes.grid": True
        })

    """modules:"""

    def add(self, fig, name : str = None) -> None:
        """
        Add a figure to the instantiated manager container with given name (dict key). 
        """
        if name == None:
            number = len(self.figures)
            name = number
        self.figures[name] = fig

    def remove(self, name : str) -> None:
        """
        Remove specified figure from container (name = dict key)
        """
        self.figures.pop(name)

    def show(self) -> None:
        plt.show()

    """SET modules:"""

    def setInteractive(self, interactive : bool = False) -> None:
        """
        Set matplotlib mode to interactive. If false (standard), figures are not shown until show() is called.
        If true, figures are shown immediately after creation.
        """
        
        if interactive == False:
            plt.ioff()
        else:
            plt.ion()

    """GET modules:"""

    def isInteractive() -> bool:
        return plt.isinteractive()
            
    """Different standard type plots: """

    def plot(self, *args, **kwargs) -> list:
        return plt.plot(*args, **kwargs)

    def plotDipolesPlane(self, lat : Lattice, ax : plt.Axes = None, plane : Literal["xy", "xz", "yz"] = "xy", title : str = None, 
                         xlim : tuple[float, float] = None, ylim : tuple[float, float] = None, ham : Hamiltonian = None, 
                         index : int = None, legend : bool = True, drawColorbar : bool = True, returnAmpls : bool = False) -> tuple[plt.Figure, plt.Axes]:
        """
        TODO: Description
        """
        from numpy import array
        firstTime = False
        if ax is None:
            plt.figure()
            ax = plt.gca()
            firstTime = True    #Raise flag for colorbar
        fig = plt.gcf()

        ampl = array([1])  #just the same color, if no ham provided
        if not ham is None and not index is None:
            ampl = ham.getAmplNorm(index=index)

        pos, dir, pola = lat.getPositions(), lat.getDisplacements(), lat.getPolarizations()

        x = pos[:,0]    #All x coordinates
        y = pos[:,1]    #All y coordinates
        z = pos[:,2]    #All z coordinates

        polax = pola[:,0]
        polay = pola[:,1]
        polaz = pola[:,2]

        ax1, ax2 = None, None
        p1, p2 = None, None
        po1, po2 = None, None
        match(plane):
            case "xy":
                ax1, ax2 = "x", "y"
                p1, p2 = x, y
                po1, po2 = polax, polay
            case "xz":
                ax1, ax2 = "x", "z"
                p1, p2 = x, z
                po1, po2 = polax, polaz
            case "yz":
                ax1, ax2 = "y", "z"
                p1, p2 = y, z
                po1, po2 = polay, polaz

        cmap = matplotlib.cm.coolwarm
        norm = matplotlib.colors.Normalize(vmin=min(ampl), vmax=max(ampl))
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])
        ax.scatter(p1, p2, c=ampl, cmap=cmap, norm=norm, label="sites")
        ax.quiver(p1, p2, po1, po2, scale=15, width=0.005, color=cmap(norm(ampl)), pivot="mid")
        if firstTime and drawColorbar:
            plt.colorbar(sm, ax=ax, label=r"Amplitude norm of $|e_j> = |c_j|$")
        if not ylim is None:
            lower, upper = ylim
            ax.set_ylim(lower, upper)
        if not xlim is None:
            lower, upper = xlim
            ax.set_xlim(lower, upper)
        ax.set_xlabel(r"$\mathbf{\hat{%c}}$" % ax1, loc="right")
        ax.set_ylabel(r"$\mathbf{\hat{%c}}$" % ax2, loc="top")
        if title is None:
            title = r"Linear lattice of $N=50$ dipoles, $\frac{d}{\lambda_0}=0.3$"
        ax.set_title(title, wrap = True, size="large")
        if legend:
            plt.legend()
        
        if returnAmpls:
            return fig, ax, ampl
        else:
            return fig, ax
    
    #Plot rates manually
    def plotRates(self, N : int, d : float, rates : ndarray, ax : plt.Axes = None, scalex : str = "linear", scaley : str = "linear", title : str = None, label : str = None) -> plt.Figure:
        """
        TODO: Description
        """
        from textwrap import wrap

        if title == None: 
            title = "\n".join(wrap(r"$N = %s$ dipoles in linear lattice, polarized in z-direction, $\frac{d}{\lambda_0} = %s$" % (N, d), 60))

        if label is None:
            label = "decay rates"

        if ax is None:
            fig = plt.figure()
            ax = plt.gca()
        else:
            fig = ax.get_figure()
        
        ax.plot(range(1, N+1), rates, 'o', label=label)
        ax.set_yscale(scaley)
        ax.set_xscale(scalex)
        ax.set_xlabel(r"$\mathbf{\xi \in [1,%s]}$" % N, loc="right")
        ax.set_ylabel(r"$\mathbf{\Gamma_\xi / \Gamma_0}$", loc="top")
        ax.set_title(title, loc="center")
        ax.legend()

        return fig

    def plotRatesLat(self, lat : Lattice, ham : Hamiltonian, ax : plt.Axes = None, scalex : str = "linear", scaley : str = "linear", title : str = None, label : str = None) -> plt.Figure:
        """
        TODO: Description

        Same as plotRates(), but automatically extracts information from Lattice and Hamiltonian objects.
        """

        #Extract information:
        N = lat.getN()
        d = lat.getd()
        rates = ham.getDecayRates()

        fig = self.plotRates(N, d, rates, ax, scalex, scaley, title, label)

        return fig
    
    def xy_scatter(self, vals):
        x = vals.real
        y = -2 * vals.imag
        return x, y

    def scatterEigenvalues(self, N : int, d : float, vals : Union[ndarray, tuple], col : tuple = None, ax : plt.Axes = None, scalex : str = "linear", scaley : str = "linear", title : str = None, **kwargs) -> plt.Figure:
        """
        Complex eigenvalues are split into real and imaginary part for scatter plot

        --- INPUT

        vals - Union[ndarray(N,), tuple(ndarray, ndarray)] - complex eigenvalues or x,y values for plot
        col - all supported types of pyplot.scatter

        --- OUTPUT

        fig - plt.Figure
        """
        from textwrap import wrap
        from numpy import real, imag, sqrt

        if title == None: 
            title = "\n".join(wrap(r"$N = %s$ eigenvalues for linear chain with $\frac{d}{\lambda_0} = %s$" % (N, d), 120))

        if ax is None:
            fig = plt.figure()
            ax = plt.gca()
        else:
            fig = ax.get_figure()
        
        if type(vals) is tuple:
            x, y = vals
        else:
            x, y = self.xy_scatter(vals)
        #magnitudes = sqrt(x**2 + y**2) #irrelevant for scatter plot ?

        #special handling
        if col is not None:
            points = ax.scatter(x, y, color=col, **kwargs)
        else:
            points = ax.scatter(x, y, **kwargs)

        ax.set_yscale(scaley)
        ax.set_xscale(scalex)
        ax.set_xlabel(r"$J_\xi = Re(\lambda_\xi)$", loc="right")
        ax.set_ylabel(r"$\Gamma_\xi = - 2 \cdot Im(\lambda_\xi)$", loc="top")
        ax.set_title(title, loc="center")

        return fig, points
    
    def scatterUpdate(self, vals : Union[ndarray, tuple], pc : matplotlib.collections.PathCollection = None, colors : list = None, sizes : list = None, ax : plt.Axes = None):
        if ax is None:
            ax = plt.gca()

        if type(vals) is tuple:
            x, y = vals
        else:
            x, y = self.xy_scatter(vals)
        
        pc.set_offsets(c_[x, y])

        if not colors is None:
            pc.set(facecolor=colors)

        if not sizes is None:
            pc.set(sizes=sizes)

        #ax.axis("auto")
    
    def scatterEigenvaluesLat(self, lat : Lattice, ham : Hamiltonian, ax : plt.Axes = None, scalex : str = "linear", scaley : str = "linear", title : str = None) -> plt.Figure:
        """
        TODO: Description

        Same as scatterEigenvalues(), but automatically extracts information from Lattice and Hamiltonian objects.
        """

        #Extract information:
        N = lat.getN()
        d = lat.getd()
        vals, vecs = ham.getEigenDecomp()

        fig = self.scatterEigenvalues(N=N, d=d, vals=vals, ax=ax, scalex=scalex, scaley=scaley, title=title)

        return fig
    
    def scatterEigenvaluesAnimate(self, parameters : ParameterSet, data : ExperimentData, frames : int, interval : int = 50, animate_instructions : Callable = None, xydata : list = None, title : str = None) -> FuncAnimation:
        """
        Matplotlib FuncAnimation of complex eigenvalues
        
        ---INPUT

         - parameters (ParameterSet) - experiment parameters
         - data (dict) - experiment data as produced by Experiment

        OPTIONAL:

         - animate_instructinos (function(i)) - animation function ... instructions to perform at every animation index besides updating x,y points
         - xy (list of complex values) - set of xpoints and ypoints to be obtained from complex (eigen)values

        ---OUTPUT

        animation (matplotlib.FuncAnimation)

        """

        #Beware, if the join function cuts off the string in middle of math field, latex will complain.
        if title is None:
            mytitle = "\n".join(wrap(r"Animated Eigenvalues", 100))
        else:
            mytitle = title

        """ unpack parameters and data """

        #parameters, for now just default to p_start
        N = parameters["N"][0]
        d = parameters["distance"][0]

        #data
        if xydata is None:
            xydata = data.getEigenvaluesTotal()
        ra = xydata.real
        ia = -2*xydata.imag

        fig, ax_data = plt.subplots()

        ax_data.set(xlim=(np.amin(ra)-1, np.amax(ra)+1), ylim=(np.amin(ia)-1, np.amax(ia)+1))

        #init scatter:
        fig, points = self.scatterEigenvalues(N, d, (ra, ia), ax = ax_data, title=mytitle)

        """animating motion"""
        def animator(i):
            #i = np.where(angles == i)
            xy = (ra[i], ia[i]) #NOTE: det har pludselig ændret sig, før behøvede jeg [i][0] for at refere til hele arrayet, men nu referer det kun til første punkt??? 
            #NOTE 2: DET er selvfølgelig fordi np.where(angles == i) returnere [index] og ikke bare index.
            self.scatterUpdate(xy, points, ax=ax_data)

            """ additional frame update instructions from animate """
            if not animate_instructions is None:
                animate_instructions(i)
            #fig.canvas.draw_idle()
            #return points, 

        animat = FuncAnimation(fig, animator, repeat=True, frames = frames-1, interval = interval)

        return animat
    
    def legend():
        return plt.legend()