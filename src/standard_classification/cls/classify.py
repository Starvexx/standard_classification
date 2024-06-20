from __future__ import annotations

import numpy as np

import astropy
from astropy import units as u

from scipy import optimize
from scipy import odr

from matplotlib import pyplot as plt
from matplotlib.colors import TABLEAU_COLORS as colors

from tqdm import tqdm


class star:
    """Star class.
    
    The star class is used to determine the infrared alpha index of a
    young stellar object (YSO).

    Attributes
    ----------
    srcID : int
        The unique source identifier of the YSO (from input catalog).
    ra : float
        Right Ascension of the YSO.
    dec : float
        Declination of the YOS
    alpha : dict
        Dictionary containing the alpha indices calculated for this YSO
    intercept : dict
        Dictionary holding the intercepts of the line fits. These can be
        used to plot the fitted line to the SED of the YSO.
    intercept_n : dict
        Dictionary holding the intercepts from the alternative alpha
        index method.
    cls : dict
        Dictionary holding the determined observational class of the
        YSO.
    fluxNames : list
        A list containing the names of the flux columns from the input
        catalog.
    fluxErrNames : list
        A list containing the names of the flux error columns from the
        input catalog.
    lambdaNames : list
        A list containing the names of the wavelength columns from the
        input catalog.
    """

    def __init__(self, source : astropy.table.Row):
        """Class constructor.
        
        Constructs the instance of a star.

        Parameters
        ----------
            source : astropy.table.Row
        """
        self.__data = source
        self.srcID = self.__data["Internal_ID"]
        self.ra = self.__data["RA"]
        self.dec = self.__data["DE"]
        self.alpha = {}
        self.n = {}
        self.intercept = {}
        self.intercept_n = {}
        self.cls = {}
        self.fluxNames = []
        self.fluxErrNames = []
        self.lambdaNames = []

        #Get the names of the different columns holding the relevant
        #data. To compute the alpha index we need all columns
        #containing flux, flux errors, as well as their respective
        #wavelengths, lambda.

        for name in self.__data.colnames:
            if ('flux' in name and 'error' not in name):
                self.fluxNames.append(name)
            elif 'flux_error' in name:
                self.fluxErrNames.append(name)
            elif 'lambda' in name:
                self.lambdaNames.append(name)


        # Extract the columns holding the flux measurements (Jy) from the table.
        self.fluxDens = np.array([col for col in self.__data[self.fluxNames].as_void()])
        self.fluxDensErrs = np.array([col for col in self.__data[self.fluxErrNames].as_void()])
        self.wlngths = np.array([col for col in self.__data[self.lambdaNames].as_void()])

        self.fluxDens[self.fluxDens == 0] = np.nan
        self.fluxDensErrs[self.fluxDensErrs == 0] = np.nan
        self.wlngths[self.wlngths == 0] = np.nan

        self.fluxDens2flux()
        self.fluxDensErr2fluxErr()



    def __line(self, x : float, k : float, d : float) -> float:
        """Just a boring old line, nothing else to see here.
        
        The line is used for the least squares fit to the SED to
        estimate the infrared spectral index of a YSO.

        Prameters
        ---------
            x : float
                The wavelength (x coordinate) at which the line is
                evaluated.
            k : float
                The slope of the line.
            d : float
                The intercept of the line.

        Returns
        -------
            float
                The flux (y-value) at the respective wavelength (x)
                of the line.
        """
        return (k * x) + d


    def __powerlaw(self, x : float, k : float, d : float) -> float:
        """Powerlaw representing a line in log-log space.
        
        This method is the powerlaw representation of a line in double
        logarithmic space. This workaround is needed to determine the
        alpha index with the ODR method which includes errors in the
        line fit.
        
        Parameters
        ----------
            param : list of floats
                This list contains slope and intercept parameters of
                the line in log/log space.
            x : float
                The wavelength at which the power law is evaluated.

        Returns
        -------
            float
                The powerlaw evaluated at wavelength x.
        """
        if (x <= 0).any():
            x_copy = x.copy()
            x_copy[x <= 0] = np.nan
            return 10 ** (k * np.log10(x_copy) + d)
        else:
            return 10 ** (k * np.log10(x) + d)
    

    def __wlClose(self, wl1: float, wl2: float, threshold: float) -> bool:
        """Tests if the two wavelengths are very close to each other.
        
        If the two wavelengths are within a threshold in % of each
        other, if they are the only two wavelengths in a range of
        wavelengths for which the alpha index is to be determined,
        they retrieved spectral index may be very far off from the 
        true value. Therefore, the two measurements need to be
        separated by a minimum distance wavelength wise to give a
        viable result. This method checks if the two measurements
        lie within the % threshold window of each other.

        Parameters
        ----------
            wl1 : float
                The first wavelength.
            wl2 : float
                The second wavelength.
            threshold : 
                The threshold in percent that defines wether the two
                measurements are too close to each other.
            
        Returns
        -------
            tooClose : bool
                Returns True if the two measurements are too close
                to each other and False otherwise.
        """
        if 0 < threshold < 100:
            rel_wl_diff = 1 - (wl1 / wl2)
            if rel_wl_diff < (threshold / 100):
                return True
            else:
                return False
        else:
            raise ValueError("threshold must be a real number between 0 and 100")

   
    def fluxDens2flux(self):
        """Converts flux densities to fluxes.
        
        Does the conversion from Jansky (Jy) to erg / s / cm^2 / um

        Returns
        -------
            fluxes : float
                The converted fluxes. 
        """

        self.fluxes = self.fluxDens * 1e-23 * 2.99792458e14 / (self.wlngths**2) * u.erg / u.s / u.cm**2 / u.um
        
        
        #self.fluxes = (self.fluxDens * u.Jy).to(u.erg / u.s / u.cm**2 / u.um,
        #    equivalencies=u.spectral_density(self.wlngths * u.um))

        #return self.fluxes


    def fluxDensErr2fluxErr(self):
        """Converts flux density errors to flux errors.
        
        Does the conversion from Jansky (Jy) to erg / s / cm^2 / um
        
        Returns
        -------
            fluxErrs : float
                Returns the converted flux errors.
        """
        
        self.fluxErrs = self.fluxDensErrs * 1e-23 * 2.99792458e14 / (self.wlngths**2) * u.erg / u.s / u.cm**2 / u.um
        
        #self.fluxErrs = (self.fluxDensErrs * u.Jy).to(u.erg / u.s / u.cm**2 / u.um,
        #        equivalencies=u.spectral_density(self.wlngths * u.um))

        #return self.fluxErrs


    def getAlphaWithErrors(self, lower : float = None, upper : float = None, Mask : np.ndarray = None, key : str = '') -> tuple:
        """Computes the alpha index considering errors.
        
        Computes the infrared spectral index also regarding measurement
        errors. The method used is orthogonal distance regression (ODR)
        which fits a line to the SED within a defined wavelength range.

        The ODR algorithm used here only works with symmetrical errors.
        Therefore, fitting is done in non logarithmic space here the
        measurement errors are symmetrical unlike to log-log space where
        they are asymmetrical. This however requires a non linear
        fitting function which is parametrized to represent a line in
        double logarithmic space needed to determine the alpha index.

        Parameters
        ----------
            lower : float
                The lower boundary of the wavelength range within
                which the alpha index is computed.
            upper : float
                The upper boundary of the wavelength range within
                which the alpha index is computed.
        
        Returns
        -------
            alpha : float
                The slope of the line in log-log space determining the
                infrared spectral index for YSO classification.
            intercept : float
                The intercept of the line in log-log space. This can be
                used to over plot the fitted function on the SED.
        """
        if Mask is None and (lower is not None or upper is not None):
            # Get the selection masks for the selected wavelength range.
            wlRangeMask = (self.wlngths > lower) & (self.wlngths < upper)
            # Get the selection mask for all measurements that have errors.
            hasError = ~np.isnan(self.fluxErrs)
            # Get the selection mask for the data used in the ODR fit.
            mask = wlRangeMask & hasError
            # Store the selection mask for the data in the selected wavelength range.
            self.wlMask = wlRangeMask
            key = f'{lower}-{upper}'
        else:
            print("using mask")
            mask = Mask

        print(key, np.sum(mask))
        
        # Check if there are enough data points to compute the alpha index.
        # If there are less than two return NaN values.
        # If there are only two measurements but there is at least one
        # measurement error missing, use the estimated vale from scipy optimize.
        # If there are at least two measurements with errors, use only those
        # that provide measurement errors to compute the alpha index with ODR.

        self.log_wl = np.log10(self.wlngths)
        
        popt, pcov = optimize.curve_fit(self.__powerlaw,
                                        self.wlngths[mask],
                                        self.fluxes[mask].value,
                                        sigma=self.fluxErrs[mask].value)
        
        self.alpha[f"{key}"] = popt[0]
        self.intercept[f"{key}"] = popt[1]

        #if (np.sum(mask) == 2):
        #    if self.__wlClose(self.wlngths[mask][0],
        #                      self.wlngths[mask][1],
        #                      5):
        #        self.alpha[f"{lower}-{upper}"] = np.nan
        #        self.intercept[f"{lower}-{upper}"] = np.nan
        #if (np.sum(wlRangeMask) <= 1) or (np.sum(mask) <= 1):
        #    # Less than two measurements with errors. Using estimate.
        #    self.alpha[f"{lower}-{upper}"] = np.nan
        #    self.intercept[f"{lower}-{upper}"] = np.nan
        #elif np.any(self.fluxErrs[wlRangeMask].value) and (np.sum(fullDataMask[wlRangeMask]) > 2):
        #    # Has at least two measurements with errors. Using LSR on data
        #    # with errors.
        #    #mask = wlRangeMask & fullDataMask
        #    self.log_wl = np.log10(self.wlngths[mask])

        #    popt, pcov = optimize.curve_fit(self.__powerlaw,
        #                                    self.wlngths[mask],
        #                                    self.fluxes[mask],
        #                                    sigma=self.fluxErrs[mask])
        #    self.alpha[f"{lower}-{upper}"] = popt[0]
        #    self.intercept[f"{lower}-{upper}"] = popt[1]
        #else:
        #    # All measurements have errors. Using full LSR.
        #    self.log_wl = np.log10(self.wlngths[wlRangeMask])

        #    popt, pcov = optimize.curve_fit(self.__powerlaw,
        #                                    self.wlngths[wlRangeMask],
        #                                    self.fluxes[wlRangeMask],
        #                                    sigma=self.fluxErrs[wlRangeMask])
        #    self.alpha[f"{lower}-{upper}"] = popt[0]
        #    self.intercept[f"{lower}-{upper}"] = popt[1]
        
        self.cls[key] = self.classify(self.alpha[key])
        return self.alpha, self.intercept


    def altAlpha(self, lower : list, upper : list) -> tuple:
        """Work in progress...
        
        This functionality is not yet implemented.

        Parameters
        ----------
            lower : list
                The lower boundary of the wavelength range within
                which the alpha index is computed.
            upper : list
                The upper boundary of the wavelength range within
                which the alpha index is computed.

        Returns
        -------
            alpha : float
                The slope of the line in log-log space determining the
                infrared spectral index for YSO classification.
            intercept : float
                The intercept of the line in log-log space. This can be
                used to over plot the fitted function on the SED.
            
        """
        lower_mask = (self.wlngths > lower[0]) & (self.wlngths < lower[1])
        upper_mask = (self.wlngths > upper[0]) & (self.wlngths < upper[1])
        
        lower_wlngths = self.wlngths[lower_mask]
        upper_wlngths = self.wlngths[upper_mask]

        lower_flux = self.fluxes[lower_mask]
        upper_flux = self.fluxes[upper_mask]
        
        if (len(lower_flux) < 1) & (len(upper_flux) < 1):
            self.n[f'{lower[0]}-{upper[-1]}'] = np.nan
            self.intercept_n[f'{lower[0]}-{upper[-1]}'] = np.nan
        else:
            n = np.log10(upper_flux[0] / lower_flux[0]) / np.log10(lower_wlngths[0] / upper_wlngths[0])
            self.n[f'{lower[0]}-{upper[-1]}'] = n
            
            print(n)


    def classify(self, alpha : float) -> str:
        """Classification method.
        
        From the computed alpha index, the observational class is
        inferred as per Grossschedl et.al., (2016).

        Parameters
        ----------
            alpha : float
                The infrared spectral index.
        
        Returns
        -------
            YSO class : string
                The observational YSO class.
        """
        if (0.3 < alpha):
            return "0/I"
        elif ((-0.3 < alpha) & (alpha < 0.3)):
            return "flat"
        elif ((-1.6 < alpha) & (alpha < -0.3)):
            return "II"
        elif ((-2.5 < alpha) & (alpha < -1.6)):
            return "III thin disk"
        elif (alpha < -2.5):
            return "III no disk / MS"
        else:
            return "not classified"
    

    def plot(self, savepath, keys=None, lower=2, upper=24):
        """Plots the data.

        This methods plots the SED of the source including any the
        lines of the fit for the alpha indices. 

        Parameters
        ----------
            savepath : str
                The path to the directory where the plots should be
                stored to.
            lower : array_like
                The lower limits of the ranges for the alpha indices.
            upper : array_like
                The upper limits of the ranges for the alpha indices.
        """
        # If no values for the lower and upper boundary are provided, 
        if (keys == None):
            keys = list(self.alpha.keys())
                    
        savepath = savepath + f'/{self.srcID:05d}.png'


        fig = plt.figure(figsize=(4,4), dpi=150)

        font = {'size'  : 7,
                'family': 'serif'}
        plt.rc('font', **font)
        plt.rc('text', usetex=True)

        ax = fig.add_subplot(111)

        for i, key in enumerate(keys):
            if i == 0:
                c = 'k'
                a = 1
            else:
                c = colors[list(colors.keys())[i]]
                a = 0.5
            
            try:
                if ~np.isnan(self.alpha[key]):
                    wl_range = np.linspace(0.5, 1000, 100)

                    ax.plot(wl_range,
                            self.__powerlaw(wl_range,
                                            self.alpha[key],
                                            self.intercept[key]),
                            lw=0.5,
                            color=c,
                            alpha=a,
                            ls='--',
                            label="$\\alpha_{"+key+"}"+f" = {self.alpha[key]:.2f}"+"$")

                    #print(self.intercept[f'key'])
                    #ax.plot(wl_range,
                    #        self.__powerlaw(wl_range,
                    #                        self.alpha[f'josefa'],
                    #                        self.intercept[f'josefa']),
                    #        lw=0.5,
                    #        color='red',
                    #        alpha=a,
                    #        ls='--',
                    #        label="$\\alpha_{"+"josefa"+"}"+f" = {self.alpha['josefa']}"+"$")
                    

            except KeyError as e:
                tqdm.write(f'{e} not found. Skipping this range.')
                #missing_handle = Line2D([0], [0], label='manual line', color=c, lw=0.5, alpha=a)
                ## access legend objects automatically created from data
                #handles, labels = plt.gca().get_legend_handles_labels()
                #handles.extend([missing_handle])


        wlRangeMask = (self.wlngths > lower) & (self.wlngths < upper)
        fullDataMask = ~np.isnan(self.fluxErrs)
        mask = wlRangeMask & fullDataMask

        if np.sum(np.isnan(self.fluxes)) == len(self.fluxes):
            tqdm.write(f"Empty data: skipping plot for source nr.: {int(self.srcID):4d}")
            return
        ax.scatter(self.wlngths[mask],
                   self.fluxes[mask],
                   marker=".",
                   c='r')
        ax.errorbar(self.wlngths,
                    self.fluxes,
                    yerr=self.fluxErrs*1,
                    fmt='.',
                    ecolor='k',
                    elinewidth=0.75,
                    barsabove=False,
                    capsize=2,
                    capthick=0.75,
                    ms=3.5)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(loc="upper right")

        c = self.cls[list(self.cls.keys())[0]]

        fig.suptitle(f"Source ID: {int(self.srcID):04d}"+"\n"+f"Class: {c} | Josefa: {self.cls['josefa']}")
        ax.set_xlabel('$\log_{10}\lambda [\mathrm{\mu m}]$')
        ax.set_ylabel('$\log_{10}\\left(\lambda F_\lambda \\left[\\frac{\mathrm{erg}}{\mathrm{s\,cm}^2}\\right]\\right)$')
        ax.set_xlim(10**(-0.5), 10**3)
        try:
            ax.set_ylim(10**(np.nanmean(np.log10(self.fluxes.value))-3),
                        10**(np.nanmean(np.log10(self.fluxes.value))+3))
        except:
            ax.set_ylim(10**(np.nanmin(np.log10(self.fluxes.value))),
                        10**(np.nanmax(np.log10(self.fluxes.value))))

        ax.axvspan(10**-0.5, lower, alpha=0.1, color='gray')
        ax.axvspan(upper, 10**3, alpha=0.1, color='gray')
        
        print(self.lambdaNames)
        
        t = np.array([10**-14, 3*10**-13, 10**-14, 3*10**-13, 10**-14, 3*10**-13, 10**-14, 3*10**-13, 10**-14, 3*10**-13, 10**-14, 10**-13])
        t = self.fluxes.value + (self.fluxes.value * 20)
        labels = [self.lambdaNames[i].replace('_lambda', '') for i in range(len(self.wlngths))]
        
        for i in range(len(self.wlngths)):
            ax.vlines(self.wlngths[i], t[i], self.fluxes[i].value, color='k', alpha=0.4, lw=0.5)
            ax.text(self.wlngths[i], t[i], labels[i], fontsize=5, rotation=0)

        plt.tight_layout()
        fmt = savepath.split('.')[-1]
        plt.savefig(savepath, format=fmt, dpi=300)
        plt.close(fig=fig)