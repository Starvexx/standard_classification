from __future__ import annotations

import numpy as np

import astropy
from astropy.table import Row
from astropy import units as u

from scipy import optimize
from scipy import odr

from matplotlib import pyplot as plt

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

        # Get the names of the different columns holding the relevant
        # data. To compute the alpha index we need all columns
        # containing flux, flux errors, as well as their respective
        # wavelengths, lambda.

        for name in self.__data.colnames[4::]:
            if (name[-4::] == 'flux'):
                self.fluxNames.append(name)
            elif name[-10::] == 'flux_error':
                self.fluxErrNames.append(name)
            elif name[-6::] == 'lambda':
                self.lambdaNames.append(name)
                
        

        # Extract the columns holding the flux measurements (Jy) from the table.
        self.fluxDens = np.array([col for col in self.__data[self.fluxNames].as_void()])
        self.fluxDensErrs = np.array([col for col in self.__data[self.fluxErrNames].as_void()])
        self.wlngths = np.array([col for col in self.__data[self.lambdaNames].as_void()])

        self.fluxDens[self.fluxDens == 0] = np.nan
        self.fluxDensErrs[self.fluxDensErrs == 0] = np.nan
        self.wlngths[self.wlngths == 0] = np.nan


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


    def __powerlaw(self, param : ArrayLike, x : float) -> float:
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
        k, d = param[0], param[1]
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
        self.fluxes = (self.fluxDens * u.Jy).to(u.erg / u.s / u.cm**2 / u.um,
                equivalencies=u.spectral_density(self.wlngths * u.um))

        return self.fluxes


    def fluxDensErr2fluxErr(self):
        """Converts flux density errors to flux errors.
        
        Does the conversion from Jansky (Jy) to erg / s / cm^2 / um
        
        Returns
        -------
            fluxErrs : float
                Returns the converted flux errors.
        """
        self.fluxErrs = (self.fluxDensErrs * u.Jy).to(u.erg / u.s / u.cm**2 / u.um,
                equivalencies=u.spectral_density(self.wlngths * u.um))

        return self.fluxErrs


    def estimateAlpha(self, lower : int, upper : int) -> tuple[float, float]:
        """Estimates the alpha index.
         
        This method computes the infrared alpha index in the chosen
        wavelength range by least squares fitting a line to the SED as
        a first estimate of the alpha index.

        Parameters
        ----------
            lower : int
                The lower boundary of the wavelength range for the 
                least squares fit.
            upper : int
                The upper boundary of the wavelength range for the
                least squares fit.

        Returns
        -------
            alpha : float
                The estimated slope of the least squares line fit
                which determines the alpha index.
            intercept : float
                The estimated intercept from the least squares line
                fit. can be used to over plot the line on a SED plot.
        
        """
        # Convert the flux density and flux density errors
        self.fluxDens2flux()
        self.fluxDensErr2fluxErr()

        # Get the mask to select data in the selected wavelength range.
        wlRangeMask = (self.wlngths > lower) & (self.wlngths < upper)

        if np.sum(wlRangeMask) <= 1:
            # If there is at most one viable flux measurement, the source
            # is not classifiable; store NaN values.
            self.alpha[f"{lower}-{upper}_est"] = np.nan
            self.intercept[f"{lower}-{upper}_est"] = np.nan
        elif (np.sum(wlRangeMask) == 2):
            # compute the relative difference between the two wavelengths.
            rel_wl_diff = self.wlngths[wlRangeMask][0] / self.wlngths[wlRangeMask][1]
            
            # If they are within 5% of the shorter wavelength the source is
            # not classifiable. Do the least squares fit otherwise.
            #if (1 - rel_wl_diff) < 0.05:
            if self.__wlClose(self.wlngths[wlRangeMask][0],
                              self.wlngths[wlRangeMask][1],
                              5):
                # If the measurements are too close to each other wavelength
                # wise, the source is not classifiable; store NaN values.
                self.alpha[f"{lower}-{upper}_est"] = np.nan
                self.intercept[f"{lower}-{upper}_est"] = np.nan
            else:
                # get the log wavelengths in the selected wavelength range.
                self.log_wl = np.log10(self.wlngths[wlRangeMask])
                # get the log fluxes in the selected wavelength range.
                self.log_fluxes = np.log10(self.fluxes[wlRangeMask].value)

                # Do the least squares line fit and store the ruslts.
                opt, _ = optimize.curve_fit(self.__line, self.log_wl, self.log_fluxes)
                self.alpha[f"{lower}-{upper}_est"] = opt[0]
                self.intercept[f"{lower}-{upper}_est"] = opt[1]
                # Free memory for garbage collection.
                del opt
        else:
            # get the log wavelengths in the selected wavelength range.
            self.log_wl = np.log10(self.wlngths[wlRangeMask])
            # get the log fluxes in the selected wavelength range.
            self.log_fluxes = np.log10(self.fluxes[wlRangeMask].value)

            # Do the least squares line fit and store the results.
            opt, _ = optimize.curve_fit(self.__line, self.log_wl, self.log_fluxes)
            self.alpha[f"{lower}-{upper}_est"] = opt[0]
            self.intercept[f"{lower}-{upper}_est"] = opt[1]
            # Free memory for garbage collection.
            del opt

        return self.alpha, self.intercept


    def getAlphaWithErrors(self, lower : float, upper : float) -> tuple:
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
        # Get the initial estimate for the ODR fit from a least squares fit.
        self.estimateAlpha(lower, upper)

        # Get the selection masks for the selected wavelength range.
        wlRangeMask = (self.wlngths > lower) & (self.wlngths < upper)
        # Get the selection mask for all measurements that have errors.
        fullDataMask = ~np.isnan(self.fluxErrs)
        # Get the selection mask for the data used in the ODR fit.
        odrMask = wlRangeMask & fullDataMask
        # Store the selection mask for the data in the selected wavelength range.
        self.wlMask = wlRangeMask

        # Check if there are enough data points to compute the alpha index.
        # If there are less than two return NaN values.
        # If there are only two measurements but there is at least one
        # measurement error missing, use the estimated vale from scipy optimize.
        # If there are at least two measurements with errors, use only those
        # that provide measurement errors to compute the alpha index with ODR.

        if (np.sum(odrMask) == 2):
            if self.__wlClose(self.wlngths[odrMask][0],
                              self.wlngths[odrMask][1],
                              5):
                self.alpha[f"{lower}-{upper}"] = self.alpha[f"{lower}-{upper}_est"]
                self.intercept[f"{lower}-{upper}"] = self.intercept[f"{lower}-{upper}_est"]
        if (np.sum(wlRangeMask) <= 1) or (np.sum(odrMask) <= 1):
            # Less than two measurements with errors. Using estimate.
            self.alpha[f"{lower}-{upper}"] = self.alpha[f"{lower}-{upper}_est"]
            self.intercept[f"{lower}-{upper}"] = self.intercept[f"{lower}-{upper}_est"]
        elif np.any(self.fluxErrs[wlRangeMask].value) and (np.sum(fullDataMask[wlRangeMask]) > 2):
            # Has at least two measurements with errors. Using ODR on data
            # with errors.
            mask = wlRangeMask & fullDataMask
            self.log_wl = np.log10(self.wlngths[mask])

            b_0 = [self.alpha[f"{lower}-{upper}_est"],
                   self.intercept[f"{lower}-{upper}_est"]]

            pl = odr.Model(self.__powerlaw)
            myData = odr.RealData(
                    self.wlngths[mask],
                    self.fluxes[mask],
                    sy=self.fluxErrs[mask])
            myODR = odr.ODR(myData, pl, beta0=b_0)
            results = myODR.run()

            self.alpha[f"{lower}-{upper}"] = results.beta[0]
            self.intercept[f"{lower}-{upper}"] = results.beta[1]
            del myODR
            del results
        else:
            # All measurements have errors. Using full ODR.
            self.log_wl = np.log10(self.wlngths[wlRangeMask])

            b_0 = [self.alpha[f"{lower}-{upper}_est"],
                   self.intercept[f"{lower}-{upper}_est"]]

            pl = odr.Model(self.__powerlaw)
            myData = odr.RealData(
                    self.wlngths[wlRangeMask],
                    self.fluxes[wlRangeMask],
                    sy=self.fluxErrs[wlRangeMask])
            myODR = odr.ODR(myData, pl, beta0=b_0)
            results = myODR.run()

            self.alpha[f"{lower}-{upper}"] = results.beta[0]
            self.intercept[f"{lower}-{upper}"] = results.beta[1]
            del myODR
            del results
        
        self.cls[f'{lower}-{upper}_est'] = self.classify(self.alpha[f'{lower}-{upper}_est'])
        self.cls[f'{lower}-{upper}'] = self.classify(self.alpha[f'{lower}-{upper}'])
        return self.alpha, self.intercept


    def altAlpha(self, lower : list, upper : list) -> tuple:
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
            return "III thin disc"
        elif (alpha < -2.5):
            return "III no disc / MS"
        else:
            return "not classified"
    

    def plot(self, savepath, lower=None, upper=None):
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

        # Check if the lower and upper boundaries are of the correct
        # type and convert to numpy array if not.
        try:
            assert isinstance(lower, (tuple, np.ndarray, list))
            assert isinstance(upper, (tuple, np.ndarray, list))
        except AssertionError:
            lower = np.array([lower])
            upper = np.array([upper])

        # If no values for the lower and upper boundary are provided, 
        if (lower == None) or (upper == None):
            key = list(self.alpha.keys())[0]
            if "_est" in key:
                lower = [int(key.split("_")[0].split('-')[0])]
                upper = [int(key.split("_")[0].split('-')[1])]
            else:
                lower = [int(key.split('-')[0])]
                upper = [int(key.split('-')[1])]
                    
        savepath = savepath + f'/{self.srcID:05d}.png'


        fig = plt.figure(figsize=(4,4), dpi=150)

        font = {'size'  : 7,
                'family': 'serif'}
        plt.rc('font', **font)
        plt.rc('text', usetex=True)

        ax = fig.add_subplot(111)

        for l, u in zip(lower, upper):
            try:
                if ~np.isnan(self.alpha[f'{l}-{u}']):
                    wl_range = np.linspace(0.5, 1000, 100)

                    ax.plot(wl_range,
                            self.__powerlaw(
                                [self.alpha[f'{l}-{u}'],
                                 self.intercept[f'{l}-{u}']],
                                wl_range),
                            lw=0.5,
                            color='k',
                            ls='--',
                            label="$\\alpha_{"+str(l)+"-"+str(u)+"}$")

                    ax.plot(wl_range,
                            self.__powerlaw(
                                [self.alpha[f'{l}-{u}_est'],
                                 self.intercept[f'{l}-{u}_est']],
                                wl_range),
                            lw=0.5,
                            color='k',
                            ls='-.',
                            label="$\\alpha_{"+str(l)+"-"+str(u)+"}$ est")

            except KeyError as e:
                tqdm.write(f'{e} not found. Skipping this range.')


        if np.sum(np.isnan(self.fluxes[self.wlMask])) == len(self.fluxes[self.wlMask]):
            tqdm.write("Empty data: skipping plot...")
            return
        ax.scatter(self.wlngths[self.wlMask],
                   self.fluxes[self.wlMask],
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


        fig.suptitle(f"Source ID: {int(self.srcID):04d}") #\nClass {args[0]['class_kiwm']}\t$\\alpha = {slopes[-1]:.2f}$")
        ax.set_xlabel('$\lambda [\mathrm{\mu m}]$')
        ax.set_ylabel('$\log\\left(\lambda F_\lambda \\left[\\frac{\mathrm{erg}}{\mathrm{s\,cm}^2}\\right]\\right)$')
        ax.set_xlim(10**(-0.5), 10**3)
        try:
            ax.set_ylim(10**(np.nanmean(np.log10(self.fluxes.value))-2),
                        10**(np.nanmean(np.log10(self.fluxes.value))+2))
        except:
            ax.set_ylim(10**(np.nanmin(np.log10(self.fluxes.value))),
                        10**(np.nanmax(np.log10(self.fluxes.value))))

        plt.tight_layout()
        fmt = savepath.split('.')[-1]
        plt.savefig(savepath, format=fmt, dpi=300)
        plt.close(fig=fig)
