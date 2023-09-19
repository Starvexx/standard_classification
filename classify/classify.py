import os

import argparse

import numpy as np

from astropy.table import Table
from astropy import units as u

from scipy import optimize
from scipy import odr

from matplotlib import pyplot as plt

from tqdm import tqdm

class star:
    """Star class."""

    def __init__(self, source):
        """Class constructor."""
        self.__data = source
        self.srcID = self.__data["Internal_ID"]
        self.alpha = {}
        self.intercept = {}
        self.cls = {}
        self.fluxNames = []
        self.fluxErrNames = []
        self.lambdaNames = []

        # Get the names of the different columns holding the relevant data.
        # To compute the alpha index we need all coumns containing flux,
        # flux errors, as well as their respective wavelengths lambda
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

   
    def fluxDens2flux(self):
        """Converts fluxdensities to fluxes"""
        self.fluxes = (self.fluxDens * u.Jy).to(u.erg / u.s / u.cm**2 / u.um,
                equivalencies=u.spectral_density(self.wlngths * u.um))

        return self.fluxes


    def fluxDensErr2fluxErr(self):
        self.fluxErrs = (self.fluxDensErrs * u.Jy).to(u.erg / u.s / u.cm**2 / u.um,
                equivalencies=u.spectral_density(self.wlngths * u.um))

        return self.fluxErrs


    def __line(self, x, k, d):
        """Just a boring old line, nothin else to see here."""
        return (k * x) + d


    def __powerlaw(self, param, x):
        """ Powerlaw representation of the line in log log space."""
        k, d = param[0], param[1]
        return 10 ** (k * np.log10(x) + d)


    def estimateAlpha(self, lower, upper):
        """Computes the alpha index in the chosen wavelength range in log-log space"""
        self.fluxDens2flux()
        wlRangeMask = (self.wlngths > lower) & (self.wlngths < upper)

        if np.sum(wlRangeMask) <= 1:
            self.alpha[f"{lower}-{upper}_est"] = np.nan
            self.intercept[f"{lower}-{upper}_est"] = np.nan
        else:
            self.log_wl = np.log10(self.wlngths[wlRangeMask])
            self.log_fluxes = np.log10(self.fluxes[wlRangeMask].value)

            opt, _ = optimize.curve_fit(self.__line, self.log_wl, self.log_fluxes)
            self.alpha[f"{lower}-{upper}_est"] = opt[0]
            self.intercept[f"{lower}-{upper}_est"] = opt[1]
            del opt

        return self.alpha, self.intercept


    def getAlphaWithErrors(self, lower, upper):
        """Computes the alpha index considering errors."""

        self.fluxDens2flux()
        self.fluxDensErr2fluxErr()
        self.estimateAlpha(lower, upper)

        wlRangeMask = (self.wlngths > lower) & (self.wlngths < upper)
        fullDataMask = ~np.isnan(self.fluxErrs)

        odrMask = wlRangeMask & fullDataMask

        self.wlMask = wlRangeMask

        #if str(self.srcID) == '11915':
        #    print(self.wlngths[wlRangeMask])
        #    print(self.fluxes[wlRangeMask])
        #    print(self.fluxErrs[wlRangeMask])
        #    print(np.sum(fullDataMask[wlRangeMask]))
        #    print(fullDataMask & wlRangeMask)
        #    results.pprint()
        #    exit(99)

        # Check if there are enough datapoints to compute the alpha index.
        # If there are less than two return NaN values.
        # If there are only two measurements but there is at least one
        # measurement error missing, use the estimated vale from scipy optimize.
        # If there are at least two measurements with errors, use only those
        # that provide measurement errors to compute the alpha index with ODR.

        if (np.sum(wlRangeMask) <= 1) or (np.sum(odrMask) <= 1):
            # Less than two measurements with errors. Using estimate.
            self.alpha[f"{lower}-{upper}"] = self.alpha[f"{lower}-{upper}_est"]
            self.intercept[f"{lower}-{upper}"] = self.intercept[f"{lower}-{upper}_est"]
        elif np.any(self.fluxErrs[wlRangeMask].value) and (np.sum(fullDataMask[wlRangeMask]) >= 2):
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


    def classify(self, alpha):
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

        # Check if the lower and upper boundaries are of the correct type
        # and convert to numpy array if not.
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


        ax.scatter(self.wlngths[self.wlMask], self.fluxes[self.wlMask], marker=".", c='r')
        ax.errorbar(self.wlngths, self.fluxes, yerr=self.fluxErrs*1, fmt='.', ecolor='k', elinewidth=0.75, barsabove=False, capsize=2, capthick=0.75, ms=3.5)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.legend(loc="upper right")


        fig.suptitle(f"Source ID: {int(self.srcID):04d}") #\nClass {args[0]['class_kiwm']}\t$\\alpha = {slopes[-1]:.2f}$")
        ax.set_xlabel('$\lambda [\mathrm{\mu m}]$')
        ax.set_ylabel('$\log\\left(\lambda F_\lambda \\left[\\frac{\mathrm{erg}}{\mathrm{s\,cm}^2}\\right]\\right)$')
        ax.set_xlim(10**(-0.5), 10**3)
        try:
            pass
            ax.set_ylim(10**(np.nanmean(np.log10(self.fluxes.value))-2), 10**(np.nanmean(np.log10(self.fluxes.value))+2))
        except:
            ax.set_ylim(10**(np.nanmin(np.log10(self.fluxes.value))), 10**(np.nanmax(np.log10(self.fluxes.value))))
            pass

        plt.tight_layout()
        fmt = savepath.split('.')[-1]
        plt.savefig(savepath, format=fmt, dpi=300)
        plt.close(fig=fig)
