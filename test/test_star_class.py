import sys

sys.path.append('/home/starvexx/Nemesis/standard_classification/src/standard_classification')

from astropy.table import Table
import numpy as np

from cls.classify import star


def load_test_data():
    #data = Table.read('./test/data/src_1453.csv')
    fluxes = [4, 3, 2, 1]
    flux_err = [0.1, 0.2, 0.3, 0.4]
    wl = [2, 4, 8, 20]
    
    data_list = [0, 1, 0, 0]
    
    for i in range(4):
        data_list.append(fluxes[i])
        data_list.append(flux_err[i])
        data_list.append(wl[i])

    data = Table(
        names=("col1", "Internal_ID", "RA", "DE",
               "f2_flux", "f2_flux_error", "f2_lambda",
               "f4_flux", "f4_flux_error", "f4_lambda",
               "f8_flux", "f8_flux_error", "f8_lambda",
               "f20_flux", "f20_flux_error", "f20_lambda"
               )
    )
    data.add_row(data_list)

    return data[0]


def test_convert_flux():
    results = [2.99792458e-09,
               5.6211085875e-10,
               9.3685143125e-11,
               7.494811450000002e-12]
    
    s = star(load_test_data())

    fluxes = s.fluxDens2flux().value
    
    for flux, result in zip(fluxes, results):
        assert flux == result
    

def test_convert_flux_errors():
    results = [7.49481145e-11,
               3.747405725e-11,
               1.405277146875e-11,
               2.9979245800000003e-12]
    
    s = star(load_test_data())
    
    flux_errs = s.fluxDensErr2fluxErr().value
    
    for err, result in zip(flux_errs, results):
        assert err == result

        
def test_alpha_estimate():
    s = star(load_test_data())
    s.fluxDens2flux()
    s.fluxDensErr2fluxErr()
    
    alpha_est, intercept_est = s.estimateAlpha(1, 24)
    
    assert alpha_est["1-24_est"] == -2.605854342722085
    assert intercept_est["1-24_est"] == -7.707497323179037
    

def test_alpha_with_errors():
    s = star(load_test_data())
    s.fluxDens2flux()
    s.fluxDensErr2fluxErr()
    
    alpha, intercept = s.getAlphaWithErrors(1, 24)
    
    assert alpha["1-24"] == -2.6710142298228625
    assert intercept["1-24"] == -7.636196229482067
    
    
def test_classify():
    s = star(load_test_data())
    
    assert s.classify(1) == "0/I"
    assert s.classify(0) == "flat"
    assert s.classify(-1) == "II"
    assert s.classify(-2) == "III thin disk"
    assert s.classify(-3) == "III no disk / MS"
    assert s.classify(np.nan) == "not classified"