from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord as sc

from multiprocessing import Pool
from multiprocessing import cpu_count
from functools import partial
from tqdm import tqdm
from tqdm.contrib import tzip

import pnicer
from pnicer import ApparentMagnitudes
from pnicer.extinction import DiscreteExtinction

def process_coord(coord_i, ext_pnicer_discrete):
    coord, i = coord_i
    results = ext_pnicer_discrete.query_position(coord, bandwidth=4/60, use_fwhm=True, nicest=True, alpha=0.25)
    return i, results

yso_catalog = Table.read("/home/starvexx/Nemesis/catalogs/NEMESIS_Orion_Nov23.csv")
extinction_table = Table.read("/home/starvexx/Dropbox/Nemesis/standard_classification/discrete_extinction.tbl.fits")
extinction_table.info()
features = None
glon = extinction_table['GLON']
glat = extinction_table['GLAT']
ext = extinction_table['Extinction']
var = extinction_table['Variance']

#Get the features from the colours class. They contain colors, and coordinates

ext_pnicer_discrete = DiscreteExtinction(features=features,
                                         extinction=ext,
                                         variance=var)

ids, ext_K, ext_V, ext_err, nr_src, src_dens = [], [], [], [], [], []

coords = sc(ra = yso_catalog["RA"] * u.deg,
            dec = yso_catalog["DE"] * u.deg,
            frame='icrs')

partial_process_coord = partial(process_coord, ext_pnicer_discrete=ext_pnicer_discrete)
max_cores = cpu_count()

with Pool(processes=max_cores) as pool:
    results = list(tqdm(pool.map(partial_process_coord, zip(coords, yso_catalog['Internal_ID'])), total=len(coords)))
    
ids, ext_K, ext_V, ext_err, nr_src, src_dens = zip(*[(i, result[0], result[0] * 8.8, result[1], result[2], result[3]) for i, result in results])

#for coord, i in tzip(coords, yso_catalog['Internal_ID']):
#    results = ext_pnicer_discrete.query_position(coord, bandwidth=4/60, use_fwhm=True, nicest=True, alpha=0.25)
#    
#    ids.append(i)
#    ext_K.append(results[0])
#    ext_V.append(results[0] * 8.8)
#    ext_err.append(results[1])
#    nr_src.append(results[2])
#    src_dens.append(results[3])
    
ext_table = Table(data=[ids, ext_K * u.mag, ext_err * u.mag, ext_V * u.mag, nr_src, src_dens],
                  names=['Internal_ID', "A_Ks", "e_A_Ks", "A_V", "NR_of_sources", "source_density"])

ext_table.write("Orion_Extinction_per_source_new.fits", format='fits', overwrite=True)