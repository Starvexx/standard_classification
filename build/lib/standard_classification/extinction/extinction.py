from pnicer import ApparentMagnitudes

from tqdm import tqdm
from tqdm.contrib import tzip

from astropy.coordinates import SkyCoord as sc
from astropy import units as u
from astropy.table import Table

from multiprocessing import Pool
from functools import partial

def process_coord(coord_i, ext_pnicer_discrete):
    coord, i = coord_i
    results = ext_pnicer_discrete.query_position(coord, bandwidth=4/60, use_fwhm=True, nicest=True, alpha=0.25)
    return i, results

yso_catalog = Table.read("/home/starvexx/Nemesis/catalogs/NEMESIS_Orion_Nov23.csv")

ids, ext_K, ext_V, ext_err, nr_src, src_dens = [], [], [], [], [], []

coords = sc(ra = yso_catalog["RA"] * u.deg,
            dec = yso_catalog["DE"] * u.deg,
            frame='icrs')

scifield = "/home/starvexx/Nemesis/catalogs/NEMESIS_FOV_clean.cat.fits"
ctrlfield = "/home/starvexx/Nemesis/catalogs/NEMESIS_CTRL_clean.cat.fits"

mag_names = ["Jmag", "Hmag", "Ksmag", "W1mag", "W2mag"]
err_names = ["e_Jmag", "e_Hmag", "e_Ksmag", "e_W1mag", "e_W2mag"]
extvec = [2.5, 1.55, 1.0, 0.97, 0.55]

science = ApparentMagnitudes.from_fits(path=scifield, extvec=extvec,
                                       mag_names=mag_names, err_names=err_names,
                                       lon_name="_Glon", lat_name="_Glat",
                                       frame="galactic", coo_unit="deg")

control = ApparentMagnitudes.from_fits(path=ctrlfield, extvec=extvec,
                                       mag_names=mag_names, err_names=err_names, 
                                       lon_name="_Glon", lat_name="_Glat",
                                       frame="galactic", coo_unit="deg")

science_color = science.mag2color()
control_color = control.mag2color()

ext_pnicer = science_color.pnicer(control=control_color)

ext_pnicer_discrete = ext_pnicer.discretize()

ext_pnicer_discrete.save_fits(path = './discrete_extinction_update_2023-12-21.tbl.fits',
                              overwrite = True)

pnicer_emap = ext_pnicer_discrete.build_map(bandwidth=4/60, metric="gaussian", sampling=2, use_fwhm=True, nicest=True, alpha=0.25)
print("finished extinction estimation.")

#partial_process_coord = partial(process_coord, ext_pnicer_discrete=ext_pnicer_discrete)
#
#with Pool(processes=22) as pool:
#    results = list(tqdm(pool.map(partial_process_coord, zip(coords, yso_catalog['Internal_ID'])), total=len(coords)))
#    
#ids, ext_K, ext_V, ext_err, nr_src, src_dens = zip(*[(i, result[0], result[0] * 8.8, result[1], result[2], result[3]) for i, result in results])
#
##for coord, i in tzip(coords, yso_catalog['Internal_ID']):
##    results = ext_pnicer_discrete.query_position(coord, bandwidth=4/60, use_fwhm=True, nicest=True, alpha=0.25)
##    
##    ids.append(i)
##    ext_K.append(results[0])
##    ext_V.append(results[0] * 8.8)
##    ext_err.append(results[1])
##    nr_src.append(results[2])
##    src_dens.append(results[3])
#
#ext_table = Table(data=[ids, ext_K * u.mag, ext_err * u.mag, ext_V * u.mag, nr_src, src_dens],
#                  names=['Internal_ID', "A_Ks", "e_A_Ks", "A_V", "NR_of_sources", "source_density"])
#
#ext_table.write("Orion_Extinction_per_source_new.fits", format='fits', overwrite=True)

pnicer_emap.save_fits(path="./ext_map_complete_update_2023-12-21.fits")