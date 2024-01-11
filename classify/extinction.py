from pnicer import ApparentMagnitudes
#import pnicer

#from pnicer.tests import orion

#orion()

scifield = "/home/starvexx/Nemesis/standard_classification/test/NEMESIS_2MASS_VISION_SPITZER_WISE_no.cat.fits"
ctrlfield = "/home/starvexx/Nemesis/standard_classification/test/Ori-A_VISIONS_WISE_ctrl.cat.fits"
#scifield = "/home/starvexx/Nemesis/standard_classification/test/Ori-A_VISIONS_WISE_sci.cat.fits"
#ctrlfield = "/home/starvexx/Nemesis/standard_classification/test/Ori-A_VISIONS_WISE_ctrl.cat.fits"
#scifield = "/home/starvexx/Nemesis/standard_classification/test/Ori-A_VISIONS_sci.cat.fits"
#ctrlfield = "/home/starvexx/Nemesis/standard_classification/test/Ori-A_VISIONS_ctrl.cat.fits"

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

#science.mag2color().plot_combinations_kde()
#control.mag2color().plot_combinations_kde()

ext_pnicer = science_color.pnicer(control=control_color)

ext_pnicer_discrete = ext_pnicer.discretize()

#ext_pnicer_discrete.save_fits(path="./temp.fits")

pnicer_emap = ext_pnicer_discrete.build_map(bandwidth=2/60, metric="gaussian", sampling=2, use_fwhm=True)
print("finished extinction estimation.")
#exit(99)

pnicer_emap.save_fits(path="./ext_map_north.fits")

#pnicer_emap.plot_map(figsize=10)