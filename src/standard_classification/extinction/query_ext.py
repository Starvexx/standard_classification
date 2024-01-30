from pnicer import ApparentMagnitudes
from pnicer.utils.auxiliary import get_resource_path
from astropy.coordinates import SkyCoord

science_path = get_resource_path(
    package="pnicer.tests_resources", resource="Orion_A_2mass.fits"
)
control_path = get_resource_path(
    package="pnicer.tests_resources", resource="CF_2mass.fits"
)

mag_names = ["Jmag", "Hmag", "Kmag"]
err_names = ["e_Jmag", "e_Hmag", "e_Kmag"]
extvec = [2.5, 1.55, 1.0]  # Indebetouw et al. 2015

science = ApparentMagnitudes.from_fits(
    path=science_path,
    extvec=extvec,
    mag_names=mag_names,
    err_names=err_names,
    lon_name="GLON",
    lat_name="GLAT",
    frame="galactic",
    coo_unit="deg",
)
# %%
control = ApparentMagnitudes.from_fits(
    path=control_path,
    extvec=extvec,
    mag_names=mag_names,
    err_names=err_names,
    lon_name="GLON",
    lat_name="GLAT",
    frame="galactic",
    coo_unit="deg",
)

science_color = science.mag2color()
control_color = control.mag2color()

ext_nicer = science.pnicer(control=control).discretize()
coo_test = SkyCoord(210.10845, -19.60590, unit="deg", frame="galactic")

ext = ext_nicer.query_position(skycoord=coo_test, bandwidth=5 / 60, nicest=True)

print(ext)