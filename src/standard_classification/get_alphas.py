#! /usr/bin/env python

import argparse
import os

import numpy as np

from astropy.table import Table

from cls.classify import star

from tqdm import tqdm

def main():
    """Main function"""

    parser = argparse.ArgumentParser(
            prog='classify',
            description='Classifies YSOs using the standard classification (α-index).',
            epilog='Enjoy :)')

    parser.add_argument(
            '-i', '--input-table',
            type=str)
    parser.add_argument(
            '-o', '--output-table',
            type=str)
    parser.add_argument(
            '-p', '--store-plots',
            action='store_true',
            default=False)
    parser.add_argument(
            '-d', '--plot-dir',
            type=str,
            default='')
    #args = parser.parse_args()
    args = parser.parse_args([
        '-i', '/home/starvexx/Nemesis/Data/NEMESIS_SEDs/Orion_YSO_fluxes_Jy_19Jun24count_update.csv',
        #'-i', '/home/starvexx/Nemesis/SOM_final/data/SEDs/Josefa_fluxes.csv',
        '-o', '/home/starvexx/Nemesis/SOM_final/data/SEDs/NEMESIS_YSOs_classified.csv',
        #'-o', '/home/starvexx/Nemesis/SOM_final/data/SEDs/Josefa_YSOs_classified.csv',
        #'-p',
        '-d', '/home/starvexx/Nemesis/SED_plots'])

    in_path = os.path.abspath(args.input_table)
    out_path = os.path.abspath(args.output_table)
    store_plots = args.store_plots
    plot_dir = os.path.abspath(args.plot_dir)

    if store_plots:
        if args.plot_dir == '':
            plot_dir = os.getcwd() + '/SED_plots'
        else:
            plot_dir = os.path.abspath(plot_dir)

        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

    data = Table.read(in_path)

    src_id = []
    ra = []
    dec = []

    alpha_2_24 = []
    intercept_2_24 = []
    cls_2_24 = []

    alpha_km = []
    intercept_km = []
    cls_km = []

    alpha_irac = []
    intercept_irac = []
    cls_irac = []

    
    wl_names = [key for key in data.colnames if 'lambda' in key]
    km = np.zeros(len(wl_names), dtype=bool)
    for i, wl in enumerate(wl_names):
        if 'irac' in wl or 'mips' in wl or 'Ks' in wl:
            km[i] = True

    irac = np.zeros(len(wl_names), dtype=bool)
    for i, wl in enumerate(wl_names):
        if 'irac' in wl:
            irac[i] = True


    print('\nDone reading!\n\nCrunching numbers ...\n')
    for source in tqdm(data):
        s = star(source)
        if len(s.fluxDens) != 0:
            s.getAlphaWithErrors(2, 24)
            s.getAlphaWithErrors(Mask=km, key='KM')
            s.getAlphaWithErrors(Mask=irac, key='irac')
            #s.getAlphaWithErrors(2, 10)
            #s.getAlphaWithErrors(10, 24)
            #s.alpha['josefa'] = source['alphaKM']
            #s.intercept['josefa'] = s.intercept['2-24']
            #s.cls['josefa'] = source['Class'].split("'")[1]

            if store_plots:
                s.plot(plot_dir)

            tqdm.write('Source ID: {}'.format(s.srcID))
            
            src_id.append(s.srcID)
            ra.append(s.ra)
            dec.append(s.dec)

            alpha_2_24.append(s.alpha['2-24'])
            intercept_2_24.append(s.intercept['2-24'])
            cls_2_24.append(s.cls['2-24'])

            alpha_km.append(s.alpha['KM'])
            intercept_km.append(s.intercept['KM'])
            cls_km.append(s.cls['KM'])

            alpha_irac.append(s.alpha['irac'])
            intercept_irac.append(s.intercept['irac'])
            cls_irac.append(s.cls['irac'])


    t = Table(
            data=[src_id,
                  ra,
                  dec,
                  alpha_2_24,
                  intercept_2_24,
                  cls_2_24,
                  alpha_km,
                  intercept_km,
                  cls_km,
                  alpha_irac,
                  intercept_irac,
                  cls_irac],
            names=('Internal_ID',
                   'RA',
                   'DE',
                   'alpha_2_24',
                   'intercept_2_24',
                   'class_2_24',
                   'alpha_km',
                   'intercept_km',
                   'class_km',
                   'alpha_irac',
                   'intercept_irac',
                   'class_irac'))#,
    t.write(out_path, format='csv', overwrite=True)
    
    del s


if __name__ == "__main__":
    main()

    exit(0)
