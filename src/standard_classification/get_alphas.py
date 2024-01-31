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
        '-i', '../data/raw_data/Orion_YSO_fluxes_24Jul23.csv',
        '-o', '../data/processed_data/Orion_YSO_fluxes_24Jul23_classified.csv',
        '-d', '../SED_plots'])

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

    stars = []
    print('\nreading Data...\n')
    #for source in tqdm(data[11085:11086]):
    for source in tqdm(data):
        stars.append(star(source))

    src_id = []
    ra = []
    dec = []
    alpha = []
    intercept = []
    alpha_est = []
    intercept_est = []
    cls = []
    cls_est = []

    s = stars[1]
    s.fluxDens2flux()
    
    #print(s.fluxes)
    #print(s.wlngths)
    #s.altAlpha([1,2], [90, 200])
    #exit(9)

    print('\nDone reading!\n\nCrunching numbers ...\n')
    for s in tqdm(stars):
        if len(s.fluxDens) != 0:
            #s.fluxDens2flux()
            s.getAlphaWithErrors(2, 20)
            s.getAlphaWithErrors(2, 10)
            s.getAlphaWithErrors(10, 25)
            #print(s.getAlpha(2, 20))

            if store_plots:
                s.plot(plot_dir, 2, 20)

            src_id.append(s.srcID)
            ra.append(s.ra)
            dec.append(s.dec)
            alpha.append(s.alpha['2-20'])
            alpha_est.append(s.alpha['2-20_est'])
            intercept.append(s.intercept['2-20'])
            intercept_est.append(s.intercept['2-20_est'])
            cls.append(s.cls['2-20'])
            cls_est.append(s.cls['2-20_est'])

    results = [src_id,
               ra]

    t = Table(
            data=[src_id,
                  ra,
                  dec,
                  alpha_est,
                  alpha,
                  intercept_est,
                  intercept,
                  cls,
                  cls_est],
            names=('Internal_ID',
                   'RA',
                   'DE',
                   'alpha_est',
                   'alpha',
                   'intercept_est',
                   'intercept',
                   'class',
                   'class_est'))
    t.write(out_path, format='csv', overwrite=True)


if __name__ == "__main__":
    main()

    exit(0)
