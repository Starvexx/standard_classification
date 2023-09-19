#! /usr/bin/env python

import argparse
import os

from astropy.table import Table

from classify.classify import star

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
    args = parser.parse_args(
        ['-i', 'Orion_YSO_fluxes_24Jul23.csv', '-o', 'Orion_YSO_fluxes_24Jul23_classified.csv']
    )

    in_path = os.path.abspath(args.input_table)
    out_path = os.path.abspath(args.output_table)
    store_plots = args.store_plots
    plot_dir = os.path.abspath(args.plot_dir)

    if store_plots:
        if args.plot_dir == '':
            plot_dir = '/'.join(in_path.split('/')[0:-1]) + '/SED_plots' #{s.srcID:05d}.png'

        if ~os.path.exists(plot_dir):
            os.makedirs(plot_dir)

    data = Table.read(in_path)

    stars = []
    print('reading Data...')
    for source in data:
        stars.append(star(source))

    src_id = []
    alpha = []
    intercept = []
    alpha_est = []
    intercept_est = []
    cls = []
    cls_est = []

    print('Done reading!\n\nCrunching numbers ...')
    for s in tqdm(stars):
        if len(s.fluxDens) != 0:
            #s.fluxDens2flux()
            s.getAlphaWithErrors(2, 20)
            #print(s.getAlpha(2, 20))

            if store_plots:
                s.plot(plot_path)

            src_id.append(s.srcID)
            alpha.append(s.alpha['2-20'])
            alpha_est.append(s.alpha['2-20_est'])
            intercept.append(s.intercept['2-20'])
            intercept_est.append(s.intercept['2-20_est'])
            cls.append(s.cls['2-20'])
            cls_est.append(s.cls['2-20_est'])

    t = Table(
            data=[src_id,
                alpha_est,
                alpha,
                intercept_est,
                intercept,
                cls,
                cls_est],
            names=('InternalID',
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
