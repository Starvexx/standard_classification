#! /usr/bin/env python

from astropy.table import Table

import numpy as np

from matplotlib import pyplot as plt


def main():
    file_path = "/home/starvexx/Nemesis/standard_classification/YSO_selection/yso_mags.fits"

    magnitudes = Table.read(file_path)

    Jmag = magnitudes["Jmag"]
    Hmag = magnitudes["Hmag"]
    Kmag = magnitudes["Ksmag"]
    e_Jmag = magnitudes["e_Jmag"]
    e_Hmag = magnitudes["e_Hmag"]
    e_Kmag = magnitudes["e_Ksmag"]
    
    color_JH = Jmag - Hmag
    color_HK = Hmag - Kmag

    e_color_JH = np.sqrt(e_Jmag**2 + e_Hmag**2)
    e_color_HK = np.sqrt(e_Hmag**2 + e_Kmag**2)

    plot(color_HK, color_JH, e_color_HK, e_color_JH)

    mask = remove(Hmag-Kmag, Jmag-Hmag)

    embedded = magnitudes[mask]['Internal_ID', 'RA', 'DE']

    embedded.write('embedded.fits', format='fits', overwrite=True)


def plot(x, y, dx, dy):
    fig = plt.figure()

    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ttl = np.array([(0.148, 0.549984), (2.888, 5.309279999999999), (3.782, 5.852831999999999), (1.054, 1.100832), (0.148, 0.549984)])
    exnir = np.array([(1.054, 1.100832), (3.782, 5.852831999999999), (5.728, 7.036), (3.0, 2.284), (1.054, 1.100832)])

    mask = remove(x-dx, y-dy)

    ax1.scatter(x, y, marker='.', color='k', alpha=0.3)
    ax1.plot(ttl[:,0], ttl[:,1])
    ax1.plot(exnir[:,0], exnir[:,1])
    ax1.vlines(0.148, -2, 0.549984)

    ax2.scatter(x-dx, y-dy, marker='.', color='k', alpha=0.3)
    ax2.scatter((x-dx)[mask], (y-dy)[mask], marker='.', color='r', alpha=0.3)
    ax2.plot(ttl[:,0], ttl[:,1])
    ax2.plot(exnir[:,0], exnir[:,1])
    xx = np.arange(0.148, 4, 0.1)
    ax2.plot(xx, line(1.7369693430656932, xx, 0.29291253722627747), color='r')
    ax2.vlines(0.148, -2, 0.549984)

    plt.show()


def line(k, x, d):
    return k * x + d


def remove(x, y):
    mask = (x > 0.148) & (y < line(1.7369693430656932, x, 0.29291253722627747))
    return mask


if __name__ == "__main__":
    main()

    exit(0)
