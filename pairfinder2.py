from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import math
import time

# The index of each of the needed values in the file

ID_INDEX =          1 - 1
RA_INDEX =          2 - 1
DEC_INDEX =         3 - 1
TYPE_INDEX =        6 - 1
REDSHIFT_INDEX =    7 - 1
LEFT_RED_INDEX =    8 - 1
RIGHT_RED_INDEX =   9 - 1
MASS_INDEX =        24 - 1

GALAXY_TYPE =       '0'

# the file catalog with the galaxies

file_name = 'galaxies.dat'

# The program's representation of a galaxy

class Galaxy:
    def __init__(self, idnum, ra, dec, redshift, lowred, highred, mass):
        self.idnum = idnum
        self.right_ascension = ra
        self.declination = dec
        self.mass = mass
        self.redshift = redshift
        self.lowred = lowred # lower bound of redshift (of 68% confidence interval)
        self.highred = highred # upper bound of redshift (of 68% confidence interval)

# locates the index where the value 'val' occurs in the list. If 'val occurs
# multiple times, return the smallest index. Finds the index in log(n) time

def findlowest(g, val):
    return lowhelper(g, 0, len(g) - 1, val)

def lowhelper(g, low, high, val):
    mid = (high + low) // 2
    if low == high:
        return low
    elif g[mid].lowred >= val:
        return lowhelper(g, low, mid, val)
    else:
        return lowhelper(g, mid + 1, high, val)

# locates the index where the value 'val' occurs in the list. If 'val occurs
# multiple times, return the largest index. Finds the index in log(n) time

def findhighest(g, val):
    return highhelper(g, 0, len(g) - 1, val)

def highhelper(g, low, high, val):
    mid = (high + low) // 2 + 1
    if low == high:
        return low
    elif g[mid].lowred <= val:
        return highhelper(g, mid, high, val)
    else:
        return highhelper(g, low, mid - 1, val)

# finds the angular separation between two vectors in spherical coordinates
# works much faster than the astropy implementation
# t0, p0 = ra0, dec0
# t1, p1 = ra1, dec1

def sep(t0, p0, t1, p1):
	theta0 = math.radians(t0)
	theta1 = math.radians(t1)
	phi0 = math.radians(p0)
	phi1 = math.radians(p1)
	return math.degrees(math.acos( math.cos(phi0) * math.cos(phi1) * math.cos(theta0 - theta1) + \
                                       math.sin(phi0) * math.sin(phi1)))

# reads the galaxy catalog (file_name) and produces an array of galaxies
# possible pairs between galaxies are tested. Pairs are put into another array
# pairs are written to a file (named paircatalog.txt). If paircatalog.txt already
# exists, it is replaced with the new file.

def main():
    print('reading...')
    f = open(file_name)
    galaxies = []
    counter = 0
    print('<' + '---i---|' * 10 + '>\n<', end = '')
    print('starting timer')
    start_time = time.time()
    for line in f:
        if counter % 4242 == 0:
            print('*', end = '')
        a = line.split()
        # adds the object to the array if it is a galaxy and is above the mass cut
        if a[TYPE_INDEX] == GALAXY_TYPE and float(a[MASS_INDEX]) >= 9.522878745280336:
            galaxies.append(\
                Galaxy(int(  a[ID_INDEX]),\
                       float(a[RA_INDEX]),\
                       float(a[DEC_INDEX]),\
                       float(a[REDSHIFT_INDEX]),\
                       float(a[LEFT_RED_INDEX]),\
                       float(a[RIGHT_RED_INDEX]),\
                       float(a[MASS_INDEX])))
        counter = counter + 1
    print('>\nfound', len(galaxies), 'galaxies')

    print('sorting galaxies by redshift...')
    galaxies = sorted(galaxies, key=lambda g: g.lowred)
    print('testing pairs')
    total = 0
    galaxypairs = []

    # tests the pairs between galaxies
    
    for i in range(0, len(galaxies)):
        g = galaxies[i]
        if g.redshift < 0:
            continue
        # the approximate degree separation of 100 kpc at the redshift of the galaxy
        # intentionally overshoots the field so pairs will not be incorrectly eliminated
        narrowed_field = (1/3600) * 8.99201355134812 ** (g.lowred**-0.5)
        for j in range( max(i + 1, findlowest(galaxies, g.lowred)), findhighest(galaxies, g.highred)):
            g2 = galaxies[j]
            # if the second galaxy is not in the narrowed field, eliminate it
            if abs(g.right_ascension - g2.right_ascension) > narrowed_field or \
               abs(g.declination - g2.declination) > narrowed_field:
                continue
            # checks that:
            # the ratio of the galaxies is less than 1:3
            # one of the galaxies has a mass greater than 10^10 solar masses
            # the comoving distance between the two galaxies is less than 100kpc
            if abs(g.mass - g2.mass) <= 0.47712125471966244 and \
               max(g.mass, g2.mass) >= 10.0 and \
               sep(g.right_ascension, g.declination, g2.right_ascension, g2.declination) * 60.0 * \
               cosmo.kpc_comoving_per_arcmin((g.redshift + g2.redshift) / 2).value <= 100.0:
                galaxypairs.append((g, g2))
    print(time.time() - start_time, 'seconds')
    print(total)

    print('writing to file')

    # writes to the paircatalog file

    outcatalog = open('paircatalog.txt', 'w')
    
    for p in galaxypairs:
        outcatalog.write(str(p[0].idnum) + '\t' +\
                         str(p[0].right_ascension) + '\t' +\
                         str(p[0].declination) + '\t' +\
                         str(p[0].redshift) + '\t' +\
                         str(p[0].lowred) + '\t' +\
                         str(p[0].highred) + '\t' +\
                         str(p[0].mass) + '\t' +\
                         str(p[1].idnum) + '\t' +\
                         str(p[1].right_ascension) + '\t' +\
                         str(p[1].declination) + '\t' +\
                         str(p[1].redshift) + '\t' +\
                         str(p[1].lowred) + '\t' +\
                         str(p[1].highred) + '\t' +\
                         str(p[1].mass) + '\n')
        
main()
