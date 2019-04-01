import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import astropy.coordinates as coord
import astropy.units as u
import ephem
import pandas
import psycopg2
from astropy.table import Table
import healpy as hp

def teffcut(band):
    if band == 'u' or band == 'g':
        return 0.2
    else:
        return 0.3

def read_file(data_file):
    df = pandas.read_csv(data_file,index_col=[0])
    return df

def make_exp_map(data_frame=None,
                 data_file=None,
                 title=None,
                 color=plt.cm.Purples,
                 band=None,
                 out_file='test.pdf',
                 show_ra_labels_in_hours=False,
                 width_of_mw_band=10.,
                 color_of_mw_band='k',
                 minteff=0.,
                 minexptime=0.,
                 mintilings=1):

    if data_file is None:
        data = Table.from_pandas(data_frame)  
    else:
        data = Table.from_pandas(pandas.read_csv(data_file,index_col=[0]))

    if band is None:
        idx = (data['teff']>=minteff) & (data['exptime']>=minexptime)  
        label = 'total exposure [log (exptime/sec)]'
    else:
        idx = (data['teff']>=minteff) & (data['exptime']>=minexptime) & (data['band']==band) 
        label = 'total '+band+'-band exposure [log (exptime/sec)]'
        
    ra = coord.Angle(data[idx]['radeg']*u.degree)
    dec = coord.Angle(data[idx]['decdeg']*u.degree)
    exptime = data[idx]['exptime']
    
    ra = ra.wrap_at(180*u.degree)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    npix_ra=360/6 
    hb=ax.hexbin(ra.radian, dec.radian, C=exptime, gridsize=npix_ra, bins='log',
                 reduce_C_function=np.sum,cmap=color)
    cb=fig.colorbar(hb,orientation='horizontal',cax=fig.add_axes([0.2, 0.17, 0.6, 0.03]))
    cb.set_label(label)
    if title is not None: ax.set_title(title+'\n')
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')
    if show_ra_labels_in_hours:
        ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    ax.grid(True)

    ra_MW, dec_MW = the_galaxy_in_equatorial_coords(width_of_mw_band)
    for i in range(len(ra_MW)):
        ra = coord.Angle(ra_MW[i]*u.degree)
        ra = ra.wrap_at(180*u.degree)        
        dec = coord.Angle(dec_MW[i]*u.degree)
        slices = unlink_wrap(ra.radian)
        for slc in slices:
            ax.plot(ra.radian[slc],dec.radian[slc],color=color_of_mw_band,lw=2)        

#    ra_fields_hours={'F1': 10, 'F2': 16, 'F3': 22}
#    ra_fields_degrees={'F1': 150, 'F2': -120, 'F3': -30}
#    dec_fields_degrees={'F1': -15, 'F2': -15, 'F3': -15}
#    ra_fields = np.array([150.,-120.,-30.])
#    dec_fields = np.array([-15.,-15.,-15.])
#    ra = coord.Angle(ra_fields*u.degree)
#    ra = ra.wrap_at(180*u.degree)
#    dec = coord.Angle(dec_fields*u.degree)
#    ax.scatter(ra.radian,dec.radian,s=100,color='red',marker='D')

    fig.savefig(out_file)
#    fig.show()

    return coverage_fraction(data[idx]['radeg'],data[idx]['decdeg'],data[idx]['exptime'],mintilings*minexptime)


def the_galaxy_in_equatorial_coords(width=0.):
    lon_array = np.arange(0,360)
    ra, dec = [] , []
    lat_list = [ -0.5*width, 0. , 0.5*width]
    for lat in lat_list:
        eq_array = np.zeros((360,2))
        for lon in lon_array:
            ga = ephem.Galactic(np.radians(lon), np.radians(lat))
            eq = ephem.Equatorial(ga)
            eq_array[lon] = np.degrees(eq.get())
        ra.append(eq_array[:,0])
        dec.append(eq_array[:,1])
    return ra,dec

def unlink_wrap(dat, lims=[-np.pi, np.pi], thresh = 0.95):
    """
    Iterate over contiguous regions of `dat` (i.e. where it does not
    jump from near one limit to the other).

    This function returns an iterator object that yields slice
    objects, which index the contiguous portions of `dat`.

    This function implicitly assumes that all points in `dat` fall
    within `lims`.

    """    
    jump = np.nonzero(np.abs(np.diff(dat)) > ((lims[1] - lims[0]) * thresh))[0]
    lasti = 0
    for ind in jump:
        yield slice(lasti, ind + 1)
        lasti = ind + 1
    yield slice(lasti, len(dat))

def bands():
    filters=['u','g','r','i','z','Y','VR']

    colormaps={'u': plt.cm.Blues,
               'g': plt.cm.Greens,
               'r': plt.cm.Oranges,
               'i': plt.cm.Reds,
               'z': plt.cm.Greys,
               'Y': plt.cm.Purples,
               'VR': plt.cm.Purples}

    coverage={'u': 0.,
              'g': 0.,
              'r': 0.,
              'i': 0.,
              'z': 0.,
              'Y': 0.,
              'VR': 0.}

    return filters, colormaps, coverage


def get_exp_info(outfile=None,minexptime=29.999,minteff=0.,test=False,limit=10,verbose=False):

    columns = """id as EXPNUM,
       TO_CHAR(date - '12 hours'::INTERVAL, 'YYYYMMDD') AS NITE,
       EXTRACT(EPOCH FROM date - '1858-11-17T00:00:00Z')/(24*60*60) AS MJD_OBS,
       ra AS RADEG,
       declination AS DECDEG,
       filter AS BAND,
       exptime AS EXPTIME,
       propid AS PROPID,
       qc_teff AS TEFF,
       object as OBJECT"""

    table = "exposure.exposure"

    criteria = "exptime>"+str(minexptime)+" and qc_teff>"+str(minteff)+" and flavor='object' "

    query = "SELECT "+columns+"\nFROM "+table+"\nWHERE "+criteria 

    if test : query = query + "\nLIMIT "+str(limit)

    if verbose : print query

    conn =  psycopg2.connect(database='decam_prd',
                             user='decam_reader',
                             host='des61.fnal.gov',
                             password='reader',
                             port=5443)
    exposures = pandas.read_sql(query, conn)
    conn.close()

    if outfile is not None: exposures.to_csv(outfile)

    return exposures



def coverage_fraction(radeg,decdeg,exptime,minexptime):

    nside = 32
    phi = np.radians(radeg)
    theta = np.radians(decdeg+90)
    npix = hp.nside2npix(nside)
    sky = np.zeros(npix)

    for idx in range(exptime.size):
        pix = hp.ang2pix(nside, theta[idx], phi[idx])
        sky[pix] += exptime[idx]


    visible_fraction = 0.65

    pix_fraction = 1. #3./hp.nside2pixarea(32, degrees=True)

    npix_covered = np.where(sky>=minexptime)[0].size

    return pix_fraction * npix_covered / (visible_fraction * npix)




