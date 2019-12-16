"""The Satellite class."""

from sgp4.ext import jday
from sgp4.propagation import sgp4, sgp4init
import datetime

minutes_per_day = 1440.


class Satellite(object):
    """An earth-orbiting satellite as represented by the SGP4 model.

    Most of this class's hundred-plus attributes are intermediate values
    of interest only to the propagation algorithm itself.  Here are the
    attributes set by ``sgp4.io.twoline2rv()`` in which users are likely
    to be interested:

    ``satnum``
        Unique satellite number given in the TLE file.
    ``epochyr``
        Full four-digit year of this element set's epoch moment.
    ``epochdays``
        Fractional days into the year of the epoch moment.
    ``jdsatepoch``
        Julian date of the epoch (computed from ``epochyr`` and ``epochdays``).
    ``ndot``
        First time derivative of the mean motion (ignored by SGP4).
    ``nddot``
        Second time derivative of the mean motion (ignored by SGP4).
    ``bstar``
        Ballistic drag coefficient B* in inverse earth radii.
    ``inclo``
        Inclination in radians.
    ``nodeo``
        Right ascension of ascending node in radians.
    ``ecco``
        Eccentricity.
    ``argpo``
        Argument of perigee in radians.
    ``mo``
        Mean anomaly in radians.
    ``no``
        Mean motion in radians per minute.
    """

    def __init__(self, whichconst, satnum, epoch, bstar, inclo, nodeo, ecco, argpo, mo, no_kozai,
                 ndot=0, nddot=0, classification='U', intldesg=None, elnum=1, revnum=1, afspc_mode=False):
        """
        :param whichconst: standard set of gravitational constants.
            `sgp4.earth_gravity.wgs72` - Standard WGS 72 model
            `sgp4.earth_gravity.wgs84` - More recent WGS 84 model
            `sgp4.earth_gravity.wgs72old` - Legacy support for old SGP4 behavior
        :param satnum: Unique satellite number given in the TLE file.
        :param epoch: naive datetime objects with TLE epoch
        :param bstar: Ballistic drag coefficient B* in inverse earth radii.
        :param inclo: Inclination in radians.
        :param nodeo: Right ascension of ascending node in radians.
        :param ecco: Eccentricity.
        :param argpo: Argument of perigee in radians.
        :param mo: Mean anomaly in radians.
        :param no_kozai: Mean motion in radians per minute (rad/min).
        :param ndot: First time derivative of the mean motion (ignored by SGP4).
        :param nddot: Second time derivative of the mean motion (ignored by SGP4).
        :param classification:  (U=Unclassified, C=Classified, S=Secret)
        :param intldesg: International designator (last two digits of launch year,
        launch number of the year, piece of the launch)
        :param elnum: Element set number. Incremented when a new TLE is generated for this object.
        :param revnum: Revolution number at epoch (revolutions)
        :param afspc_mode: Normally, computations are made using various recent improvements to the
         algorithm.  If you want to turn some of these off and go back into "afspc" mode, then set
          `afspc_mode` to `True`.
        """
        self.whichconst = whichconst
        self.satnum = satnum
        self.epoch = epoch
        self.epochyr = epoch.year
        self.jdsatepoch = jday(epoch.year, epoch.month, epoch.day, epoch.hour, epoch.minute,
                               epoch.second+1e-6*epoch.microsecond)
        days_in_the_year = (epoch.date() - datetime.date(epoch.year, 1, 1)).days + 1
        self.epochdays = days_in_the_year + (epoch.hour + epoch.minute / 60 +
                                             (epoch.second + 1e-6 * epoch.microsecond) / 60 / 60) / 24.0
        self.bstar = bstar
        self.inclo = inclo
        self.nodeo = nodeo
        self.ecco = ecco
        self.argpo = argpo
        self.mo = mo
        self.no_kozai = no_kozai
        self.a = pow(self.no_kozai * whichconst.tumin, (-2.0 / 3.0))
        self.alta = self.a * (1.0 + self.ecco) - 1.0
        self.altp = self.a * (1.0 - self.ecco) - 1.0
        self.ndot = ndot
        self.nddot = nddot
        self.classification = classification
        self.intldesg = intldesg
        self.elnum = elnum
        self.revnum = revnum
        self.error = 0

        sgp4init(whichconst, afspc_mode, self.satnum, self.jdsatepoch - 2433281.5, self.bstar,
                 self.ecco, self.argpo, self.inclo, self.mo, self.no_kozai,
                 self.nodeo, self)

    def propagate(self, year, month=1, day=1, hour=0, minute=0, second=0.0):
        """Return a position and velocity vector for a given date and time."""

        j = jday(year, month, day, hour, minute, second)
        m = (j - self.jdsatepoch) * minutes_per_day
        r, v = sgp4(self, m)
        return r, v
