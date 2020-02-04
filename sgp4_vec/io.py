"""Read the TLE earth satellite file format.

This is a minimally-edited copy of "sgp4io.cpp".

"""
import re
from datetime import datetime
from math import pi, pow, log10
from sgp4_vec.ext import days2mdhms
from sgp4_vec.model import Satellite

INT_RE = re.compile(r'[+-]?\d*')
FLOAT_RE = re.compile(r'[+-]?\d*(\.\d*)?')

LINE1 = '1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN'
LINE2 = '2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN'

error_message = """TLE format error

The Two-Line Element (TLE) format was designed for punch cards, and so
is very strict about the position of every period, space, and digit.
Your line does not quite match.  Here is the official format for line {0}
with an N where each digit should go, followed by the line you provided:

{1}
{2}"""

"""
/*     ----------------------------------------------------------------
*
*                               sgp4io.cpp
*
*    this file contains a function to read two line element sets. while 
*    not formerly part of the sgp4 mathematical theory, it is 
*    required for practical implemenation.
*
*                            companion code for
*               fundamentals of astrodynamics and applications
*                                    2007
*                              by david vallado
*
*       (w) 719-573-2600, email dvallado@agi.com
*
*    current :
*              27 Aug 10  david vallado
*                           fix input format and delete unused variables in twoline2rv
*    changes :
*               3 sep 08  david vallado
*                           add operationmode for afspc (a) or improved (i)
*               9 may 07  david vallado
*                           fix year correction to 57
*              27 mar 07  david vallado
*                           misc fixes to manual inputs
*              14 aug 06  david vallado
*                           original baseline
*       ----------------------------------------------------------------      */
"""

"""
/* -----------------------------------------------------------------------------
*
*                           function twoline2rv
*
*  this function converts the two line element set character string data to
*    variables and initializes the sgp4 variables. several intermediate varaibles
*    and quantities are determined. note that the result is a structure so multiple
*    satellites can be processed simultaneously without having to reinitialize. the
*    verification mode is an important option that permits quick checks of any
*    changes to the underlying technical theory. this option works using a
*    modified tle file in which the start, stop, and delta time values are
*    included at the end of the second line of data. this only works with the
*    verification mode. the catalog mode simply propagates from -1440 to 1440 min
*    from epoch and is useful when performing entire catalog runs.
*
*  author        : david vallado                  719-573-2600    1 mar 2001
*
*  inputs        :
*    longstr1    - first line of the tle
*    longstr2    - second line of the tle
*    typerun     - type of run                    verification 'v', catalog 'c', 
*                                                 manual 'm'
*    typeinput   - type of manual input           mfe 'm', epoch 'e', dayofyr 'd'
*    afspc_mode  - True for afspc calculations, False for 'improved' mode
*    whichconst  - which set of constants to use  72, 84
*
*  outputs       :
*    satrec      - structure containing all the sgp4 satellite information
*
*  coupling      :
*    getgravconst-
*    days2mdhms  - conversion of days to month, day, hour, minute, second
*    jday        - convert day month year hour minute second into julian date
*    sgp4init    - initialize the sgp4 variables
*
*  references    :
*    norad spacetrack report #3
*    vallado, crawford, hujsak, kelso  2006
  --------------------------------------------------------------------------- */
"""

def twoline2rv(longstr1, longstr2, whichconst, afspc_mode=False):
    """Return a Satellite imported from two lines of TLE data.

    Provide the two TLE lines as strings `longstr1` and `longstr2`,
    and select which standard set of gravitational constants you want
    by providing `gravity_constants`:

    `sgp4_vec.earth_gravity.wgs72` - Standard WGS 72 model
    `sgp4_vec.earth_gravity.wgs84` - More recent WGS 84 model
    `sgp4_vec.earth_gravity.wgs72old` - Legacy support for old SGP4 behavior

    Normally, computations are made using various recent improvements
    to the algorithm.  If you want to turn some of these off and go
    back into "afspc" mode, then set `afspc_mode` to `True`.

    """

    deg2rad  =   pi / 180.0;         #    0.0174532925199433
    xpdotp   =  1440.0 / (2.0 *pi);  #  229.1831180523293

    line = longstr1.rstrip()

    if (len(line) >= 64 and
        line.startswith('1 ') and
        line[8] == ' ' and
        line[23] == '.' and
        line[32] == ' ' and
        line[34] == '.' and
        line[43] == ' ' and
        line[52] == ' ' and
        line[61] == ' ' and
        line[63] == ' '):

        _saved_satnum = int(line[2:7])
        classification = line[7] or 'U'
        intldesg = line[9:17]
        two_digit_year = int(line[18:20])
        epochdays = float(line[20:32])
        ndot = float(line[33:43])
        nddot = float(line[44] + '.' + line[45:50])
        nexp = int(line[50:52])
        bstar = float(line[53] + '.' + line[54:59])
        ibexp = int(line[59:61])
        # numb = int(line[62])
        elnum = int(line[64:68])
    else:
        raise ValueError(error_message.format(1, LINE1, line))

    line = longstr2.rstrip()

    if (len(line) >= 69 and
        line.startswith('2 ') and
        line[7] == ' ' and
        line[11] == '.' and
        line[16] == ' ' and
        line[20] == '.' and
        line[25] == ' ' and
        line[33] == ' ' and
        line[37] == '.' and
        line[42] == ' ' and
        line[46] == '.' and
        line[51] == ' '):

        satnum = int(line[2:7])
        if _saved_satnum != satnum:
            raise ValueError('Object numbers in lines 1 and 2 do not match')

        inclo = float(line[8:16])
        nodeo = float(line[17:25])
        ecco = float('0.' + line[26:33].replace(' ', '0'))
        argpo = float(line[34:42])
        mo = float(line[43:51])
        no_kozai = float(line[52:63])
        revnum = int(line[63:68])
    #except (AssertionError, IndexError, ValueError):
    else:
        raise ValueError(error_message.format(2, LINE2, line))

    #  ---- find no, ndot, nddot ----
    no_kozai = no_kozai / xpdotp; #   rad/min
    nddot= nddot * pow(10.0, nexp);
    bstar= bstar * pow(10.0, ibexp);

    #  ---- convert to sgp4 units ----
    ndot = ndot  / (xpdotp*1440.0);  #   ? * minperday
    nddot= nddot / (xpdotp*1440.0*1440);

    #  ---- find standard orbital elements ----
    inclo = inclo  * deg2rad;
    nodeo = nodeo  * deg2rad;
    argpo = argpo  * deg2rad;
    mo    = mo     * deg2rad;

    """
    // ----------------------------------------------------------------
    // find sgp4epoch time of element set
    // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
    // and minutes from the epoch (time)
    // ----------------------------------------------------------------

    // ---------------- temp fix for years from 1957-2056 -------------------
    // --------- correct fix will occur when year is 4-digit in tle ---------
    """
    if two_digit_year < 57:
        year = two_digit_year + 2000;
    else:
        year = two_digit_year + 1900;

    mon,day,hr,minute,sec = days2mdhms(year, epochdays);
    sec_whole, sec_fraction = divmod(sec, 1.0)
    epoch = datetime(year, mon, day, hr, minute, int(sec_whole),
                            int(sec_fraction * 1000000.0 // 1.0))

    satrec = Satellite(whichconst=whichconst, satnum=satnum, epoch=epoch, bstar=bstar,
                       inclo=inclo, nodeo=nodeo, ecco=ecco, argpo=argpo, mo=mo, no_kozai=no_kozai,
                       ndot=ndot, nddot=nddot, classification=classification, intldesg=intldesg,
                       elnum=elnum, revnum=revnum, afspc_mode=afspc_mode)
    return satrec


def verify_checksum(*lines):
    """Verify the checksum of one or more TLE lines.

    Raises `ValueError` if any of the lines fails its checksum, and
    includes the failing line in the error message.

    """
    for line in lines:
        checksum = line[68:69]
        if not checksum.isdigit():
            continue
        checksum = int(checksum)
        computed = compute_checksum(line)
        if checksum != computed:
            complaint = ('TLE line gives its checksum as {}'
                         ' but in fact tallies to {}:\n{}')
            raise ValueError(complaint.format(checksum, computed, line))


def fix_checksum(line):
    """Return a new copy of the TLE `line`, with the correct checksum appended.

    This discards any existing checksum at the end of the line, if a
    checksum is already present.

    """
    return line[:68].ljust(68) + str(compute_checksum(line))


def compute_checksum(line):
    """Compute the TLE checksum for the given line."""
    return sum((int(c) if c.isdigit() else c == '-') for c in line[0:68]) % 10


def rv2twoline(satrec):
    deg2rad = pi / 180.0
    xpdotp = 1440.0 / (2.0 * pi)
    line1_out_list = list(" " * 69)
    line1_out_list[0] = str(1)
    line1_out_list[2:7] = "{:05d}".format(satrec.satnum)
    line1_out_list[7] = satrec.classification
    line1_out_list[9:15] = satrec.intldesg
    line1_out_list[18:20] = str(satrec.epochyr)[-2:]
    line1_out_list[20:32] = "{:>012.8f}".format(satrec.epochdays)
    ndot_min = satrec.ndot * (xpdotp * 1440.0)
    if ndot_min < 0:
        line1_out_list[33:43] = "-."+"{:=8.8f}".format(ndot_min).split(".")[-1]
    else:
        line1_out_list[33:43] = " ." + "{:=8.8f}".format(ndot_min).split(".")[-1]
    if satrec.nddot == 0.0:
        line1_out_list[45:52] = "00000-0"
    else:
        nddot = satrec.nddot * (xpdotp * 1440.0 * 1440)
        nddot_module = abs(nddot)
        nexp = int(log10(nddot_module))
        ndot_fraction = pow(10.0, log10(nddot_module) - nexp)
        line1_out_list[44] = "-" if nddot < 0 else " "
        line1_out_list[45:50] = "{:1.5f}".format(ndot_fraction)[2:]
        line1_out_list[50:52] = str(nexp) if nexp < 0 else "+" + str(nexp)
    if satrec.bstar == 0.0:
        line1_out_list[54:61] = "00000+0"
    else:
        bstar_module = abs(satrec.bstar)
        ibexp = int(log10(bstar_module))
        if ibexp - log10(bstar_module) == 0:
            ibexp = ibexp + 1
        bstar_frac = pow(10.0, log10(bstar_module) - ibexp)
        line1_out_list[53] = "-" if satrec.bstar < 0.0 else " "
        line1_out_list[54:59] = "{:1.5f}".format(bstar_frac)[2:]
        line1_out_list[59:61] = "-"+str(ibexp) if ibexp == 0 else str(ibexp)
    line1_out_list[62] = "0"  # Its always 0 (originally this should have been "Ephemeris type")
    line1_out_list[64:68] = "{:>4d}".format(satrec.elnum)
    line1_out_list[68] = str(compute_checksum("".join(line1_out_list)))
    line1 = "".join(line1_out_list)

    line2_out_list = list(" " * 69)
    line2_out_list[0] = str(2)
    line2_out_list[2:7] = "{:05d}".format(satrec.satnum)
    if int(log10(satrec.inclo / deg2rad)) == 2:
        line2_out_list[8:16] = "{:>8.4f}".format(satrec.inclo / deg2rad)
    elif int(log10(satrec.inclo / deg2rad)) == 1:
        line2_out_list[9:16] = "{:>7.4f}".format(satrec.inclo / deg2rad)
    else:
        line2_out_list[10:16] = "{:>6.4f}".format(satrec.inclo / deg2rad)
    line2_out_list[17:25] = "{:>8.4f}".format(satrec.nodeo / deg2rad)
    line2_out_list[26:33] = "{:>8.7f}".format(satrec.ecco)[2:]
    line2_out_list[34:42] = "{:>8.4f}".format(satrec.argpo / deg2rad)
    line2_out_list[43:51] = "{:>8.4f}".format(satrec.mo / deg2rad)
    #func = lambda no: unkozai(no, satrec.ecco, satrec.inclo, whichconst) - satrec.no
    #sol = scipy.optimize.root_scalar(func, x0=satrec.no, x1=satrec.no - 1e-4, xtol=1e-13)
    #line2_out_list[52:63] = "{:3.8f}".format(sol.root * xpdotp)
    line2_out_list[52:63] = "{:>11.8f}".format(satrec.no_kozai * xpdotp)
    line2_out_list[63:68] = "{:>5d}".format(satrec.revnum)
    line2_out_list[68] = str(compute_checksum("".join(line2_out_list)))
    line2 = "".join(line2_out_list)

    return line1, line2
