The "signal-something.dat" files hold information regarding fusion events, necessary to run firstMacro.The naming convention goes:
    signal-<projectile name><projectile isotope><target name><target isotope>to<daughter isotope><daughter name>.dat

The information within them is stored in the following way:
    energy[keV], charge[% of Beam Energy], peak Height[fraction of beam bragg peak], ~distance travelled [mm]
    
    

The 'similarRuns.txt' contains a list of all the runs for the fusion experiment, grouped by similar parameters (e.g.: same pressure, same voltage, etc.). Refer to elog for details. It is strictly for humanoïd usage.


8Li-crossSec.txt is a misnomer, as it contains list of cross sections for several fusion events, not all involving 8Li. Also strictly for humanoïd use.
