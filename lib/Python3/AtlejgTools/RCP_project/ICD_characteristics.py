# -- coding: utf-8     # handle filenames with norwegian special characters
from pylab import *
import AtlejgTools.SimulationTools.UnitConversion as u
import AtlejgTools.Utils as ut
import AtlejgTools.RCP_project.Valve_characteristics as vc

def reynolds_number(rho, v, d, mu):
    mu *= u.CP
    v /= u.HOUR # m/hour -> m/s
    return rho*v*d/mu

def _icd_geom(h, b, w, unit=1.):
    h *= unit; b *= unit; w *= unit;
    A = h*(b + w) / 2.
    perimeter = b + w + 2*(sqrt(((w-b)/2)**2 + h**2))
    hydr_diam = 4*A / perimeter
    return A, hydr_diam

def icd_geometry(type):
    ''' see mail from sigurd 9/4-2010, 08:43 '''
    if type == 1.6 or type == 3.2:
        return _icd_geom(0.209, 0.180, 0.288, u.INCH)
    else:
        print("no such type")

class Struct : pass

def plot_icd_curves(file, col1, col2, line_from, line_to, sheetnm, rho, mu, type, sec_order_only, label):
    data = [ut.read_excel_column(sheetnm, col1, line_from, line_to, file=file),
            ut.read_excel_column(sheetnm, col2, line_from, line_to, file=file)]
    q = array(data[0])
    dp = array(data[1])
    # will fit dp = polynom(q) wrt polynom-coeff so must do some tricks
    if sec_order_only:
        A = vstack([q**2, zeros(len(q))]).T
        c = lstsq(A, dp)[0][0]
        dp_fit = c*q**2
    else:
        A = vstack([q**2, q]).T
        c1 = lstsq(A, dp)[0][0]; c2 = lstsq(A, dp)[0][1]
        dp_fit = c1*q**2 + c2*q
    #
    # laminar / turbulent transition
    A_channel, hydr_diam = icd_geometry(type);
    re = reynolds_number(rho, q/A_channel, hydr_diam, mu)
    q_lamin = 2300*mu*u.CP/hydr_diam/rho * A_channel * u.HOUR
    q_turb  = 4000*mu*u.CP/hydr_diam/rho * A_channel * u.HOUR
    print("fluid = %s q_lamin = %.1f q_turb = %1f" % (label, q_lamin, q_turb))
    #
    # eclipse vals
    a_icd = 0.0034600 * type/3.2
    dp_ecl = vc.icd_dp(rho, mu, a_icd, q*1000) # ref values are default eclipse values
    # compare
    print('fluid= %s: best fit: ratio = %.3e' % (label, dp_fit[0]/dp_ecl[0]))
    # plotting
    figure()
    title('ICD Baker %.1fbar fluid = %s' % (type, label))
    plot(q,dp, '*-', label='lab')
    plot(q, dp_ecl, '-*', label='eclipse')
    plot(q, dp_fit, '-*', label='best fit')
    plot([q_lamin, q_lamin], [0,max(dp)], '--k', label='laminar below')
    plot([q_turb, q_turb], [0,max(dp)], '--r', label='turbulent above')
    xlabel('flowrate [m3/h]')
    ylabel('dp [bar]')
    legend(loc='best')
    axis([0, max(q)*1.1, 0, max(dp_fit)*1.05])
    gca().grid()
    savefig('icd_%s.png' % label)
    #
    s = Struct()
    s.q = q; s.dp = dp; s.dp_fit = dp_fit; s.dp_ecl = dp_ecl; s.nm = label
    return s

if __name__ == '__main__':

    #if not sys.platform == 'win32':
        #print "gotta run on windows"
        #sys.exit(1)

    type = 3.2   # bar
    file = '/private/agy/Misc/Baker_3,2bar_ICD-tests_50_09-Final.xls'
    #file = u'G:\oile\Porsgrunn\part_1\otk\Projects-2009\Backup Kjøreplaner\Kjøreplaner bacup for Saudi Aramco ICD\Baker_3,2bar_ICD tests 50_09-Final.xls'
    #
    # gas viscosity found from Perry 3-249/250
    # gas density should be function of dp, but using an average value (from spreadsheet)
    oil12 = plot_icd_curves(file, 'l', 'k',    9, 16, 'Baker 3,2 bar', 878.,  12.,    type, False, 'oil12')
    gas   = plot_icd_curves(file, 'p', 'k',   22, 35, 'Baker 3,2 bar', 65.,    0.011, type, False, 'gas')
    wat   = plot_icd_curves(file, 'n', 'k',   46, 53, 'Baker 3,2 bar', 1005.,  0.5,   type, False, 'water')
    oil60 = plot_icd_curves(file, 'l', 'k', 100, 107, 'Baker 3,2 bar', 900.,  60.,    type, False, 'oil60')
    #stop
    figure(); hold('on')
    plot(oil12.q, oil12.dp, '*-k', label=oil12.nm + ' lab')
    plot(oil12.q, oil12.dp_ecl,'*--k',  label=oil12.nm + ' ecl')
    plot(oil60.q, oil60.dp, '*-r', label=oil60.nm + ' lab')
    plot(oil60.q, oil60.dp_ecl,'*--r',  label=oil60.nm + ' ecl')
    plot(gas.q, gas.dp, '*-b', label=gas.nm + ' lab')
    plot(gas.q, gas.dp_ecl, '*--b', label=gas.nm + ' ecl')
    plot(wat.q, wat.dp, '*-m', label=wat.nm + ' lab')
    plot(wat.q, wat.dp_ecl, '*--m', label=wat.nm + ' ecl')
    xlabel('flowrate [m3/h]')
    ylabel('dp [bar]')
    legend(loc='best')
    gca().grid()
    savefig('comp.png')
    show()
