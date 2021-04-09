import numpy as np
from py_wake.site._site import UniformWeibullSite
from py_wake.wind_turbines import OneTypeWindTurbines, WindTurbines

# Turbine positions
fnm='DOW_WTGpos.csv'
wt_x,wt_y=np.loadtxt(fnm, delimiter=';', unpack=True)

# Turbine names (optional)
turbine_names=[ "A01", \
                "A02", \
                "A03", \
                "A04", \
                "A05", \
                "B01", \
                "B02", \
                "B03", \
                "B04", \
                "B05", \
                "B06", \
                "C01", \
                "C02", \
                "C03", \
                "C04", \
                "C05", \
                "C06", \
                "D01", \
                "D02", \
                "D03", \
                "D04", \
                "D05", \
                "D06", \
                "E01", \
                "E02", \
                "E03", \
                "E04", \
                "E05", \
                "E06", \
                "F01", \
                "F02", \
                "F03", \
                "F04", \
                "F05", \
                "F06", \
                "G01", \
                "G02", \
                "G03", \
                "G04", \
                "G05", \
                "G06", \
                "H01", \
                "H02", \
                "H03", \
                "H04", \
                "H05", \
                "H06", \
                "J01", \
                "J02", \
                "J03", \
                "J04", \
                "J05", \
                "K01", \
                "K02", \
                "K03", \
                "K04", \
                "K05", \
                "L01", \
                "L02", \
                "L03", \
                "L04", \
                "L05", \
                "T01", \
                "T02", \
                "T03", \
                "T04", \
                "T05" ]

# Turbine data
D = 154 # Turbine diameter
hh = 102 # Turbine hub height
power_curve = np.array([[ 2., 0.00000000e+00],
                        [ 3., 1.04732169e+01],
                        [ 4., 2.81111181e+02],
                        [ 5., 5.80419291e+02],
                        [ 6., 1.04732169e+03],
                        [ 7., 1.67876083e+03],
                        [ 8., 2.52927047e+03],
                        [ 9., 3.55966160e+03],
                        [10., 4.70039799e+03],
                        [11., 5.60327371e+03],
                        [12., 5.93400147e+03],
                        [13., 5.99053175e+03],
                        [14., 6.01066816e+03],
                        [15., 6.00669792e+03],
                        [16., 5.98146180e+03],
                        [17., 5.99746825e+03],
                        [18., 5.98821577e+03],
                        [19., 6.02544247e+03],
                        [20., 6.02381102e+03],
                        [21., 6.02240775e+03],
                        [22., 5.95252926e+03],
                        [23., 5.96883234e+03],
                        [24., 5.99314435e+03],
                        [25., 6.06088940e+03],
                        [30., 0.00000000e+00]])

ct_curve = np.array([[ 2., 0.05 ],
                     [ 3., 0.768],
                     [ 4., 0.768],
                     [ 5., 0.762],
                     [ 6., 0.764],
                     [ 7., 0.768],
                     [ 8., 0.758],
                     [ 9., 0.763],
                     [10., 0.688],
                     [11., 0.608],
                     [12., 0.434],
                     [13., 0.328],
                     [14., 0.259],
                     [15., 0.209],
                     [16., 0.171],
                     [17., 0.143],
                     [18., 0.122],
                     [19., 0.104],
                     [20., 0.089],
                     [21., 0.078],
                     [22., 0.07 ],
                     [23., 0.062],
                     [24., 0.055],
                     [25., 0.05 ],
                     [30., 0.05 ]])


class Equinor6MW(OneTypeWindTurbines):
    def __init__(self):
        OneTypeWindTurbines.__init__(self, name='Equinor6MW', diameter=D, hub_height=hh,
                                     ct_func=self._ct, power_func=self._power, power_unit='kW')

    def _ct(self, u):
        return np.interp(u, ct_curve[:, 0], ct_curve[:, 1])

    def _ct(self, u):
        return np.interp(u, ct_curve[:, 0], ct_curve[:, 1])


    def _power(self, u):
        return np.interp(u, power_curve[:, 0], power_curve[:, 1])

class TwoWTGs(WindTurbines):
#   wt = Dudgeon2.TwoWTGs(n1, len(wt_x)-n1, 1e-3, Dudgeon2.D, 1., Dudgeon2.hh)
    def __init__(self, n1, n2, d1, d2, h1, h2, scaler=0.5):
        #names = ['WTG1']*n1 + ['WTG2']*n2
        #diameters = [d1]*n1 + [d2]*n2
        #hub_heights = [h1]*n1 + [h2]*n2
        #ct_funcs = [self._ct]*(n1+n2)
        #power_funcs = [self._power]*(n1+n2)
        #
        self.scaler = scaler
        #
        names       = ['WTG-old', 'WTG-new']
        diameters   = [d1, d2]
        hub_heights = [h1, h2]
        ct_funcs    = [self._ct_old, self._ct_new]
        power_funcs = [self._power_old, self._power_new]
        #
        WindTurbines.__init__(self, names=names, diameters=diameters, hub_heights=hub_heights,
                                     ct_funcs=ct_funcs, power_funcs=power_funcs, power_unit='kW')
#
    def _ct_new(self, u):
        return np.interp(u, ct_curve[:, 0], ct_curve[:, 1])
#
    def _ct_old(self, u):
        return self.scaler * self._ct_new(u)
#
    def _power_new(self, u):
        return np.interp(u, power_curve[:, 0], power_curve[:, 1])
#
    def _power_old(self, u):
        return  self.scaler * self._power_new(u)



# Site data (12 sector Weibull distribution in this case)
class DudgeonSite(UniformWeibullSite):
    def __init__(self):
        f = [0.06358372, 0.05154005, 0.05598169, 0.05601585, 0.06062832,
             0.07878778, 0.09757931, 0.12542494, 0.15489349, 0.11049422,
             0.07620821, 0.06886243]
        a = [ 9.3142286 ,  8.73501422,  9.13492833,  9.4372641 ,  9.21722703,
             10.09130741, 10.9808964 , 12.38269603, 12.52082084, 11.7641739 ,
             10.29622427, 10.37081206]
        k = [2.12890285, 2.22624468, 2.2883152 , 2.20179583, 2.24467622,
             2.43758064, 2.36221578, 2.50003135, 2.54711259, 2.3353032 ,
             2.25678044, 2.09631299]
        ti = 0.06
        UniformWeibullSite.__init__(self, f, a, k, ti)
        self.initial_position = np.array([wt_x, wt_y]).T
