import numpy as np
import scipy.constants as constants
from scipy.integrate import odeint
from scipy.optimize import root

class Material:
    def __init__(self, element, verbose=True):
        self.verbose = verbose
        if self.verbose:
            print 'initializing Material class using', element
        self.set_properties(element)
        self.sink_strength = 0.0
        self.sink_values = dict()  # holds individual contributions from different sink types

    def add_dislocation_line_sink(self, density, name='line'):
        line_sink_strength = 2.0*constants.pi*density/np.log(1.0/self.Rd/np.sqrt(constants.pi*density))  # sink strength for dislocation line in cm^-2
        self.sink_strength += line_sink_strength  # add to total sink strength
        self.sink_values[name] = line_sink_strength  # save contribution from this sink type
        self.sink_values['all'] = line_sink_strength + self.sink_values.get('all', 0.0)  # add to total

    def add_dislocation_loop_sink(self, density, radius, name='loop'):
        loop_sink_strength = 2.0*constants.pi*density*radius  # sink strength for dislocation loops in cm^-2
        self.sink_strength += loop_sink_strength  # add to total sink strength
        self.sink_values[name] = loop_sink_strength  # save contribution from this sink type
        self.sink_values['all_loops'] = loop_sink_strength + self.sink_values.get('all_loops', 0.0)  # track contribution from all loops
        self.sink_values['all'] = loop_sink_strength + self.sink_values.get('all', 0.0)  # add to total

    def set_temperature(self, T):
        k_boltz = constants.value('Boltzmann constant in eV/K')
        self.Di = self.Di0*np.exp(-self.Eim/k_boltz/T)          # interstitial diffusion coefficient in cm^2/s
        self.Dv = self.Dv0*np.exp(-self.Evm/k_boltz/T)          # vacancy diffusion coefficient in cm^2/s
        # self.K_recom = self.Ziv*self.Di/self.a/self.a/self.N  # recombination coefficient in cm^3/s
        self.K_recom = 1.0e-8 * self.Ziv * self.Di              # recombination coefficient in cm^3/s (for backwards compatibility)

    def set_flux(self, flux):
        self.flux = flux
        self.P = self.dpa_norm*self.N*self.casc_eff*flux  # defect production rate in defects/cm^3/s

    def ode_system(self, y, t):
        """RHS of ODE"""
        Ci, Cv = y
        return np.array([self.P - self.K_recom*Ci*Cv - self.sink_strength*self.Di*Ci,
                         self.P - self.K_recom*Ci*Cv - self.sink_strength*self.Dv*Cv])

    def transient(self, initial_conditions, ss_tol=1e-9):
        t = np.logspace(-8, 6, num=100)
        t = np.insert(t, 0, np.zeros(1))
        solution = odeint(self.ode_system, initial_conditions, t)
        Ci = solution[:,0]
        Cv = solution[:,1]

        # trim off extra solution times past steady state
        max_i = t.size
        for i in reversed(range(0, t.size)):
            Ci_slope, Cv_slope = self.ode_system([Ci[i], Cv[i]], t[i])
            if np.abs(Ci_slope/Ci[-1]) > ss_tol or np.abs(Cv_slope/Cv[-1]) > ss_tol:
                if i == t.size - 1:
                    print '*** Warning: simulation was possibly not run to steady state! ***'
                else:
                    max_i = i + 1
                    if self.verbose:
                        print 'steady state detected at:', t[i], 'seconds'
                break

        return [Ci[0:max_i], Cv[0:max_i]], t[0:max_i]

    def steady_state(self):
        # use no recombination case as initial guess
        initial_guess = [self.P/self.sink_strength/self.Di,
                         self.P/self.sink_strength/self.Dv]
        sol = root(self.ode_system, initial_guess, args=0)  # self.ode_system expects more arguments than passed
        return sol.x

    def get_sink_strength(self, name='all'):
        return self.sink_values[name]

    def set_properties(self, element):
        if element == 'Zr':
            self.Di0 = 3.5e-4        # interstitial diffusion pre-exponential factor in cm^2/s
            self.Dv0 = 2.2e-2        # vacancy diffusion pre-exponential factor in cm^2/s
            self.Eim = 0.06          # interstitial migration energy in eV
            self.Evm = 0.93          # vacancy migration energy in eV
            self.Evf = 1.79          # vacancy formation energy in eV
            self.Ziv = 500.0         # combinatorial factor for recombination (unitless)
            self.a = 3.234e-8        # lattice constant in cm
            self.Rd = 7.0e-8         # dislocation core radius in cm
            self.N = 4.29e22         # atomic density in atoms/cm^3
            self.dpa_norm = 1.6e-21  # specific dpa rate to fluence ratio
            self.casc_eff = 1.0/3.0  # cascade survival efficiency

if __name__ == '__main__':
    # test this script
    import matplotlib.pyplot as plt

    # create an instance
    m = Material('Zr')

    # add a sink
    m.add_dislocation_line_sink(1.0e8)

    # set test conditions
    m.set_temperature(700)
    m.set_flux(1.0e14)

    # check parameters
    print 'Di: {:.2e} (should be about 1.29E-04)'.format(m.Di)
    print 'Dv: {:.2e} (should be about 4.43E-09)'.format(m.Dv)
    print 'P: {:.2e} (should be about 2.29E+15)'.format(m.P)
    print 'K: {:.2e} (should be about 6.47E-10)'.format(m.K_recom)
    print 'Q: {:.2e} (should be about 9.39E+07)'.format(m.get_sink_strength())
    print 'Q from lines: {:.2e} (should also be about 9.39E+07)'.format(m.get_sink_strength('line'))

    Ci_steady, Cv_steady = m.steady_state()
    print 'steady state values:'
    print 'Ci: {:.2e}'.format(Ci_steady)
    print 'Cv: {:.2e}'.format(Cv_steady)

    print 'now compare those with the plot:'
    initial_conditions = np.zeros(2)
    sol, t = m.transient(initial_conditions)
    Ci, Cv = sol
    plt.loglog(t, Ci, label='Ci')
    plt.loglog(t, Cv, label='Cv')
    plt.legend(loc='best')
    plt.grid()
    plt.xlabel('time (seconds)')
    plt.ylabel('defect concentration (defects/cm^3)')

    plt.show()
