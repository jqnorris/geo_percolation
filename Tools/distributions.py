import numpy as np
from scipy.special import erf
from scipy.special import gamma
from scipy.special import gammainc

# Functions Commonly used in distributions
def heaviside(x):
    return np.floor(np.sign(x))

def Theta(x):
    return 0.5 * (1 + erf(x/np.sqrt(2.)))


def format_params(params, errors = None):
    error_precision = (np.floor(np.log10(abs(errors)))).astype(int)
    rounded_params = [round(par, -error_precision[i]) for i, par in enumerate(params)]
    error_size = [('{0:e}'.format(err))[0] for err in errors]
    return ['{0}'.format(par) +'({0})'.format(error_size[i]) for i, par in enumerate(rounded_params)]

class normal:
    class pdf:
        @staticmethod
        def calc(x, m, s):
            return 1./(s * np.sqrt(2. * np.pi)) * np.exp(-(x-m)**2./(2. * s**2.))
        @staticmethod
        def latex(m = None, s = None, short = False, error = None):
            latex_string = r'\begin{gather*}f(x, \mu, \sigma)'
            if short == False:
                latex_string += (r' =\displaystyle\frac{1}{\sigma \sqrt{2 \pi}}'
                                + r' e^{-\displaystyle\frac{(x - \mu)^2}{2 \sigma^2}}')
            else:
                latex_string += r'\sim e^{-\displaystyle\frac{(x - \mu)^2}{2 \sigma^2}}'      
            if m != None and s != None:
                if len(error) == 2:
                    latex_string += r'\\ \mu={0},\; \sigma={1}'.format(*format_params([m, s], error))
                else:
                    latex_string += r'\\ \mu={0},\; \sigma={1}'.format(m, s)
            latex_string += r'\end{gather*}'
            return latex_string
    class cdf:
        @staticmethod
        def calc(x, m, s):
            return 0.5*(1. + erf((x-m)/(s * np.sqrt(2))))
        @staticmethod
        def latex(m = None, s = None, short = False, error = None):
            latex_string = (r'\begin{gather*}F(x, \mu, \sigma)'
                            + r'=\displaystyle\frac{1}{2} \left[1'
                            + r'+ \mathrm{erf}\left(\frac{x - \mu}{\sigma \sqrt{2}}\right)\right]')
            if m != None and s != None:            
                if len(error) == 2:
                    latex_string += r'\\ \mu={0},\; \sigma={1}'.format(*format_params([m, s], error))
                else:
                    latex_string += r'\\ \mu={0},\; \sigma={1}'.format(m, s)
            latex_string += r'\end{gather*}'
            return latex_string
  
class exponential:
    class pdf:
        @staticmethod
        def calc(x, l):
            return l * np.exp(-l * x) * heaviside(x)
        @staticmethod
        def latex(l = None, short = False, error = None):
            latex_string = r'\begin{gather*}f(x, \lambda)'
            if short == False:
                latex_string += r'= \lambda e^{-\lambda x} H(x)'
            else:
                latex_string += r'\sim e^{-\lambda x}'
            if l != None:
                if len(error) == 1:
                    latex_string += r'\\ \lambda={0}'.format(*format_params([l,], error))
                else:
                    latex_string += r'\\ \lambda={0}'.format(l)                
            latex_string += '\end{gather*}'
            return latex_string
    class cdf:
        @staticmethod
        def calc(x, l):
            return (1. - np.exp(-l * x)) * heaviside(x)
        @staticmethod
        def latex(l = None, short = False, error = None):
            latex_string = r'\begin{gather*}F(x, \lambda)'
            if short == False:
                latex_string += r' = (1 - e^{-\lambda x}) H(x)'
            else:
                latex_string += ' = 1 - e^{-\lambda x}'
            if l != None:
                if len(error) == 1:
                    latex_string += r'\\ \lambda={0}'.format(*format_params([l,], error))
                else:
                    latex_string += r'\\ \lambda={0}'.format(l) 
            latex_string += r'\end{gather*}'
            return latex_string

class lognormal:
    class pdf:
        @staticmethod
        def calc(x, m, s):
            return (1./(x * s * np.sqrt(2. * np.pi))* np.exp( - (np.log(x) - m)**2/(2. * s**2))) * heaviside(x)
        @staticmethod
        def latex(m = None, s = None, short = False, error = None):
            latex_string = r'\begin{gather*}f(x, \mu, \sigma)'
            if short == False:
                latex_string += (r' = \displaystyle\frac{1}{x \sigma \sqrt{2 \pi}}'
                                + r'e^{-\displaystyle\frac{(\ln{x} - \mu)^2}{2 \sigma^2}}H(x)')
            else:
                latex_string += r'\sim e^{-\displaystyle\frac{(\ln{x} - \mu)^2}{2 \sigma^2}}'
            if m != None and s != None:
                if len(error) == 2:
                    latex_string += r'\\ \mu={0},\; \sigma={1}'.format(*format_params([m, s], error))
                else:
                    latex_string += r'\\ \mu={0},\; \sigma={1}'.format(m, s)
            latex_string += r'\end{gather*}'
            return latex_string
    class cdf:
        @staticmethod
        def calc(x, m, s):
            return 0.5 * (1 + erf((np.log(x) - m)/(s * np.sqrt(2))))
        @staticmethod
        def latex(m = None, s = None, short = False, error = None):
            latex_string = (r'\begin{gather*}f(x, \mu, \sigma)'
                            + r'= \displaystyle\frac{1}{2}'
                            + r'\left[1 + \mathrm{erf}\left('
                            + r'\displaystyle\frac{\ln{x} - \mu}{\sigma \sqrt{2}}\right)\right]')
            if m != None and s != None:
                if len(error) == 2:
                    latex_string += r'\\ \mu={0},\; \sigma={1}'.format(*format_params([m, s], error))
                else:
                    latex_string += r'\\ \mu={0},\; \sigma={1}'.format(m, s)
            latex_string += r'\end{gather*}'
            return latex_string

class weibull:
    class pdf:
        @staticmethod
        def calc(x, l, k):
            return (k/float(l) * (x/float(l))**(k-1) * np.exp(-(x/float(l))**k)) * heaviside(x)
        @staticmethod
        def latex(l = None, k = None, short = False, error = None):
            latex_string = r'\begin{gather*}f(x, \lambda, k)'
            if short == False:
                latex_string += (r'= \displaystyle\frac{k}{\lambda}'
                                + r'\left(\displaystyle\frac{x}{\lambda}\right)^{k-1}'
                                + r'e^{-\left(x / \lambda\right)^k} H(x)')
            else:
                latex_string += (r'\sim \left(\displaystyle\frac{x}{\lambda}\right)^{k-1}'
                                + r'e^{-\left(x / \lambda\right)^k}')
            if l != None and k != None:
                if len(error) == 2:                    
                    latex_string += r'\\ \lambda={0},\; k={1},\; \alpha={2}'.format(*format_params([l, k], error))
                else:
                    latex_string += r'\\ \lambda={0},\; k={1},\; \alpha={2}'.format(l, k)
            latex_string += r'\end{gather*}'
            return latex_string

    class cdf:
        @staticmethod
        def calc(x, l, k):
            return (1. - np.exp(-((x/float(l))**k))) * heaviside(x)
        @staticmethod
        def latex(l = None, k = None, short = False, error = None):
            latex_string = r'\begin{gather*}F(x, \lambda, k)'
            if short != False:
                latex_string += r' = \left[1 - e^{-\left(x / \lambda\right)^k}\right] H(x)'
            else:
                latex_string += r' = 1 - e^{-\left(x / \lambda\right)^k}'
            if l != None and k != None:
                if len(error) == 2:                    
                    latex_string += r'\\ \lambda={0},\; k={1},\; \alpha={2}'.format(*format_params([l, k], error))
                else:
                    latex_string += r'\\ \lambda={0},\; k={1},\; \alpha={2}'.format(l, k)
            latex_string += r'\end{gather*}'
            return latex_string

class exponentiated_weibull:
    class pdf:
        @staticmethod
        def calc(x, l, k, a):
            return (a * k/float(l) * (x/float(l))**(k-1) * (1 - np.exp(-(x/float(l))**k))**(a-1) * np.exp(- (x/float(l))**k)) * heaviside(x)
        @staticmethod
        def latex(l = None, k = None, a = None, short = False, error = None):
            latex_string = r'\begin{gather*}f(x, \lambda, k)'
            if short == False:
                latex_string += (r' = \alpha \displaystyle\frac{k}{\lambda}'
                                + r'\left[\displaystyle\frac{x}{\lambda}\right]^{k-1}'
                                + r'\left[1 - e^{-\left(x / \lambda\right)^k}\right]^{\alpha -1}'
                                + r'e^{-\left(x / \lambda\right)^k} H(x)')
            else:
                latex_string += (r'\sim \left[\displaystyle\frac{x}{\lambda}\right]^{k-1}'
                                + r'\left[1 - e^{-\left(x / \lambda\right)^k}\right]^{\alpha -1}'
                                + r'e^{-\left(x / \lambda\right)^k}')
            if l != None and k != None and a != None:
                if len(error) == 3:                    
                    latex_string += r'\\ \lambda={0},\; k={1},\; \alpha={2}'.format(*format_params([l, k, a], error))
                else:
                    latex_string += r'\\ \lambda={0},\; k={1},\; \alpha={2}'.format(l, k, a)
            latex_string += r'\end{gather*}'
            return latex_string
    class cdf:
        @staticmethod
        def calc(x, l, k, a):
            return (1. - np.exp(-((x/float(l))**k)))**a * heaviside(x)
        @staticmethod
        def latex(l = None, k = None, a = None, short = False, error = None):
            latex_string = r'\begin{gather*}F(x, \lambda, k)'
            if short == False:
                latex_string += r' = \left[1 - e^{- \displaystyle \left(x / \lambda\right)^k}\right]^{\alpha} H(x)'
            else:
                latex_string += r' = \left[1 - e^{- \displaystyle \left(x / \lambda\right)^k}\right]^{\alpha}'
            if l != None and k != None and a != None:
                if len(error) == 3:                    
                    latex_string += r'\\ \lambda={0},\; k={1},\; \alpha={2}'.format(*format_params([l, k, a], error))
                else:
                    latex_string += r'\\ \lambda={0},\; k={1},\; \alpha={2}'.format(l, k, a)
            latex_string += r'\end{gather*}'
            return latex_string

class gamma_dist:
    class pdf:
        @staticmethod
        def calc(x, a, b):
            return ((b**a * x**(a-1.) * np.exp(-x * b))/gamma(a)) * heaviside(x)
        @staticmethod
        def latex(a = None, b = None, short = False, error = None):
            latex_string = r'\begin{gather*}f(x, \alpha, \beta)'
            if short == False:
                latex_string += r'= \displaystyle\frac{\beta^\alpha x^{\alpha-1} e^{-x \beta}}{\Gamma(\alpha)} H(x)'
            else:
                latex_string += r'\sim x^{\alpha-1} e^{-x \beta}'
            if a != None and b != None:
                if len(error) == 2:
                    latex_string += r'\\ \alpha={0},\; \beta={1}'.format(*format_params([a, b], error))
                else:
                    latex_string += r'\\ \alpha={0},\; \beta={1}'.format(a, b)
            latex_string += r'\end{gather*}'
            return latex_string
    class cdf:
        @staticmethod
        def calc(x, a, b):
            return gammainc(a, b*x)/gamma(a) * heaviside(x)
        @staticmethod
        def latex(a = None, b = None, short = False, error = None):
            latex_string = r'\begin{gather*}F(x, \alpha, \beta)'
            if short == False:
                latex_string += r'= \displaystyle\frac{\gamma(\alpha, \beta x)}{\Gamma(\alpha)}'
            else:
                latex_string += r'= \displaystyle\frac{\gamma(\alpha, \beta x)}{\Gamma(\alpha)}'
            if a != None and b != None:
                if len(error) == 2:
                    latex_string += r'\\ \alpha={0},\; \beta={1}'.format(*format_params([a, b], error))
                else:
                    latex_string += r'\\ \alpha={0},\; \beta={1}'.format(a, b)
            latex_string += r'\end{gather*}'
            return latex_string

class inverse_gaussian:
    class pdf:
        @staticmethod
        def calc(x, m, l):
            return np.sqrt(l / (2. * np.pi * x**3.)) * np.exp((-l * (x - m)**2.)/(2. * m**2 * x))
        @staticmethod
        def latex(m = None, l = None, short = False, error = None):
            latex_string = r'\begin{gather*}f(x, \mu, \lambda)'
            if short == False:
                latex_string += (r' = \left[\displaystyle\frac{\lambda}{2 \pi \sigma x^3}\right]^{1 / 2}'
                                + r' e^{\displaystyle\frac{-\lambda (x - \mu)^2 }{2 \mu^2 x}}')
            else:
                latex_string += (r'\sim \displaystyle\frac{1}{\sqrt{x^3}}'
                                + r' e^{\displaystyle\frac{-\lambda (x - \mu)^2 }{2 \mu^2 x}}')
            if m != None and s != None:
                if len(error) == 2:
                    latex_string += r'\\ \mu={0},\; \lambda={1}'.format(*format_params([m, l], error))
                else:
                    latex_string += r'\\ \mu={0},\; \lambda={1}'.format(m, s)
            latex_string += r'\end{gather*}'
            return latex_string
    class cdf:
        @staticmethod
        def calc(x, m, l):
            return Theta(np.sqrt(float(l)/x) * (x/float(m) - 1)) + np.exp(2. * l /float(m)) * Theta(-np.sqrt(float(l)/x) * (x/float(m) + 1.))
        @staticmethod
        def latex(m = None, s = None, short = False, error = None):
            latex_string = (r'\begin{gather*}\begin{align*}F(x, \mu, \lambda)'
                            + r'=&\Phi\left(\sqrt{\frac{\lambda}{x}}'
                            + r'\left(\frac{x}{\mu}-1 \right)\right)'
                            + r'\\ &+\exp\left(\frac{2 \lambda}{\mu}\right)'
                            + r'\Phi\left(-\sqrt{\frac{\lambda}{x}}'
                            + r'\left(\frac{x}{\mu}+1 \right)\right)\end{align*}')
            if m != None and s != None:            
                if len(error) == 2:
                    latex_string += r'\\ \mu={0},\; \lambda={1}'.format(*format_params([m, s], error))
                else:
                    latex_string += r'\\ \mu={0},\; \lambda={1}'.format(m, s)
            latex_string += r'\end{gather*}'
            return latex_string
                                                            
    

# Maintain a list of names that maps to classes
names = {'Normal':normal, 'Exponential':exponential, 'Lognormal':lognormal,
         'Weibull':weibull, 'Exponentiated Weibull':exponentiated_weibull,
         'Gamma':gamma_dist, 'Inverse Gaussian':inverse_gaussian}