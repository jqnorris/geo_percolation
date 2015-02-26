import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.optimize import curve_fit
from distributions import *


def plot_with_fits(raw_data):
    # Separate data into x and y arrays
    if len(raw_data) == 2:
        x = raw_data[0]
        y = raw_data[1]
        data = (np.array(raw_data)).T
    else:
        try:
            x = data[:, 0]
            y = data[:, 1]
            data = np.array(raw_data)
        except:
            print "Data is not a 2D array, cannot plot or fit."
    
    # Plot data
    plt.close()
    pdf_fig = plt.figure('pdf', figsize=(6,4))
    pdf_ax = pdf_fig.add_axes([0.1, 0.1, 0.6, 0.8])
    pdf_ax.plot(x, y, label='Data', color='black', ls='--', lw=2)
    
    cdf_fig = plt.figure('cdf', figsize=(6,4))
    cdf_ax = cdf_fig.add_axes([0.1, 0.1, 0.6, 0.8])
    cdf_ax.plot(x, np.cumsum(y), label='Data', color='black', ls='--', lw=2)


    # Normal fit and plot
    if True:
        try:
            fit_params, pcov = curve_fit(normal.cdf.calc, x, np.cumsum(y))
            pdf_ax.plot(x, normal.pdf.calc(x, *fit_params), label='Normal')
            cdf_ax.plot(x, normal.cdf.calc(x, *fit_params), label='Normal')
            KS = np.max(np.abs(normal.cdf.calc(x, *fit_params) - np.cumsum(y)))
            print KS
        except:
            pass

    # Exponential fit and plot
    if True:
        try:
            fit_params, pcov = curve_fit(exponential.cdf.calc, x, np.cumsum(y))
            pdf_ax.plot(x, exponential.pdf.calc(x, *fit_params), label='Exponetial')
            cdf_ax.plot(x, exponential.cdf.calc(x, *fit_params), label='Exponential')
            KS = np.max(np.abs(exponential.cdf(x, *fit_params) - np.cumsum(y)))
            print KS
        except:
            pass

    # Lognormal fit and plot
    if True:
        try:
            fit_params, pcov = curve_fit(lognormal.cdf.calc, x, np.cumsum(y))
            pdf_ax.plot(x, lognormal.pdf.calc(x, *fit_params), label='Lognormal')
            cdf_ax.plot(x, lognormal.cdf.calc(x, *fit_params), label='Lognormal')
            KS = np.max(np.abs(lognormal.cdf.calc(x, *fit_params) - np.cumsum(y)))
            print KS
        except:
            pass

    weibull_fit = [1,]
    # Weibull fit and plot
    if True:
        try:
            fit_params, pcov = curve_fit(weibull.cdf.calc, x, np.cumsum(y))            
            pdf_ax.plot(x, weibull.pdf.calc(x, *fit_params), label='Weibull')
            cdf_ax.plot(x, weibull.cdf.calc(x, *fit_params), label='Weibull')
            KS = np.max(np.abs(weibull.cdf.calc(x, *fit_params) - np.cumsum(y)))
            print KS
        except:
            pass
    
    # Exponentiated Weibull fit and plot
    if True:
        try:            
            fit_params, pcov = curve_fit(exponentiated_weibull.cdf.calc, x, np.cumsum(y))
            weibull_fit = fit_params
            pdf_ax.plot(x, exponentiated_weibull.pdf.calc(x, *fit_params), label='Exponentiated Weibull')
            cdf_ax.plot(x, exponentiated_weibull.cdf.calc(x, *fit_params), label='Exponentiated Weibull')
            KS = np.max(np.abs(exponentiated_weibull.cdf.calc(x, *fit_params) - np.cumsum(y)))
            print KS
        except:
            pass
    
    # Gamma fit and plot
    if True:
        try:
            fit_params, pcov = curve_fit(gamma_dist.cdf.calc, x, np.cumsum(y))
            pdf_ax.plot(x, gamma_dist.pdf.calc(x, *fit_params), label='Gamma')
            cdf_ax.plot(x, gamma_dist.cdf.calc(x, *fit_params), label='Gamma')
            KS = np.max(np.abs(gamma_dist.cdf.calc(x, *fit_params) - np.cumsum(y)))
            print KS
        except:
            pass
        

    plt.figure('pdf')
    pdf_ax.annotate(exponentiated_weibull.pdf.latex(*weibull_fit, short=True), xy=(0.5, 0.5), xycoords='axes fraction', ha='center')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('maybe_1.pdf')
    plt.close('pdf')

    plt.figure('cdf')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('maybe_2.pdf')
    plt.close('cdf')

    return
