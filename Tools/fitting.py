import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.optimize import curve_fit
from distributions import *


def plot_with_fits(raw_data, type='cdf'):
    # Separate data into bin_edge and count arrays
    if len(raw_data) == 2:
        bin_edges = raw_data[0]
        bin_centers = bin_edges[:-1] + 0.5 * (bin_edges[1:] - bin_edges[:-1])
        data = raw_data[1]
    else:
        try:
            bin_edges = raw_data[:, 0]
            bin_centers = bin_edges[:-1] + 0.5 * (bin_edges[1:] - bin_edges[:-1])
            data = raw_data[:, 1]
        except:
            print "Data is not a 2D array, cannot plot or fit."
    
    # Calculate pdf and cdf data
    if(type == 'cdf'):        
        data_pdf = data / (bin_edges[1:] - bin_edges[:-1])
        data_cdf = data
    
    if(type == 'pdf'):
        data_pdf = data
        data_cdf = np.cumsum(data * (bin_edges[1:] - bin_edges[:-1]))
    
    # Plot data
    plt.close()
    pdf_fig = plt.figure('pdf', figsize=(6,4))
    pdf_ax = pdf_fig.add_axes([0.1, 0.1, 0.6, 0.8])
    pdf_ax.plot(bin_centers, data_pdf, label='Data', color='black', ls='None', marker='.')
    
    
    cdf_fig = plt.figure('cdf', figsize=(6,4))
    cdf_ax = cdf_fig.add_axes([0.1, 0.1, 0.6, 0.8])
    cdf_ax.plot(bin_edges[1:], data_cdf, label='Data', color='black', ls='--', lw=2)

    fits = {}

    for dist in names:
        #print eval(dist).string
        #fit_params, pcov = curve_fit(eval(dist).cdf.calc, x, data_cdf, p0 = eval(dist).cdf.initial_guess)
        #KS = np.max(np.abs(eval(dist).cdf.calc(x, *fit_params) - data_cdf))
        #print eval(dist).string, KS
        try:
            fit_params, pcov = curve_fit(eval(dist).cdf.calc, bin_edges[1:], data_cdf, p0 = eval(dist).cdf.initial_guess)
            pdf_ax.plot(bin_centers, eval(dist).pdf.calc(bin_centers, *fit_params), label=eval(dist).string)
            cdf_ax.plot(bin_edges[1:], eval(dist).cdf.calc(bin_edges[1:], *fit_params), label=eval(dist).string)
            KS = np.max(np.abs(eval(dist).cdf.calc(bin_edges[1:], *fit_params) - data_cdf))
            print eval(dist).string, KS, fit_params
            fits[KS] = eval(dist).pdf.latex(*fit_params, short=True)
        except:
            pass
        
    best_fit = [fits[key] for key in sorted(fits.keys())][0]

    plt.figure('pdf')
    pdf_ax.annotate(best_fit, xy=(0.5, 0.5), xycoords='axes fraction', ha='center')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=12)
    #plt.xscale('log')
    plt.savefig('maybe_1.pdf')
    plt.close('pdf')

    plt.figure('cdf')
    cdf_ax.annotate(best_fit, xy=(0.5, 0.5), xycoords='axes fraction', ha='center')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=12)
    #plt.xscale('log')
    plt.savefig('maybe_2.pdf')
    plt.close('cdf')

    return
