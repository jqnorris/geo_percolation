import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import distributions


def plot_distribution(raw_data, cumulative = False, fit = 'All', file_name = None, title = None, x_label = None, y_label= None):
    # Process data into x, y, and y_sum arrays
    if len(raw_data) == 2:
        x = np.array(raw_data[0])
        if cumulative == True:
            y_sum = raw_data[1]
            y = np.insert(np.diff(y_sum), 0, y_sum[0])
        else:
            y = raw_data[1]
            y_sum = np.cumsum(y)
    else:
        try:            
            x = raw_data[:, 0]
            if cumulative == True:
                y_sum = raw_data[:, 1]
                y = np.insert(np.diff(y_sum), 0, y_sum[0])
            else:
                y = raw_data[:, 1]
                y_sum = np.cumsum(y)
        except:
            print "Data is not a 2D array, cannot plot or fit."
    
    # Create pdf and cdf plots without any fits
    if fit == None:
        0
        
    
    # Generate plots for pdf and cdf
    pdf_fig = plt.figure('pdf', figsize=(8,4))
    pdf_ax = pdf_fig.add_axes([0.1, 0.1, 0.6, 0.8])
    
    cdf_fig = plt.figure('cdf', figsize=(8,4))
    cdf_ax = cdf_fig.add_axes([0.1, 0.1, 0.6, 0.8])
    
    
    # Plot the data
    pdf_ax.plot(x, y, label='Data', color='black', ls='--', lw=2)
    cdf_ax.plot(x, y_sum, label='Data', color='black', ls='--', lw=2)
    
     # Create a dictionary for choosing which distribution is the best fit
     # KS:[name, pdf_artist, cdf_artist]
    fits = {}
    
    # Loop over all distributions
    for name, dist in distributions.names.iteritems():
    #try:
        print name
        fit_params, pcov = curve_fit(dist.cdf.calc, x, y_sum)
        print 'Fitting Done!'
        print fit_params
        plt_1, = pdf_ax.plot(x, dist.pdf.calc(x, *fit_params), label=name)
        y_sum_fit = dist.cdf.calc(x, *fit_params)
        plt_2, = cdf_ax.plot(x, y_sum_fit,  label=name)
        KS = np.max(np.abs(y_sum_fit - y_sum))
        fits[KS] =  [name, plt_1, plt_2]
    #except:
    #    pass
       
    # Sort the fits by KS distance
    sorted_keys = sorted(fits.keys())
    
    print sorted_keys
    
    

    plt.figure('pdf')
    pdf_ax.annotate(distributions.inverse_gaussian.cdf.latex(), xy=(0.5, 0.5), xycoords='axes fraction', ha='center')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('maybe_1.pdf')
    plt.close('pdf')

    plt.figure('cdf')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.savefig('maybe_2.pdf')
    plt.close('cdf')

    return
