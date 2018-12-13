import pylab as plt


def plot_init():

    fs = 12
    fig, ax = plt.subplots()
    #fig = plt.figure()

    ax = plt.gca()
    ax.yaxis.label.set_size(fs)
    ax.xaxis.label.set_size(fs)

    plt.tick_params(axis='both', which='major', labelsize=fs)

    plt.rcParams.update({'figure.autolayout': True})

    #return fig, ax
    return ax


