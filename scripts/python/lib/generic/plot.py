import pylab as plt
from itertools import cycle


class color:
    def __init__(self, palette):

        if palette == 'tableau20':

            # These are the "Tableau 20" colors as RGB
            colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
                      (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
                      (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
                      (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
                      (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

        # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
        for i in range(len(colors)):    
            r, g, b = colors[i]    
            colors[i] = (r / 255., g / 255., b / 255.)

        # cycle
        self.colors_cyc = cycle(colors)


    def next(self):
        return next(self.colors_cyc)


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


