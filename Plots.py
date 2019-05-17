import os,sys
import matplotlib
r = os.system('{0} -c "import matplotlib.pyplot as plt; plt.figure()"'.format(sys.executable))
if r != 0:
    matplotlib.use("Agg")
from matplotlib import pyplot
from matplotlib import animation
import numpy
from nbody6tools.Reader import read_snapshot,get_number_of_snapshots

def make_animation(folder,output=None,**kw):
    """
    Quick animation of the simulation
    If output is given animation is saved in that file.
    Else shows animation in screen
    """
    ##TODDO - add time to title
    ##      - make it look better, but make sure is fast.

    def update_line(num,line, folder):
            sn = read_snapshot(folder,num+1)
            x = sn.stars["x"]
            y = sn.stars["y"]
            line.set_data(x,y)
            return line,

    nsnap = get_number_of_snapshots(folder)
    fig1 = pyplot.figure()

    data = numpy.random.rand(2, 25)
    l, = pyplot.plot([], [], 'k.')
    pyplot.xlim(-10, 10)
    pyplot.ylim(-10, 10)
    pyplot.title('test')
    movie = animation.FuncAnimation(fig1, update_line, nsnap, fargs=(l,folder),
                                       interval=10, blit=True)
    if output is None :
        pyplot.show()
    else:
        movie.save(output)

def evol(folder,output=None,parameter="fbin",**kw):
    if parameter == "fbin":
        p_fbin(folder,**kw)

    if output is None:
        pyplot.show()
    else:
        pyplot.savefig(output)
    return

def plot(folder,output=None,snapshot=1,projection="xy",space="position",ax=pyplot.gca(),**kw):

    st = read_snapshot(folder,snapshot)
    st.to_physical()
    if space == "position" : 
        fmt = "%s" 
        xlabel = "$%s$ (pc)"%projection[0]
        ylabel = "$%s$ (pc)"%projection[1]
    else: 
        fmt = "v%s"
        xlabel = "$v_%s$ (km/s)"%projection[0]
        ylabel = "$v_%s$ (km/s)"%projection[1]

    x = st.stars[fmt%projection[0]]
    y = st.stars[fmt%projection[1]]
    m = st.stars["mass"]
    #s = (numpy.log10(m)+1)**2
    s = m/500

    ax.axis("equal")

    ax.scatter(x,y,s = s,c="k",alpha=0.7,edgecolors="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if output is None:
        pyplot.show()
    else:
        pyplot.savefig(output)
    return


def p_fbin(folder,ax = pyplot.gca(),**kw):
    """ Plot binary fractions using file global.30
        This means shows only the regularized binaries
        Assumes the first line of global.30 is correct.
    """
    globalfile = "%s/global.30"%folder
    header = numpy.array(open(globalfile,"r").readline().split())
    data = numpy.loadtxt(globalfile,skiprows=1).T
    np = data[ numpy.where(header == "NPAIRS")[0][0] -1 ]
    ns = data[ numpy.where(header == "NS")[0][0] -1 ]
    t = data[1]  

    ax.set_ylabel("regularized binary fraction")
    ax.set_xlabel("time [Myr]")
    ax.plot(t, np/(ns+np), "-")


