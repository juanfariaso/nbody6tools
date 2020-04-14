import os,sys
import matplotlib
r = os.system('{0} -c "import matplotlib.pyplot as plt; plt.figure()"'.format(sys.executable))
if r != 0:
    matplotlib.use("Agg")
from matplotlib import pyplot
from matplotlib import animation
import numpy
from nbody6tools.Reader import read_snapshot,get_number_of_snapshots,parse_inputfile,Options

def make_animation(folder,output=None,xy="xy",fps=10,dpi=None,boxsize=None,show_bound=False,**kw):
    """
    Quick animation of the simulation
    If output is given animation is saved in that file.
    Else shows animation in screen
    """
    ##TODDO - add time to title
    ##      - make it look better, but make sure is fast.
    print("showbound",show_bound)

    def update_line(num,line,bline,folder,sn0=None):

            if sn0 is None:
                sn = read_snapshot(folder,num)
            else:
                if num>0:
                    sn0.step()
                sn = sn0

            sn.to_physical()

            x = sn.stars[xy[0]]
            y = sn.stars[xy[1]]
            s = sn.stars["mass"]/1
            line.set_offsets(numpy.c_[x,y])
            line.set_sizes(s)
            if show_bound :
                bound_set = sn.stars.bound_set()
                xb = bound_set[xy[0]]
                yb = bound_set[xy[1]]
                bline.set_offsets(numpy.c_[xb,yb])
                bline.set_sizes(s)
            print("snapshot: ", num, sn.time) #TODO: Put nice progress info
            #title.set_text('Time %.2f Myr'%(sn.parameters["time"]*sn.parameters["tscale"] ) )
            return line,bline
    opt = parse_inputfile(folder+"input")
    nsnap = get_number_of_snapshots(folder)
    sn0 = read_snapshot(folder,0) if Options.getboolean("CONFIG","singleFile") else None
    fig1 = pyplot.gcf()
    l = pyplot.scatter([],[],[],c="k",alpha=0.8)
    bl = pyplot.scatter([],[],[],c="r",alpha=0.8)
    if boxsize is None:
        size = 10*opt["RBAR"]
    else:
        size = boxsize/2.0
    pyplot.xlim(-size, size)
    pyplot.ylim(-size, size)
    pyplot.xlabel( "%s (pc)"%xy[0] )
    pyplot.ylabel( "%s (pc)"%xy[1] )
    ax = pyplot.gca()
    ax.set_aspect("equal")
    #title = ax.set_title("Time : %.2f"%0.0) #not working
    movie = animation.FuncAnimation(fig1, update_line, nsnap, fargs=(l,bl,folder,sn0),
                                       interval=10, blit=True)
    if output is None :
        pyplot.show()
    else:
        movie.save(output,fps=fps,dpi=dpi)

def evol(folder,output=None,parameter="fbin",**kw):
    if parameter == "fbin":
        p_fbin(folder,**kw)
    elif parameter == "lagr":
        p_lagr(folder,**kw)

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

def p_lagr(folder,ax = pyplot.gca(),kind="lagr",inputfile="input",
        #mass_fractions = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        mass_fractions = [0.01,0.1,0.5,1.0]
        ,**kw):
    """ Plots the lagrangian radii from lagr.7 file. Requires KZ(7) >= 3 """

    #opt = parse_inputfile(folder+inputfile)
    sn = read_snapshot(folder)
    opt = sn.inputfile
    tscale = sn.parameters["tscale"]

    if opt["KZ"][7] < 3:
        raise("ERROR: No lagr.7 file. KZ(7)<3 in inputfile")
    lind = []
    if kind == "lagr" :
        lind = list(range(1,19)) 
    data = numpy.loadtxt(folder+"lagr.7",skiprows=2).T
    header = open(folder+"lagr.7","r").readlines()[1].split()
    mfs = [ "%.2E"%m for m in mass_fractions ] 
    ylabel = "log(r) for "
    time = data[0]

    j = 0
    for i in  lind :
        if header[i] in mfs :
            if j > 0 : ylabel+=", "
            ylabel+=" %.0f%%"%(mass_fractions[j]*100)
            ax.plot(time*tscale,numpy.log10(data[i]*opt["RBAR"]),"-k",lw=0.5)
            j+=1
    ax.set_ylabel(ylabel)
    ax.set_xlabel("time (Myr)")
    #ax.set_xscale("log")
