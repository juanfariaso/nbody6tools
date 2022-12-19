import os,sys
import matplotlib
r = os.system('{0} -c "import matplotlib.pyplot as plt; plt.figure()"'.format(sys.executable))
if r != 0:
    matplotlib.use("Agg")
from matplotlib import pyplot
from matplotlib import animation
import numpy
from nbody6tools.Reader import (read_snapshot,get_number_of_snapshots,
                                parse_inputfile,Options,get_globals,
                                get_orbit_interpolator,get_times)

def make_animation(folder,output=None,xy="xy",fps=10,dpi=None,boxsize=None,show_bound=False,**kw):
    """
    Quick animation of the simulation
    If output is given animation is saved in that file.
    Else shows animation in screen
    """
    ##TODDO - add time to title
    ##      - make it look better, but make sure is fast.
    #print("showbound",show_bound)

    def update_line(num,line,bline,cline,folder,title=None,sn0=None,draw_orbit=False):
            orbit = None
            if sn0 is None:
                sn = read_snapshot(folder,num)
            else:
                if num>0:
                    sn0.step()
                sn = sn0

            sn.to_physical()

            if (0 < sn.inputfile['KZ'][14] <= 3) and draw_orbit :
                if orbit is None:
                    orbit = get_orbit_interpolator(folder)
                #orbit.times = sn.time + 0.000001
                orbit.times = numpy.linspace(0.00001,sn.time,1000)
                xg,yg,zg = orbit.RG
            else: 
                xg,yg,zg = numpy.array([[0],[0],[0]])
            XG = dict(x = xg, y = yg, z = zg)
            #print(XG[xy[0]] - XG[xy[0]][-1])

            x = sn.stars[xy[0]] 
            y = sn.stars[xy[1]] 
            s = sn.stars["mass"]/1
            line.set_offsets(numpy.c_[x,y])
            line.set_sizes(s)
            if show_bound :
                bound_set = sn.bound_stars
                xb = bound_set[xy[0]]
                yb = bound_set[xy[1]]
                sb = bound_set["mass"]/0.5

                bline.set_offsets(numpy.c_[xb,yb])
                bline.set_sizes(sb)
            cline.set_xdata( XG[xy[0]] - XG[xy[0]][-1] )
            cline.set_ydata( XG[xy[1]] - XG[xy[1]][-1] )
            print("snapshot: ", num, sn.time) #TODO: Put nice progress info
            #title.set_text('Time %.2f Myr'%(sn.parameters["time"]*sn.parameters["tscale"] ) )
            return line,bline,cline
    opt = parse_inputfile(folder+"input")
    nsnap = get_number_of_snapshots(folder)
    if Options.has_section("CONFIG"):
        sn0 = read_snapshot(folder,0) if Options.getboolean("CONFIG","singleFile") else None
    else:
        sn0 = None

    O = None
    if sn0 is not None:
        if sn0.inputfile['KZ'][16] <= 3:
                O = get_orbit_interpolator(folder)
                O.times = get_times('./') 
                xg,yg,zg = O.RG
                XG = dict(x = xg, y = yg, z = zg)

    fig1 = pyplot.gcf()
    l = pyplot.scatter([],[],[],c="k",alpha=0.8)
    bl = pyplot.scatter([],[],[],c="r",alpha=0.8)
    if O is not None:
        pyplot.plot( XG ,c='c',zorder=-1 )

    CG,  = pyplot.plot([],[],c='c',zorder=-1)
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
    movie = animation.FuncAnimation(fig1, update_line, nsnap, fargs=(l,bl,CG,folder,sn0),
                                       interval=10, blit=True)
    if output is None :
        pyplot.show()
    else:
        movie.save(output,fps=int(fps),dpi=dpi)


def plot(folder,output=None,snapshot=1,time=None,projection="xy",
        space="position",ax=pyplot.gca(),
        point_scale = 1,
        **kw):
    if time is None:
        st = read_snapshot(folder,snapshot)
    else:
        st  = read_snapshot(folder,time=time)
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
    cmx = numpy.median(x)
    cmy = numpy.median(y)
    #s = (numpy.log10(m)+1)**2
    s = m*point_scale*10
    rh = st.stars.half_mass_radius()
    rbox = 3*rh

    #ax.axis("equal")

    ax.scatter(x,y,s = s,c="k",alpha=0.7,edgecolors="none")
    ax.set_aspect('equal')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim( cmx - rbox, cmx+rbox  )
    ax.set_ylim( cmy - rbox, cmy+rbox  )

    if output is None:
        pyplot.show()
    else:
        pyplot.savefig(output)
    return

def plot_evol(folder,output=None,parameter="fbin",
        show_parameters=False,
        yscale = "linear",
        xscale = 'linear',**kw):
    if show_parameters:
        p_globals(folder,None,show_parameters=True)
        return

    if parameter == "fbin":
        p_fbin(folder,**kw)
    elif parameter == "lagr":
        p_lagr(folder,**kw)
    elif parameter == "total_energy":
        p_total_energy(folder,**kw)
    else:
        p_globals(folder,parameter)

    pyplot.xscale(xscale)
    pyplot.yscale(yscale)

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
        mass_fractions = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        #mass_fractions = [0.01,0.1,0.5,1.0]
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

def p_total_energy(folder,ax= pyplot.gca(), **kw):
    gl = get_globals(folder)
    ax.plot(gl["Time[Myr]"], gl["BE(3)"],"-k" )
    ax.set_xlabel("Time [Myr]")
    ax.set_ylabel("Total Energy (Nbunits)")

def p_globals(folder, parameter,
        show_parameters=False,ax = pyplot.gca(),**kw):
    """Plot parameters from global.30.
    Call with  show_parameters=True for a list of available parameters.
    For definitions we refer to Nbody6 manual.
    """
    gl  =  get_globals(folder)
    if show_parameters:
        print("Available parameters from %sglboal.30: "%folder)
        print([ k   for k in gl.keys() if "TIME" not in k])
        return
    t = gl['TIME[Myr]']
    val = gl[parameter]

    ax.plot(t,val, "-")
    ax.set_ylabel(parameter)
    ax.set_xlabel('time (Myr)')


