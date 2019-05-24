from optparse import OptionParser
from argparse import ArgumentParser
from nbody6tools import Plots

def executable_options():
    parser = ArgumentParser()
    subparsers = parser.add_subparsers(help='Available Commands')
    options = parser.add_argument_group("options")

    parser.add_argument("folder",help="simulation folder", default="./") 
    ##common options
    options.add_argument("-o","--output",help="output file. Default shows without saving", default=None) 

    #animate
    parse_animate = subparsers.add_parser("animate",help="shows animation of simulation, may be slow in big clusters.")
    parse_animate.set_defaults(fun = Plots.make_animation )
    parse_animate.add_argument("--fps",help = "frames per second if an output file is given",default=10)
    parse_animate.add_argument("--xy",help="Projection of animation. Default 'xy'",default = "xy",
            choices=["xy","yx","xz","zx","yz","zy"])

    #evol
    parse_evol = subparsers.add_parser("evol",help="shows the evolution of a parameter")
    parse_evol.set_defaults(fun = Plots.evol )
    parse_evol.add_argument("-p","--parameter",help="paramter to show. Default fbin.",default = "fbin",
            choices=["fbin","lagr"])
    
    #plot
    parse_evol = subparsers.add_parser("plot",help="plot position[velocity] of stars")
    parse_evol.set_defaults(fun = Plots.plot )
    parse_evol.add_argument("-s","--space",help="space to plot. Default position.",default = "position",
            choices=["position","velocity"])
    return parser
   
def execute_actions(args) :
    args.fun(**vars(args))




