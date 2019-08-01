

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Plotting program for Tristan-MP files.')
    #        parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                        help='The maximum file number to consider')
    #        parser.add_argument('--foo', nargs='?', const='c', default='d')
    #        parser.add_argument('bar', nargs='?', default='d')
    parser.add_argument('-n', nargs = '?',# dest='accumulate', action='store_const',
                            const=-1, default=-1,
                            help='Maximum file # to consider')
    parser.add_argument('-framerate', nargs = '?',# dest='accumulate', action='store_const',
                        const=10, default=10,
                        help='FPS for the movie')

    parser.add_argument('-outmovie', nargs = '?',# dest='accumulate', action='store_const',
                        const='out.mov', default='out.mov',
                        help='FPS for the movie')

    parser.add_argument('-O', nargs = '+',# dest='accumulate', action='store_const',
                            default=[''],
                            help='Directory Iseult will open. Default is output')

    parser.add_argument('-p', nargs = '?',# dest='accumulate', action='store_const',
                            const='Default', default='Default',
                            help='''Open Iseult with the given saved view.
                                  If the name of view contains whitespace,
                                  either it must be enclosed in quotation marks or given
                                  with whitespace removed. Name is case sensitive.''')
    parser.add_argument("-b", help="Run Iseult from bash script. Makes a movie.",
                            action="store_true")

    parser.add_argument("-name", nargs = '+',# dest='accumulate', action='store_const',
                         default=[''],
                            help='Plot Title')


    #parser.add_argument("--wait", help="Wait until current simulation is finished before making movie.",
    #                    action="store_true")

    cmd_args = parser.parse_args()
    import sys
    sys.path.append('./src')

    if not cmd_args.b:
        from main_app import runMe
        runMe(cmd_args)
    else:
        from oengus import runMe
        print(cmd_args.name, cmd_args.O)
        runMe(cmd_args)
