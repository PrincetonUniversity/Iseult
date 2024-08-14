#! /usr/bin/env python

if __name__ == '__main__':

    import argparse
    import warnings

    parser = argparse.ArgumentParser(description='Plotting program for Tristan-MP files.')

    parser.add_argument('-n', nargs = '?',# dest='accumulate', action='store_const',
                            const=-1, default=-1, type=int,
                            help='Maximum file # to consider')
    parser.add_argument('-framerate', nargs = '?',# dest='accumulate', action='store_const',
                        const=10, default=10, type=int,
                        help='FPS for the movie')

    parser.add_argument('-outmovie', nargs = '?',# dest='accumulate', action='store_const',
                        const='out.mov', default='out.mov', type=str,
                        help='FPS for the movie')

    parser.add_argument('-O', nargs = '+',# dest='accumulate', action='store_const',
                            default=[''], type=str,
                            help='Directory Iseult will open. Default is output')

    parser.add_argument('-p', nargs = '?',# dest='accumulate', action='store_const',
                            const='Default', default='Default', type=str,
                            help='''Open Iseult with the given saved view.
                                  If the name of view contains whitespace,
                                  either it must be enclosed in quotation marks or given
                                  with whitespace removed. Name is case sensitive.''')
    parser.add_argument("-b", help="Run Iseult from bash script. Makes a movie.", type=bool,
                            action="store_true")

    parser.add_argument("-name", nargs = '+',# dest='accumulate', action='store_const',
                         default=[''], type=str,
                            help='Plot Title')

    parser.add_argument("--wait", help="Wait until current simulation is finished before making movie.",
                        action="store_true", type=bool)

    parser.add_argument("-e",
                        "--electron-spectra",
                        default=None,
                        type=str,
                        choices=['n'+str(i) for i in range(1,11)],
                        help="Which data set to use for the electron spectra. Should be of the format 'nX` where X is the integer number of the dataset")
    parser.add_argument("-i",
                        "--ion-spectra",
                        default=None,
                        type=str,
                        choices=['n'+str(i) for i in range(1,11)],
                        help="Which data set to use for the ion spectra. Should be of the format 'nX` where X is the integer number of the dataset")

    cmd_args = parser.parse_args()

    import sys, os
    sys.path.append(os.path.join(os.path.dirname(__file__), 'src'))

    if not cmd_args.b:
        from main_app import runMe
        runMe(cmd_args)
    else:
        warnings.warn('The `-b` option has not been verified to be fully functional.', RuntimeWarning)
        if cmd_args.wait:
            import subprocess, time
            slurm_num = sys.stdin.read().split[-1]
            print(slurm_num)
            num = 0
            done = False
            while num < 2000 and not done:
                    slurm_queue = subprocess.check_output(["squeue"])
                    if slurm_queue.find(slurm_num) != -1:
                        num += 1
                        time.sleep(300)

        from oengus import runMe
        print(cmd_args.name, cmd_args.O)
        runMe(cmd_args)
