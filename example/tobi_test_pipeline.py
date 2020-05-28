import sys
sys.path.append('/Users/tliaudat/Documents/PhD/codes/venv_p3/ShapePipe/')
sys.argv = ['/Users/tliaudat/Documents/PhD/codes/venv_p3/ShapePipe/shapepipe_run.py',
            '-c',
            '/Users/tliaudat/Documents/PhD/codes/venv_p3/ShapePipe/my_config/mccd_complete.ini']

import shapepipe_run as shpipe

shpipe.main()