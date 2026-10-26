# canfar_run.py

# Scripts of high-level functions that call canfar classes. The python
# executables are created by pyproject.toml.

import sys
from shapepipe.canfar import canfar_monitor, canfar_log_monitor

def run_log(args=None):
    """Run Log.

    Handles job logging with the canfar library

    Parameters
    ----------
    args : list, optional
        command line arguments, default is ``None``

    """
    # Create instance
    obj = canfar_monitor.Log()

    try:
        # Run instance and exit
        sys.exit(obj.run(args=args))
    except ValueError as e:
        print(e)
        sys.exit(1)


def run_monitor_log(args=None):
    """Run Monitor Log.

    Handles job monitoring with the canfar library

    Parameters
    ----------
    args : list, optional
        command line arguments, default is ``None``

    """
    # Create instance
    obj = canfar_log_monitor.LogMonitor()

    # Run instance
    obj.run(args=args)
