import sys
from shapepipe.canfar import canfar_submit, canfar_monitor, canfar_log_monitor

def run_job(args=None):

    # Create instance
    obj = canfar_submit.Job()

    obj.run(args=args)


def run_log(args=None):

    # Create instance
    obj = canfar_monitor.Log()
    try:
        sys.exit(obj.run(args=args))
    except ValueError as e:
        print(e)
        sys.exit(1)


def run_monitor_log(args=None):

    # Create instance
    obj = canfar_log_monitor.LogMonitor()
    obj.run(args=args)


def main(argv=None):
    """Main

    Main program
    """
    # Scripts to call canfar classes are created by pyproject.toml
    return 0
