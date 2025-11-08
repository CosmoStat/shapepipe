import sys
from shapepipe.canfar import canfar_submit, canfar_monitor

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

def main(argv=None):
    """Main

    Main program
    """
    # Scripts to call canfar classes are created by pyproject.toml
    return 0
