"""COVERAGE_RUN

Console entry point for plotting a coverage map.

Building the map is the Snakemake ``coverage_map`` rule's job; plotting a
finished ``.hsp`` is a human act on a durable product, so it stays a
hand-run command, on the same argument that keeps ``run_report.py`` out of the
DAG. The plot windows for the UNIONS SGC and NGC fields are in
``workflow/config.yaml``'s ``coverage:`` block.

Author: Martin Kilbinger <martin.kilbinger@cea.fr>

"""

from shapepipe.utilities.coverage_plotter import CoveragePlotter


def run_plot_coverage(args=None):
    """Run Plot Coverage.

    Plot HealSparse coverage map.

    Parameters
    ----------
    args : list, optional
        command line arguments

    Returns
    -------
    int
        exit code

    """
    obj = CoveragePlotter()
    return obj.run(args=args)
