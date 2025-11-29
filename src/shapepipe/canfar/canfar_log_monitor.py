# canfar_log_monitor.py

# Monitor CANFAR job logs and save non-empty outputs

import os
import sys
import time
import subprocess
from datetime import datetime
from pathlib import Path

from canfar.sessions import Session

from cs_util import args as cs_args
from cs_util import logging


class LogMonitor:
    """Monitor CANFAR job logs continuously."""

    def __init__(self):
        self.params_default()
        self._tracked_jobs = {}  # Track which jobs we've seen

    def params_default(self):
        """Set default parameters."""
        self._params = {
            "status": "Running",
            "interval": 180,
            "output_dir": "canfar_logs",
            "kind": "headless",
            "verbose": False,
        }

        self._short_options = {
            "status": "-s",
            "interval": "-i",
            "output_dir": "-o",
            "kind": "-k",
        }

        self._types = {
            "interval": "int",
        }

        self._help_strings = {
            "status": "job status to monitor (Running, Completed, or both); default is {}",
            "interval": "polling interval in seconds; default is {}",
            "output_dir": "directory to save log files; default is {}",
            "kind": "session kind (headless, notebook, etc.); default is {}",
        }

    def set_params_from_command_line(self, args):
        """Set parameters from command line."""
        options = cs_args.parse_options(
            self._params,
            self._short_options,
            self._types,
            self._help_strings,
        )
        self._params = options

        # Save calling command
        logging.log_command(args)

    def get_job_ids(self):
        """Fetch job IDs based on status filter."""
        session = Session()

        kind = self._params["kind"] if self._params["kind"] != "all" else None

        try:
            info = session.fetch(kind=kind)
        except Exception as e:
            print(f"Error fetching sessions: {e}")
            return []

        if not info:
            return []

        # Filter by status
        status_filter = self._params["status"]
        job_ids = []

        for job in info:
            job_status = job.get("status", "")
            job_id = job.get("id", "")
            job_name = job.get("name", "")

            # Handle "both" or comma-separated statuses
            if status_filter.lower() == "both":
                if job_status in ["Running", "Completed"]:
                    job_ids.append((job_id, job_name, job_status))
            elif "," in status_filter:
                statuses = [s.strip() for s in status_filter.split(",")]
                if job_status in statuses:
                    job_ids.append((job_id, job_name, job_status))
            else:
                if job_status == status_filter:
                    job_ids.append((job_id, job_name, job_status))

        return job_ids

    def fetch_logs(self, job_id):
        """Fetch logs for a specific job ID using canfar logs command."""
        try:
            result = subprocess.run(
                ["canfar", "logs", job_id],
                capture_output=True,
                text=True,
                timeout=30,
            )
            return result.stdout
        except subprocess.TimeoutExpired:
            print(f"  Timeout fetching logs for {job_id}")
            return None
        except Exception as e:
            if self._params["verbose"]:
                print(f"  Error fetching logs for {job_id}: {e}")
            return None

    def is_log_empty(self, log_content):
        """Check if log content is just the header (3 lines with no real content)."""
        if not log_content:
            return True

        lines = log_content.strip().split('\n')

        # If only 3 or fewer lines, and they're mostly empty, consider it empty
        if len(lines) <= 3:
            non_empty_lines = [line for line in lines if line.strip() and "Logs for session" not in line]
            return len(non_empty_lines) == 0

        return False

    def save_log(self, job_id, job_name, log_content):
        """Save log content to a file."""
        output_dir = Path(self._params["output_dir"])
        output_dir.mkdir(parents=True, exist_ok=True)

        # Use job name if available, otherwise use ID
        filename = f"{job_name}_{job_id}.log" if job_name else f"{job_id}.log"
        filepath = output_dir / filename

        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        with open(filepath, 'w') as f:
            f.write(f"# Captured at {timestamp}\n")
            f.write(f"# Job ID: {job_id}\n")
            f.write(f"# Job Name: {job_name}\n\n")
            f.write(log_content)

        return filepath

    def monitor_loop(self):
        """Main monitoring loop."""
        print(f"Starting CANFAR log monitor...")
        print(f"  Status filter: {self._params['status']}")
        print(f"  Polling interval: {self._params['interval']} seconds")
        print(f"  Output directory: {self._params['output_dir']}")
        print(f"  Kind: {self._params['kind']}")
        print()

        iteration = 0

        try:
            while True:
                iteration += 1
                timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                print(f"[{timestamp}] Poll #{iteration}")

                # Get current job IDs
                job_ids = self.get_job_ids()

                if not job_ids:
                    print(f"  No jobs found with status '{self._params['status']}'")
                else:
                    print(f"  Found {len(job_ids)} job(s)")

                    for job_id, job_name, job_status in job_ids:
                        print(f"  Checking {job_id} ({job_name}, {job_status})...", end=" ", flush=True)

                        # Fetch logs
                        log_content = self.fetch_logs(job_id)

                        if log_content is None:
                            print("failed to fetch")
                            continue

                        # Check if logs have content
                        if self.is_log_empty(log_content):
                            print("empty")
                            # Track that we've seen this job
                            if job_id not in self._tracked_jobs:
                                self._tracked_jobs[job_id] = {"name": job_name, "has_logs": False, "num_lines": 0}
                        else:
                            # Count number of lines in current log
                            current_num_lines = len(log_content.strip().split('\n'))

                            # Get previous number of lines
                            if job_id in self._tracked_jobs:
                                previous_num_lines = self._tracked_jobs[job_id].get("num_lines", 0)
                            else:
                                previous_num_lines = 0

                            # Only save if there's new content
                            if current_num_lines > previous_num_lines:
                                print(f"has content! ({current_num_lines} lines, +{current_num_lines - previous_num_lines})", end=" ")

                                # Save the logs
                                filepath = self.save_log(job_id, job_name, log_content)
                                print(f"saved to {filepath}")

                                # Update tracking
                                self._tracked_jobs[job_id] = {
                                    "name": job_name,
                                    "has_logs": True,
                                    "num_lines": current_num_lines
                                }
                            else:
                                print(f"no new content ({current_num_lines} lines)")
                                # Update tracking even if no new content
                                if job_id not in self._tracked_jobs:
                                    self._tracked_jobs[job_id] = {
                                        "name": job_name,
                                        "has_logs": True,
                                        "num_lines": current_num_lines
                                    }

                print()
                time.sleep(self._params['interval'])

        except KeyboardInterrupt:
            print("\n\nMonitoring stopped by user")
            print(f"Total jobs tracked: {len(self._tracked_jobs)}")
            jobs_with_logs = sum(1 for j in self._tracked_jobs.values() if j["has_logs"])
            print(f"Jobs with logs saved: {jobs_with_logs}")

    def run(self, args=None):
        """Main entry point."""
        if args is None:
            args = sys.argv

        self.set_params_from_command_line(args)
        self.monitor_loop()
