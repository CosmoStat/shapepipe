                # Get JobHandler objects
                timeout = jh.timeout
                # Get file handler objects
                run_dirs = jh.filehd.module_run_dirs
                module_runner = jh.filehd.module_runners[module]
                worker_log = jh.filehd.get_worker_log_name
                # Define process list
                process_list = jh.filehd.process_list
                # Define job list
                jobs = split_mpi_jobs(process_list, comm.size)
                del process_list
        else:
            job_type, run_dirs, module_runner, worker_log, timeout, jobs = \
             (None, None, None, None, None, None)

        job_type = comm.bcast(job_type, root=0)

        if job_type == 'parallel':

            # Broadcast objects to all nodes
            run_dirs = comm.bcast(run_dirs, root=0)
            module_runner = comm.bcast(module_runner, root=0)
            worker_log = comm.bcast(worker_log, root=0)
            timeout = comm.bcast(timeout, root=0)
            jobs = comm.scatter(jobs, root=0)

            # Submit the MPI jobs and gather results
            results = comm.gather(submit_mpi_jobs(jobs, config, timeout,
                                  run_dirs, module_runner, worker_log,
                                  verbose), root=0)

            del run_dirs, module_runner, timeout, jobs

        if master:
            # Assign worker dictionaries
