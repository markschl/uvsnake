# Memory and time constraints for HPC:
# Approximations based on benchmarks and non-linear regression
# (could always be inproved)
# Some safety margin is added by multiplying the memory and time constraints with a factor
# In addition, constraints are doubled at every attempt of re-running after a fail
_mem_f = 2
_time_f = 4


def mem_func(mem=5, f=0, max_mem=50000):
    def _mem_func(wildcards, input, attempt):
        _mem = mem + f * input.size_mb
        return round(min(max_mem, _mem * 2 ** (attempt - 1)))

    return _mem_func


def time_func(time=1, f=0, max_time=24 * 60 * 20):
    def _time_func(wildcards, input, attempt):
        _time = time + f * input.size_mb
        return round(min(max_time, _time * 2 ** (attempt - 1)))

    return _time_func


def unoise3_memfunc(program, min_size, threads):
    if program == "usearch":
        return lambda _, input, attempt: (
            _mem_f * (180 * input.size_mb * min_size**-0.9 + 50)
        ) * 2 ** (attempt - 1)
    return lambda _, input, attempt: (
        _mem_f * (90 * input.size_mb * min_size**-1.4 * (threads - 0.46) + 50)
    ) * 2 ** (attempt - 1)


def unoise3_timefunc(program, min_size, maxrejects, threads):
    if program == "usearch":
        return lambda _, input, attempt: (
            _time_f * (30 * input.size_mb**1.6 * min_size**-2 + 5)
        ) * 2 ** (attempt - 1)
    return lambda _, input, attempt: (
        _time_f
        * (
            0.4
            * input.size_mb**1.52
            * min_size**-2
            * (threads**-0.63 - 0.58)
            * (maxrejects + 4)
            + 5
        )
    ) * 2 ** (attempt - 1)


def uparse_memfunc(min_size):
    return lambda _, input, attempt: (
        _mem_f * (30 * input.size_mb**0.5 * min_size**0.2 + 50)
    ) * 2 ** (attempt - 1)


def uparse_timefunc(min_size, maxaccepts):
    return lambda _, input, attempt: (
        _time_f * (2 * input.size_mb * min_size**-2 * maxaccepts**0.75 + 5)
    ) * 2 ** (attempt - 1)


from os.path import getsize
from math import log


def cluster_cores(program, n_amplicons, max_cores):
    if program == "usearch":
        # USEARCH v11/v12 do not appear to use more than 1 thread
        return 1
    else:
        # VSEARCH uses all cores.
        # But in case of multiple amplicons, leave some cores free for the
        # merging/trimming processes that may still be ongoing
        # TODO: may still be inefficient
        other_amplicons = n_amplicons - 1
        return max(1, max_cores//2, max_cores * (1 - other_amplicons*0.2))


def otutab_cores(program, n_amplicons, max_cores):
    if program == "usearch":
        # Only one core used for clustering:
        # don't use all available cores, but leave some for the clustering
        # processes of other primer combinations, which may take longer.
        # If all cores were used, the OTU table creation would have to wait
        # until ready other clustering/preparation processes are finished.
        return max(1, max_cores // 2, max_cores - n_amplicons)
    else:
        # VSEARCH uses multiple cores for clustering, and for OTU table
        # creation we just use everything available
        return max_cores

def otutab_memfunc(program, threads):
    if program == "usearch":
        return lambda _, input, attempt: (
            _mem_f
            * (
                4e-4
                * (log(getsize(input.uniques)) ** 8.3e-3 - 1)
                * (getsize(input.otus) + 2e6)
                * (threads + 5.1)
                + 50
            )
        ) * 2 ** (attempt - 1)
    return lambda _, input, attempt: (
        _mem_f * (2e-5 * getsize(input.otus) * (threads - 0.8) + 20)
    ) * 2 ** (attempt - 1)


def otutab_timefunc(program, maxaccepts, maxrejects, threads):
    if program == "usearch":
        return lambda _, input, attempt: (
            _time_f
            * (
                1.5e-13
                * threads**-1
                * getsize(input.uniques)
                * (getsize(input.otus) + 4.2e5)
                * maxrejects**0.65
                + 5
            )
        ) * 2 ** (attempt - 1)
    return lambda _, input, attempt: (
        _time_f
        * (
            4e-15
            * threads**-1
            * getsize(input.uniques)
            * (getsize(input.otus) + 8.4e6)
            * (maxrejects + 20)
            + 5
        )
    ) * 2 ** (attempt - 1)
