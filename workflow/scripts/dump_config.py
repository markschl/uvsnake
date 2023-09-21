from contextlib import contextmanager
import sys
import yaml


def dump_config(config, output):
    config = {k: v for k, v in config.items() if not k.startswith("_")}
    with open(output, "w") as o:
        yaml.safe_dump(config, o)


@contextmanager
def file_logging(f):
    with open(f, "w") as handle:
        sys.stderr = sys.stdout = handle
        yield


with file_logging(snakemake.log[0]):
    dump_config(
        snakemake.params.config,
        snakemake.output[0]
    )
