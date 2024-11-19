from contextlib import contextmanager
import sys
import yaml


def dump_config(config, output):
    config = {k: v for k, v in config.items() if not k.startswith("_")}
    with open(output, "w") as o:
        yaml.safe_dump(config, o)


# same code in every Python script due to https://github.com/snakemake/snakemake/issues/1632
@contextmanager
def file_logging(f):
    with open(f, "w") as log:
        sys.stderr = log
        try:
            yield
        except Exception:
            import traceback
            log.write(traceback.format_exc())
            sys.exit(1)


with file_logging(snakemake.log[0]):
    dump_config(
        snakemake.params.config,
        snakemake.output[0]
    )
