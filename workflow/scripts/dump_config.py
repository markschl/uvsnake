import yaml

from utils import file_logging


def dump_config(config, output):
    config = {k: v for k, v in config.items() if not k.startswith("_")}
    with open(output, "w") as o:
        yaml.safe_dump(config, o)


with file_logging(snakemake.log[0]):
    dump_config(
        snakemake.params.config,
        snakemake.output[0]
    )
