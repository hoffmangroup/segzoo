import snakemake
from os import path
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')

def main():
    """Entry point for the application script"""
    here = path.abspath(path.dirname(__file__))
    snakemake.snakemake(path.join(here, "Snakefile")) #dryrun=True, printshellcmds=True)
