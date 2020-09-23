#!/home-02/eanderson/Virtual_Envs/SnakeMake/bin/python 
import os, sys
import subprocess
import yaml
from pathlib import Path

def main():
    test_frag_len_snakefile()
        

def test_frag_len_snakefile():
    # Make sure we're in the appropriate directory to run this script
    cwd = Path.cwd()
    print(f"Current Directory: {cwd}", file=sys.stderr)
    target_dir = Path("Rule_Testing_Scripts/Calc_Frag_Len_Tests")
    try:
        # expectation is that someone would attempt to run
        # this from the main directory
        # if so attempt to move to the Calc_Frag_Len_Tests dir
        if cwd.name != "Calc_Frag_Len_Tests":
            if target_dir.exists():
                try:
                    os.chdir(target_dir)
                except:
                    print(f"{target_dir} exists but can't be set as current working directory", file=sys.stderr)
                    sys.exit(1)

            # try another route
            elif target_dir.name.exists():
                try:
                    os.chdir(target_dir.name)
                except:
                    print(f"{target_dir.name} exists but can't be set as current working directory", file=sys.stderr)
                    sys.exit(1)

            # give up and exit
            else:
                print("Can't find appropriate testing directory to test split_reads.snakefile", file=sys.stderr)
                sys.exit(1)
            
        call_frag_len_snakefile()
    
    finally:
        try:
            os.chdir(cwd)
        except:
            print(f"Couldn't return to directory {cwd}", file=sys.stderr)
            sys.exit(1)
    

def setup_files(sample):
    # create links to the test files
    # first get the path to the test bam
    bampath = "../Test_Files/hg001.sort.rmdup.bam"
    # resolve to it's absolute path
    bamfile = Path(bampath).resolve()
    # add .bai for the index
    bamindex = Path(bampath + ".bai").resolve()

    # Create an Align/ dir and desired bam path
    Path("Align").mkdir()
    bamlinkpath = "Align/hg001.sort.rmdup.bam"
    # link the bam and bai at the desired location using absolute paths
    bamlink = Path(bamlinkpath).resolve()
    bamindexlink = Path(bamlinkpath + ".bai").resolve()
    bamlink.symlink_to(bamfile)
    bamindexlink.symlink_to(bamindex)


def remove_files():
    # iterate through Align/ and unlink all files
    for path in Path("./").iterdir():
        if path.is_dir() and not path.name.startswith("."):
            for sub_path in path.iterdir():
                sub_path.unlink()
            print(f"Cleaning dir {path}")
            # remove the Align/ directory
            path.rmdir()

    
def call_frag_len_snakefile():
    try:
        # load in config file with pyyaml
        with open("../../config.yaml", 'r') as config_path:
            config = yaml.load(config_path, yaml.SafeLoader)

        # pass config ID to setup files to ensure everything is named correctly
        setup_files(config['samples']['id'])
        # try running snakemake for calc_frag_len
        try:
            subprocess.check_call(['snakemake', '-R', 'calc_frag_len'])
        # except called process error and note failure
        except CalledProcessError:
            print("Calc_Frag_Len test failed during execution", file=sys.stderr)
            sys.exit(1)
    # clean mock files
    finally:
        remove_files()
        
        
if __name__ == "__main__":
    main()

