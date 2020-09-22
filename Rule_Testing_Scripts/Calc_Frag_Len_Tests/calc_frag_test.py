#!/home-02/eanderson/Virtual_Envs/SnakeMake/bin/python 
import os, sys
import subprocess
import yaml
from pathlib import Path

def main():
    test_frag_len_snakefile()
        

def test_frag_len_snakefile():
    cwd = Path.cwd()
    print(f"Current Directory: {cwd}", file=sys.stderr)
    target_dir = Path("Rule_Testing_Scripts/Calc_Frag_Len_Tests")
    try:
        if cwd.name != "Calc_Frag_Len_Tests":
            if target_dir.exists():
                try:
                    os.chdir(target_dir)
                except:
                    print(f"{target_dir} exists but can't be set as current working directory", file=sys.stderr)
                    sys.exit(1)

            elif target_dir.name.exists():
                try:
                    os.chdir(target_dir.name)
                except:
                    print(f"{target_dir.name} exists but can't be set as current working directory", file=sys.stderr)
                    sys.exit(1)

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
    bampath = "../Test_Files/hg001.sort.rmdup.bam"
    bamfile = Path(bampath).resolve()
    bamindex = Path(bampath + ".bai").resolve()

    Path("Align").mkdir()
    bamlinkpath = "Align/hg001.sort.rmdup.bam"
    bamlink = Path(bamlinkpath).resolve()
    bamindexlink = Path(bamlinkpath + ".bai").resolve()
    bamlink.symlink_to(bamfile)
    bamindexlink.symlink_to(bamindex)


def remove_files():
    for path in Path("./").iterdir():
        if path.is_dir() and not path.name.startswith("."):
            for sub_path in path.iterdir():
                sub_path.unlink()
            print(f"Cleaning dir {path}")
            path.rmdir()

    
def call_frag_len_snakefile():
    try:
        with open("../../config.yaml", 'r') as config_path:
            config = yaml.load(config_path, yaml.SafeLoader)

        setup_files(config['samples']['id'])
        try:
            subprocess.check_call(['snakemake', '-R', 'calc_frag_len'])
        except CalledProcessError:
            print("Calc_Frag_Len test failed during execution", file=sys.stderr)
            sys.exit(1)
    finally:
        remove_files()
        
        
if __name__ == "__main__":
    main()

