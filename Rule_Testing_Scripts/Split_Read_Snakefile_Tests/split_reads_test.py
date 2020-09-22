#!/home-02/eanderson/Virtual_Envs/SnakeMake/bin/python 
import os, sys
import subprocess
from pathlib import Path
import yaml

def main():
    test_splitreads_snakefile()
        

def test_splitreads_snakefile():
    cwd = Path(os.getcwd())
    print(f"Current Directory: {cwd}", file=sys.stderr)
    target_dir = Path("Rule_Testing_Scripts/Split_Read_Snakefile_Tests")
    try:
        if cwd.name != "Split_Read_Snakefile_Tests":
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
            
        call_splitreads_snakefile()
    
    finally:
        try:
            os.chdir(cwd)
        except:
            print(f"Couldn't return to directory {cwd}", file=sys.stderr)
            sys.exit(1)
    

def setup_files(config):
    reads_path = "../Test_Files/1000count"
    reads_dir = Path(reads_path).resolve()
    Path(config['samples']['fq_path']).mkdir()
    reads_link_path = Path(config['samples']['fq_path']) / config['samples']['lanes'][0]
    reads_link = reads_link_path.resolve()
    reads_link.symlink_to(reads_dir)
    

def clean_files(config):
    reads_link_path = Path(config['samples']['fq_path']) / config['samples']['lanes'][0]
    reads_link_path.unlink()
    Path(config['samples']['fq_path']).rmdir()

    
def call_splitreads_snakefile():
    with open("../../config.yaml", 'r') as config_path:
        config = yaml.load(config_path, yaml.SafeLoader)

    try:
        setup_files(config)

        try:
            subprocess.check_call(['snakemake', '-R', 'split_reads'])
        except CalledProcessError:
            print("Splitreads test failed during execution", file=sys.stderr)
            sys.exit(1)
    finally:
        clean_files(config)
        
    file_tests = {}
    data_path = Path("data")
    print(os.path.isdir(data_path))
    if os.path.isdir(data_path):
        file_tests['data'] = True
    else:
        file_tests['data'] = False
    for i in ["1", "2"]:
        link = "read_" + i + ".fq.gz"
        link_path = data_path / link 
        print(os.path.islink(link_path))
        if os.path.islink(link_path):
            file_tests['read' + i] = True
            subprocess.call(['unlink', link_path])
        else:
            file_tests['read' + i] = False
          
        split_read = "split_read." + i + ".fq.gz"
        split_read_path = data_path / split_read
        print(os.path.isfile(split_read_path))
        if os.path.isfile(split_read_path):
            file_tests['split_read' + i] = True
            subprocess.call(['rm', split_read_path])
        else:
            file_tests['split_read' + i] = False
            
    subprocess.call(['rm', '-r', data_path])
    
    split_stat_path = Path("split_stat_read1.log") 
    if os.path.isfile(split_stat_path):
        file_tests['split_err'] = True
        subprocess.call(['rm',  split_stat_path])
    else:
        file_tests['split_err'] = False
    print("Splitreads file output:", file=sys.stderr)
    for key, value in file_tests.items():
        print(f"{key}: {value}")
        
    if not all(file_tests.values()):
        print("Splitreads file output failed", file=sys.stderr)
        sys.exit(1)
    else:
        print("Splitreads file output succeeded", file=sys.stderr)

            
            
if __name__ == "__main__":
    main()

