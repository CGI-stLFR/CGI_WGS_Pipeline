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
    # set the read directory path and resolve it to its absolute path
    reads_path = "../Test_Files/1000count"
    reads_dir = Path(reads_path).resolve()
    # Make the samples fq directory
    Path(config['samples']['fq_path']).mkdir()
    # create the symlink target and resolve it
    reads_link_path = Path(config['samples']['fq_path']) / config['samples']['lanes'][0]
    reads_link = reads_link_path.resolve()
    # symlink the directory as the lane within the fq_path dir
    reads_link.symlink_to(reads_dir)
    

def clean_files(config):
    # unlink the symlinked lane
    reads_link_path = Path(config['samples']['fq_path']) / config['samples']['lanes'][0]
    reads_link_path.unlink()
    # remove the fq_path dir
    Path(config['samples']['fq_path']).rmdir()

    
def call_splitreads_snakefile():
    with open("../../config.yaml", 'r') as config_path:
        config = yaml.load(config_path, yaml.SafeLoader)

    try:
        # attempt to setup files
        # clean if fails
        setup_files(config)

        try:
            # attempt to run snakemake with rule split_reads
            subprocess.check_call(['snakemake', '-R', 'split_reads'])
        except CalledProcessError:
            print("Splitreads test failed during execution", file=sys.stderr)
            sys.exit(1)
    finally:
        clean_files(config)
        

    # Test to make sure the correct output files 
    # are generated as expected
    file_tests = {}
    # Make sure there's a directory )
    data_path = Path("data")
    print(os.path.isdir(data_path))
    # Note success of first test
    if os.path.isdir(data_path):
        file_tests['data'] = True
    else:
        file_tests['data'] = False
    # check for the correctly linked starting files
    for i in ["1", "2"]:
        link = "read_" + i + ".fq.gz"
        link_path = data_path / link 
        print(os.path.islink(link_path))
        # Note outcome of test, clean file if True
        if os.path.islink(link_path):
            file_tests['read' + i] = True
            subprocess.call(['unlink', link_path])
        else:
            file_tests['read' + i] = False
        
        # check for the correct output files
        split_read = "split_read." + i + ".fq.gz"
        split_read_path = data_path / split_read
        print(os.path.isfile(split_read_path))
        # note outcome of test, clean if True
        if os.path.isfile(split_read_path):
            file_tests['split_read' + i] = True
            subprocess.call(['rm', split_read_path])
        else:
            file_tests['split_read' + i] = False
            
    # remove data/
    subprocess.call(['rm', '-r', data_path])
    
    # check log file exists
    split_stat_path = Path("split_stat_read1.log")
    # Note outcome and clean
    if os.path.isfile(split_stat_path):
        file_tests['split_err'] = True
        subprocess.call(['rm',  split_stat_path])
    else:
        file_tests['split_err'] = False
    # Print all test results
    print("Splitreads file output:", file=sys.stderr)
    for key, value in file_tests.items():
        print(f"{key}: {value}")
        
    # Print overall success
    if not all(file_tests.values()):
        print("Splitreads file output failed", file=sys.stderr)
        sys.exit(1)
    else:
        print("Splitreads file output succeeded", file=sys.stderr)

            
            
if __name__ == "__main__":
    main()

