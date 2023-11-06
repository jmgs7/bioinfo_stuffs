import argparse
import json


################################### FUNCTIONS #############################################
def mkdir(dir):
    """
        Function to create a directory 
    """
    import os
    if not os.path.exists(dir):
        os.mkdir(dir)

def list_dir_files(dir,pattern = "None"):
    """
        Function to list the files of a directory
    """
    import glob
    if pattern == "None":
        files = glob.glob(f"{dir}/*") #f para indicar que es un string literal y así introducir la variable "dir" en la cadena.
    else:
        files = glob.glob(f"{dir}/*{pattern}*")
    return files


def get_sample_name(file_names):
    """
        Function to list the sample names of a list of fastq files
    """
    import os
    return(list(set([os.path.basename(file).split("_R1_")[0] for file in file_names if "_R1_" in file]))) # R1????


def eval_fastqc_file(args):
    """
        Function to run fastqc
    """
    sample_dict,output,threads,run = args # Fíjate que la notación es al revés de cualquier variable.
    import subprocess # Este módulo permite correr subprocesos de bash en tu código de python usando las pipelines de I/O.
    if run == "1": # Variable de control para evitar ejecuciones erróneas, supongo. ??????
        subprocess.run(f"fastqc {sample_dict['R1']} -o {output} -t {threads}",shell=True) # Replicate 1.
        subprocess.run(f"fastqc {sample_dict['R2']} -o {output} -t {threads}",shell=True) # Replicate 2.

def eval_fastq_files(sample_dict,output,run): # Paralelización de la función anterior. Para todas las muestras, no para ir de 1 en 1.
    """
        Function to run fastqc in multiple threads
    """
    import multiprocessing                                          # 8 threads por muestra por algo?
    with multiprocessing.Pool(len(sample_dict)) as pool:
        pool.map(eval_fastqc_file,[(sample_dict[sample_name],output,8,run) for sample_name in sample_dict]) # OJO! Mapeamos a todos los procesadores lógicos.

def run_trimming(args):  # Quitamos los adaptadores de los reads.
    """
        Function to run cutadapt    
    """
    import subprocess
    sample_name, sample_dict, num_threads, run = args
    fin1 = sample_dict["R1"]
    fin2 = sample_dict["R2"]
    fout1 = f"02_trim/{sample_name}_R1.fastq.gz"
    fout2 = f"02_trim/{sample_name}_R2.fastq.gz"
    if run == "1":
        subprocess.run(f"cutadapt --quiet -j {num_threads} -u 1 -U 1 -m 20:20 -q 20 -o {fout1} -p {fout2} {fin1} {fin2}",shell=True)
        
    return({sample_name:{"R1":fout1, "R2":fout2}})

def trimming_files(sample_dict,adapter,run): # Paralelización de la anterior, igual que lo hemos hecho antes.
    """
        Function to run cutadapt in multiple threads
    """
    import multiprocessing
    import collections
    with multiprocessing.Pool(len(sample_dict)) as pool:
        sample_dict = pool.map(run_trimming,[(sample_name,sample_dict[sample_name],adapter,8,run) for sample_name in sample_dict])
    sample_dict = dict(collections.ChainMap(*sample_dict))
    return(sample_dict)

###################################################################################################################


parser = argparse.ArgumentParser()
parser.add_argument("-I", "--input-dir")
parser.add_argument("-L", "--run")
args = vars(parser.parse_args())

input_dir = args["input_dir"]
adapter = args["adapter"]
run = args["run"]

filenames = list_dir_files(input_dir,"fastq.gz")
sample_names = get_sample_name(filenames)

sample_dict = {}
for sample_name in sample_names:
    fastq_file_r1 = [x for x in filenames if sample_name in x and "_R1_" in x][0]
    sample_dict[sample_name]["R1"] = fastq_file_r1
    fastq_file_r2 = [x for x in filenames if sample_name in x and "_R2_" in x][0]
    sample_dict[sample_name]["R2"] = fastq_file_r2

mkdir("FastQC")
mkdir("FastQC/Raw")
mkdir("FastQC/Trim")
mkdir("02_trim")

eval_fastq_files(sample_dict,"FastQC/Raw",run)
sample_dict = trimming_files(sample_dict,run)
eval_fastq_files(sample_dict,"FastQC/Trim",run)

with open("00_log/1_2_fastq.json","w") as jsonfile:
    json.dump(sample_dict,jsonfile,indent=4)
