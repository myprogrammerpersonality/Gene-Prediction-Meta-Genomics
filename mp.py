import multiprocessing
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

def run_blast(sequence):
    print("Running blast for sequence {}".format(sequence))
    result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
    print(f"Finished blast for sequence {sequence} with result {result_handle}")
    return result_handle.read()

def save_result(result):
    print("Saving result", result)
    result 
    with open("blast_results.txt", "a") as f:
        f.write(result)

if __name__ == '__main__':
    sequences = ['MNLRPMTYQMN','MLNRP','MLNRP','MLNRP']
    pool = multiprocessing.Pool(processes=4)
    pool.map_async(run_blast, sequences)
    pool.close()
    pool.join()
