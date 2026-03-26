algo_list = ["hbf", "hpf-p", "hpf-p-serial"]

parallel_algo = ["pbbs"]
number_cores = [1, 2, 4, 8, 16, 32, 64]

#test_list = ["BVZ-tsukuba11.max"]
test_list = []

timeout = 3600 #1h

print_output = False

extra_data = []

import logging
import subprocess
import re
import os
import datetime

def compile_all():
    print("Compiling...")

    #hbf
    subprocess.run(["make", "clean"], cwd="hbf", stdout=None)
    subprocess.run(["make"], cwd="hbf")

    #hpf-p
    subprocess.run(["make", "clean"], cwd="hpf-p", stdout=None)
    subprocess.run(["make"], cwd="hpf-p")

    #hpf-p-serial
    subprocess.run(["make", "serial"], cwd="hpf-p")

    #dinic and push-relabel
    #subprocess.run(["cmake", ".."], cwd="../../build", stdout=None)
    #subprocess.run(["make"], cwd="../../build")

    #pbbs todo
    print("Compiling done")

def run_alg(algo, test, cores=1):

    if algo == "hbf":
        proc = subprocess.run(["./pseudo_lifo", "/data/mchiriac/" + test], cwd="hbf/bin",
                                encoding='utf-8', stdout=subprocess.PIPE,
                                timeout=timeout)
        output = proc.stdout

        time = float(re.search(r'Time to min cut\s*:\s*(\S+)', output).group(1))
        answer = int(re.search(r'Max Flow\s*:\s*(\S+)', output).group(1))
        
        for line in output.splitlines():
            if line.startswith('data'):
                extra_data.append(line)

        #print(extra_data)
        
        return (time, answer)

    if algo == "hpf-p":
        proc = subprocess.run(["./pseudo_fifo", "/data/mchiriac/" + test], cwd="hpf-p/bin",
                                encoding='utf-8', stdout=subprocess.PIPE,
                                timeout=timeout)
        output = proc.stdout

        time = float(re.search(r'Time to min cut\s*:\s*(\S+)', output).group(1))
        answer = int(re.search(r'Max Flow\s*:\s*(\S+)', output).group(1))

        for line in output.splitlines():
            if line.startswith('data'):
                extra_data.append(line)

        return (time, answer)

    
    if algo == "hpf-p-serial":
        proc = subprocess.run(["./pseudo_fifo_serial", "/data/mchiriac/" + test], cwd="hpf-p/bin",
                                encoding='utf-8', stdout=subprocess.PIPE,
                                timeout=timeout)
        output = proc.stdout

        time = float(re.search(r'Time to min cut\s*:\s*(\S+)', output).group(1))
        answer = int(re.search(r'Max Flow\s*:\s*(\S+)', output).group(1))

        for line in output.splitlines():
            if line.startswith('data'):
                extra_data.append(line)

        return (time, answer)


    if algo == "pbbs":
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = str(cores)

        proc = subprocess.run(["./maxFlow", "/data/mchiriac/" + test], cwd="pbbs",
                        encoding='utf-8', stdout=subprocess.PIPE,
                        env=env,
                        timeout=timeout)

        output = proc.stdout
        time = float(re.search(r'PBBS-time:\s*(\S+)', output).group(1))
        answer = int(int(re.search(r'flow=(\S+)', output).group(1)))

        if print_output: 
            print(output)

        return (time, answer)

    return (-1, -1)

def run_all():
    for test in test_list:
        logger.info("Running test " + test)
        print("Running test ", test)

        for algo in algo_list:
            if algo in parallel_algo:
                    for cores in number_cores:
                        try:
                            result = run_alg(algo, test, cores)

                            logger.info("Algorithm " + algo + " " +
                            str(cores) + " cores "
                                " got: time: " + str(result[0]) +
                                ", flow: " + str(result[1]))

                        except subprocess.TimeoutExpired:
                            logger.info("Algorithm " + algo + " " +
                            str(cores) + " cores timed out")

                        except:
                            logger.info("Algorithm " + algo + " " +
                            str(cores) + " cores failed")

            else:
                try:
                    result = run_alg(algo, test)

                    logger.info("Algorithm " + algo + 
                        " got: time: " + str(result[0]) +
                        ", flow: " + str(result[1]))

                    #print(extra_data)

                    for line in extra_data:
                        logger.info(line)
                    extra_data.clear()

                except subprocess.TimeoutExpired:
                    logger.info("Algorithm " + algo + " timed out")

                except:
                    logger.info("Algorithm " + algo + " failed")

        print("")



date = datetime.datetime.now()
logger = logging.getLogger(__name__)
logging.basicConfig(filename='logs/results' + date.strftime("%d%b_%H:%M") + '.log', 
                    encoding='utf-8', level=logging.DEBUG)

compile_all()

if len(test_list) == 0:
    test_list = os.listdir('/data/mchiriac')
    print(test_list)

run_all()