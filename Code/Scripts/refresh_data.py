import subprocess
from multiprocessing import Pool
from subprocess import Popen

import numpy as np


def refresh_data(cmd):
    p = Popen(
        cmd, cwd=r"D:\Work\git\GIS_tools\Code\Scripts", shell=True, stdout=subprocess.DEVNULL
    )
    return p.communicate()


def runs(cmd, file_list):
    N = len(file_list)
    if N == 0:
        raise NameError("Can't do zero runs")

    # Prepare individual workers
    cpus = 8
    pool = Pool(processes=cpus)

    # Split number of model runs in cycles with the same size as the number of workers
    cycles = np.ceil(N / cpus)
    subNrange = range(cpus)

    # Start model runs
    for c in range(cycles.astype(int)):
        # Create jobs
        results = []
        for n in subNrange:
            ix = n + c * cpus
            if ix < N:
                _cmd = cmd + file_list[ix]
                results.append(pool.apply_async(refresh_data, kwds={"cmd": _cmd}))
                print(file_list[ix])

    for r in results:
        _, stderr = r.get()

    pool.close()
    pool.join()
    return True


if __name__ == "__main__":

    cmd = r"C:\Users\Koen.Reef\Miniconda3\envs\P1414test\python.exe -m"
    file_list = [
        "ark_nzk",
        "hdsr",
        "hhd",
        "hhr",
        "hhsk",
        "markermeer",
        "noordzee",
        "noordzee_hoog",
        "ontbrekende_stuwen",
        "rijnmaasmonding_open",
        "rijnmaasmonding_closed",
        "rijnmaasmonding_hij_closed",
        "rijntakken",
        "wagv",
    ]
    runs(cmd=cmd, file_list=file_list)
