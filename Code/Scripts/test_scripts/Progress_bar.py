# Test progress bar
from tqdm import tqdm
import time

tot = 6
pbar = tqdm(total=tot)

for i in range(tot):
    pbar.set_description(f"Progress {round(i/(tot-1),1)}")
    pbar.update(1)
    time.sleep(2)


