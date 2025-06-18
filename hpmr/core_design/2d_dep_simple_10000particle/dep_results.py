import openmc.deplete
import matplotlib.pyplot as plt
import numpy as np 

results = openmc.deplete.Results("depletion_results.h5")
time, keff = results.get_keff()
time_days = [t / 86400 for t in time] 

plt.plot(time_days, keff[:, 0])
plt.xlabel("Time [days]")
plt.ylabel("k-effective")
plt.title("k-effective vs. Depletion Time")
plt.grid(True)
plt.savefig("keff_vs_time.png")

# Finding the cycle length by interpolation
cycle_length = None
for i in range(1, len(keff[:, 0])):
    if keff[i, 0] < 1.0 and keff[i-1, 0] >= 1.0:
        # Linear interpolation to find the exact time when k=1.0
        k1, k2 = keff[i-1, 0], keff[i, 0]
        t1, t2 = time_days[i-1], time_days[i]
        cycle_length = t1 + (1.0 - k1) * (t2 - t1) / (k2 - k1)
        break

cycle_length = round(cycle_length, 0)
print("Cycle length when k=1.0:", cycle_length)


