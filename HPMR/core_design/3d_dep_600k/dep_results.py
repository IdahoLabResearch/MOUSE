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

keff_mean = keff[:, 0]
k1, k2 = keff_mean[-2] , keff_mean[-1]
t1, t2 = time_days[-2] , time_days[-1] 

slope = (k2-k1)/(t2-t1)
cycle_lenght = t1 + (1.0-k1)/slope
round_cycle_lenght = round(cycle_lenght,0)
print(round_cycle_lenght)