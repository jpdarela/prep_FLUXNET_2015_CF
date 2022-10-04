# test_mods.py
import matplotlib.pyplot as plt
import pandas as pd
from netCDF4 import Dataset

idx = pd.date_range("19890101", "20141231", freq="D")


# PR
pr = Dataset("pr_FLUXNET2015.nc")
pr_dec = Dataset("pr_SSP3-7.0_NearTerm-_FLUXNET2015.nc")

s1 = []
s2 = []
for j in range(12):
    s1.append(pd.Series(pr.variables['pr'][j, :]))
    s2.append(pd.Series(pr_dec.variables['pr'][j, :]))

plt.figure(figsize=(14,5))

# for x in range(12):
#     s1[x].rolling(365).mean().plot(color='b', alpha=0.5) 
#     s2[x].rolling(365).mean().plot(color='r', alpha=0.5)

s1[0].rolling(100).mean().plot(color='b', alpha=0.4) 
s2[0].rolling(100).mean().plot(color='r', alpha=0.4)


plt.savefig("test_pr.png", dpi=300)
plt.clf()
pr.close()
pr_dec.close()

# TAS
ts = Dataset("tas_FLUXNET2015.nc")
ts_dec = Dataset("tas_SSP3-7.0_NearTerm-_FLUXNET2015.nc")
s1 = []
s2 = []
for j in range(12):
    s1.append(pd.Series(ts.variables['tas'][j, :]))
    s2.append(pd.Series(ts_dec.variables['tas'][j, :]))

plt.figure(figsize=(14,5))

# for x in range(12):
#     s1[x].rolling(365).mean().plot(color='b', alpha=0.5) 
#     s2[x].rolling(365).mean().plot(color='r', alpha=0.5)

s1[0].rolling(100).mean().plot(color='b', alpha=0.4) 
s2[0].rolling(100).mean().plot(color='r', alpha=0.4)


plt.savefig("test_tas.png", dpi=300)
plt.clf()
ts.close()
ts_dec.close()