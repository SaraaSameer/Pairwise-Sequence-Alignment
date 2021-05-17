import matplotlib.pyplot as plt
import numpy as np


Methods=['Test_Case01','Test_Case02','Test_Case03']      #Lists
SeriesTime=[0.6,0.8,0.9]
PThreadsTime=[0.4,0.3,0.66]
OpenMPTime=[0.3,0.3,0.2]

xAxis=np.arange(len(Methods))   #Converting the Methods into array of numbers 0 , 1 , 2
plt.xticks(xAxis,Methods)   #Replacing 0 1 2 with Methods strings

plt.title("Sequence Alignment Problem")
plt.ylabel("Time")

plt.bar(xAxis-0.2,SeriesTime,width=0.2,color = "blue",label='Series')
plt.bar(xAxis,PThreadsTime,width=0.2,color = "red",label='PThreads')
plt.bar(xAxis+0.2,OpenMPTime,width=0.2,color = "green",label='OpenMP')
plt.legend()
plt.show()