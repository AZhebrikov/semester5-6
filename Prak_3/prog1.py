from scipy import signal
import matplotlib.pyplot as plt
import numpy as np
from math import log10

with open('ZhAV_EO') as f:
    lines = f.readlines()
lines.pop(0)
new_lines = list()
first_column = []
second_column = []
third_column = []
for line in lines:
    l = line.split()
    first_column.append(np.float64(l[0]))
    second_column.append(np.float64(l[1]))
    third_column.append(np.float64(l[2]))
f, Pxx_den = signal.welch(np.array(first_column),
                          fs=50,
                          nfft=1024,
                          noverlap=524,
                          nperseg=1024,
                          # window='hamming',
                          # detrend=False
                          )
plt.semilogy(f, Pxx_den)
plt.xlim([0, 2]) # я ниче не настроил норм по иксу но вроде картинка
plt.ylim([1e-1, 25])  # менять для себя
plt.xlabel('F')
plt.ylabel('P')
plt.yscale('linear')
# plt.gca()
plt.savefig('OE_frontal.png')

# plt.show()


plt.cla()
f, Pxx_den = signal.welch(np.array(second_column),
                          fs=50,
                          nfft=1024,
                          noverlap=524,
                          nperseg=1024
                          )
plt.semilogy(f, Pxx_den)
plt.xlim([0, 2])
plt.ylim([1e-1, 60])  # менять для себя
plt.xlabel('F')
plt.ylabel('P')
plt.yscale('linear')
plt.savefig('OE_saggital.png')

plt.cla()
f, Pxx_den = signal.welch(np.array(third_column),
                          fs=50/7,
                          nfft=1024,
                          noverlap=524,
                          nperseg=1024
                          )
plt.semilogy(f, Pxx_den)
plt.xlim([0, 12])
plt.ylim([1.5e-20, 0.0025])  # менять для себя
plt.xlabel('F')
plt.ylabel('P')
plt.yscale('linear')
plt.savefig('OE_weight.png')
plt.cla()


# тут все тоже самое для другого файла
with open('ZhAV_EC') as f:
    new_lines = f.readlines()
new_lines.pop(0)
# new_lines = list()
first_column = []
second_column = []
third_column = []
for line in new_lines:
    l = line.split()
    first_column.append(np.float64(l[0]))
    second_column.append(np.float64(l[1]))
    third_column.append(np.float64(l[2]))
f, Pxx_den = signal.welch(np.array(first_column),
                          fs=50,
                          nfft=1024,
                          noverlap=524,
                          nperseg=1024
                          )
plt.semilogy(f, Pxx_den)
plt.xlim([0, 2]) # я ниче не настроил норм по иксу но вроде картинка
plt.ylim([1.5e-20, 45])
plt.xlabel('F')
plt.ylabel('P')
plt.yscale('linear')
plt.savefig('CE_frontal.png')

# plt.show()


plt.cla()
f, Pxx_den = signal.welch(np.array(second_column),
                          fs=50,
                          nfft=1024,
                          noverlap=524,
                          nperseg=1024
                          )
plt.semilogy(f, Pxx_den)
plt.xlim([0, 2])
plt.ylim([1.5e-20, 75])  # менять для себя
plt.xlabel('F')
plt.ylabel('P')
plt.yscale('linear')
plt.savefig('CE_saggital.png')

plt.cla()
f, Pxx_den = signal.welch(np.array(third_column),
                          fs=50,
                          nfft=1024,
                          noverlap=524,
                          nperseg=1024
                          )
plt.semilogy(f, Pxx_den)
plt.xlim([0, 12])
plt.ylim([1.5e-20, 0.005])
plt.xlabel('F')
plt.ylabel('P')
plt.yscale('linear')
plt.savefig('CE_weight.png')