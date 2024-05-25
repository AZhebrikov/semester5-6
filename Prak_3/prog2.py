from scipy import signal
import matplotlib.pyplot as plt
import numpy as np


with open('ZhAV_EO') as f:
    lines = f.readlines()
lines.pop(0)
new_lines = list()
first_column = []
second_column = []
third_column = []
s = 1
for line in lines:
    l = line.split()
    first_column.append(np.float64(l[0]))
    second_column.append(np.float64(l[1]))
    third_column.append(np.float64(l[2]))
    s = 0
f, Pxx_den = signal.welch(np.array(first_column),
                          fs=50,
                          nfft=1024,
                          noverlap=524,
                          nperseg=1024,
                          # window='hamming',
                          # detrend=False
                          )
plt.plot(f, Pxx_den)
plt.xlim([0, 2]) # я ниче не настроил норм по иксу но вроде картинка
plt.ylim([1e-1, 15])  # менять для себя
plt.xlabel('F')
plt.ylabel('P')
# plt.yscale('linear')
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
plt.ylim([1e-1, 45])  # менять для себя
plt.xlabel('F')
plt.ylabel('P')
plt.yscale('linear')
plt.savefig('OE_saggital.png')

plt.cla()
third_column_new = []
for i in range(0, len(third_column), 7):
    third_column_new.append(third_column[i])
print(third_column)
print(third_column_new)


f, Pxx_den = signal.welch(np.array(third_column_new),
                          fs=50,
                          nfft=1024,
                          noverlap=234,
                          nperseg=468
                          )
plt.semilogy(f, Pxx_den)
plt.xlim([0, 12])
plt.ylim([1.5e-4, 0.11125])  # менять для себя
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
plt.ylim([0.5e-1, 50])
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
plt.ylim([1e-1, 100])  # менять для себя
plt.xlabel('F')
plt.ylabel('P')
plt.yscale('linear')
plt.savefig('CE_saggital.png')

plt.cla()
third_column_new = []
for i in range(0, len(third_column), 7):
    third_column_new.append(third_column[i])
f, Pxx_den = signal.welch(np.array(third_column_new),
                          fs=50,
                          nfft=1024,
                          noverlap=233,
                          nperseg=466
                          )
plt.semilogy(f, Pxx_den)
plt.xlim([0, 12])
plt.ylim([1e-4, 0.025])
plt.xlabel('F')
plt.ylabel('P')
plt.yscale('linear')
plt.savefig('CE_weight.png')