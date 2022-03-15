import scipy.interpolate
import numpy as np, matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from reconstruction import CreateHSP
import math

nn =64
nn2 = 2*nn
x= np.array([0, 0.25, 0.5, 0.75, 1])
# y= np.array([1, 0.815, 0.4564, 0.1636, 0])
# y= np.array([1, 1.0485, 1.17, 1.2202, 0.9201])
y= np.array([1, 0.9338, 0.7441, 0.4425, 0.0531])
# xx = np.linspace(x.min(), x.max(), 100)
# fig, ax = plt.subplots(figsize=(8, 4))
# ax.scatter(x, y)

f = interp1d(x, y, kind = 'quadratic')
# ax.plot(xx, f(xx), label= 'quadratic')
# ax.legend()
# ax.set_ylabel(r"$y$", fontsize=18)
# ax.set_xlabel(r"$x$", fontsize=18)

FFT_F = np.zeros(nn2)
# for i in range(nn):
#     FFT_F[i] = f((nn-i-0.003)/nn)
#     FFT_F[nn2 - i - 1] = f((nn-i-1+0.003)/nn)
for i in range(nn):
    FFT_F[i] = f((i)/nn)*0.997 * (i+0.003) / nn2
    FFT_F[nn2 - i - 1] = f((i)/nn)* 0.997*(i + 1 + 0.003) / nn2
# FFT_F= FFT_F * complex(0,1)


# WindowType = 2
# print(WindowType)
# # nn = int(math.pow(2, (math.ceil(math.log2(abs(YL))) + 1)))
# nn =64
# HS = CreateHSP(nn, WindowType)
#
# nn2 = nn * 2
# k = int(nn / 2)
# TempF = np.zeros(nn2)
# TempF[0:k] = HS[k:nn]
# TempF[k + nn:nn2] = HS[0:k]
# HS = TempF
# P = HS
# FFT_F = np.fft.fft(HS)


plt.plot(FFT_F)
plt.show()
# nn=64
# nn2 = 128
# FFT_F = np.zeros(nn2)
# for i in range(nn):
#     FFT_F[i] = f((nn - i) / nn) * 0.997 * (i + 0.003) / nn2
#     FFT_F[nn2 - i - 1] = f((nn - i - 1 + 0.003) / nn) * 0.997 * (i + 1 + 0.003) / nn2
#
# plt.plot(FFT_F)
# plt.show()
# Length = 8
# YL = 2*Length +1
# # Hk = np.ones(2*Length+1)
# Center = int((Length) / 2)
# PI = 3.14159265358979
#
# for i in range(Length+1):
#     # HS[i] = -2 / (PI * PI * (4 * (i - Center) * (i - Center) - 1))
#     Hk[i+Length] = i
#     Hk[Length - i] = i

# WindowType = 1
# print(WindowType)
# nn = int(math.pow(2, (math.ceil(math.log2(abs(YL))) + 1)))
# HS = CreateHSP(nn, WindowType)

# nn2 = nn * 2
# k = int(nn / 2)
# TempF = np.zeros(nn2)
# TempF[0:k] = HS[k:nn]
# TempF[k + nn:nn2] = HS[0:k]
# HS = TempF
# P = HS
# FFT_F = np.fft.fft(HS)
# print(nn2)
# Hk = np.ones(nn2)
# Ht = np.ones(nn2)

# Ht[0] = 0
# Ht[nn] = 0.25
# for i in range(1, nn):
#     Ht[i] = -np.power((math.sin(PI * (i - nn) / 2)), 2) / (PI * PI * (i - nn) * (i - nn))
#
# for k in range((nn + 1), nn2):
#     Ht[k] = -np.power((math.sin(PI * (k - nn) / 2)), 2) / (PI * PI * (k - nn) * (k - nn))

# for i in range(nn):
#     Hk[i] =  i/nn2
#     Hk[nn2 - i-1] = (i+1)/nn2
#     # Hk[i] = PI * math.sin(i*PI/nn2)
#     # Hk[nn2 - i-1] = PI * math.sin(i*PI/nn2)
#
# Hkif = np.fft.ifft(Hk, nn2)
# Hkif = Hkif.real

# print(ktf.real[60:68])
# print(np.shape(Ht))
# plt.plot(FFT_F)
# plt.show()