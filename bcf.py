#==============================================================================
#  Copyright (c) 2014-2018 Joel de Guzman. All rights reserved.
#
#  Distributed under the Boost Software License, Version 1.0. (See accompanying
#  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#==============================================================================
from matplotlib.pyplot import figure, show
from numpy import arange, sin, pi

t = arange(0.0, 2.0, 0.01)
input = 0.3*sin(2*pi*t) + 0.4*sin(4*pi*t) + 0.3*sin(6*pi*t)

fig = figure(1)

ax1 = fig.add_subplot(311)
ax1.plot(t, input)
ax1.grid(True)
ax1.set_ylim((-1, 1))

class zero_cross:
   def __init__(self):
      self.y = 0

   def __call__(self, s):
      if s < -0.1:
         self.y = 0
      elif s > 0.1:
         self.y = 1
      return self.y

zc = zero_cross()
trig = [zc(s) for s in input]

ax2 = fig.add_subplot(312)
ax2.plot(t, trig)
ax2.grid(True)
ax2.set_ylim((-0.1, 1.1))

def count_ones(l):
   r = 0
   for e in l:
      if e:
         r += 1
   return r

cross = trig
results = []
len = len(trig)/2
for i in range(len):
   x = [a ^ b for a, b in zip(trig[0:len], cross[i:i+len])]
   results.append(count_ones(x))
results.extend(results)

ax3 = fig.add_subplot(313)
ax3.plot(results)
ax3.grid(True)
ax3.set_ylim((-5, 80))

show()



