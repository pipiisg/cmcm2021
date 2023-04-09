import numpy as np
from collections.abc import Iterable

def toArray(str):
    l = str.split("\n")
    print(list(map(float, l)))

toArray("""270
264
267
265
270
260
250
240
230
226
211
196
195
189
180
172
165
152
145
152
126
130
112
98
94
83""")