# *-* coding: utf-8 *-*

import wave
import os.path
import codecs
import matplotlib.pyplot as plt

import Spherical_Codebook
import coding_scheme

import struct

# VARIABLES:
Lv = 5
N_sp = 5
M_r = 8
p = 10
window_size = 320


# Create the codebook

codebook = Spherical_Codebook.S_Codebook(Lv, N_sp, M_r)
print 'Number of centroids: '+str(codebook.centroids_count)

audio_file_name = 'Audio_cav.wav'
audio_path = os.path.join('./audio', audio_file_name)
audio_output = os.path.join('./audio', 'output.wav')
audio_reader = wave.open(audio_path, 'r')
audio_writer = wave.open(audio_output, 'w')

sample_width = audio_reader.getsampwidth()
print sample_width
framerate = audio_reader.getframerate()
print framerate
frame_number = audio_reader.getnframes()
print frame_number

audio_writer.setnchannels(1)
audio_writer.setsampwidth(sample_width)
audio_writer.setframerate(framerate)

print 'Reading the frames'
frames = audio_reader.readframes(frame_number)
x = []
for b in frames:
    x.append(int(codecs.encode(b, 'hex'), 16)) #int.from_bytes(frames)
print 'x: ______________________________________________________'
print x


print 'Coding frames'
codewords_indxs, lpcs = coding_scheme.encode_simple(x, codebook, window_size, p)

print 'Decoding frames'
x_q = coding_scheme.decode(codewords_indxs, lpcs, codebook)

for i in x_q:
    packed = struct.pack('>B',i)
    audio_writer.writeframesraw(packed)


plt.plot(x,'r-o', label='Original sequence')
plt.plot(x_q,'g-o',label='Synthetic sequence')
plt.legend()
plt.title('Simple encoding simulation')
plt.show()
audio_writer.writeframes('')
audio_reader.close()
audio_writer.close()
print 'DONE :)'
