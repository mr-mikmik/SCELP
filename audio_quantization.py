# *-* coding: utf-8 *-*

import wave
import os.path
import codecs

import Spherical_Codebook
import coding_scheme

import struct

# VARIABLES:
Lv = 3
N_sp = 3
M_r = 3
p = 10
window_size = 320


# Create the codebook

codebook = Spherical_Codebook.S_Codebook(Lv, N_sp, M_r)

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

print 'Reading the frames'
frames = audio_reader.readframes(frame_number)
x = []
for b in frames:
    x.append(int(codecs.encode(b, 'hex'), 16)) #int.from_bytes(frames)

print 'Coding frames'
codewords_indxs, lpcs = coding_scheme.encode_simple(x, codebook, window_size, p)

print 'Decoding frames'
x_q = coding_scheme.decode(codewords_indxs, lpcs, codebook)
for i in x_q:
    packed = struct.pack(i)
    audio_writer.writeframesraw(packed)


audio_reader.close()
audio_writer.close()
print 'DONE :)'
