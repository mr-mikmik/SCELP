# *-* coding: utf-8 *-*

import wave
import os.path
import codecs

audio_file_name = 'Audio_cav.wav'
audio_path = os.path.join('./audio', audio_file_name)

audio_reader = wave.open(audio_path, 'r')

sample_width = audio_reader.getsampwidth()
print sample_width
framerate = audio_reader.getframerate()
print framerate
frame_number = audio_reader.getnframes()
print frame_number

for i in range(1000):
    frames = audio_reader.readframes(20)
    b = int(codecs.encode(frames, 'hex'), 16)#int.from_bytes(frames)
    print b

audio_reader.close()


