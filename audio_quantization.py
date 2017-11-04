import wave
import os.path

audio_file_name = 'Audio_cav.wav'
audio_path = os.path.join('./audio', audio_file_name)

with wave.open(audio_path, 'r') as audio_reader:
    sample_width = audio_reader.getsampwidth()
    print sample_width

    audio_reader.close()


