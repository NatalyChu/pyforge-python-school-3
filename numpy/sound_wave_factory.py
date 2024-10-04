import numpy as np
from scipy.io import wavfile

SAMPLING_RATE = 44100  # 44.1 KHz
DURATION_SECONDS = 5
# SOUND_ARRAY_LEN = SAMPLING_RATE * DURATION_SECONDS
MAX_AMPLITUDE = 2 ** 13

# Piano key frequencies
NOTES = {
    '0': 0, 'e0': 20.60172, 'f0': 21.82676, 'f#0': 23.12465, 'g0': 24.49971, 'g#0': 25.95654, 'a0': 27.50000, 'a#0': 29.13524,
    'b0': 30.86771, 'c0': 32.70320, 'c#0': 34.64783, 'd0': 36.70810, 'd#0': 38.89087,
    'e1': 41.20344, 'f1': 43.65353, 'f#1': 46.24930, 'g1': 48.99943, 'g#1': 51.91309, 'a1': 55.00000, 'a#1': 58.27047,
    'b1': 61.73541, 'c1': 65.40639, 'c#1': 69.29566, 'd1': 73.41619, 'd#1': 77.78175,
    'e2': 82.40689, 'f2': 87.30706, 'f#2': 92.49861, 'g2': 97.99886, 'g#2': 103.8262, 'a2': 110.0000, 'a#2': 116.5409,
    'b2': 123.4708, 'c2': 130.8128, 'c#2': 138.5913, 'd2': 146.8324, 'd#2': 155.5635,
    'e3': 164.8138, 'f3': 174.6141, 'f#3': 184.9972, 'g3': 195.9977, 'g#3': 207.6523, 'a3': 220.0000, 'a#3': 233.0819,
    'b3': 246.9417, 'c3': 261.6256, 'c#3': 277.1826, 'd3': 293.6648, 'd#3': 311.1270,
    'e4': 329.6276, 'f4': 349.2282, 'f#4': 369.9944, 'g4': 391.9954, 'g#4': 415.3047, 'a4': 440.0000, 'a#4': 466.1638,
    'b4': 493.8833, 'c4': 523.2511, 'c#4': 554.3653, 'd4': 587.3295, 'd#4': 622.2540,
    'e5': 659.2551, 'f5': 698.4565, 'f#5': 739.9888, 'g5': 783.9909, 'g#5': 830.6094, 'a5': 880.0000, 'a#5': 932.3275,
    'b5': 987.7666, 'c5': 1046.502, 'c#5': 1108.731, 'd5': 1174.659, 'd#5': 1244.508,
    'e6': 1318.510, 'f6': 1396.913, 'f#6': 1479.978, 'g6': 1567.982, 'g#6': 1661.219, 'a6': 1760.000, 'a#6': 1864.655,
    'b6': 1975.533, 'c6': 2093.005, 'c#6': 2217.461, 'd6': 2349.318, 'd#6': 2489.016,
    'e7': 2637.020, 'f7': 2793.826, 'f#7': 2959.955, 'g7': 3135.963, 'g#7': 3322.438, 'a7': 3520.000, 'a#7': 3729.310,
    'b7': 3951.066, 'c7': 4186.009, 'c#7': 4434.922, 'd7': 4698.636, 'd#7': 4978.032,
}

class SoundWaveFactory:
    def __init__(self, sampling_rate=SAMPLING_RATE, duration_seconds=DURATION_SECONDS, max_amplitude=MAX_AMPLITUDE):
        self.sampling_rate = sampling_rate
        self.duration_seconds = duration_seconds
        self.sound_array_len = sampling_rate * duration_seconds
        self.max_amplitude = max_amplitude
        self.common_timeline = np.linspace(0, duration_seconds, num=self.sound_array_len)

    def get_normed_sin(self, timeline, frequency):
        return self.max_amplitude * np.sin(2 * np.pi * frequency * timeline)

    def get_soundwave(self, timeline, note):
        return self.get_normed_sin(timeline, NOTES[note])

    # Generates and saves a sound wave for a given note.
    def create_note(self, note="a4", name=None):
        sound_wave = self.get_soundwave(self.common_timeline, note).astype(np.int16)
        if name is None:
            file_name = f"{note}_sin.wav".replace("#", "s")
        else:
            file_name = f"{name}.wav"
        wavfile.write(file_name, self.sampling_rate, sound_wave)
        return sound_wave

    # Reads wave data from a txt file
    def read_wave_from_txt(self, file_name):
        return np.loadtxt(file_name, dtype=np.int16)

    # Prints details important about the wave
    def print_wave_details(self, wave):
        print(f"Wave Shape: {wave.shape}")
        print(f"Max Amplitude: {np.max(wave)}, Min Amplitude: {np.min(wave)}")
        print(f"Wave Duration: {wave.shape[0] / self.sampling_rate} seconds")

    def normalize_sound_waves(self, waves):
        # Normalize length to the shortest wave
        min_len = min(wave.shape[0] for wave in waves)
        waves = [wave[:min_len] for wave in waves]
        # Normalize amplitude
        max_amplitude = max(np.max(np.abs(wave)) for wave in waves)
        normalized_waves = [(wave / max_amplitude * self.max_amplitude).astype(np.int16) for wave in waves]
        return normalized_waves

    def save_wave(self, wave, name="output", file_type='txt'):
        if file_type == 'WAV':
            wavfile.write(f"{name}.wav", self.sampling_rate, wave)
        else:
            np.savetxt(f"{name}.txt", wave)

    def create_triangle_wave(self, timeline, frequency):
        return self.max_amplitude * 2 * np.abs(2 * (timeline * frequency - np.floor(0.5 + timeline * frequency)))

    def create_square_wave(self, timeline, frequency):
        return self.max_amplitude * np.sign(np.sin(2 * np.pi * frequency * timeline))

    def apply_adsr(self, wave, attack, decay, sustain, release):
        # ADSR envelope modifications to the amplitude of the wave
        total_length = wave.shape[0]
        attack_len = int(total_length * attack)
        decay_len = int(total_length * decay)
        release_len = int(total_length * release)
        sustain_len = total_length - attack_len - decay_len - release_len

        attack_curve = np.linspace(0, 1, attack_len)
        decay_curve = np.linspace(1, sustain, decay_len)
        sustain_curve = np.full(sustain_len, sustain)
        release_curve = np.linspace(sustain, 0, release_len)

        envelope = np.concatenate([attack_curve, decay_curve, sustain_curve, release_curve])
        return (wave * envelope).astype(np.int16)

    def combine_waves(self, waves):
        return np.concatenate(waves)

    def read_melody(self, melody_text):
        melody = []
        for segment in melody_text.split():
            note = segment[:-1]
            duration = float(segment[-1])  # Assuming the last character is duration
            wave = self.get_soundwave(self.common_timeline[:int(duration * self.sampling_rate)], note)
            melody.append(wave)
        return self.combine_waves(melody)

if __name__ == "__main__":
    factory = SoundWaveFactory()
    wave_a4 = factory.create_note()
    factory.save_wave(wave_a4, 'a4', 'WAV')
