pip freezeerfgf"""
Ripple (.nf3/.nf6/.nev) -> EDF+ converter, EMG muscle channels only.
=====================================================================

Why this exists
----------------
`.ns2`/`.ns5` are standard Blackrock "NEURALCD" files: neo.BlackrockIO /
brpylib read them fine, but on this rig they only contain the NANO
stimulation-electrode traces (channels labelled "raw 129".."raw 160" /
"lfp 129".."lfp 160") plus 30 auxiliary analog inputs -- NOT the EMG.

The real, correctly-labelled EMG (e.g. "RGLU hifreq", "LTA hifreq", ...)
lives in `.nf3` (2 kHz, hardware "hi-res" stream) and `.nf6` (7.5 kHz,
hardware "hifreq" stream). Both declare the file-type ID "NEUCDFLT",
which is a Ripple-specific extension that neither `neo` nor `brpylib`
recognizes -- their NSx readers refuse anything that isn't "NEURALCD"/
"NEURALSG", so these files get silently skipped by both libraries.

Byte-for-byte, `.nf3`/`.nf6` use the *same* basic-header and per-channel
extended-header layout as standard "NEURALCD" 2.2 files (confirmed by
hand-parsing real files: the recovered high/low filter corner
frequencies match the device's own JSON config exactly, e.g. 30-950 Hz).
The only real differences are:
  1. the 8-byte magic string ("NEUCDFLT" instead of "NEURALCD"), and
  2. samples are stored as 4-byte float32 (already-scaled physical
     values), not 2-byte int16 raw ADC codes.
  3. the per-channel digital_min/digital_max header fields are not
     populated meaningfully for these derived streams (a Ripple
     export quirk, not a parsing bug) -- physical min/max for the EDF
     header are computed from the actual data instead, same as the
     original script already did.

So this script reads .nf3/.nf6 itself (no neo/brpylib dependency for
that part) and only calls `neo` for the `.nev` file (a standard
Blackrock format that neo handles natively) to pull event timestamps
for annotations / _events.tsv.

Output: one BIDS-ish EDF+ per run, containing ONLY the 32 EMG muscle
channels (no NANO stim-electrode channels, no auxiliary analog
channels) -- plus the matching _channels.tsv, _events.tsv and .json
sidecars.
"""

import os
import re
import json
import struct
import csv
import numpy as np
import pyedflib
from neo.io import BlackrockIO

# --------------------------------------------------------------------
# CONFIG -- edit these for your setup
# --------------------------------------------------------------------
RAW_DIR = "D:/Data_Pour_Arnaud/DATA_ARKEMA/AMP_2026-06-17/EMG/" #PATH OF RAW DATA
EDF_DIR = "D:/Data_Pour_Arnaud/DATA_ARKEMA/AMP_2026-06-17/edf/" #PATH WEHRE .EDF AND .JSON ARE SAED
SUBJECT = "sub-01"          # BIDS subject label
STREAM = "nf6"              # "nf6" (7.5 kHz, hifreq) or "nf3" (2 kHz, hi-res)
POWERLINE_HZ = 60           # 60 Hz in Canada

# Matches Ripple's native naming, e.g.:
#   ARKEMA_AMP_2026-06-17_(1-circ)#(2-circ)_0#1_0001.nf6
# and also tolerates a sanitized variant (parentheses/# replaced with "_"),
# e.g. from files re-exported through a filesystem that strips those chars:
#   ARKEMA_AMP_2026-06-17__1-circ___2-circ__0_1_0001.nf6
# Adjust this if your naming convention differs from both.
FNAME_RE = re.compile(
    r'^(?P<subject>[A-Za-z0-9]+)_(?P<task>[A-Za-z]+)_(?P<date>\d{4}-\d{2}-\d{2})[_(]+'
    r'(?P<el1>\d+)-circ[_)#(]+(?P<el2>\d+)-circ[_)#]+'
    r'(?P<p1>\d)[_#]+(?P<p2>\d)_(?P<run>\d+)$'
)


# --------------------------------------------------------------------
# Reader for Ripple's "NEUCDFLT" filtered EMG streams (.nf3 / .nf6)
# --------------------------------------------------------------------
def read_filtered_nsx(path):
    """Read a Ripple .nf3/.nf6 file. Returns fs, channel labels, and a
    (n_samples, n_channels) float32 array of already-filtered EMG data
    in the units declared by the header (normally uV)."""
    with open(path, 'rb') as f:
        magic = f.read(8).decode('ascii', errors='replace').strip('\x00')
        if magic not in ('NEUCDFLT', 'NEURALCD'):
            raise ValueError(f"Unexpected file signature {magic!r} in {path}")
        f.read(2)  # file spec major/minor -- not needed
        bytes_in_headers = struct.unpack('<I', f.read(4))[0]
        f.read(16)   # label
        f.read(256)  # comment
        period = struct.unpack('<I', f.read(4))[0]
        timeres = struct.unpack('<I', f.read(4))[0]
        fs = timeres / period
        f.read(16)   # time origin
        channel_count = struct.unpack('<I', f.read(4))[0]

        labels, units, hp_hz, lp_hz = [], [], [], []
        for _ in range(channel_count):
            hdr = f.read(66)
            label_ch = hdr[4:20].decode('ascii', errors='replace').strip('\x00')
            unit = hdr[30:46].decode('ascii', errors='replace').strip('\x00')
            hi_fc = struct.unpack('<I', hdr[46:50])[0] / 1000.0
            lo_fc = struct.unpack('<I', hdr[56:60])[0] / 1000.0
            labels.append(label_ch)
            units.append(unit)
            hp_hz.append(hi_fc)
            lp_hz.append(lo_fc)

        f.seek(bytes_in_headers)
        marker = f.read(1)
        if marker != b'\x01':
            raise ValueError(f"Unexpected data marker {marker!r} in {path}")
        f.read(4)  # timestamp (uint32) -- assumes one continuous segment
        npoints = struct.unpack('<I', f.read(4))[0]

        raw = np.fromfile(f, dtype='<f4', count=npoints * channel_count)
        if raw.size != npoints * channel_count:
            raise ValueError(
                f"{path}: expected {npoints * channel_count} samples, "
                f"got {raw.size} -- file may be truncated or have a second "
                f"data segment (paused recording) not handled by this reader."
            )
        data = raw.reshape(npoints, channel_count)

    return dict(fs=fs, labels=labels, units=units, hp_hz=hp_hz, lp_hz=lp_hz, data=data)


def clean_muscle_label(raw_label):
    """'RGLU hifreq' / 'LTA hi-res' -> 'RGLU' / 'LTA'"""
    return raw_label.replace('hifreq', '').replace('hi-res', '').strip()


# --------------------------------------------------------------------
# Best-effort event extraction from .nev (standard Blackrock format,
# so neo can read it directly -- only the nf3/nf6 EMG needed a custom
# reader).
# --------------------------------------------------------------------
def read_nev_events(nev_path, t_start_offset=0.0):
    """Returns a list of (onset_seconds, label) tuples. Never raises --
    prints a warning and returns [] if anything goes wrong, so a missing
    or unreadable .nev never blocks the EMG conversion."""
    events = []
    try:
        reader = BlackrockIO(filename=nev_path.rsplit('.nev', 1)[0], load_nev=True)
        reader.parse_header()
        ev_dict = reader.get_event_timestamps()
        if ev_dict and len(ev_dict[0]) > 0:
            fs_events = getattr(reader, 'sample_frequency_event', 30000.0)
            for ts, label in zip(ev_dict[0], ev_dict[2]):
                t = float(ts) / fs_events - t_start_offset
                events.append((t, str(label)))
    except Exception as e:
        print(f"  [events] Could not read {nev_path}: {e}")
    return events


# --------------------------------------------------------------------
# Main conversion for a single run
# --------------------------------------------------------------------
def convert_run(nfx_path, nev_path, edf_dir, subject=SUBJECT):
    base = os.path.basename(nfx_path)
    stem = os.path.splitext(base)[0]
    m = FNAME_RE.match(stem)
    if not m:
        raise ValueError(
            f"Filename '{stem}' doesn't match the expected pattern -- "
            f"edit FNAME_RE to match your naming convention."
        )
    info = m.groupdict()
    session = "ses-" + info['date'].replace('-', '')
    task = "task-" + info['task']
    run = f"run-{int(info['run']):03d}"
    acq = f"acq-circ{info['el1']}{info['el2']}set{info['p1']}"

    print(f"Converting {base}  ->  {subject}_{session}_{task}_{acq}_{run}_emg.edf")

    result = read_filtered_nsx(nfx_path)
    fs = result['fs']
    data = result['data']            # (n_samples, n_channels)
    n_samples, n_channels = data.shape

    channel_headers = []
    channel_signals = []
    for ch in range(n_channels):
        sig = data[:, ch]
        label = clean_muscle_label(result['labels'][ch])

        p_max = float(np.max(sig)) if len(sig) else 1000.0
        p_min = float(np.min(sig)) if len(sig) else -1000.0
        if p_max == p_min:
            p_max += 1.0
            p_min -= 1.0

        channel_headers.append({
            'label': label,
            'dimension': result['units'][ch] or 'uV',
            'sample_frequency': fs,
            'physical_max': p_max,
            'physical_min': p_min,
            'digital_max': 32767,
            'digital_min': -32768,
            'transducer': 'Ripple System (hardware-filtered EMG)',
            'prefilter': f"HP:{result['hp_hz'][ch]:.0f}Hz|LP:{result['lp_hz'][ch]:.0f}Hz|Notch:{POWERLINE_HZ}Hz",
        })
        channel_signals.append(sig)

    os.makedirs(edf_dir, exist_ok=True)
    base_out = f"{subject}_{session}_{task}_{acq}_{run}"
    edf_path = os.path.join(edf_dir, base_out + "_emg.edf")
    json_path = os.path.join(edf_dir, base_out + "_emg.json")
    channels_tsv_path = os.path.join(edf_dir, base_out + "_channels.tsv")
    events_tsv_path = os.path.join(edf_dir, base_out + "_events.tsv")

    duration_sec = n_samples / fs

    # --- events (best effort) ---
    events = read_nev_events(nev_path)
    events = [(t, label) for t, label in events if 0 <= t <= duration_sec]
    events.sort(key=lambda x: x[0])

    # --- EDF ---
    with pyedflib.EdfWriter(edf_path, n_channels, file_type=pyedflib.FILETYPE_EDFPLUS) as w:
        w.setSignalHeaders(channel_headers)
        w.writeSamples(channel_signals)
        for t, label in events:
            w.writeAnnotation(t, -1, label)

    # --- BIDS json sidecar ---
    bids_json = {
        "TaskName": info['task'],
        "SamplingFrequency": float(fs),
        "EMGChannelCount": int(n_channels),
        "RecordingDuration": round(duration_sec, 3),
        "PowerLineFrequency": POWERLINE_HZ,
        "Manufacturer": "Ripple Neuro",
        "SoftwareVersions": "custom NEUCDFLT reader / pyedflib",
        "RecordingType": "continuous",
        "SourceStream": STREAM,
    }
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(bids_json, f, indent=4, ensure_ascii=False)

    # --- channels.tsv ---
    with open(channels_tsv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['name', 'type', 'units', 'sampling_frequency', 'description'])
        for ch in channel_headers:
            writer.writerow([ch['label'], 'EMG', ch['dimension'], float(fs), 'Surface EMG channel'])

    # --- events.tsv (now actually populated) ---
    with open(events_tsv_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['onset', 'duration', 'trial_type'])
        for t, label in events:
            writer.writerow([round(t, 6), 'n/a', label])

    print(f"  -> {edf_path}  ({n_channels} EMG channels, {duration_sec:.1f} s, {len(events)} events)")
    return edf_path


# --------------------------------------------------------------------
# Batch over every matching stream file in RAW_DIR
# --------------------------------------------------------------------
if __name__ == "__main__":
    ext = "." + STREAM
    runs = sorted(f for f in os.listdir(RAW_DIR) if f.endswith(ext))
    if not runs:
        print(f"No .{STREAM} files found in {RAW_DIR}")
    for fname in runs:
        nfx_path = os.path.join(RAW_DIR, fname)
        nev_path = os.path.join(RAW_DIR, os.path.splitext(fname)[0] + ".nev")
        convert_run(nfx_path, nev_path, EDF_DIR)
