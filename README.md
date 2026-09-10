# ripple_conversion_raw2edf
# A short script aiming to recover the EMG data from ripples Nano2+ RAW files (nf3 and nf6) and rewritte them into a comprehensible .edf scaled to 7.5 kHz

1. Clone the repository:
   ```bash
   git clone https://github.com/OBernardReMo/ripple_conversion_raw2edf.git
   cd ripple-conversion-raw2edf
   ```
 
2. Install dependencies:
   ```bash
   uv sync
   ```
 
   <details>
   <summary>Alternative: pip</summary>
 
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   ```
   </details>
