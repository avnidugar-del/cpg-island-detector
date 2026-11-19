
import streamlit as st
import subprocess
import sys
from pathlib import Path

st.set_page_config(
    page_title="CpG Island Detector",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 CpG Island Detector")
st.write(
    "Analyse DNA sequences and identify potential CpG islands."
)

st.subheader("Enter your DNA sequence")

sequence = st.text_area(
    "Paste your DNA sequence below:",
    height=180,
    placeholder="ATGCGCGCGATCG..."
)

if st.button("Analyse DNA Sequence"):
    if not sequence.strip():
        st.warning("Please enter a DNA sequence.")
    else:
        script = Path(__file__).with_name("cpg_detector.py")

        if not script.exists():
            st.error("Could not find cpg_detector.py")
        else:
            with st.spinner("Analysing your DNA sequence..."):
                try:
                    result = subprocess.run(
                        [
                            sys.executable,
                            str(script),
                            "-s",
                            sequence.strip()
                        ],
                        capture_output=True,
                        text=True,
                        timeout=120
                    )

                    if result.returncode == 0:
                        st.success("Analysis completed!")

                        if result.stdout.strip():
                            st.subheader("Analysis Results")
                            st.code(result.stdout)

                        if result.stderr.strip():
                            st.text(result.stderr)

                        if not result.stdout.strip():
                            st.info(
                                "The program completed but did not "
                                "print results. Check whether it "
                                "saves results to a file."
                            )
                    else:
                        st.error("The analysis encountered an error.")
                        st.code(result.stderr or result.stdout)

                except subprocess.TimeoutExpired:
                    st.error("The analysis took too long.")
                except Exception as e:
                    st.error(str(e))
