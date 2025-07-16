import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
from io import StringIO

# Konfigurasi halaman
st.set_page_config(
    page_title="ChemPathFinder",
    page_icon="🧪",
    layout="wide"
)

# CSS untuk styling
st.markdown("""
<style>
    .stSlider>div>div>div>div {
        background: #4CAF50;
    }
    .st-b7 {
        color: #FF5722;
    }
</style>
""", unsafe_allow_html=True)