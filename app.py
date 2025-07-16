import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
from io import StringIO
st.title('🪄 Reaction Pathway Predictor')
st.subheader('Prediksi Jalur Reaksi Kimia Berdasarkan Parameter Analisis')
# Sidebar untuk parameter input
with st.sidebar:
    st.header('Parameter Reaksi')
    temperature = st.slider('Suhu (°C)', 0, 200, 25)
    ph = st.slider('pH', 0, 14, 7)
    concentration = st.slider('Konsentrasi (mol/L)', 0.1, 5.0, 1.0)
    solvent = st.selectbox('Pelarut', ['Air', 'Etanol', 'Aseton', 'Dietil Eter', 'DMSO'])
    catalyst = st.checkbox('Ada Katalis?')
    if catalyst:
        catalyst_type = st.selectbox('Jenis Katalis', [
            'Asam', 'Basa', 'Logam Transisi', 'Enzim'
        ])
# Input senyawa awal
st.header('Masukkan Reaktan')
col1, col2 = st.columns(2)
with col1:
    reactant1 = st.text_input('Reaktan 1 (SMILES)', 'C=O')
with col2:
    reactant2 = st.text_input('Reaktan 2 (SMILES)', 'O')
# Fungsi untuk validasi SMILES
def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return True
# Fungsi untuk memprediksi jalur reaksi (simulasi)
def predict_reaction_pathway(r1, r2, temp, ph, conc, solvent, catalyst=None):
    # Ini adalah simulasi - dalam aplikasi nyata akan menggunakan library kimia
    pathways = [] 
    # Contoh sederhana untuk reaksi aldehid + air
    if r1 == 'C=O' and r2 == 'O':
        pathways.append({
            'name': 'Hidrasi Aldehid',
            'steps': [
                ('Pembentukan intermediate tetrahedral', 50),
                ('Proton transfer', 30),
                ('Pembentukan produk hidrat', 20)
            ],
            'products': ['OC(O)C'],
            'activation_energy': 75 - (ph * 2) - (temp * 0.1),
            'thermodynamics': 'Eksotermik'
        })
    # Contoh lain untuk reaksi esterifikasi
    elif r1 == 'CC(=O)O' and r2 == 'CO':
        pathways.append({
            'name': 'Esterifikasi',
            'steps': [
                ('Protonasi karbonil', 60),
                ('Serangan nukleofilik', 40),
                ('Dehidrasi', 30)
            ],
            'products': ['CC(=O)OC', 'O'],
            'activation_energy': 90 - (ph * 3) - (temp * 0.2),
            'thermodynamics': 'Endotermik'
        })
    else:
        # Default pathway untuk senyawa acak
        pathways.append({
            'name': 'Reaksi Adisi',
            'steps': [
                ('Inisiasi', 40),
                ('Propagasi', 30),
                ('Terminasi', 20)
            ],
            'products': [f'{r1}{r2}'],
            'activation_energy': 85 - (ph * 1.5) - (temp * 0.15),
            'thermodynamics': 'Eksotermik'
        })
    return pathways
# Tombol prediksi
if st.button('Prediksi Jalur Reaksi'):
    if not validate_smiles(reactant1) or not validate_smiles(reactant2):
        st.error('Format SMILES tidak valid untuk satu atau lebih reaktan!')
    else:
        with st.spinner('Menganalisis kemungkinan jalur reaksi...'):
            # Dapatkan prediksi
            pathways = predict_reaction_pathway(
                reactant1, reactant2, temperature, ph, 
                concentration, solvent, catalyst if catalyst else None
            )            
            # Tampilkan hasil
            st.success('Prediksi berhasil!')
            for pathway in pathways:
                st.subheader(f"Jalur: {pathway['name']}")
                # Diagram energi
                st.markdown("### Diagram Energi Potensial")
                fig, ax = plt.subplots()
                steps = [step[0] for step in pathway['steps']]
                energies = [step[1] for step in pathway['steps']]
                # Tambahkan reaktan dan produk
                steps = ['Reaktan'] + steps + ['Produk']
                energies = [0] + energies + [energies[-1] - 20]
                ax.plot(steps, energies, marker='o')
                ax.set_ylabel('Energi (kJ/mol)')
                ax.set_xlabel('Tahapan Reaksi')
                ax.set_title('Profil Energi Reaksi')
                plt.xticks(rotation=45)
                plt.tight_layout()
                st.pyplot(fig)
                # Tampilkan parameter kunci
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Energi Aktivasi", f"{pathway['activation_energy']:.1f} kJ/mol")
                with col2:
                    st.metric("Termodinamika", pathway['thermodynamics'])
                # Tampilkan tahapan
                st.markdown("### Tahapan Reaksi")
                for i, step in enumerate(pathway['steps'], 1):
                    st.markdown(f"{i}. **{step[0]}** - Energi: {step[1]} kJ/mol")
                
                # Tampilkan produk prediksi
                st.markdown("### Produk Prediksi")
                products = pathway['products']
                
                mols = []
                legends = []
                for p in products:
                    mol = Chem.MolFromSmiles(p)
                    if mol:
                        mols.append(mol)
                        legends.append(p)
                
                if mols:
                    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(300,300), legends=legends)
                    st.image(img)
                
                st.divider()

# Bagian analisis tambahan
st.header('Analisis Lanjutan')
with st.expander('Optimasi Kondisi Reaksi'):
    st.write("""
    Berdasarkan parameter yang dimasukkan, berikut rekomendasi optimasi:
    """)
    
    data = {
        'Parameter': ['Suhu', 'pH', 'Konsentrasi', 'Pelarut'],
        'Nilai Saat Ini': [f"{temperature}°C", ph, f"{concentration} M", solvent],
        'Rekomendasi': [
            f"{temperature + 10}°C" if temperature < 100 else "Cukup",
            "Asam (pH 3-6)" if ph > 6 else "Basa (pH 8-11)" if ph < 8 else "Netral (baik)",
            f"{concentration * 1.2:.1f} M" if concentration < 3 else "Cukup",
            "Pertahankan" if solvent in ['Air', 'Etanol'] else "Coba ganti ke Air"]
    }
    
    st.table(pd.DataFrame(data))

with st.expander('Simpan Hasil Analisis'):
    st.download_button(
        label="Unduh Laporan PDF",
        data="Ini adalah simulasi laporan PDF. Dalam aplikasi nyata akan berisi data aktual.",
        file_name="reaction_analysis.pdf",
        mime="application/pdf"
    )

# Catatan kaki
st.caption("""
Aplikasi ini menggunakan model prediksi sederhana untuk tujuan demonstrasi. 
Hasil aktual mungkin berbeda tergantung kondisi eksperimen yang sebenarnya.
""")
