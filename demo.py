import streamlit as st
from stmol import showmol
import py3Dmol

import mols2grid
import pandas as pd
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors
import requests

@st.cache(allow_output_mutation=True)
def download_dataset():
    """Loads once then cached for subsequent runs"""
    df = pd.read_csv(
        "https://raw.githubusercontent.com/aspuru-guzik-group/chemical_vae/master/models/zinc_properties/250k_rndm_zinc_drugs_clean_3.csv"
    ).dropna()
    return df[['smiles']].sample(50)

url = 'https://passer.smu.edu/api'

def pocket_detection(protein):
    data = {"pdb": protein, "chain": "A"}
    results = requests.post(url, data=data)
    pocket_residues = results.json()["1"]["residues"].split(" ")[4:]
    pocket_residues = [eval(i) for i in pocket_residues]
    return pocket_residues

# Calculate descriptors
def calc_mw(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the molecular weight"""
    mol = Chem.MolFromSmiles(smiles_string)
    return ExactMolWt(mol)


def calc_logp(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the LogP"""
    mol = Chem.MolFromSmiles(smiles_string)
    return MolLogP(mol)


def calc_NumHDonors(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the NumHDonors"""
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHDonors(mol)


def calc_NumHAcceptors(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the NumHAcceptors"""
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHAcceptors(mol)

def generate_ligands():

    st.markdown("""
    Generated small molecules for allosteric binding site
    """)

    # Copy the dataset so any changes are not applied to the original cached version
    df = download_dataset().copy()
    df["MW"] = df.apply(lambda x: calc_mw(x["smiles"]), axis=1)
    df["LogP"] = df.apply(lambda x: calc_logp(x["smiles"]), axis=1)
    df["NumHDonors"] = df.apply(lambda x: calc_NumHDonors(x["smiles"]), axis=1)
    df["NumHAcceptors"] = df.apply(lambda x: calc_NumHAcceptors(x["smiles"]), axis=1)

    # Sidebar panel
    st.sidebar.header('Set parameters')
    st.sidebar.write('*Note: Display compounds having values less than the following thresholds*')
    weight_cutoff = st.sidebar.slider(
        label="Molecular weight",
        min_value=0,
        max_value=1000,
        value=500,
        step=10,
    )
    logp_cutoff = st.sidebar.slider(
        label="LogP",
        min_value=-10,
        max_value=10,
        value=5,
        step=1,
    )
    NumHDonors_cutoff = st.sidebar.slider(
        label="NumHDonors",
        min_value=0,
        max_value=15,
        value=5,
        step=1,
    )
    NumHAcceptors_cutoff = st.sidebar.slider(
        label="NumHAcceptors",
        min_value=0,
        max_value=20,
        value=10,
        step=1,
    )

    df_result = df[df["MW"] < weight_cutoff]
    df_result2 = df_result[df_result["LogP"] < logp_cutoff]
    df_result3 = df_result2[df_result2["NumHDonors"] < NumHDonors_cutoff]
    df_result4 = df_result3[df_result3["NumHAcceptors"] < NumHAcceptors_cutoff]

    st.write(df_result4.shape)
    st.write(df_result4)

    raw_html = mols2grid.display(df_result4,
                                 # subset=["Name", "img"],
                                 subset=["img", "MW", "LogP", "NumHDonors", "NumHAcceptors"],
                                 mapping={"smiles": "SMILES"})._repr_html_()
    components.html(raw_html, width=900, height=1100, scrolling=False)


def main():
    st.sidebar.title('AlloKey.ai')
    prot_str='4JSV,1BML,1D5M,1D5X,1D5Z,1D6E,1DEE,1E9F,1FC2,1FCC,1G4U,1GZS,1HE1,1HEZ,1HQR,1HXY,1IBX,1JBU,1JWM,1JWS'
    prot_list=prot_str.split(',')
    bcolor = st.sidebar.color_picker('Pick A Color', '#00f900')
    protein=st.sidebar.selectbox('select protein',prot_list)
    style = st.sidebar.selectbox('style', ['cartoon', 'cross', 'stick', 'sphere', 'line', 'clicksphere'])
    residues = pocket_detection(protein)
    xyzview = py3Dmol.view(query='pdb:'+protein)
    xyzview.setStyle({style:{'color':'spectrum'}})
    xyzview.addStyle({'within': {'distance': 3, 'sel': {'chain': "A", 'resi': residues}}}, {'sphere': {'color': 'red'}})
    xyzview.setBackgroundColor(bcolor)
    xyzview.spin(True)
    button = st.sidebar.button('run')

    if button:
        showmol(xyzview, height=500, width=800)
        generate_ligands()


if __name__ == "__main__":
    main()
