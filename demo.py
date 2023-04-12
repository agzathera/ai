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
    df = pd.read_csv("1eea_smiles.csv ")
    return df

url = 'https://passer.smu.edu/api'

def pocket_detection(protein):
    data = {"pdb": '3ahr', "chain": "A"}
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
    Generated small molecules for binding site
    """)

    # Copy the dataset so any changes are not applied to the original cached version
    df = pd.read_csv("./SMILES1.csv")

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

    raw_html = mols2grid.display(df_result4,
                                 # subset=["Name", "img"],
                                 subset=["img"],
                                 mapping={"smiles": "SMILES", "MW": "MW", "LogP": "LogP", "NumHDonors": "NumHDonors", "NumHAcceptors": "NumHAcceptors"},
                                 tooltip = ["MW", "LogP", "NumHDonors", "NumHAcceptors", "SMILES"])._repr_html_()

    components.html(raw_html, width=900, height=1100, scrolling=True)


def main():
    st.sidebar.title('AlloKey.ai')
    prot_str="1EEA"
    prot_list=prot_str.split(',')
    bcolor = st.sidebar.color_picker('Pick A Color', '#00f900')
    protein=st.sidebar.selectbox('select protein',prot_list)
    style = st.sidebar.selectbox('style', ['cartoon', 'line', 'cross', 'stick', 'sphere', 'clicksphere'])
#     residues = pocket_detection(protein)
    residues = [70, 72, 84, 85, 117, 118, 119, 121, 122, 123, 130, 199, 200, 287, 288, 290, 330, 331, 334, 335, 440, 444]
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
