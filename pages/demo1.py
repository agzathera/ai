import requests
import streamlit as st

url = 'https://passer.smu.edu/api'

def pocket_detection(protein):
    data = {"pdb": protein, "chain": "A"}
    results = requests.post(url, data=data)
    pocket_residues = results.json()["1"]["residues"].split(" ")[4:]
    pocket_residues = [eval(i) for i in pocket_residues]
    return pocket_residues

def main():
    st.sidebar.title('AlloKey.ai')
    prot_str='1A2C,1BML,1D5M,1D5X,1D5Z,1D6E,1DEE,1E9F,1FC2,1FCC,1G4U,1GZS,1HE1,1HEZ,1HQR,1HXY,1IBX,1JBU,1JWM,1JWS'
    prot_list=prot_str.split(',')
    protein=st.sidebar.selectbox('select protein',prot_list)
    button = st.sidebar.button('run')

    if button:
        st.text(pocket_detection(protein))



if __name__ == "__main__":
    main()