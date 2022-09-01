import streamlit as st
from streamlit_option_menu import option_menu

import base64


def body_bg(side_bg):
    '''
    A function to unpack an image from url and set as bg.
    Returns
    -------
    The background.
    '''
    side_bg_ext = 'png'

    st.markdown(
        f"""
         <style>
         .stApp {{
            background: url(data:image/{side_bg_ext};base64,{base64.b64encode(open(side_bg, "rb").read()).decode()});
            background-size: cover
         }}
         </style>
         """,
        unsafe_allow_html=True
    )


def sidebar_bg(side_bg):

   side_bg_ext = 'png'

   st.markdown(
      f"""
      <style>
      [data-testid="stBody"] > div:first-child {{
          background: url(data:image/{side_bg_ext};base64,{base64.b64encode(open(side_bg, "rb").read()).decode()});
      }}
      </style>
      """,
      unsafe_allow_html=True,
      )

body_bg('./images/background2.jpeg')
sidebar_bg('./images/background2.jpeg')


selected = option_menu(
    menu_title=None,
    options = ["Home", "Technology", "Team"],
    # icons = []
    default_index = 0,
    orientation = "horizontal",
)

if selected == "Home":
    st.title(f"You have selected {selected}")
if selected == "Technology":
    st.title(f"You have selected {selected}")
if selected == "Team":
    col1, col2 = st.columns([1,3])

    with col1:
        st.image("./images/dauken.jpeg")

    with col2:
        st.markdown("Dauken Seitkali, CEO")
        st.markdown("Dauken is a data scientist and entrepreneur. Prior to founding AgzaThera, he was a core contributor at Cerebra.ai, a leading stroke detection medtech company in Asia, where he worked as a product owner and data scientist. Dauken believes that success in such interdisciplinary project like AgzaThera depends on people that possess ability to dive deep in both fields and can apply creativity to find points of contact. His primary due as a CEO is to be that person and attract such people. ")

    col3, col4 = st.columns([1,3])

    with col3:
        st.image("./images/mohamad.png")

    with col4:
        st.markdown("Dr. Mohamad Aljofan, Drug Discovery expert")
        st.markdown("Dr. Mohamad Aljofan graduated with a BSc (Pharmaceutical Sciences) (Honors) from the Royal Melbourne Institute of Technology Melbourne (Australia). He has a Master’s Degree in Clinical Pharmacy and a PhD in drug discovery that was conducted at Australia’s highest research institute “The Commonwealth Scientific and Industrial Research Organisation”. Dr. Aljofan did his Postdoctoral training at Prince Henry’s Institute of Medical Research and the School of Medicine, Nursing and Health Sciences at Monash University- Australia. After, he became the Vice Dean and Department Chair of Clinical Pharmacy at the University of Hail in Saudi Arabia where he established his mixed Clinical and Basic research group. His Clinical research focuses on medication adherence and public health and his Basic research interests revolve around drug discovery and development, particularly early drug discovery. Dr. Aljofan received several research grants and funding from different sources including the National Health and Medical Research Council of Australia, Monash University, University of Hail and KACST in Saudi Arabia. He is a recipient of many awards and recognitions including the prestigious Australian NHMRC- Postdoctoral Fellowship, CSIRO-Postgraduate Scholarship and RMIT Undergraduate Scholarship.")
if selected == "Careers":
    st.title(f"You have selected {selected}")

if selected == "Contacts":
    st.title(f"You have selected {selected}")

if selected == "Demo":
    st.title(f"You have selected {selected}")