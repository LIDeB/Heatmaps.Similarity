# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 19:23:28 2021

@author: Lucas
"""

import streamlit as st
import pandas as pd
import base64
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs # para calcular Tanimoto similarity
import plotly.graph_objects as go
from pathlib import Path
import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt


#---------------------------------#
# Page layout
## Page expands to full width
st.set_page_config(page_title='LIDEB Tools - Molecular Similarity Heatmap',
    layout='wide')

######
# Funcion para poner una imagen    
def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded

from PIL import Image
image = Image.open('cropped-header-heatmaps.png')
st.image(image)

st.write("&nbsp[![Website](https://img.shields.io/badge/website-LIDeB-blue)](https://lideb.biol.unlp.edu.ar)&nbsp[![Twitter Follow](https://img.shields.io/twitter/follow/LIDeB_UNLP?style=social)](https://twitter.com/intent/follow?screen_name=LIDeB_UNLP)")
st.subheader(":pushpin:" "About Us")
st.markdown("We are a drug discovery team with an interest in the development of publicly available open-source customizable cheminformatics tools to be used in computer-assisted drug discovery. We belong to the Laboratory of Bioactive Research and Development (LIDeB) of the National University of La Plata (UNLP), Argentina. Our research group is focused on computer-guided drug repurposing and rational discovery of new drug candidates to treat epilepsy and neglected tropical diseases.")



#---------------------------------#
st.write("""
# LIDeB Tools - Similarity Heatmaps

This WebApp builds a heatmap of molecular similarity. These plots of inter-molecular similarity (computed as Tanimoto similarity coefficient using Morgan fingerprints or other molecular fingerprinting systems) allow for a fast, visual inspection of the molecular diversity of the datasets, and also preliminary detection of clusters within a dataset. The resulting plots are downloadable as .png files through a simple right click on your mouse!.
""")


#---------------------------------#
# Sidebar - Collects user input features into dataframe
st.sidebar.header('Upload your set of compounds (SMILES)')

st.sidebar.markdown("Upload your set of compounds using SMILES chemical notation, one compound per line in a .txt file")
st.sidebar.markdown("""
[Example TXT input file](https://raw.githubusercontent.com/Capigol/molecular_similarity_heatmaps/main/molecules_1.txt)
""")
uploaded_file_1 = st.sidebar.file_uploader("Upload your first set of compounds", type=["txt"])
uploaded_file_2 = st.sidebar.file_uploader("Upload your second set of compounds", type=["txt"])

# Sidebar - Specify parameter settings
st.sidebar.header('Morgan FP Radio')
split_size = st.sidebar.slider('Morgan fingerprint Radio', 2, 4, 2, 1)

st.sidebar.subheader('Morgan FP Lengh')
parameter_n_estimators = st.sidebar.slider('Set the fingerprint lenght', 512, 2048, 1024, 512)

similarity_metric = st.sidebar.selectbox("Select the similarity metric", ("TanimotoSimilarity", "DiceSimilarity", "CosineSimilarity", "SokalSimilarity", "RusselSimilarity", "KulczynskiSimilarity", "McConnaugheySimilarity"),0)


st.sidebar.header('Type of Plot')
type_plot = st.sidebar.checkbox('Interactive Plot')
if type_plot == True:
    plotly_color = st.sidebar.selectbox("Select the heatmap color", 
                         ('Blackbody','Bluered','Blues','Earth','Electric','Greens',
                          'Greys','Hot','Jet','Picnic','Portland','Rainbow','RdBu','Reds','Viridis','YlGnBu','YlOrRd'),
                         15)
else:
    sns_color = st.sidebar.selectbox("Select the heatmap color", ("rocket", "mako", "flare","crest","magma","viridis"),5)

st.sidebar.title(":speech_balloon: Contact Us")
st.sidebar.info(
"""
If you are looking to contact us, please
[:e-mail:](mailto:lideb@biol.unlp.edu.ar) or [Twitter](https://twitter.com/LIDeB_UNLP)
""")

#---------------------------------#
# Main panel

# Displays the dataset
st.subheader('Dataset')

#%%

# GENERAMOS LAS FINGERPRINTS DE CADA DATASET POR SEPARADO

def similarity(df_1,df_2):
    df_1 = df_1[0].tolist()
    df_2 = df_2[0].tolist()
    # df_2 = df_1.copy()
    lenght_dataset_1 = len(df_1)
    lenght_dataset_2 = len(df_2)
    st.markdown('Dataset 1 have: ' + str(lenght_dataset_1) + " molecules")
    st.markdown('Dataset 2 have: ' + str(lenght_dataset_2) + " molecules")
    
    fps_1 = []    
    for m in df_1:
        mol = Chem.MolFromSmiles(m)
        fp_1 = AllChem.GetMorganFingerprintAsBitVect(mol,split_size,nBits=parameter_n_estimators,useFeatures=False)
        fps_1.append(fp_1)
    
    fps_2 = []
    for m1 in df_2:
        mole = Chem.MolFromSmiles(m1)
        fp_2 = AllChem.GetMorganFingerprintAsBitVect(mole,split_size,nBits=parameter_n_estimators,useFeatures=False)
        fps_2.append(fp_2)
    
    # COMPARAMOS LAS FINGERPRINTS POR TANIMOTO Y GENERAMOS LA MATRIZ
    
    matriz_tanimoto = []
    for finger1 in fps_1:
        lista = []
        for finger2 in fps_2:
            similarity_metric_ok = getattr(DataStructs, similarity_metric)
            tan_sim_ac_in= similarity_metric_ok(finger1, finger2)
            lista.append(tan_sim_ac_in)
        matriz_tanimoto.append(lista)
    
    df_ok = pd.DataFrame(matriz_tanimoto)
    filas = list(range(1,len(fps_1)+1,1))
    columnas = list(range(1,len(fps_2)+1,1))
    df_ok.index = filas
    df_ok.columns = columnas
    return df_ok


def heatmap(df_ok):
    #-----Plot-----#
    if type_plot == True:
        # color = "YlGnBu"
        fig = go.Figure(go.Heatmap(z=df_ok,x0=1,dx=1, y0=1,dy=1, hoverongaps = False,showscale=True, colorscale=plotly_color,zmax=1,zmin=0))
        fig.update_xaxes(title_text='Dataset 2')
        fig.update_yaxes(title_text='Dataset 1')
        fig.update_layout(margin = dict(t=60,r=20,b=20,l=20),
        width = 800, height = 800,
        autosize = False )
        st.plotly_chart(fig)
    else:
        df_ok.sort_index(axis=0, ascending=False,inplace=True)
        ax = sns.heatmap(df_ok, xticklabels=False, yticklabels=False,cmap=sns_color)
        plt.xlabel("Dataset 2")
        plt.ylabel("Dataset 1")
        plt.savefig('image.png', dpi = 300)       
        st.set_option('deprecation.showPyplotGlobalUse', False)
        st.pyplot()
        return ax

# PARA GENERAR LA MATRIX Y EL HEATMAP

# ---------------------------------#

if uploaded_file_1 is not None and uploaded_file_2 is not None:
    df_1 = pd.read_csv(uploaded_file_1,sep="\t",header=None)
    df_2 = pd.read_csv(uploaded_file_2,sep="\t",header=None)
    df_ok = similarity(df_1,df_2)
    plot = heatmap(df_ok)
    if type_plot == False:
        st.markdown("You can download the heatmap by Right Click in the image and then **'save image as'** :blush: ")
    else:
        st.markdown("You can download the heatmap by clicking on the camera icon at the top of the plot :blush: ")

    
# Example file
else:
    st.info('Awaiting for TXT file to be uploaded.')
    if st.button('Press to use Example Dataset'):
        df_1 = pd.read_csv("molecules_1.txt",sep="\t",header=None)
        df_2 = pd.read_csv("molecules_1.txt",sep="\t",header=None)
        df_ok = similarity(df_1,df_2)    
        plot = heatmap(df_ok)
        st.markdown('A dataset of **40 SMILES** has been used as example.')
        if type_plot == False:
            st.markdown("You can download the heatmap by Right Click in the image and then **'save image as'** :blush: ")
        else:
            st.markdown("You can download the heatmap by clicking on the camera icon at the top of the plot :blush: ")


        
#Footer edit

footer="""<style>
a:link , a:visited{
color: blue;
background-color: transparent;
text-decoration: underline;
}
a:hover,  a:active {
color: red;
background-color: transparent;
text-decoration: underline;
}
.footer {
position: fixed;
left: 0;
bottom: 0;
width: 100%;
background-color: white;
color: black;
text-align: center;
}
</style>
<div class="footer">
<p>Made in  🐍 and <img style='display: ; ' href="https://streamlit.io" src="https://i.imgur.com/iIOA6kU.png" target="_blank"></img> Developed with ❤️ by <a style='display: ; text-align: center' href="https://twitter.com/capigol" target="_blank">Lucas Alberca</a> for <a style='display:; text-align: center;' href="https://lideb.biol.unlp.edu.ar/" target="_blank">LIDeB</a></p>
</div>
"""
st.markdown(footer,unsafe_allow_html=True)


